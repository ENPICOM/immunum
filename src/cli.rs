//! The `immunum` command-line tool, shared by the native binary and the `immunum` command the
//! Python package installs.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufWriter, IsTerminal, Write};
use std::str::FromStr;

use clap::{Args, Parser, Subcommand};

use crate::{
    per_domain, read_input, Annotator, NumberedRecord, OutputFormat, Record, SegmentedRecord,
    DEFAULT_MIN_CONFIDENCE,
};

#[derive(Parser)]
#[command(name = "immunum", version, about = "Immune receptor sequence numbering")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Number sequences using a specified scheme and chain type
    Number(AnnotateArgs),
    /// Split sequences into FR/CDR regions, one row per sequence
    Segment(AnnotateArgs),
}

#[derive(Args)]
struct AnnotateArgs {
    /// Input sequence string, FASTA file path, or "-" for stdin (default)
    input: Option<String>,
    /// Output file path (default: stdout)
    output: Option<String>,
    /// Numbering scheme: imgt (i), kabat (k), chothia (c), martin (m), aho (a)
    #[arg(short, long, default_value = "imgt")]
    scheme: String,
    /// Chain filter: h,k,l,a,b,g,d or groups: ig, tcr, all
    #[arg(short, long, default_value = "ig")]
    chain: String,
    /// Output format: tsv, json, jsonl
    #[arg(short, long, default_value = "tsv")]
    format: String,
    /// Minimum confidence threshold (0.0-1.0). A sequence below it gets an error record with
    /// `error_kind` `low_confidence`.
    #[arg(long, default_value_t = DEFAULT_MIN_CONFIDENCE)]
    min_confidence: f32,
    /// Process every variable domain in each sequence (e.g. both domains of an scFv), one record
    /// per domain with its 0-based `domain` index. A sequence without a domain gets one error record.
    #[arg(long)]
    all_domains: bool,
}

/// Runs the command line `args`, whose first item is the program name, and returns the exit
/// status: 0 on success, even when some sequences couldn't be numbered (they get error records);
/// 1 for a setup or input error, printed to stderr; 2 for a usage error.
pub fn run<I, T>(args: I) -> i32
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let cli = match Cli::try_parse_from(args) {
        Ok(cli) => cli,
        Err(e) => {
            // Help goes to stdout with status 0; a usage error to stderr with status 2
            let _ = e.print();
            return e.exit_code();
        }
    };
    let result = match &cli.command {
        Commands::Number(args) => run_number(args),
        Commands::Segment(args) => run_segment(args),
    };
    match result {
        Ok(()) => 0,
        // The reader stopped early, as `immunum number seqs.fasta | head` does: nothing is wrong
        Err(Failure::BrokenPipe) => 0,
        Err(Failure::Message(e)) => {
            eprintln!("error: {}", e);
            1
        }
    }
}

enum Failure {
    Message(String),
    BrokenPipe,
}

impl From<String> for Failure {
    fn from(message: String) -> Self {
        Failure::Message(message)
    }
}

fn write_error(e: io::Error) -> Failure {
    match e.kind() {
        io::ErrorKind::BrokenPipe => Failure::BrokenPipe,
        _ => Failure::Message(format!("write error: {}", e)),
    }
}

// What every command needs: the annotator, the output format, the input records and the output
struct Run {
    annotator: Annotator,
    format: OutputFormat,
    records: Vec<Record>,
    writer: BufWriter<Box<dyn Write>>,
}

impl Run {
    fn new(args: &AnnotateArgs, command: &str) -> Result<Self, String> {
        if args.input.is_none() && io::stdin().is_terminal() {
            return Err(format!(
                "no input provided. Usage: immunum {command} <sequence|file.fasta> or pipe via stdin"
            ));
        }

        let annotator = Annotator::from_names(
            args.chain.split(','),
            &args.scheme,
            Some(args.min_confidence),
        )
        .map_err(|e| e.to_string())?;
        let format = OutputFormat::from_str(&args.format)?;
        let records = read_input(args.input.as_deref())?;
        let writer: BufWriter<Box<dyn Write>> = match &args.output {
            Some(path) => {
                let file =
                    File::create(path).map_err(|e| format!("cannot create '{}': {}", path, e))?;
                BufWriter::new(Box::new(file))
            }
            None => BufWriter::new(Box::new(io::stdout().lock())),
        };
        Ok(Self {
            annotator,
            format,
            records,
            writer,
        })
    }
}

fn run_number(args: &AnnotateArgs) -> Result<(), Failure> {
    let Run {
        annotator,
        format,
        records,
        mut writer,
    } = Run::new(args, "number")?;
    format
        .write_header(&mut writer, args.all_domains)
        .map_err(write_error)?;

    let total = records.len();
    let mut written = 0;
    let mut failed = 0;
    for rec in records {
        let numbered: Vec<_> = if args.all_domains {
            per_domain(annotator.number_domains(&rec.sequence))
                .into_iter()
                .enumerate()
                .map(|(domain, result)| {
                    let numbered =
                        NumberedRecord::new(rec.id.clone(), rec.sequence.clone(), result);
                    match numbered.result {
                        Ok(_) => numbered.in_domain(domain),
                        Err(_) => numbered,
                    }
                })
                .collect()
        } else {
            let result = annotator.number(&rec.sequence);
            vec![NumberedRecord::new(rec.id, rec.sequence, result)]
        };
        if numbered.iter().any(|r| r.result.is_err()) {
            failed += 1;
        }
        for numbered in &numbered {
            format
                .write_record(&mut writer, numbered, written, args.all_domains)
                .map_err(write_error)?;
            written += 1;
        }
    }

    format.write_footer(&mut writer).map_err(write_error)?;
    writer.flush().map_err(write_error)?;
    report_failures(failed, total);
    Ok(())
}

fn run_segment(args: &AnnotateArgs) -> Result<(), Failure> {
    let Run {
        annotator,
        format,
        records,
        mut writer,
    } = Run::new(args, "segment")?;
    format
        .write_segment_header(&mut writer, args.all_domains)
        .map_err(write_error)?;

    let total = records.len();
    let mut written = 0;
    let mut failed = 0;
    for rec in records {
        let segmented: Vec<_> = if args.all_domains {
            per_domain(annotator.segment_domains(&rec.sequence))
                .into_iter()
                .enumerate()
                .map(|(domain, result)| {
                    let segmented = SegmentedRecord::new(rec.id.clone(), result);
                    match segmented.result {
                        Ok(_) => segmented.in_domain(domain),
                        Err(_) => segmented,
                    }
                })
                .collect()
        } else {
            vec![SegmentedRecord::new(
                rec.id,
                annotator.segment(&rec.sequence),
            )]
        };
        if segmented.iter().any(|r| r.result.is_err()) {
            failed += 1;
        }
        for segmented in &segmented {
            format
                .write_segment_record(&mut writer, segmented, written, args.all_domains)
                .map_err(write_error)?;
            written += 1;
        }
    }

    format.write_footer(&mut writer).map_err(write_error)?;
    writer.flush().map_err(write_error)?;
    report_failures(failed, total);
    Ok(())
}

// Per-sequence errors don't stop the batch or change the exit code, but shouldn't go unnoticed
fn report_failures(failed: usize, total: usize) {
    if failed > 0 {
        eprintln!(
            "warning: {failed} of {total} sequences had errors; see the error and error_kind columns"
        );
    }
}
