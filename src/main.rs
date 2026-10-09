use std::fs::File;
use std::io::{self, BufWriter, IsTerminal, Write};
use std::str::FromStr;

use clap::{Args, Parser, Subcommand};
use immunum::{
    Annotator, NumberedRecord, OutputFormat, Record, SegmentedRecord, DEFAULT_MIN_CONFIDENCE,
};

#[derive(Parser)]
#[command(name = "immunum", about = "Immune receptor sequence numbering")]
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
    /// Minimum confidence threshold (0.0-1.0). Sequences below this are skipped.
    #[arg(long, default_value_t = DEFAULT_MIN_CONFIDENCE)]
    min_confidence: f32,
    /// Process every variable domain in each sequence (e.g. both domains of an scFv), one record
    /// per domain with its 0-based `domain` index. A sequence without a domain gets no record.
    #[arg(long)]
    all_domains: bool,
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
        let records = immunum::read_input(args.input.as_deref())?;
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

fn write_error(e: io::Error) -> String {
    format!("write error: {}", e)
}

fn run_number(args: &AnnotateArgs) -> Result<(), String> {
    let Run {
        annotator,
        format,
        records,
        mut writer,
    } = Run::new(args, "number")?;
    format
        .write_header(&mut writer, args.all_domains)
        .map_err(write_error)?;

    let mut written = 0;
    for rec in records {
        let numbered = if args.all_domains {
            match annotator.number_domains(&rec.sequence) {
                Ok(results) => results
                    .into_iter()
                    .enumerate()
                    .map(|(domain, result)| {
                        NumberedRecord::success(rec.id.clone(), rec.sequence.clone(), result)
                            .in_domain(domain)
                    })
                    .collect(),
                Err(e) => vec![NumberedRecord::failure(rec.id, rec.sequence, e.to_string())],
            }
        } else {
            vec![match annotator.number(&rec.sequence) {
                Ok(result) => NumberedRecord::success(rec.id, rec.sequence, result),
                Err(e) => NumberedRecord::failure(rec.id, rec.sequence, e.to_string()),
            }]
        };
        for numbered in &numbered {
            format
                .write_record(&mut writer, numbered, written, args.all_domains)
                .map_err(write_error)?;
            written += 1;
        }
    }

    format.write_footer(&mut writer).map_err(write_error)
}

fn run_segment(args: &AnnotateArgs) -> Result<(), String> {
    let Run {
        annotator,
        format,
        records,
        mut writer,
    } = Run::new(args, "segment")?;
    format
        .write_segment_header(&mut writer, args.all_domains)
        .map_err(write_error)?;

    let mut written = 0;
    for rec in records {
        let segmented = if args.all_domains {
            match annotator.segment_domains(&rec.sequence) {
                Ok(results) => results
                    .into_iter()
                    .enumerate()
                    .map(|(domain, result)| {
                        SegmentedRecord::success(rec.id.clone(), result).in_domain(domain)
                    })
                    .collect(),
                Err(e) => vec![SegmentedRecord::failure(rec.id, e.to_string())],
            }
        } else {
            vec![match annotator.segment(&rec.sequence) {
                Ok(result) => SegmentedRecord::success(rec.id, result),
                Err(e) => SegmentedRecord::failure(rec.id, e.to_string()),
            }]
        };
        for segmented in &segmented {
            format
                .write_segment_record(&mut writer, segmented, written, args.all_domains)
                .map_err(write_error)?;
            written += 1;
        }
    }

    format.write_footer(&mut writer).map_err(write_error)
}

fn main() {
    let cli = Cli::parse();
    let result = match &cli.command {
        Commands::Number(args) => run_number(args),
        Commands::Segment(args) => run_segment(args),
    };
    if let Err(e) = result {
        eprintln!("error: {}", e);
        std::process::exit(1);
    }
}
