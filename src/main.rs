use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

use clap::{Args, Parser, Subcommand};
use immunum::{Annotator, Chain, Scheme};

#[derive(Parser)]
#[command(name = "immunum", about = "Immune receptor sequence numbering")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Number immune receptor sequences
    Number(NumberArgs),
}

#[derive(Args)]
struct NumberArgs {
    /// Input: sequence string, FASTA file path, or "-" for stdin (default)
    input: Option<String>,
    /// Output file path (default: stdout)
    output: Option<String>,
    /// Numbering scheme: imgt (i), kabat (k)
    #[arg(short, long, default_value = "imgt")]
    scheme: String,
    /// Chain filter: h,k,l,a,b,g,d or groups: ig, tcr, all
    #[arg(short, long, default_value = "all")]
    chain: String,
    /// Output format: tsv, json, jsonl
    #[arg(short, long, default_value = "tsv")]
    format: String,
}

// -- Input handling --

enum InputKind {
    Stdin,
    File(String),
    Sequence(String),
}

struct Record {
    id: String,
    sequence: String,
}

fn detect_input(input: Option<&str>) -> InputKind {
    match input {
        None | Some("-") => InputKind::Stdin,
        Some(s) => {
            let path = Path::new(s);
            if path.extension().is_some() || path.exists() {
                InputKind::File(s.to_string())
            } else {
                InputKind::Sequence(s.to_string())
            }
        }
    }
}

fn read_fasta(reader: impl BufRead) -> Vec<Record> {
    let mut records = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };
        let line = line.trim_end();
        if line.starts_with('>') {
            if !current_id.is_empty() && !current_seq.is_empty() {
                records.push(Record {
                    id: current_id,
                    sequence: current_seq,
                });
                current_seq = String::new();
            }
            current_id = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
        } else if !line.is_empty() {
            current_seq.push_str(line);
        }
    }
    if !current_id.is_empty() && !current_seq.is_empty() {
        records.push(Record {
            id: current_id,
            sequence: current_seq,
        });
    }
    records
}

fn read_input(input: Option<&str>) -> Result<Vec<Record>, String> {
    match detect_input(input) {
        InputKind::Stdin => {
            let stdin = io::stdin();
            let reader = BufReader::new(stdin.lock());
            Ok(read_fasta(reader))
        }
        InputKind::File(path) => {
            let file =
                File::open(&path).map_err(|e| format!("cannot open '{}': {}", path, e))?;
            let reader = BufReader::new(file);
            Ok(read_fasta(reader))
        }
        InputKind::Sequence(seq) => Ok(vec![Record {
            id: "seq_1".to_string(),
            sequence: seq,
        }]),
    }
}

// -- Parsing helpers --

fn parse_scheme(s: &str) -> Result<Scheme, String> {
    Scheme::from_str(s).map_err(|_| format!("unknown scheme '{}' (options: imgt, kabat)", s))
}

fn parse_chains(s: &str) -> Result<Vec<Chain>, String> {
    let s = s.to_lowercase();
    match s.as_str() {
        "all" => Ok(vec![
            Chain::IGH,
            Chain::IGK,
            Chain::IGL,
            Chain::TRA,
            Chain::TRB,
            Chain::TRG,
            Chain::TRD,
        ]),
        "ig" => Ok(vec![Chain::IGH, Chain::IGK, Chain::IGL]),
        "tcr" => Ok(vec![Chain::TRA, Chain::TRB, Chain::TRG, Chain::TRD]),
        _ => s
            .split(',')
            .map(|c| {
                Chain::from_str(c.trim())
                    .map_err(|_| format!("unknown chain '{}' (options: h,k,l,a,b,g,d,ig,tcr,all)", c.trim()))
            })
            .collect(),
    }
}

#[derive(Clone, Copy)]
enum OutputFormat {
    Tsv,
    Json,
    Jsonl,
}

fn parse_format(s: &str) -> Result<OutputFormat, String> {
    match s.to_lowercase().as_str() {
        "tsv" => Ok(OutputFormat::Tsv),
        "json" => Ok(OutputFormat::Json),
        "jsonl" => Ok(OutputFormat::Jsonl),
        _ => Err(format!("unknown format '{}' (options: tsv, json, jsonl)", s)),
    }
}

// -- Output formatting --

struct NumberedRecord {
    id: String,
    chain: Chain,
    scheme: Scheme,
    numbering: Vec<(String, char)>, // (position_label, residue)
}

fn write_tsv(writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
    writeln!(writer, "sequence_id\tchain\tscheme\tposition\tresidue")?;
    for rec in records {
        for (pos, res) in &rec.numbering {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                rec.id, rec.chain, rec.scheme, pos, res
            )?;
        }
    }
    Ok(())
}

fn write_json(writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
    let json_records: Vec<serde_json::Value> = records.iter().map(record_to_json).collect();
    serde_json::to_writer_pretty(&mut *writer, &json_records)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    writeln!(writer)?;
    Ok(())
}

fn write_jsonl(writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
    for rec in records {
        let json = record_to_json(rec);
        serde_json::to_writer(&mut *writer, &json)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        writeln!(writer)?;
    }
    Ok(())
}

fn record_to_json(rec: &NumberedRecord) -> serde_json::Value {
    let numbering: serde_json::Map<String, serde_json::Value> = rec
        .numbering
        .iter()
        .map(|(pos, res)| (pos.clone(), serde_json::Value::String(res.to_string())))
        .collect();

    serde_json::json!({
        "sequence_id": rec.id,
        "chain": rec.chain.to_string(),
        "scheme": rec.scheme.to_string(),
        "numbering": numbering,
    })
}

// -- Main logic --

fn run_number(args: &NumberArgs) -> Result<(), String> {
    let scheme = parse_scheme(&args.scheme)?;
    let chains = parse_chains(&args.chain)?;
    let format = parse_format(&args.format)?;

    let annotator =
        Annotator::new(&chains, scheme).map_err(|e| format!("failed to create annotator: {}", e))?;

    let input_records = read_input(args.input.as_deref())?;

    let mut numbered: Vec<NumberedRecord> = Vec::new();
    for rec in &input_records {
        let result = annotator
            .annotate(&rec.sequence)
            .map_err(|e| format!("failed to annotate '{}': {}", rec.id, e))?;

        let positions = result.numbering(scheme);
        let seq_chars: Vec<char> = rec.sequence.chars().collect();

        let numbering: Vec<(String, char)> = positions
            .iter()
            .zip(seq_chars.iter())
            .map(|(pos, &ch)| (pos.to_string(), ch))
            .collect();

        numbered.push(NumberedRecord {
            id: rec.id.clone(),
            chain: result.chain,
            scheme,
            numbering,
        });
    }

    let stdout = io::stdout();
    let mut writer: BufWriter<Box<dyn Write>> = match &args.output {
        Some(path) => {
            let file =
                File::create(path).map_err(|e| format!("cannot create '{}': {}", path, e))?;
            BufWriter::new(Box::new(file))
        }
        None => BufWriter::new(Box::new(stdout.lock())),
    };

    match format {
        OutputFormat::Tsv => write_tsv(&mut writer, &numbered),
        OutputFormat::Json => write_json(&mut writer, &numbered),
        OutputFormat::Jsonl => write_jsonl(&mut writer, &numbered),
    }
    .map_err(|e| format!("write error: {}", e))?;

    Ok(())
}

fn main() {
    let cli = Cli::parse();

    let result = match &cli.command {
        Commands::Number(args) => run_number(args),
    };

    if let Err(e) = result {
        eprintln!("error: {}", e);
        std::process::exit(1);
    }
}
