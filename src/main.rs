use std::fs::File;
use std::io::{self, BufWriter, IsTerminal, Write};
use std::str::FromStr;

use clap::{Args, Parser, Subcommand};
use immunum::{Annotator, Chain, NumberedRecord, OutputFormat, Scheme};

#[derive(Parser)]
#[command(name = "immunum", about = "Immune receptor sequence numbering")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Number sequences using a specified scheme and chain type
    Number(NumberArgs),
}

#[derive(Args)]
struct NumberArgs {
    /// Input sequence string, FASTA file path, or "-" for stdin (default)
    input: Option<String>,
    /// Output file path (default: stdout)
    output: Option<String>,
    /// Numbering scheme: imgt (i), kabat (k)
    #[arg(short, long, default_value = "imgt")]
    scheme: String,
    /// Chain filter: h,k,l,a,b,g,d or groups: ig, tcr, all
    #[arg(short, long, default_value = "ig")]
    chain: String,
    /// Output format: tsv, json, jsonl
    #[arg(short, long, default_value = "tsv")]
    format: String,
}

fn run_number(args: &NumberArgs) -> Result<(), String> {
    if args.input.is_none() && io::stdin().is_terminal() {
        return Err(
            "no input provided. Usage: immunum number <sequence|file.fasta> or pipe via stdin"
                .to_string(),
        );
    }

    let scheme =
        Scheme::from_str(&args.scheme).map_err(|_| format!("unknown scheme '{}'", args.scheme))?;
    let chains = Chain::parse_chain_spec(&args.chain).map_err(|e| e.to_string())?;
    let format = OutputFormat::from_str(&args.format)?;
    let mut annotator = Annotator::new(&chains, scheme)
        .map_err(|e| format!("failed to create annotator: {}", e))?;

    let records = immunum::read_input(args.input.as_deref())?;

    let stdout = io::stdout();
    let mut writer: BufWriter<Box<dyn Write>> = match &args.output {
        Some(path) => {
            let file =
                File::create(path).map_err(|e| format!("cannot create '{}': {}", path, e))?;
            BufWriter::new(Box::new(file))
        }
        None => BufWriter::new(Box::new(stdout.lock())),
    };

    format
        .write_header(&mut writer)
        .map_err(|e| format!("write error: {}", e))?;

    for (i, rec) in records.into_iter().enumerate() {
        let result = annotator
            .number(&rec.sequence)
            .map_err(|e| format!("failed to number '{}': {}", rec.id, e))?;
        let numbered = NumberedRecord {
            id: rec.id,
            sequence: rec.sequence,
            result,
        };
        format
            .write_record(&mut writer, &numbered, i)
            .map_err(|e| format!("write error: {}", e))?;
    }

    format
        .write_footer(&mut writer)
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
