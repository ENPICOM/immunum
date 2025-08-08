mod annotator;
mod cli;
mod consensus_scoring;
mod constants;
mod sequence_io;
mod gap_penalty;
mod insertion_naming;
mod needleman_wunsch;
mod numbering_scheme;
mod prefiltering;
mod result;
mod schemes;
mod scoring_matrix;
mod types;

use clap::Parser;
use cli::Cli;
use sequence_io::SequenceStream;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::annotator::Annotator;
use crate::constants::{get_scoring_params, ScoringParams};
use crate::result::OutputFormat;
use crate::types::Chain;

const DEFAULT_CHAINS: [Chain; 7] = [
    Chain::IGH,
    Chain::IGK,
    Chain::IGL,
    Chain::TRA,
    Chain::TRB,
    Chain::TRG,
    Chain::TRD,
];

fn main() {
    let cli = Cli::parse();
    // Determine chains and log args
    let chains: Vec<Chain> = cli
        .chains
        .clone()
        .unwrap_or_else(|| DEFAULT_CHAINS.to_vec());
    log_cli(&cli, &chains);

    // Build annotator
    let scoring_params = build_scoring_params_from_cli(&cli);
    let annotator = match build_annotator(&cli, chains, scoring_params) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error creating annotator: {}", e);
            std::process::exit(1);
        }
    };

    // Create output writer (either file or stdout)
    let mut output_writer = open_output_writer(&cli.output);

    // If parallel requested and input is a file, use batch processing
    if cli.parallel && Path::new(&cli.input).exists() {
        if let Err(e) = process_file(&mut *output_writer, &annotator, &cli.input, cli.all_chains, cli.format.clone()) {
            eprintln!("Error processing file in parallel: {}", e);
            std::process::exit(1);
        }
        return;
    }

    // Fallback: streaming (sequential) processing
    if let Err(e) = process_stream(&mut *output_writer, &annotator, &cli.input, cli.all_chains, cli.format.clone()) {
        eprintln!("{}", e);
        std::process::exit(1);
    }
}

fn log_cli(cli: &Cli, chains: &Vec<Chain>) {
    eprintln!("Scheme: {:?}", cli.scheme);
    eprintln!("Chains: {:?}", chains);
    eprintln!("Pre-filtering: {}", cli.prefilter);
    eprintln!("All-chains mode: {}", cli.all_chains);
    eprintln!("Parallel: {}", cli.parallel);
    eprintln!("Output format: {:?}", cli.format);
}

fn build_scoring_params_from_cli(cli: &Cli) -> Option<ScoringParams> {
    if cli.gap_pen_cp.is_some()
        || cli.gap_pen_fr.is_some()
        || cli.gap_pen_ip.is_some()
        || cli.gap_pen_op.is_some()
        || cli.gap_pen_cdr.is_some()
        || cli.gap_pen_other.is_some()
        || cli.cdr_increase.is_some()
        || cli.pen_leap_insertion_point_imgt.is_some()
        || cli.pen_leap_insertion_point_kabat.is_some()
    {
        let mut params = get_scoring_params();
        if let Some(val) = cli.gap_pen_cp {
            params.gap_pen_cp = val;
        }
        if let Some(val) = cli.gap_pen_fr {
            params.gap_pen_fr = val;
        }
        if let Some(val) = cli.gap_pen_ip {
            params.gap_pen_ip = val;
        }
        if let Some(val) = cli.gap_pen_op {
            params.gap_pen_op = val;
        }
        if let Some(val) = cli.gap_pen_cdr {
            params.gap_pen_cdr = val;
        }
        if let Some(val) = cli.gap_pen_other {
            params.gap_pen_other = val;
        }
        if let Some(val) = cli.cdr_increase {
            params.cdr_increase = val;
        }
        if let Some(val) = cli.pen_leap_insertion_point_imgt {
            params.pen_leap_insertion_point_imgt = val;
        }
        if let Some(val) = cli.pen_leap_insertion_point_kabat {
            params.pen_leap_insertion_point_kabat = val;
        }
        Some(params)
    } else {
        None
    }
}

fn build_annotator(cli: &Cli, chains: Vec<Chain>, scoring_params: Option<ScoringParams>) -> Result<Annotator, String> {
    Annotator::new(cli.scheme, chains, scoring_params, Some(cli.prefilter))
}

fn open_output_writer(output: &Option<String>) -> Box<dyn Write> {
    match output {
        Some(output_path) => Box::new(BufWriter::new(File::create(output_path).unwrap())),
        None => Box::new(BufWriter::new(io::stdout())),
    }
}

fn emit_formatted(
    writer: &mut dyn Write,
    base_name: &str,
    idx: usize,
    all_chains: bool,
    formatted: String,
) -> io::Result<()> {
    let header = if all_chains && idx > 0 {
        format!("{}_chain_{}", base_name, idx + 1)
    } else {
        base_name.to_string()
    };
    writeln!(writer, "# {}", header)?;
    writeln!(writer, "{}", formatted)?;
    Ok(())
}

fn process_file(
    writer: &mut dyn Write,
    annotator: &Annotator,
    input_path: &str,
    all_chains: bool,
    format: OutputFormat,
) -> Result<(), String> {
    let results = annotator.number_file(input_path, all_chains, true)?;
    for (name, multi) in results {
        for (idx, res) in multi.into_iter().enumerate() {
            match res {
                Ok(ar) => {
                    let s = ar.to_string(format.clone());
                    emit_formatted(writer, &name, idx, all_chains, s).map_err(|e| e.to_string())?;
                }
                Err(e) => eprintln!("Error numbering sequence '{}': {}", name, e),
            }
        }
    }
    Ok(())
}

fn process_stream(
    writer: &mut dyn Write,
    annotator: &Annotator,
    input: &str,
    all_chains: bool,
    format: OutputFormat,
) -> Result<(), String> {
    let sequence_stream = SequenceStream::new(input)
        .map_err(|e| format!("Error creating sequence stream: {}", e))?;

    for record_result in sequence_stream {
        let record = match record_result {
            Ok(r) => r,
            Err(e) => return Err(format!("Error processing sequence: {}", e)),
        };

        let results = annotator.number_sequence(&record.sequence, all_chains);
        if results.is_empty() {
            eprintln!("No chains found in sequence '{}'", record._name);
            continue;
        }

        for (idx, res) in results.into_iter().enumerate() {
            match res {
                Ok(ar) => {
                    let s = ar.to_string(format.clone());
                    emit_formatted(writer, &record._name, idx, all_chains, s)
                        .map_err(|e| e.to_string())?;
                }
                Err(e) => eprintln!(
                    "Error numbering chain {} in sequence '{}': {}",
                    idx + 1,
                    record._name,
                    e
                ),
            }
        }
    }

    Ok(())
}
