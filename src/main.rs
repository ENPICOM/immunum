mod annotation;
mod annotator;
mod cli;
mod consensus_scoring;
mod constants;
mod sequence;
mod insertion_numbering;
mod needleman_wunsch;
mod numbering_scheme_type;
mod prefiltering;
mod result;
mod schemes;
mod scoring_matrix;
mod types;

use clap::Parser;
use cli::Cli;
use rayon::prelude::*;
use sequence::SequenceStream;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::{Arc, Mutex};

use crate::annotator::Annotator;
use crate::constants::get_scoring_params;
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

    // Populate the `chains` field if it's empty
    let default_chains_vec = DEFAULT_CHAINS.to_vec();
    let chains = cli.chains.as_ref().unwrap_or(&default_chains_vec);

    // Build custom scoring parameters if any are provided
    let scoring_params = if cli.gap_pen_cp.is_some()
        || cli.gap_pen_fr.is_some()
        || cli.gap_pen_ip.is_some()
        || cli.gap_pen_op.is_some()
        || cli.gap_pen_cdr.is_some()
        || cli.gap_pen_other.is_some()
        || cli.cdr_increase.is_some()
        || cli.pen_leap_insertion_point_imgt.is_some()
        || cli.pen_leap_insertion_point_kabat.is_some()
    {
        // Start with defaults and override with provided values
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
    };

    // Create annotator with custom parameters
    let annotator = match Annotator::new(
        cli.scheme,
        chains.clone(),
        scoring_params,
        Some(cli.prefilter),
    ) {
        Ok(annotator) => annotator,
        Err(e) => {
            eprintln!("Error creating annotator: {}", e);
            std::process::exit(1);
        }
    };

    // Create sequence stream
    let sequence_stream = match SequenceStream::new(&cli.input) {
        Ok(stream) => stream,
        Err(e) => {
            eprintln!("Error creating sequence stream: {e}");
            std::process::exit(1);
        }
    };

    // Set thread pool size if specified
    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
    }

    // Create output writer (either file or stdout)
    let output_writer: Box<dyn Write + Send> = match &cli.output {
        Some(output_path) => Box::new(BufWriter::new(File::create(output_path).unwrap())),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let output_writer = Arc::new(Mutex::new(output_writer));

    // Collect all sequences first for parallel processing
    let sequences: Result<Vec<_>, _> = sequence_stream.collect();
    let sequences = match sequences {
        Ok(sequences) => sequences,
        Err(e) => {
            eprintln!("Error collecting sequences: {e}");
            std::process::exit(1);
        }
    };

    // Process sequences in parallel
    sequences.par_iter().for_each(|record| {
        let output_writer = Arc::clone(&output_writer);
        if cli.paired {
            // Use paired sequence numbering to find multiple chains
            let results = annotator.number_paired_sequence_with_id(&record.sequence, record.name.clone());
            if results.is_empty() {
                eprintln!("No chains found in sequence '{}'", record.name);
                return;
            }

            for (chain_index, result) in results.into_iter().enumerate() {
                match result {
                    Ok(annotation_result) => {
                        // Add chain index to sequence name for multiple results
                        let _chain_name = if chain_index > 0 {
                            format!("{}_chain_{}", record.name, chain_index + 1)
                        } else {
                            record.name.clone()
                        };
                        // Write to output with thread-safe access
                        let mut writer = output_writer.lock().unwrap();
                        writeln!(
                            writer,
                            "{}",
                            annotation_result.to_string(cli.format.clone())
                        )
                        .unwrap();
                    }
                    Err(e) => {
                        eprintln!(
                            "Error numbering chain {} in sequence '{}': {}",
                            chain_index + 1,
                            record.name,
                            e
                        );
                    }
                }
            }
        } else {
            // Use single sequence numbering (original behavior)
            match annotator.number_sequence_with_id(&record.sequence, record.name.clone()) {
                Ok(result) => {
                    let mut writer = output_writer.lock().unwrap();
                    writeln!(writer, "{}", result.to_string(cli.format.clone()))
                        .unwrap();
                }
                Err(e) => {
                    eprintln!("Error numbering sequence '{}': {}", record.name, e);
                }
            }
        }
    });
}
