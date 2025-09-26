mod annotator;
mod buffer_pool;
mod cli;
mod consensus_scoring;
mod constants;
mod insertion_numbering;
mod needleman_wunsch;
mod numbering_scheme_type;
mod prefiltering;
mod schemes;
mod scoring_matrix;
mod sequence;
mod types;

use clap::Parser;
use cli::Cli;
use rayon::prelude::*;
use sequence::SequenceStream;
use std::fs::File;
use std::io::BufWriter;

use crate::annotator::Annotator;
use crate::types::{Chain, SequenceResult};

const DEFAULT_CHAINS: [Chain; 3] = [Chain::IGH, Chain::IGK, Chain::IGL];

fn main() {
    let cli = Cli::parse();

    // Populate the `chains` field if it's empty
    let default_chains_vec = DEFAULT_CHAINS.to_vec();
    let chains = cli.chains.as_ref().unwrap_or(&default_chains_vec);

    // Create annotator with custom parameters
    let annotator = match Annotator::new(
        cli.scheme,
        chains.clone(),
        cli.cdr_definitions,
        cli.disable_prefiltering,
        Some(cli.min_confidence),
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
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    // Create output file writer
    let mut output_writer = BufWriter::new(File::create(&cli.output).unwrap());

    // Collect all sequences first for parallel processing
    let sequences: Result<Vec<_>, _> = sequence_stream.collect();
    let sequences = match sequences {
        Ok(sequences) => {
            eprintln!("Loaded {} sequences for processing", sequences.len());
            sequences
        }
        Err(e) => {
            eprintln!("Error: Failed to collect sequences - {}", e);
            std::process::exit(1);
        }
    };

    let total_sequences = sequences.len();

    if !cli.verbose && total_sequences > 10 {
        eprintln!("Processing sequences...");
    }

    // Process sequences in parallel and collect results
    let sequence_results: Vec<SequenceResult> = sequences
        .par_iter()
        .map(|sequence| {
            let header = String::from_utf8_lossy(&sequence.name).to_string();
            let annotation_result = annotator.number_sequence(sequence, Some(cli.max_chains));

            match annotation_result {
                Ok(chains) => {
                    if cli.verbose {
                        if !chains.is_empty() {
                            eprintln!(
                                "  Processed '{}': found {} {}",
                                header,
                                chains.len(),
                                if chains.len() == 1 { "chain" } else { "chains" }
                            );
                        } else {
                            eprintln!("  Processed '{}': no chains found", header);
                        }
                    }
                    SequenceResult::new(header, chains)
                }
                Err(err) => {
                    if cli.verbose {
                        eprintln!("  Error processing '{}': {}", header, err);
                    }
                    SequenceResult::new(header, Vec::new())
                }
            }
        })
        .collect();

    // Calculate summary statistics
    let total_chains: usize = sequence_results
        .iter()
        .map(|result| result.chains.len())
        .sum();
    let successful_sequences = sequence_results
        .iter()
        .filter(|result| !result.chains.is_empty())
        .count();

    // Write clean JSON using serde
    serde_json::to_writer_pretty(&mut output_writer, &sequence_results).unwrap();

    eprintln!();
    eprintln!("Done!");
    eprintln!("* Sequences processed: {}", total_sequences);
    eprintln!("* Total chains found: {}", total_chains);
    eprintln!("* Successful sequences: {}", successful_sequences);
    eprintln!("* Results written to: {}", cli.output);
}
