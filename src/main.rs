mod annotator;
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
use sequence::SequenceStream;
use std::fs::File;
use std::io::{self, BufWriter, Write};

use crate::annotator::Annotator;
use crate::types::Chain;

const DEFAULT_CHAINS: [Chain; 3] = [Chain::IGH, Chain::IGK, Chain::IGL];

fn main() {
    let cli = Cli::parse();

    // Populate the `chains` field if it's empty
    let default_chains_vec = DEFAULT_CHAINS.to_vec();
    let chains = cli.chains.as_ref().unwrap_or(&default_chains_vec);

    // Create annotator with custom parameters
    let annotator = match Annotator::new(cli.scheme, chains.clone(), cli.prefilter) {
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
    let mut output_writer: Box<dyn Write + Send> = match &cli.output {
        Some(output_path) => Box::new(BufWriter::new(File::create(output_path).unwrap())),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    // ?NOTE Loading everything into memory first might not be optimal for huge files
    // Collect all sequences first for parallel processing
    let sequences: Result<Vec<_>, _> = sequence_stream.collect();
    let sequences = match sequences {
        Ok(sequences) => sequences,
        Err(e) => {
            eprintln!("Error collecting sequences: {e}");
            std::process::exit(1);
        }
    };

    for (header, annotation_result) in annotator.number(sequences, Some(cli.max_chains)).iter() {
        match annotation_result {
            // TODO improve the json output
            Ok(annotation_result) => {
                for chain_numbering in annotation_result.iter() {
                    // println!("{}", chain_numbering.to_json_string(header.to_owned()))
                    writeln!(output_writer, "{}", chain_numbering.to_json_string(header.to_owned())).unwrap();
                }
            }
            Err(err) => eprintln!("Error numbering sequence: {}", err),
        }
    }
}
