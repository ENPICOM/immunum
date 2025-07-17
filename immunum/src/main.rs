mod annotation;
mod cli;
mod consensus_scoring;
mod constants;
mod fastx;
mod insertion_numbering;
mod needleman_wunsch;
mod numbering;
mod numbering_scheme_type;
mod schemes;
mod sequence_stream;
mod types;

use clap::Parser;
use cli::Cli;
use sequence_stream::SequenceStream;
use std::fs::File;
use std::io::{self, BufWriter, Write};

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

    // Display parsed arguments to stderr for logging
    eprintln!("Scheme: {:?}", cli.scheme);
    eprintln!("Chains: {:?}", cli.chains);

    // Create sequence stream
    let sequence_stream = match SequenceStream::new(&cli.input) {
        Ok(stream) => stream,
        Err(e) => {
            eprintln!("Error creating sequence stream: {e}");
            std::process::exit(1);
        }
    };

    // Create output writer (either file or stdout)
    let mut output_writer: Box<dyn Write> = match &cli.output {
        Some(output_path) => Box::new(BufWriter::new(File::create(output_path).unwrap())),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    // Process each sequence and write results
    for record_result in sequence_stream {
        match record_result {
            Ok(record) => {
                let numbered_sequence =
                    numbering::number_sequence(&record.sequence, &cli.scheme, chains);
                writeln!(output_writer, "{numbered_sequence}").unwrap();
            }
            Err(e) => {
                eprintln!("Error processing sequence: {e}");
                std::process::exit(1);
            }
        }
    }
}
