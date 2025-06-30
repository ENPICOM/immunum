mod cli;
mod fastx;
mod numbering;
mod types;

use clap::Parser;
use cli::Cli;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Validates if a string looks like a biological sequence
fn is_valid_sequence(input: &str) -> bool {
    // Check if input string contains valid protein chars
    input.chars().all(|c| "ACDEFGHIKLMNPQRSTVWY".contains(c))
}

fn main() {
    let cli = Cli::parse();

    // Display parsed arguments to stderr for logging
    eprintln!("Scheme: {:?}", cli.scheme);
    eprintln!("Chains: {:?}", cli.chains);

    let input_path = Path::new(&cli.input);

    // Create sequence stream
    let sequence_stream: Box<dyn Iterator<Item = Result<fastx::FastxRecord, fastx::FastxError>>> = 
        if input_path.exists() {
            // Input is a file path - use the fastx parser
            match fastx::from_path(&cli.input) {
                Ok(reader) => Box::new(reader),
                Err(e) => {
                    eprintln!("Error reading file: {}", e);
                    std::process::exit(1);
                }
            }
        } else if is_valid_sequence(&cli.input) {
            // Input is a valid sequence string - create a single record stream
            let record = fastx::FastxRecord::new(
                "INPUT SEQUENCE".to_string(),
                cli.input.clone(),
                None, // No quality scores for direct sequence input
            );
            Box::new(std::iter::once(Ok(record)))
        } else {
            eprintln!("Error: Input is neither a valid file path nor a valid sequence");
            std::process::exit(1);
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
                let numbered_sequence = numbering::number_sequence(&record.sequence, &cli.scheme, &cli.chains);
                writeln!(output_writer, "{}", numbered_sequence).unwrap();
            }
            Err(e) => {
                eprintln!("Error processing sequence: {}", e);
                std::process::exit(1);
            }
        }
    }
}
