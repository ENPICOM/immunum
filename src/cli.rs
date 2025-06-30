use clap::Parser;
use crate::types::{Scheme, Chain};

#[derive(Parser)]
#[command(name = "immunumber")]  // !Tool name is up for discussion
#[command(about = "A CLI tool for numbering Antibody and T-cell receptor sequences")]
#[command(version = "0.1.0")]
pub struct Cli {
    /// Input: either a file path (FASTA/FASTQ) or a single sequence string
    #[arg(
        short,
        long,
        help = "File path to sequences (FASTA/FASTQ) or a single sequence string"
    )]
    pub input: String,
    /// Numbering scheme to use
    #[arg(
        short, long, value_enum, ignore_case = true, default_value_t = Scheme::IMGT, help = "Numbering scheme")]
    pub scheme: Scheme,
    /// Chain types to process (supports multiple values)
    #[arg(
        short, long, value_enum, ignore_case = true, num_args = 1.., help = "Chain types to process (multiple values allowed)")]
    pub chains: Vec<Chain>,
    /// Output file path (optional, defaults to stdout)
    #[arg(
        short,
        long,
        help = "Output file path (if omitted, output goes to stdout)"
    )]
    pub output: Option<String>,
}
