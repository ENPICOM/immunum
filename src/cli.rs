use crate::types::{Chain, Scheme};
use clap::Parser;

#[derive(Parser)]
#[command(name = "immunumber")] // !Tool name is up for discussion
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
    /// Chain types to try numbering with (supports multiple values)
    #[arg(
        short, long, value_enum, ignore_case = true, num_args = 1.., help = "Chain types to number with (multiple values allowed)")]
    pub chains: Option<Vec<Chain>>,
    /// Output file path (required)
    #[arg(short, long, help = "Output file path for results")]
    pub output: String,
    /// Enable pre-filtering to reduce chain types tested based on sequence characteristics
    #[arg(
        long,
        help = "Enable pre-filtering to speed up numbering by reducing tested chain types"
    )]
    pub prefilter: bool,
    /// Maximum number of non-overlapping chains to find in each sequence (0 = unlimited, 1 = single chain mode)
    #[arg(
        long,
        default_value_t = 1,
        help = "Maximum number of non-overlapping chains to find in each sequence (0 = unlimited, 1 = single chain mode)"
    )]
    pub max_chains: usize,
    /// Output format for results
    // #[arg(
    //     short = 'f',
    //     long = "format",
    //     value_enum,
    //     ignore_case = true,
    //     default_value_t = OutputFormat::Simple,
    //     help = "Output format for results"
    // )]
    // pub format: OutputFormat,
    /// Number of threads for parallel processing (defaults to number of CPU cores)
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads for parallel processing (defaults to number of CPU cores)"
    )]
    pub threads: Option<usize>,
    /// Enable verbose output with per-sequence progress logging
    #[arg(
        short = 'v',
        long = "verbose",
        help = "Show detailed progress for each sequence processed"
    )]
    pub verbose: bool,
}
