use crate::{
    constants::MINIMAL_CHAIN_IDENTITY,
    kmer_prefiltering::MIN_KMER_OVERLAP,
    types::{CdrDefinitions, Chain, Scheme},
};
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
    /// CDR definition scheme to use (defaults to same as numbering scheme)
    #[arg(
        long = "cdr-definitions",
        value_enum,
        ignore_case = true,
        help = "CDR definition scheme (defaults to same as numbering scheme)"
    )]
    pub cdr_definitions: Option<CdrDefinitions>,
    /// Chain types to try numbering with (supports multiple values)
    #[arg(
        short, long, value_enum, ignore_case = true, num_args = 1.., help = "Chain types to number with (multiple values allowed)")]
    pub chains: Option<Vec<Chain>>,
    /// Output file path (required)
    #[arg(short, long, help = "Output file path for results")]
    pub output: String,
    /// Disable pre-filtering (prefiltering is enabled by default to speed up numbering)
    #[arg(
        long = "disable-prefiltering",
        help = "Disable pre-filtering (prefiltering is enabled by default to speed up numbering)"
    )]
    pub disable_prefiltering: bool,
    /// Maximum number of non-overlapping chains to find in each sequence (0 = unlimited, 1 = single chain mode)
    #[arg(
        long,
        default_value_t = 1,
        help = "Maximum number of non-overlapping chains to find in each sequence (0 = unlimited, 1 = single chain mode)"
    )]
    pub max_chains: usize,
    /// Number of threads for parallel processing (defaults to number of CPU cores)
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = num_cpus::get(),
        help = "Number of threads for parallel processing (defaults to number of CPU cores)"
    )]
    pub threads: usize,
    /// Minimal confidence score of alignment to a scheme to pass
    #[arg(
        short = 'm',
        long = "min-confidence",
        default_value_t = MINIMAL_CHAIN_IDENTITY,
        help = "Minimal confidence score of alignment to a scheme to pass"

    )]
    pub min_confidence: f64,
    /// Enable verbose output with per-sequence progress logging
    #[arg(
        short = 'v',
        long = "verbose",
        help = "Show detailed progress for each sequence processed"
    )]
    pub verbose: bool,
    /// Minimum k-mer overlap threshold for prefiltering (fraction of query k-mers that must match)
    #[arg(
        long = "min-kmer-overlap",
        default_value_t = MIN_KMER_OVERLAP,
        help = "Minimum k-mer overlap threshold for prefiltering (fraction of query k-mers that must match)"
    )]
    pub min_kmer_overlap: f64,
}
