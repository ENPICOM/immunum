use crate::{
    result::OutputFormat,
    types::{Chain, Scheme},
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
    /// Chain types to try numbering with (supports multiple values)
    #[arg(
        short, long, value_enum, ignore_case = true, num_args = 1.., help = "Chain types to number with (multiple values allowed)")]
    pub chains: Option<Vec<Chain>>,
    /// Output file path (optional, defaults to stdout)
    #[arg(
        short,
        long,
        help = "Output file path (if omitted, output goes to stdout)"
    )]
    pub output: Option<String>,
    /// Enable pre-filtering to reduce chain types tested based on sequence characteristics
    #[arg(
        long,
        help = "Enable pre-filtering to speed up numbering by reducing tested chain types"
    )]
    pub prefilter: bool,
    /// Return all detected chains for each input (instead of best single)
    #[arg(
        long = "all-chains",
        help = "Return all detected chains within each sequence (not only the best match)"
    )]
    pub all_chains: bool,
    /// Enable parallel processing of sequences (file inputs only)
    #[arg(
        long,
        help = "Enable parallel processing when input is a FASTA/FASTQ(.gz) file"
    )]
    pub parallel: bool,
    /// Gap penalty for conserved position (default: 55.0)
    #[arg(long, help = "Gap penalty for conserved positions")]
    pub gap_pen_cp: Option<f64>,
    /// Gap penalty for framework regions (default: 26.0)
    #[arg(long, help = "Gap penalty for framework regions")]
    pub gap_pen_fr: Option<f64>,
    /// Gap penalty for insertion points (default: 1.5)
    #[arg(long, help = "Gap penalty for insertion points")]
    pub gap_pen_ip: Option<f64>,
    /// Gap penalty for other positions (default: 1.0)
    #[arg(long, help = "Gap penalty for other positions")]
    pub gap_pen_op: Option<f64>,
    /// Gap penalty for CDR regions (default: 2.5)
    #[arg(long, help = "Gap penalty for CDR regions")]
    pub gap_pen_cdr: Option<f64>,
    /// Gap penalty for other non-classified positions (default: 11.0)
    #[arg(long, help = "Gap penalty for other non-classified positions")]
    pub gap_pen_other: Option<f64>,
    /// CDR increase penalty per position away from insertion point (default: 0.5)
    #[arg(
        long,
        help = "CDR increase penalty per position away from insertion point"
    )]
    pub cdr_increase: Option<f64>,
    /// Penalty for leaping insertion point in IMGT scheme (default: 1.0)
    #[arg(long, help = "Penalty for leaping insertion point in IMGT scheme")]
    pub pen_leap_insertion_point_imgt: Option<f64>,
    /// Penalty for leaping insertion point in KABAT scheme (default: 10.0)
    #[arg(long, help = "Penalty for leaping insertion point in KABAT scheme")]
    pub pen_leap_insertion_point_kabat: Option<f64>,
    /// Output format for results
    #[arg(
        short = 'f',
        long = "format",
        value_enum,
        ignore_case = true,
        default_value_t = OutputFormat::Simple,
        help = "Output format for results"
    )]
    pub format: OutputFormat,
}
