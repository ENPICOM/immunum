use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::env;
use std::fs;
use std::path::Path;

// Standard amino acid order for indexing
const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

// Gap and insertion penalty constants
// These control how the alignment algorithm handles gaps and insertions
const CDR_BASE_GAP_PENALTY: f32 = -6.0;
const FR_BASE_GAP_PENALTY: f32 = -20.0;

const HIGHLY_CONSERVED_MULTIPLIER: f32 = 3.0;
const CONSERVED_MULTIPLIER: f32 = 2.0;
const VARIABLE_MULTIPLIER: f32 = 1.0;

const NO_INSERTION_PENALTY: f32 = -50.0;
const CDR_INSERTION_PENALTY: f32 = -3.0;
const FR_INSERTION_PENALTY: f32 = -15.0;

// BLOSUM62 matrix (20x20 for standard amino acids)
const BLOSUM62: [[i8; 20]; 20] = [
    [
        4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2,
    ], // A
    [
        0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2,
    ], // C
    [
        -2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3,
    ], // D
    [
        -1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2,
    ], // E
    [
        -2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3,
    ], // F
    [
        0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3,
    ], // G
    [
        -2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2,
    ], // H
    [
        -1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1,
    ], // I
    [
        -1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2,
    ], // K
    [
        -1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1,
    ], // L
    [
        -1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1,
    ], // M
    [
        -2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2,
    ], // N
    [
        -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3,
    ], // P
    [
        -1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1,
    ], // Q
    [
        -1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2,
    ], // R
    [
        1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2,
    ], // S
    [
        0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2,
    ], // T
    [
        0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1,
    ], // V
    [
        -3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2,
    ], // W
    [
        -2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7,
    ], // Y
];

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    // Create output directory for scoring matrices
    let matrix_out = Path::new(&out_dir).join("matrices");
    fs::create_dir_all(&matrix_out).unwrap();

    // Map of chain names to their TSV files
    let chains = [
        ("IGH", "ab_imgt_H_consensus.tsv"),
        ("IGK", "ab_imgt_K_consensus.tsv"),
        ("IGL", "ab_imgt_L_consensus.tsv"),
        ("TRA", "tcr_imgt_A_consensus.tsv"),
        ("TRB", "tcr_imgt_B_consensus.tsv"),
        ("TRG", "tcr_imgt_G_consensus.tsv"),
        ("TRD", "tcr_imgt_D_consensus.tsv"),
    ];

    for (chain, tsv_file) in &chains {
        let tsv_path = Path::new(&manifest_dir)
            .join("resources")
            .join("consensus")
            .join(tsv_file);

        println!("cargo:rerun-if-changed={}", tsv_path.display());

        if !tsv_path.exists() {
            eprintln!("Warning: {} not found, skipping", tsv_path.display());
            continue;
        }

        // Read and process TSV
        let content = fs::read_to_string(&tsv_path).unwrap();
        let positions = process_consensus_tsv(&content);

        // Generate JSON scoring matrix
        let output_path = matrix_out.join(format!("{}.json", chain));
        write_scoring_matrix(&output_path, chain, &positions).unwrap();

        println!(
            "Generated scoring matrix for {}: {}",
            chain,
            output_path.display()
        );
    }
}

#[derive(Serialize, Deserialize)]
struct ScoringMatrix {
    chain: String,
    positions: Vec<PositionScores>,
}

#[derive(Serialize, Deserialize)]
struct PositionScores {
    position: u32,
    scores: HashMap<char, f32>,
    gap_penalty: f32,
    insertion_penalty: f32,
}

struct PositionData {
    position: u32,
    aa_frequencies: Vec<(char, f32)>,
    occupancy: f32,
    conservation_class: String,
    region: String,
    allows_insertion: bool,
}

fn process_consensus_tsv(content: &str) -> Vec<PositionData> {
    let mut positions = Vec::new();

    let mut lines = content.lines();
    // Skip header
    lines.next();

    for line in lines {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue;
        }

        // Parse position
        let position: u32 = match parts[0].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };

        // Parse consensus amino acids and frequencies
        let aas_field = parts[1];
        let freq_field = parts[2];
        let occupancy: f32 = parts[3].parse().unwrap_or(1.0);
        let conservation_class = parts[6].to_string();
        let region = parts[7].to_string();
        let allows_insertion: bool = parts[8] == "true";

        let aas: Vec<&str> = aas_field.split(',').collect();
        let freqs: Vec<f32> = freq_field
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();

        let aa_frequencies: Vec<(char, f32)> = aas
            .iter()
            .zip(freqs.iter())
            .filter_map(|(aa, freq)| {
                if aa.len() == 1 {
                    Some((aa.chars().next().unwrap(), *freq))
                } else {
                    None
                }
            })
            .collect();

        positions.push(PositionData {
            position,
            aa_frequencies,
            occupancy,
            conservation_class,
            region,
            allows_insertion,
        });
    }

    positions
}

fn write_scoring_matrix(
    path: &Path,
    chain: &str,
    positions: &[PositionData],
) -> std::io::Result<()> {
    let mut position_scores = Vec::new();

    for pos_data in positions {
        // Calculate scores for all 20 amino acids
        let scores_array = calculate_position_scores(&pos_data.aa_frequencies);

        // Convert to HashMap for JSON
        let mut scores = HashMap::new();
        for (i, &aa_byte) in AMINO_ACIDS.iter().enumerate() {
            scores.insert(aa_byte as char, scores_array[i]);
        }

        // Calculate gap penalties with region-specific and conservation-based penalties
        let (gap_penalty, insertion_penalty) = calculate_gap_penalties(
            pos_data.occupancy,
            &pos_data.region,
            &pos_data.conservation_class,
            pos_data.allows_insertion,
        );

        position_scores.push(PositionScores {
            position: pos_data.position,
            scores,
            gap_penalty,
            insertion_penalty,
        });
    }

    let matrix = ScoringMatrix {
        chain: chain.to_string(),
        positions: position_scores,
    };

    // Write JSON to file
    let json = serde_json::to_string_pretty(&matrix)?;
    fs::write(path, json)?;

    Ok(())
}

fn calculate_position_scores(aa_frequencies: &[(char, f32)]) -> [f32; 20] {
    let mut scores = [0.0f32; 20];

    // For each of the 20 standard amino acids
    for (i, &query_aa_byte) in AMINO_ACIDS.iter().enumerate() {
        let query_aa = query_aa_byte as char;

        // Calculate score based on frequency-weighted BLOSUM62
        let mut score = 0.0f32;
        let mut total_freq = 0.0f32;

        for (consensus_aa, freq) in aa_frequencies {
            if let (Some(qi), Some(ci)) = (aa_to_index(query_aa), aa_to_index(*consensus_aa)) {
                score += (*freq) * (BLOSUM62[qi][ci] as f32);
                total_freq += freq;
            }
        }

        // Normalize by total frequency
        if total_freq > 0.0 {
            scores[i] = score / total_freq;
        } else {
            // No consensus data, use neutral score
            scores[i] = -4.0;
        }
    }

    scores
}

/// Return the calculated gap penalty and insertion penalty.
/// Gap penalty reduced in low occupancy positions, but increased in highly conserved positions.
/// Insertion penalty should be very high in places that dont allow insertions
/// In a CDR the gap penalty is always much lower than a FR region
fn calculate_gap_penalties(
    occupancy: f32,
    region: &str,
    conservation_class: &str,
    allows_insertion: bool,
) -> (f32, f32) {
    let is_cdr = matches!(region, "CDR1" | "CDR2" | "CDR3");

    // Base gap penalties: much lower (less negative) for CDR, higher (more negative) for FR
    let base_gap_penalty = if is_cdr {
        CDR_BASE_GAP_PENALTY
    } else {
        FR_BASE_GAP_PENALTY
    };

    // Conservation multiplier: increase penalty significantly for conserved positions
    // This forces the alignment to respect conserved anchor points in FR regions
    let conservation_multiplier = match conservation_class {
        "highly_conserved" => HIGHLY_CONSERVED_MULTIPLIER,
        "conserved" => CONSERVED_MULTIPLIER,
        _ => VARIABLE_MULTIPLIER,
    };

    // Gap penalty: scale by occupancy and conservation
    // Low occupancy -> less negative (more permissive)
    // High conservation -> more negative (less permissive)
    let gap_penalty = base_gap_penalty * occupancy * conservation_multiplier;

    // Insertion penalty: very high when not allowed, moderate in CDRs, high in FRs
    let insertion_penalty = if !allows_insertion {
        NO_INSERTION_PENALTY
    } else if is_cdr {
        CDR_INSERTION_PENALTY
    } else {
        FR_INSERTION_PENALTY
    };

    (gap_penalty, insertion_penalty)
}

fn aa_to_index(aa: char) -> Option<usize> {
    AMINO_ACIDS.iter().position(|&a| a == aa as u8)
}
