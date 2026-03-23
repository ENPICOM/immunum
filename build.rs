use serde::{Deserialize, Serialize};
use std::env;
use std::fs;
use std::path::Path;

// Standard amino acid order for indexing
const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

// Gap and insertion penalty constants
// These control how the alignment algorithm handles gaps and insertions
const CDR_BASE_GAP_PENALTY: f32 = -6.0;
const FR_BASE_GAP_PENALTY: f32 = -12.0;

const HIGHLY_CONSERVED_MULTIPLIER: f32 = 3.0;
const CONSERVED_MULTIPLIER: f32 = 2.0;
const VARIABLE_MULTIPLIER: f32 = 1.0;

const NO_INSERTION_PENALTY: f32 = -50.0;
const CDR_INSERTION_PENALTY: f32 = -3.0;

// IMGT CDR center positions where insertions are absorbed
const CDR_CENTER_POSITIONS: &[u32] = &[32, 33, 60, 61, 111, 112];

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
    // Set rpath to Python's lib dir so cargo test can find libpython
    // see https://github.com/astral-sh/uv/issues/11006
    #[cfg(feature = "python")]
    if cfg!(unix) {
        if let Some(lib_dir) = &pyo3_build_config::get().lib_dir {
            println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib_dir);
        }
    }

    let out_dir = env::var("OUT_DIR").unwrap();
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    // Create output directory for scoring matrices
    let matrix_out = Path::new(&out_dir).join("matrices");
    fs::create_dir_all(&matrix_out).unwrap();

    // Map of chain names to their CSV files
    let chains = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"];

    for chain in &chains {
        let csv_path = Path::new(&manifest_dir)
            .join("resources")
            .join("consensus")
            .join(format!("{}.csv", chain));

        println!("cargo:rerun-if-changed={}", csv_path.display());

        if !csv_path.exists() {
            eprintln!("Warning: {} not found, skipping", csv_path.display());
            continue;
        }

        // Read and process CSV
        let content = fs::read_to_string(&csv_path).unwrap();
        let positions = process_consensus_csv(&content);

        // Generate JSON scoring matrix
        let output_path = matrix_out.join(format!("{}.json", chain));
        write_scoring_matrix(&output_path, &positions).unwrap();

        println!(
            "Generated scoring matrix for {}: {}",
            chain,
            output_path.display()
        );
    }
}

#[derive(Serialize, Deserialize)]
struct ScoringMatrix {
    positions: Vec<PositionScores>,
}

#[derive(Serialize, Deserialize)]
struct PositionScores {
    position: u32,
    scores: [f32; 26],
    gap_penalty: f32,
    insertion_penalty: f32,
    max_score: f32,
    counts_for_confidence: bool,
}

struct PositionData {
    position: u32,
    aa_frequencies: Vec<(char, f32)>,
    occupancy: f32,
    region: String,
}

fn process_consensus_csv(content: &str) -> Vec<PositionData> {
    let mut positions = Vec::new();

    let mut lines = content.lines();
    // Skip header
    lines.next();

    for line in lines {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 5 {
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
        let region = parts[4].to_string();

        let aas: Vec<&str> = aas_field.split('|').collect();
        let freqs: Vec<f32> = freq_field
            .split('|')
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
            region,
        });
    }

    positions
}

fn write_scoring_matrix(path: &Path, positions: &[PositionData]) -> std::io::Result<()> {
    let mut position_scores = Vec::new();

    for pos_data in positions {
        // Calculate scores for all 20 amino acids into a [f32; 26] array
        let scores_20 = calculate_position_scores(&pos_data.aa_frequencies);

        // Map 20 amino acid scores into 26-slot array indexed by (byte - b'A')
        let mut scores = [-4.0f32; 26];
        for (i, &aa_byte) in AMINO_ACIDS.iter().enumerate() {
            scores[(aa_byte - b'A') as usize] = scores_20[i];
        }

        // Calculate gap penalties with region-specific and conservation-based penalties
        let max_freq = pos_data
            .aa_frequencies
            .first()
            .map(|(_, f)| *f)
            .unwrap_or(0.0);
        let (gap_penalty, insertion_penalty) = calculate_gap_penalties(
            pos_data.position,
            pos_data.occupancy,
            max_freq,
            &pos_data.region,
        );

        let max_score = scores.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        let counts_for_confidence = pos_data.occupancy > 0.9;

        position_scores.push(PositionScores {
            position: pos_data.position,
            scores,
            gap_penalty,
            insertion_penalty,
            max_score,
            counts_for_confidence,
        });
    }

    let matrix = ScoringMatrix {
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
/// Insertion penalty is permissive only at CDR center positions where length variation is absorbed.
/// In a CDR the gap penalty is always much lower than a FR region.
fn calculate_gap_penalties(
    position: u32,
    occupancy: f32,
    max_freq: f32,
    region: &str,
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
    let conservation_multiplier = if occupancy >= 0.5 && max_freq >= 0.9 {
        HIGHLY_CONSERVED_MULTIPLIER
    } else if occupancy >= 0.5 && max_freq >= 0.7 {
        CONSERVED_MULTIPLIER
    } else {
        VARIABLE_MULTIPLIER
    };

    // Gap penalty: scale by occupancy and conservation
    // Low occupancy -> less negative (more permissive)
    // High conservation -> more negative (less permissive)
    let gap_penalty = base_gap_penalty * occupancy * conservation_multiplier;

    // Insertion penalty: permissive at CDR center positions, harsh everywhere else
    let is_cdr_center = is_cdr && CDR_CENTER_POSITIONS.contains(&position);
    let insertion_penalty = if is_cdr_center {
        CDR_INSERTION_PENALTY
    } else {
        NO_INSERTION_PENALTY
    };

    (gap_penalty, insertion_penalty)
}

fn aa_to_index(aa: char) -> Option<usize> {
    AMINO_ACIDS.iter().position(|&a| a == aa as u8)
}
