//! Debug tool for validation alignment visualization

use immunum2::annotator::Annotator;
use immunum2::types::{Chain, Scheme};
use immunum2::validation::{load_validation_csv, validate_entry};
use immunum2::Position;
use std::env;
use std::path::PathBuf;
use std::process;

/// Get region boundaries for a given scheme (start, end) for each region
/// Returns: (FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4)
fn get_region_boundaries(scheme: Scheme) -> [(&'static str, u8, u8); 7] {
    match scheme {
        Scheme::IMGT => [
            ("FR1", 1, 26),
            ("CDR1", 27, 38),
            ("FR2", 39, 55),
            ("CDR2", 56, 65),
            ("FR3", 66, 104),
            ("CDR3", 105, 117),
            ("FR4", 118, 128),
        ],
        Scheme::Kabat => [
            ("FR1", 1, 25),
            ("CDR1", 26, 35),
            ("FR2", 36, 50),
            ("CDR2", 51, 57),
            ("FR3", 58, 92),
            ("CDR3", 93, 100),
            ("FR4", 101, 113),
        ],
    }
}

fn print_usage() {
    eprintln!("Usage: debug-validation <CHAIN> [SCHEME] [HEADER]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  CHAIN    Chain type (TRA, TRB, TRG, TRD, IGH, IGK, IGL)");
    eprintln!("  SCHEME   Optional: Numbering scheme (imgt, kabat). Default: imgt");
    eprintln!("  HEADER   Optional: Specific sequence header to debug");
    eprintln!("           If not provided, shows all sequences with imperfect alignment");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  debug-validation TRA");
    eprintln!("  debug-validation TRA imgt 1AC6_A");
    eprintln!("  debug-validation IGH kabat");
    eprintln!("  debug-validation IGH kabat \"5ijv_F|Heavy|F\"");
}

/// Get the region name for a position number given the scheme
fn get_region_for_position(pos_num: u8, scheme: Scheme) -> &'static str {
    let boundaries = get_region_boundaries(scheme);
    for (name, start, end) in boundaries {
        if pos_num >= start && pos_num <= end {
            return name;
        }
    }
    "?"
}

/// Print alignment comparison between expected and actual positions
pub fn print_alignment_comparison(
    sequence: &str,
    expected_positions: &[(Position, char)],
    actual_positions: &[Position],
    scheme: Scheme,
) {
    let seq_len = sequence.len();

    // Build position strings with proper spacing
    let mut seq_line = String::new();
    let mut exp_line = String::new();
    let mut act_line = String::new();
    let mut match_line = String::new();
    let mut region_line = String::new();

    let max_pos_width = actual_positions
        .iter()
        .chain(expected_positions.iter().map(|(p, _)| p))
        .map(|p| p.to_string().len())
        .max()
        .unwrap_or(3)
        .max(4); // At least 4 to fit region names like "CDR1"

    println!("\nAlignment Comparison:");
    println!("{}", "=".repeat(80));

    // Print in chunks of 20 residues for readability
    let chunk_size = 20;
    for chunk_start in (0..seq_len).step_by(chunk_size) {
        let chunk_end = (chunk_start + chunk_size).min(seq_len);

        seq_line.clear();
        exp_line.clear();
        act_line.clear();
        match_line.clear();
        region_line.clear();

        let mut prev_region = "";

        for i in chunk_start..chunk_end {
            let aa = sequence.chars().nth(i).unwrap_or('?');
            let actual_pos = actual_positions.get(i);

            // Find expected position for this sequence index
            let expected_pos = if i < expected_positions.len() {
                Some(&expected_positions[i].0)
            } else {
                None
            };

            // Determine region from expected position (or actual if no expected)
            let current_region = if let Some(pos) = expected_pos {
                get_region_for_position(pos.number, scheme)
            } else if let Some(pos) = actual_pos {
                get_region_for_position(pos.number, scheme)
            } else {
                "?"
            };

            // Show region name at boundaries, otherwise fill with dashes
            let region_str = if current_region != prev_region {
                prev_region = current_region;
                format!("<{:<width$}", current_region, width = max_pos_width - 1)
            } else {
                format!("{:-<width$}", "", width = max_pos_width)
            };

            // Format strings with proper alignment
            let aa_str = format!("{:^width$}", aa, width = max_pos_width);
            let exp_str = if let Some(pos) = expected_pos {
                format!("{:^width$}", pos.to_string(), width = max_pos_width)
            } else {
                format!("{:^width$}", "-", width = max_pos_width)
            };
            let act_str = if let Some(pos) = actual_pos {
                format!("{:^width$}", pos.to_string(), width = max_pos_width)
            } else {
                format!("{:^width$}", "-", width = max_pos_width)
            };

            // Match indicator
            let match_str = if let (Some(exp), Some(act)) = (expected_pos, actual_pos) {
                if exp == act {
                    format!("{:^width$}", "|", width = max_pos_width)
                } else {
                    format!("{:^width$}", "X", width = max_pos_width)
                }
            } else {
                format!("{:^width$}", " ", width = max_pos_width)
            };

            region_line.push_str(&region_str);
            region_line.push(' ');
            seq_line.push_str(&aa_str);
            seq_line.push(' ');
            exp_line.push_str(&exp_str);
            exp_line.push(' ');
            act_line.push_str(&act_str);
            act_line.push(' ');
            match_line.push_str(&match_str);
            match_line.push(' ');
        }

        println!("\nPositions {}-{}:", chunk_start, chunk_end - 1);
        println!("Rgn:  {}", region_line);
        println!("Seq:  {}", seq_line);
        println!("Exp:  {}", exp_line);
        println!("      {}", match_line);
        println!("Act:  {}", act_line);
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        print_usage();
        process::exit(1);
    }

    // Parse chain argument
    let chain: Chain = match args[1].parse() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error: Invalid chain '{}': {}", args[1], e);
            eprintln!();
            print_usage();
            process::exit(1);
        }
    };

    // Parse optional scheme argument (default to IMGT)
    let (scheme, header_arg_idx) = if args.len() >= 3 {
        match args[2].to_lowercase().as_str() {
            "imgt" => (Scheme::IMGT, 3),
            "kabat" => (Scheme::Kabat, 3),
            _ => (Scheme::IMGT, 2), // Not a scheme, treat as header
        }
    } else {
        (Scheme::IMGT, 2)
    };

    // Determine CSV path based on chain and scheme
    let csv_filename = match (chain, scheme) {
        (Chain::TRA, Scheme::IMGT) => "fixtures/validation/tcr_A_imgt.csv",
        (Chain::TRB, Scheme::IMGT) => "fixtures/validation/tcr_B_imgt.csv",
        (Chain::TRG, Scheme::IMGT) => "fixtures/validation/tcr_G_imgt.csv",
        (Chain::TRD, Scheme::IMGT) => "fixtures/validation/tcr_D_imgt.csv",
        (Chain::IGH, Scheme::IMGT) => "fixtures/validation/ab_H_imgt.csv",
        (Chain::IGK, Scheme::IMGT) => "fixtures/validation/ab_K_imgt.csv",
        (Chain::IGL, Scheme::IMGT) => "fixtures/validation/ab_L_imgt.csv",
        (Chain::IGH, Scheme::Kabat) => "fixtures/validation/ab_H_kabat.csv",
        (Chain::IGK, Scheme::Kabat) => "fixtures/validation/ab_K_kabat.csv",
        (Chain::IGL, Scheme::Kabat) => "fixtures/validation/ab_L_kabat.csv",
        (chain, Scheme::Kabat) => {
            eprintln!("Error: Kabat scheme not supported for chain {}", chain);
            process::exit(1);
        }
    };

    let csv_path = PathBuf::from(csv_filename);
    if !csv_path.exists() {
        eprintln!("Error: Validation file not found: {}", csv_filename);
        process::exit(1);
    }

    // Load validation data
    let entries = match load_validation_csv(&csv_path) {
        Ok(e) => e,
        Err(err) => {
            eprintln!("Error loading validation CSV: {}", err);
            process::exit(1);
        }
    };

    if entries.is_empty() {
        eprintln!("Error: No entries found in {}", csv_filename);
        process::exit(1);
    }
    println!("{}", "=".repeat(80));
    println!("DEBUG VALIDATION ALIGNMENT");

    // Create annotator with the selected scheme
    let annotator = match Annotator::new(&[chain], scheme) {
        Ok(a) => a,
        Err(err) => {
            eprintln!("Error creating annotator: {}", err);
            process::exit(1);
        }
    };

    // If a specific header is provided, debug that entry
    // Otherwise, show all imperfect entries
    if args.len() > header_arg_idx {
        let header = &args[header_arg_idx];
        let entry = match entries.iter().find(|e| e.header == *header) {
            Some(e) => e,
            None => {
                eprintln!("Error: Header '{}' not found in validation data", header);
                eprintln!();
                eprintln!("Available headers (first 10):");
                for (i, e) in entries.iter().take(10).enumerate() {
                    eprintln!("  {}: {}", i + 1, e.header);
                }
                if entries.len() > 10 {
                    eprintln!("  ... and {} more", entries.len() - 10);
                }
                process::exit(1);
            }
        };

        debug_entry(entry, &annotator, scheme, csv_filename);
    } else {
        // Show all imperfect entries
        let mut imperfect_count = 0;
        for entry in &entries {
            let result = match validate_entry(entry, &annotator, scheme) {
                Ok(r) => r,
                Err(err) => {
                    eprintln!("Error validating entry '{}': {}", entry.header, err);
                    continue;
                }
            };

            // Check if not perfect (has incorrect or missing positions)
            if result.incorrect_positions > 0 || result.missing_positions > 0 {
                imperfect_count += 1;
                debug_entry(entry, &annotator, scheme, csv_filename);
                println!("\n");
            }
        }

        if imperfect_count == 0 {
            println!("All {} entries have perfect alignment!", entries.len());
        } else {
            println!("{}", "=".repeat(80));
            println!(
                "Summary: {}/{} entries have imperfect alignment",
                imperfect_count,
                entries.len()
            );
        }
    }
}

fn debug_entry(
    entry: &immunum2::validation::ValidationEntry,
    annotator: &Annotator,
    scheme: Scheme,
    csv_filename: &str,
) {
    // Validate and get results
    let result = match validate_entry(entry, annotator, scheme) {
        Ok(r) => r,
        Err(err) => {
            eprintln!("Error validating entry: {}", err);
            return;
        }
    };

    let numbering = match annotator.annotate(&entry.sequence) {
        Ok(r) => r.numbering(scheme),
        Err(err) => {
            eprintln!("Error annotating sequence: {}", err);
            return;
        }
    };

    // Print results
    println!("{}", "=".repeat(80));
    println!("File: {}", csv_filename);
    println!("Header: {}", entry.header);
    println!("Sequence length: {} aa", entry.sequence.len());
    println!("Chain detected: {}", result.detected_chain);
    println!("Scheme: {}", scheme);
    println!("Alignment confidence: {:.2}", result.alignment_confidence);
    println!();
    println!("Validation Results:");
    println!("  Accuracy: {:.1}%", result.accuracy() * 100.0);
    println!(
        "  Correct positions: {}/{}",
        result.correct_positions, result.total_positions
    );
    println!("  Incorrect positions: {}", result.incorrect_positions);
    println!("  Missing positions: {}", result.missing_positions);

    print_alignment_comparison(
        &entry.sequence,
        &entry.expected_positions,
        &numbering,
        scheme,
    );
}
