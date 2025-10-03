use immunum::schemes::get_scheme_with_cdr_definition;
use immunum::types::{CdrDefinitions, Chain, Scheme};
use std::collections::{HashMap, HashSet};

/// Extract consensus sequence from a NumberingScheme as a string
fn extract_consensus_sequence(scheme: &immunum::numbering_scheme_type::NumberingScheme) -> String {
    let mut positions: Vec<u32> = scheme.consensus_amino_acids.keys().cloned().collect();
    positions.sort();

    let mut consensus = String::new();
    for pos in positions {
        if let Some(amino_acids) = scheme.consensus_amino_acids.get(&pos) {
            // Take the first amino acid if multiple are present
            if !amino_acids.is_empty() {
                consensus.push(amino_acids[0] as char);
            }
        }
    }
    consensus
}

/// Analyze the diversity of amino acids at each position
fn analyze_consensus_diversity(
    scheme: &immunum::numbering_scheme_type::NumberingScheme,
    name: &str,
) {
    let mut positions: Vec<u32> = scheme.consensus_amino_acids.keys().cloned().collect();
    positions.sort();

    let mut multiple_aa_positions = 0;
    let mut total_alternatives = 0;

    println!("\n🔬 Consensus Diversity Analysis for {}", name);
    println!("Position -> Amino Acids");

    for pos in positions.iter().take(20) {
        // Show first 20 positions as example
        if let Some(amino_acids) = scheme.consensus_amino_acids.get(pos) {
            let aa_chars: Vec<char> = amino_acids.iter().map(|&aa| aa as char).collect();
            println!("  {}: {:?}", pos, aa_chars);

            if amino_acids.len() > 1 {
                multiple_aa_positions += 1;
                total_alternatives += amino_acids.len() - 1;
            }
        }
    }

    println!("  ... (showing first 20 positions)");

    // Count total positions with multiple options
    for pos in &positions {
        if let Some(amino_acids) = scheme.consensus_amino_acids.get(pos) {
            if amino_acids.len() > 1 {
                multiple_aa_positions += 1;
                total_alternatives += amino_acids.len() - 1;
            }
        }
    }

    println!("📊 Summary:");
    println!("  Total positions: {}", positions.len());
    println!(
        "  Positions with multiple amino acids: {}",
        multiple_aa_positions
    );
    println!("  Total alternative amino acids: {}", total_alternatives);
    println!(
        "  Diversity percentage: {:.1}%",
        (multiple_aa_positions as f64 / positions.len() as f64) * 100.0
    );
}

/// Generate all possible k-mers considering amino acid alternatives at each position
fn generate_comprehensive_kmers(
    scheme: &immunum::numbering_scheme_type::NumberingScheme,
    k: usize,
) -> HashSet<String> {
    let mut positions: Vec<u32> = scheme.consensus_amino_acids.keys().cloned().collect();
    positions.sort();

    // Build all possible sequences considering alternatives
    let mut all_possible_sequences = vec![String::new()];

    for pos in positions {
        if let Some(amino_acids) = scheme.consensus_amino_acids.get(&pos) {
            let mut new_sequences = Vec::new();

            for existing_seq in &all_possible_sequences {
                for &aa in amino_acids {
                    let mut new_seq = existing_seq.clone();
                    new_seq.push(aa as char);
                    new_sequences.push(new_seq);
                }
            }

            all_possible_sequences = new_sequences;

            // Prevent exponential explosion - if we have too many sequences, sample them
            if all_possible_sequences.len() > 1000 {
                all_possible_sequences.truncate(1000);
                println!("  ⚠️  Truncated to 1000 sequences to prevent memory explosion");
                break;
            }
        }
    }

    // Generate k-mers from all possible sequences
    let mut all_kmers = HashSet::new();
    for sequence in &all_possible_sequences {
        let kmers = generate_kmers(sequence, k);
        all_kmers.extend(kmers);
    }

    all_kmers
}

/// Generate k-mers from a sequence
fn generate_kmers(sequence: &str, k: usize) -> HashSet<String> {
    let mut kmers = HashSet::new();
    if sequence.len() >= k {
        for i in 0..=sequence.len() - k {
            kmers.insert(sequence[i..i + k].to_string());
        }
    }
    kmers
}

/// Calculate Jaccard similarity between two sets
fn jaccard_similarity(set1: &HashSet<String>, set2: &HashSet<String>) -> f64 {
    if set1.is_empty() && set2.is_empty() {
        return 1.0;
    }

    let intersection: HashSet<_> = set1.intersection(set2).collect();
    let union: HashSet<_> = set1.union(set2).collect();

    intersection.len() as f64 / union.len() as f64
}

/// Analyze k-mer divergence for different schemes
fn analyze_kmer_divergence() {
    println!("🧬 K-mer Divergence Analysis for Immunoglobulin Schemes");
    println!("======================================================");

    // Create schemes for different chain types
    let heavy_scheme =
        get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT);
    let kappa_scheme =
        get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGK, CdrDefinitions::IMGT);
    let lambda_scheme =
        get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGL, CdrDefinitions::IMGT);

    let schemes = vec![
        ("IGH", &heavy_scheme),
        ("IGK", &kappa_scheme),
        ("IGL", &lambda_scheme),
    ];

    // Extract consensus sequences
    println!("\n📝 Consensus Sequences (first amino acid only):");
    println!("-----------------------------------------------");
    let mut consensus_sequences = HashMap::new();
    for (name, scheme) in &schemes {
        let consensus = extract_consensus_sequence(scheme);
        println!("{}: {} amino acids", name, consensus.len());
        println!("   {}", consensus);
        consensus_sequences.insert(*name, consensus);
    }

    // Analyze diversity at each position
    for (name, scheme) in &schemes {
        analyze_consensus_diversity(scheme, name);
    }

    // Test different k-mer sizes
    for k in 3..=6 {
        println!("\n🔍 K-mer Analysis (k={}) - Simple Approach", k);
        println!("------------------------------------------");

        // Generate k-mer sets for each scheme (simple approach)
        let mut kmer_sets_simple = HashMap::new();
        for (name, sequence) in &consensus_sequences {
            let kmers = generate_kmers(sequence, k);
            println!("{}: {} unique {}-mers (simple)", name, kmers.len(), k);
            kmer_sets_simple.insert(*name, kmers);
        }

        println!("\n🔍 K-mer Analysis (k={}) - Comprehensive Approach", k);
        println!("--------------------------------------------------");

        // Generate k-mer sets considering all amino acid alternatives
        let mut kmer_sets_comprehensive = HashMap::new();
        for (name, scheme) in &schemes {
            let kmers = generate_comprehensive_kmers(scheme, k);
            println!(
                "{}: {} unique {}-mers (comprehensive)",
                name,
                kmers.len(),
                k
            );
            kmer_sets_comprehensive.insert(*name, kmers);
        }

        // Calculate pairwise similarities for both approaches
        let names = ["IGH", "IGK", "IGL"];

        println!("\n📊 Pairwise Jaccard Similarities (Simple):");
        for i in 0..names.len() {
            for j in i + 1..names.len() {
                let similarity =
                    jaccard_similarity(&kmer_sets_simple[names[i]], &kmer_sets_simple[names[j]]);
                println!("  {} vs {}: {:.3}", names[i], names[j], similarity);
            }
        }

        println!("\n📊 Pairwise Jaccard Similarities (Comprehensive):");
        for i in 0..names.len() {
            for j in i + 1..names.len() {
                let similarity = jaccard_similarity(
                    &kmer_sets_comprehensive[names[i]],
                    &kmer_sets_comprehensive[names[j]],
                );
                println!("  {} vs {}: {:.3}", names[i], names[j], similarity);
            }
        }

        // Calculate unique k-mers for simple approach
        println!("\n🎯 Unique K-mers - Simple Approach:");
        for name in &names {
            let current_set = &kmer_sets_simple[name];
            let mut unique_count = 0;

            for kmer in current_set {
                let mut is_unique = true;
                for other_name in &names {
                    if other_name != name && kmer_sets_simple[other_name].contains(kmer) {
                        is_unique = false;
                        break;
                    }
                }
                if is_unique {
                    unique_count += 1;
                }
            }

            let unique_percentage = (unique_count as f64 / current_set.len() as f64) * 100.0;
            println!(
                "  {}: {} unique k-mers ({:.1}%)",
                name, unique_count, unique_percentage
            );
        }

        // Calculate unique k-mers for comprehensive approach
        println!("\n🎯 Unique K-mers - Comprehensive Approach:");
        for name in &names {
            let current_set = &kmer_sets_comprehensive[name];
            let mut unique_count = 0;

            for kmer in current_set {
                let mut is_unique = true;
                for other_name in &names {
                    if other_name != name && kmer_sets_comprehensive[other_name].contains(kmer) {
                        is_unique = false;
                        break;
                    }
                }
                if is_unique {
                    unique_count += 1;
                }
            }

            let unique_percentage = (unique_count as f64 / current_set.len() as f64) * 100.0;
            println!(
                "  {}: {} unique k-mers ({:.1}%)",
                name, unique_count, unique_percentage
            );
        }
    }

    // Discrimination analysis
    println!("\n🎛️  Discrimination Analysis Comparison");
    println!("======================================");

    for k in 3..=6 {
        // Simple approach
        let mut kmer_sets_simple = HashMap::new();
        for (name, sequence) in &consensus_sequences {
            let kmers = generate_kmers(sequence, k);
            kmer_sets_simple.insert(*name, kmers);
        }

        // Comprehensive approach
        let mut kmer_sets_comprehensive = HashMap::new();
        for (name, scheme) in &schemes {
            let kmers = generate_comprehensive_kmers(scheme, k);
            kmer_sets_comprehensive.insert(*name, kmers);
        }

        // Calculate average inter-scheme similarity for both approaches
        let names = ["IGH", "IGK", "IGL"];

        // Simple approach similarities
        let mut similarities_simple = Vec::new();
        for i in 0..names.len() {
            for j in i + 1..names.len() {
                let similarity =
                    jaccard_similarity(&kmer_sets_simple[names[i]], &kmer_sets_simple[names[j]]);
                similarities_simple.push(similarity);
            }
        }

        // Comprehensive approach similarities
        let mut similarities_comprehensive = Vec::new();
        for i in 0..names.len() {
            for j in i + 1..names.len() {
                let similarity = jaccard_similarity(
                    &kmer_sets_comprehensive[names[i]],
                    &kmer_sets_comprehensive[names[j]],
                );
                similarities_comprehensive.push(similarity);
            }
        }

        let avg_similarity_simple =
            similarities_simple.iter().sum::<f64>() / similarities_simple.len() as f64;
        let discrimination_potential_simple = 1.0 - avg_similarity_simple;

        let avg_similarity_comprehensive = similarities_comprehensive.iter().sum::<f64>()
            / similarities_comprehensive.len() as f64;
        let discrimination_potential_comprehensive = 1.0 - avg_similarity_comprehensive;

        println!(
            "k={}: Simple: avg={:.3}, disc={:.3} | Comprehensive: avg={:.3}, disc={:.3}",
            k,
            avg_similarity_simple,
            discrimination_potential_simple,
            avg_similarity_comprehensive,
            discrimination_potential_comprehensive
        );
    }

    println!("\n✅ Analysis Complete!");
    println!("\nRecommendations:");
    println!("- Lower similarity values indicate better discrimination");
    println!("- Higher unique k-mer percentages suggest good prefiltering potential");
    println!("- Optimal k-mer size balances discrimination vs. noise resistance");
}

fn main() {
    analyze_kmer_divergence();
}
