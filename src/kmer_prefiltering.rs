use crate::numbering_scheme_type::NumberingScheme;
use std::collections::{HashMap, HashSet};

/// K-mer size used for prefiltering (based on analysis, k=4 provides good balance)
pub const KMER_SIZE: usize = 4;

/// Minimum k-mer overlap percentage required to pass prefiltering
pub const MIN_KMER_OVERLAP: f64 = 0.2;

/// Generate comprehensive k-mers from a NumberingScheme considering all amino acid alternatives
/// This approach generates k-mers directly from consensus positions without creating intermediate sequences
pub fn generate_scheme_kmers(consensus_amino_acids: HashMap<u32, Vec<u8>>) -> HashSet<String> {
    let mut positions: Vec<u32> = consensus_amino_acids.keys().cloned().collect();
    positions.sort();

    if positions.len() < KMER_SIZE {
        return HashSet::new();
    }

    let mut all_kmers = HashSet::new();

    // Generate k-mers by sliding window over consensus positions
    for window_start in 0..=positions.len() - KMER_SIZE {
        let window_positions = &positions[window_start..window_start + KMER_SIZE];

        // Generate all possible k-mers for this window
        generate_kmers_from_positions(&consensus_amino_acids, window_positions, &mut all_kmers);
    }

    all_kmers
}

/// Generate all possible k-mers from a specific window of consensus positions
fn generate_kmers_from_positions(
    consensus_amino_acids: &HashMap<u32, Vec<u8>>,
    positions: &[u32],
    kmers: &mut HashSet<String>,
) {
    // Get amino acid options for each position in the window
    let mut amino_acid_options: Vec<&Vec<u8>> = Vec::new();
    for &pos in positions {
        if let Some(amino_acids) = consensus_amino_acids.get(&pos) {
            amino_acid_options.push(amino_acids);
        } else {
            return; // Skip if any position missing
        }
    }

    // Generate all combinations recursively
    let mut current_kmer = Vec::new();
    generate_kmer_combinations(&amino_acid_options, 0, &mut current_kmer, kmers);
}

/// Recursively generate all amino acid combinations for a k-mer
fn generate_kmer_combinations(
    amino_acid_options: &[&Vec<u8>],
    position_index: usize,
    current_kmer: &mut Vec<u8>,
    kmers: &mut HashSet<String>,
) {
    if position_index == amino_acid_options.len() {
        // We've filled all positions, convert to string and add to set
        let kmer_string = String::from_utf8_lossy(current_kmer).to_string();
        kmers.insert(kmer_string);
        return;
    }

    // Try each amino acid option at the current position
    for &amino_acid in amino_acid_options[position_index] {
        current_kmer.push(amino_acid);
        generate_kmer_combinations(amino_acid_options, position_index + 1, current_kmer, kmers);
        current_kmer.pop();
    }
}

/// Generate k-mers from a query sequence (amino acid bytes)
pub fn generate_kmers_from_query(query_sequence: &[u8], k: usize) -> HashSet<String> {
    let mut kmers = HashSet::new();
    if query_sequence.len() >= k {
        for i in 0..=query_sequence.len() - k {
            let kmer = String::from_utf8_lossy(&query_sequence[i..i + k]).to_string();
            kmers.insert(kmer);
        }
    }
    kmers
}

/// Calculate k-mer overlap score between query and scheme k-mers
pub fn calculate_kmer_overlap(
    query_kmers: &HashSet<String>,
    scheme_kmers: &HashSet<String>,
) -> f64 {
    if query_kmers.is_empty() || scheme_kmers.is_empty() {
        return 0.0;
    }

    let intersection: HashSet<_> = query_kmers.intersection(scheme_kmers).collect();
    // println!("Intersection size: {}", intersection.len());

    // Use percentage of query k-mers that match (more forgiving than Jaccard)
    intersection.len() as f64 / query_kmers.len() as f64
}

/// Apply k-mer based prefiltering to select compatible schemes
pub fn apply_kmer_prefiltering<'a>(
    query_sequence: &'a [u8],
    schemes: &'a [NumberingScheme],
    min_kmer_overlap: f64,
) -> Vec<&'a NumberingScheme> {
    // Generate k-mers from query sequence
    let query_kmers = generate_kmers_from_query(query_sequence, KMER_SIZE);

    if query_kmers.is_empty() {
        // If query too short for k-mers, return all schemes
        return schemes.iter().collect();
    }

    let mut compatible_schemes = Vec::new();

    for scheme in schemes {
        // Calculate overlap score
        let overlap_score = calculate_kmer_overlap(&query_kmers, &scheme.kmer_set);

        // Keep scheme if it meets minimum overlap threshold
        if overlap_score >= min_kmer_overlap {
            compatible_schemes.push(scheme);
        }
    }

    // If no schemes pass k-mer filtering, return all schemes to avoid false negatives
    if compatible_schemes.is_empty() {
        schemes.iter().collect()
    } else {
        compatible_schemes
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schemes::get_scheme_with_cdr_definition;
    use crate::types::{CdrDefinitions, Chain, Scheme};

    #[test]
    fn test_kmer_generation() {
        let sequence = "QVQLV";
        let kmers = generate_kmers_from_query(sequence.as_bytes(), 3);

        let expected: HashSet<String> = ["QVQ", "VQL", "QLV"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(kmers, expected);
    }

    #[test]
    fn test_kmer_overlap_calculation() {
        let set1: HashSet<String> = ["ABC", "BCD", "CDE"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let set2: HashSet<String> = ["BCD", "CDE", "DEF"]
            .iter()
            .map(|s| s.to_string())
            .collect();

        let overlap = calculate_kmer_overlap(&set1, &set2);
        // Intersection: {BCD, CDE} = 2, Query set: {ABC, BCD, CDE} = 3
        // Query overlap: 2/3 ≈ 0.667
        assert!((overlap - 0.6667).abs() < 0.001);
    }

    #[test]
    fn test_scheme_kmer_generation() {
        let scheme = get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT);
        let kmers = generate_scheme_kmers(scheme.consensus_amino_acids);

        // Should generate k-mers from the consensus sequence
        assert!(!kmers.is_empty());
        println!("Generated {} k-mers for IGH scheme", kmers.len());
    }

    #[test]
    fn test_kmer_prefiltering_real_sequences() {
        // Representative immunoglobulin sequences for testing k-mer prefiltering
        // Note: These are curated test sequences, not directly from abpdseq_agreed.fasta
        let test_cases = vec![
            // Heavy chain sequences (IGH)
            (">test_IGH_1|IMGT_reference", "SEVQLVESGGGLVQPGGSLRLSCAASGFNLYYYSIHWVRQAPGKGLEWVASISPYSSSTSYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARGRWYRRALDYWGQGTLVTVSSASTKGPSVFPLAPSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSC".as_bytes(), Chain::IGH),
            (">test_IGH_2|human_heavy", "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS".as_bytes(), Chain::IGH),
            (">test_IGH_3|variable_region", "QVQLQQWGAGLLKPSETLSLTCAVYGGSFSGYYWSWIRQPPGKGLEWIGEINHSGVNTNYNEKFKSRVTISADTSKNQFSLKLSSVTAADTAVYYCAREGTTGWGWLGKPIGAFAHWGQGTLVTVSS".as_bytes(), Chain::IGH),

            // Kappa chain sequences (IGK)
            (">test_IGK_1|kappa_light", "DIQMTQSPSSLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYKASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYPFTFGQGTKVEIK".as_bytes(), Chain::IGK),
            (">test_IGK_2|kappa_variant", "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSDAKTYLAWYQQKPGQAPKLLIFWASTRESGVPDRFSGSGSGTDFTLTIRSLEPEDFAVYYCQQYYIYLTFGGGTKVEIK".as_bytes(), Chain::IGK),
            (">test_IGK_3|kappa_diverse", "DIQMTQSPSSLSASVGDRVTITCKASQDVSIGVAWYQQKPGKAPKLLIYSASYRYTGVPDRFSGSGSGTDFTLTISSLQPEDFATYYCQQHYSTATFGQGTKLEIK".as_bytes(), Chain::IGK),

            // Lambda chain sequences (IGL)
            (">test_IGL_1|lambda_light", "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTRVFGTGTKVTVL".as_bytes(), Chain::IGL),
            (">test_IGL_2|lambda_variant", "SYELTQPPSVSVSPGQTARITCSGDALPKQYAYWYQQKPGQAPVLVIYGKDDRPSGIPERFSGSNSGNTATLTISGAQADDESLYFCSSYTSSSTFVFGTGTKVTVL".as_bytes(), Chain::IGL),
            (">test_IGL_3|lambda_diverse", "QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAGVFGTGTKVTVL".as_bytes(), Chain::IGL),
        ];

        // Create all schemes
        let schemes = vec![
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT),
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGK, CdrDefinitions::IMGT),
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGL, CdrDefinitions::IMGT),
        ];

        for (name, sequence, expected_chain) in test_cases {
            println!(
                "\n=== Testing {} (expected: {:?}) ===",
                name, expected_chain
            );

            // Apply k-mer prefiltering
            let filtered_schemes = apply_kmer_prefiltering(sequence, &schemes, MIN_KMER_OVERLAP);

            // Generate query k-mers for detailed analysis
            let query_kmers = generate_kmers_from_query(sequence, KMER_SIZE);

            // Check overlap scores for each scheme
            let mut best_overlap = 0.0;
            let mut best_chain = None;

            for scheme in &schemes {
                let overlap_score = calculate_kmer_overlap(&query_kmers, &scheme.kmer_set);
                println!("  {:?}: overlap = {:.4}", scheme.chain_type, overlap_score);

                if overlap_score > best_overlap {
                    best_overlap = overlap_score;
                    best_chain = Some(scheme.chain_type);
                }
            }

            // Verify the expected chain has the highest overlap
            assert_eq!(
                best_chain,
                Some(expected_chain),
                "Expected {:?} to have highest overlap for sequence {}",
                expected_chain,
                name
            );

            // Verify the correct scheme was included in filtered results
            let expected_scheme_included = filtered_schemes
                .iter()
                .any(|s| s.chain_type == expected_chain);
            assert!(
                expected_scheme_included,
                "Expected scheme {:?} should be included in filtered results for {}",
                expected_chain, name
            );

            // Verify that at least the expected chain meets minimum threshold
            let expected_scheme = schemes
                .iter()
                .find(|s| s.chain_type == expected_chain)
                .unwrap();
            let expected_overlap = calculate_kmer_overlap(&query_kmers, &expected_scheme.kmer_set);
            assert!(expected_overlap >= MIN_KMER_OVERLAP,
                "Expected chain {:?} should meet minimum overlap threshold ({:.3}) for sequence {}, got {:.3}",
                expected_chain, MIN_KMER_OVERLAP, name, expected_overlap);

            println!(
                "  ✓ Best match: {:?} (overlap: {:.4})",
                best_chain.unwrap(),
                best_overlap
            );
            println!(
                "  ✓ Schemes passing filter: {}/{}",
                filtered_schemes.len(),
                schemes.len()
            );
        }
    }

    #[test]
    fn test_kmer_prefiltering_discrimination() {
        // Test that different chain types have significantly different overlap scores
        let heavy_seq = b"QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS";
        let kappa_seq = b"DIQMTQSPSSLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYKASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYPFTFGQGTKVEIK";
        let lambda_seq = b"QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTRVFGTGTKVTVL";

        let heavy_scheme =
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT);
        let kappa_scheme =
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGK, CdrDefinitions::IMGT);
        let lambda_scheme =
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGL, CdrDefinitions::IMGT);

        // Test heavy sequence against all schemes
        let heavy_query_kmers = generate_kmers_from_query(heavy_seq, KMER_SIZE);
        let heavy_vs_heavy = calculate_kmer_overlap(&heavy_query_kmers, &heavy_scheme.kmer_set);
        let heavy_vs_kappa = calculate_kmer_overlap(&heavy_query_kmers, &kappa_scheme.kmer_set);
        let heavy_vs_lambda = calculate_kmer_overlap(&heavy_query_kmers, &lambda_scheme.kmer_set);

        println!("Heavy sequence overlaps:");
        println!("  vs Heavy: {:.4}", heavy_vs_heavy);
        println!("  vs Kappa: {:.4}", heavy_vs_kappa);
        println!("  vs Lambda: {:.4}", heavy_vs_lambda);

        assert!(
            heavy_vs_heavy > heavy_vs_kappa,
            "Heavy sequence should have higher overlap with heavy scheme than kappa"
        );
        assert!(
            heavy_vs_heavy > heavy_vs_lambda,
            "Heavy sequence should have higher overlap with heavy scheme than lambda"
        );

        // Test kappa sequence against all schemes
        let kappa_query_kmers = generate_kmers_from_query(kappa_seq, KMER_SIZE);
        let kappa_vs_heavy = calculate_kmer_overlap(&kappa_query_kmers, &heavy_scheme.kmer_set);
        let kappa_vs_kappa = calculate_kmer_overlap(&kappa_query_kmers, &kappa_scheme.kmer_set);
        let kappa_vs_lambda = calculate_kmer_overlap(&kappa_query_kmers, &lambda_scheme.kmer_set);

        println!("Kappa sequence overlaps:");
        println!("  vs Heavy: {:.4}", kappa_vs_heavy);
        println!("  vs Kappa: {:.4}", kappa_vs_kappa);
        println!("  vs Lambda: {:.4}", kappa_vs_lambda);

        assert!(
            kappa_vs_kappa > kappa_vs_heavy,
            "Kappa sequence should have higher overlap with kappa scheme than heavy"
        );
        assert!(
            kappa_vs_kappa > kappa_vs_lambda,
            "Kappa sequence should have higher overlap with kappa scheme than lambda"
        );

        // Test lambda sequence against all schemes
        let lambda_query_kmers = generate_kmers_from_query(lambda_seq, KMER_SIZE);
        let lambda_vs_heavy = calculate_kmer_overlap(&lambda_query_kmers, &heavy_scheme.kmer_set);
        let lambda_vs_kappa = calculate_kmer_overlap(&lambda_query_kmers, &kappa_scheme.kmer_set);
        let lambda_vs_lambda = calculate_kmer_overlap(&lambda_query_kmers, &lambda_scheme.kmer_set);

        println!("Lambda sequence overlaps:");
        println!("  vs Heavy: {:.4}", lambda_vs_heavy);
        println!("  vs Kappa: {:.4}", lambda_vs_kappa);
        println!("  vs Lambda: {:.4}", lambda_vs_lambda);

        assert!(
            lambda_vs_lambda > lambda_vs_heavy,
            "Lambda sequence should have higher overlap with lambda scheme than heavy"
        );
        assert!(
            lambda_vs_lambda > lambda_vs_kappa,
            "Lambda sequence should have higher overlap with lambda scheme than kappa"
        );
    }

    #[test]
    fn test_min_overlap_threshold_tuning() {
        // Test various overlap thresholds to ensure our current setting is appropriate
        let test_sequences = vec![
            (b"QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS".as_slice(), Chain::IGH),
            (b"DIQMTQSPSSLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYKASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYPFTFGQGTKVEIK".as_slice(), Chain::IGK),
            (b"QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTRVFGTGTKVTVL".as_slice(), Chain::IGL),
        ];

        let schemes = vec![
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT),
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGK, CdrDefinitions::IMGT),
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGL, CdrDefinitions::IMGT),
        ];

        let thresholds = vec![0.05, 0.08, 0.10, 0.12, 0.15, 0.20];

        println!("\nOverlap threshold analysis:");
        println!("Sequence -> Expected | Threshold: Schemes passing filter");

        for (sequence, expected_chain) in &test_sequences {
            let query_kmers = generate_kmers_from_query(sequence, KMER_SIZE);

            // Find the expected scheme's overlap score
            let expected_scheme = schemes
                .iter()
                .find(|s| s.chain_type == *expected_chain)
                .unwrap();
            let expected_overlap = calculate_kmer_overlap(&query_kmers, &expected_scheme.kmer_set);

            print!(
                "{:?} chain (overlap: {:.3}) | ",
                expected_chain, expected_overlap
            );

            for &threshold in &thresholds {
                let passing_schemes = schemes
                    .iter()
                    .filter(|s| calculate_kmer_overlap(&query_kmers, &s.kmer_set) >= threshold)
                    .count();

                print!("{:.2}: {}/3  ", threshold, passing_schemes);
            }
            println!();

            // Assert that our current threshold works for the expected scheme
            assert!(
                expected_overlap >= MIN_KMER_OVERLAP,
                "Current threshold {:.3} should allow {:?} chain (overlap: {:.3})",
                MIN_KMER_OVERLAP,
                expected_chain,
                expected_overlap
            );
        }
    }

    #[test]
    fn test_kmer_prefiltering_paired_sequences() {
        // Test k-mer prefiltering with paired chain constructs (heavy+light)
        println!("\n🧬 Testing K-mer Prefiltering with Paired Chain Sequences");
        println!("========================================================");

        // Paired sequences from annotator tests - these should contain multiple chains
        let paired_sequences = vec![
            (
                ">paired_IGH_IGK|heavy+kappa_construct",
                // Heavy chain followed by kappa light chain
                "QVQLVQSGAVIKTPGSSVKISCRASGYNFRDYSIHWVRLIPDKGFEWIGWIKPLWGAVSYARQLQGRVSMTRQLSQDPDDPDWGVAYMEFSGLTPADTAEYFCVRRGSCDYCGDFPWQYWCQGTVVVVSSASTKGPSVFPLAPSSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKDIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSATDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC".as_bytes(),
                vec![Chain::IGH, Chain::IGK]
            ),
            (
                ">paired_IGH_IGL|heavy+lambda_construct",
                // Heavy chain followed by lambda light chain
                "QVQLVQSGAVIKTPGSSVKISCRASGYNFRDYSIHWVRLIPDKGFEWIGWIKPLWGAVSYARQLQGRVSMTRQLSQDPDDPDWGVAYMEFSGLTPADTAEYFCVRRGSCDYCGDFPWQYWCQGTVVVVSSASTKGPSVFPLAPSSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSALTQPPSASGSLGQSVTISCTGTSSDVGGYNYVSWYQQHAGKAPKVIIYEVNKRPSGVPDRFSGSKSGNTASLTVSGLQAEDEADYYCSSYEGSDNFVFGTGTKVTVLGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS".as_bytes(),
                vec![Chain::IGH, Chain::IGL]
            ),
            (
                ">paired_IGK_IGH|kappa+heavy_construct",
                // Kappa light chain followed by heavy chain (reverse order)
                "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSATDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNECQVQLVQSGAVIKTPGSSVKISCRASGYNFRDYSIHWVRLIPDKGFEWIGWIKPLWGAVSYARQLQGRVSMTRQLSQDPDDPDWGVAYMEFSGLTPADTAEYFCVRRGSCDYCGDFPWQYWCQGTVVVVSSASTKGPSVFPLAPSSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK".as_bytes(),
                vec![Chain::IGK, Chain::IGH]
            ),
        ];

        let schemes = vec![
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT),
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGK, CdrDefinitions::IMGT),
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGL, CdrDefinitions::IMGT),
        ];

        for (name, sequence, expected_chains) in paired_sequences {
            println!(
                "\n=== Testing {} (expected: {:?}) ===",
                name, expected_chains
            );

            // Apply k-mer prefiltering
            let filtered_schemes = apply_kmer_prefiltering(sequence, &schemes, MIN_KMER_OVERLAP);

            // Generate query k-mers for analysis
            let query_kmers = generate_kmers_from_query(sequence, KMER_SIZE);
            println!("  Query sequence length: {} amino acids", sequence.len());
            println!("  Query k-mers generated: {}", query_kmers.len());

            // Analyze overlap scores for each scheme
            let mut overlap_scores = Vec::new();
            for scheme in &schemes {
                let overlap_score = calculate_kmer_overlap(&query_kmers, &scheme.kmer_set);
                overlap_scores.push((scheme.chain_type, overlap_score));
                println!("  {:?}: overlap = {:.4}", scheme.chain_type, overlap_score);
            }

            // Sort by overlap score (highest first)
            overlap_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

            // Check if the expected chains have the highest overlaps
            println!("  📊 Overlap ranking:");
            for (i, (chain, score)) in overlap_scores.iter().enumerate() {
                let marker = if expected_chains.contains(chain) {
                    "✓"
                } else {
                    " "
                };
                println!("    {}. {} {:?}: {:.4}", i + 1, marker, chain, score);
            }

            // For paired sequences, we expect both chains to have higher overlap than unrelated chains
            // but the exact ranking can vary due to k-mer mixing between chains
            let expected_min_overlap = MIN_KMER_OVERLAP;

            let mut all_expected_chains_above_threshold = true;
            for expected_chain in &expected_chains {
                let (_chain, score) = overlap_scores
                    .iter()
                    .find(|(chain, _)| chain == expected_chain)
                    .unwrap();

                if *score < expected_min_overlap {
                    all_expected_chains_above_threshold = false;
                    println!(
                        "    ⚠️ Expected chain {:?} has overlap {:.4} below threshold {:.3}",
                        expected_chain, score, expected_min_overlap
                    );
                }
            }

            // The important test: all expected chains should be above threshold
            assert!(
                all_expected_chains_above_threshold,
                "All expected chains should meet minimum overlap threshold for paired sequence {}",
                name
            );

            // Check filtering results
            println!("  🎯 Filtering results:");
            println!("    Schemes before filtering: {}", schemes.len());
            println!("    Schemes after filtering: {}", filtered_schemes.len());

            // Verify that expected chains are included in filtered results
            let filtered_chain_types: Vec<Chain> =
                filtered_schemes.iter().map(|s| s.chain_type).collect();

            for expected_chain in &expected_chains {
                assert!(
                    filtered_chain_types.contains(expected_chain),
                    "Expected chain {:?} should pass k-mer filtering for sequence {}",
                    expected_chain,
                    name
                );
            }

            // Verify that expected chains meet minimum threshold
            for expected_chain in &expected_chains {
                let expected_scheme = schemes
                    .iter()
                    .find(|s| s.chain_type == *expected_chain)
                    .unwrap();
                let expected_overlap =
                    calculate_kmer_overlap(&query_kmers, &expected_scheme.kmer_set);

                assert!(expected_overlap >= MIN_KMER_OVERLAP,
                    "Expected chain {:?} should meet minimum overlap threshold ({:.3}) for sequence {}, got {:.3}",
                    expected_chain, MIN_KMER_OVERLAP, name, expected_overlap);
            }

            println!("  ✅ All expected chains correctly identified and filtered");
        }

        println!("\n🎉 Paired sequence k-mer prefiltering test completed successfully!");
        println!("📋 Key findings:");
        println!(
            "   ✅ All expected chains from paired constructs passed minimum overlap threshold"
        );
        println!("   ✅ K-mer prefiltering successfully retained relevant schemes for multi-chain analysis");
        println!(
            "   ⚠️  Note: K-mer ranking can be affected by chain mixing in concatenated sequences"
        );
        println!(
            "   💡 K-mer prefiltering is best used as initial filter before segment-based analysis"
        );
        println!("\n📊 Limitation identified:");
        println!(
            "   • K-mer prefiltering works at whole-sequence level, not individual chain segments"
        );
        println!("   • For paired sequences, subsequent alignment-based analysis is needed for precise identification");
        println!("   • However, k-mer prefiltering still effectively reduces computational load by filtering out irrelevant schemes");
    }
}
