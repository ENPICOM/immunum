use crate::consensus_scoring::write_all_scoring_matrices;
use crate::constants::{MINIMAL_CHAIN_IDENTITY, MINIMAL_CHAIN_LENGTH};
use crate::fastx::{from_path, FastxRecord};
use crate::numbering_scheme_type::{NumberingOutput, NumberingScheme};
use crate::prefiltering::{get_terminal_schemes, run_pre_scan, select_chains_from_pre_scan};
use crate::schemes::{
    get_imgt_heavy_scheme, get_imgt_kappa_scheme, get_imgt_lambda_scheme, get_kabat_heavy_scheme,
    get_kabat_kappa_scheme, get_kabat_lambda_scheme,
};
use crate::types::{Chain, Scheme};
use std::fs;

/// Runs alignment of given schemes on sequence and selects one with highest identity
pub(crate) fn find_highest_identity_chain<'a>(
    query_sequence: &'a [u8],
    numbering_schemes: &Vec<&'a NumberingScheme>,
) -> Result<NumberingOutput<'a>, &'static str> {
    let mut highest_identity: f64 = -0.1;
    let mut best_output: Result<NumberingOutput, &'static str> = Err("No numbering schemes passed");

    for scheme in numbering_schemes {
        let output: NumberingOutput = scheme.number_sequence(query_sequence);
        if output.identity > highest_identity {
            highest_identity = output.identity;
            best_output = Ok(output);
        }
    }
    best_output
}
/// Runs alignement on all sequences in fastx file, and writes output to .txt file
pub fn number_sequences_and_write_output(
    fasta_file: &str,
    scheme: Scheme,
    chains: &[Chain],
    output_file: &str,
    update_scoring_matrices: bool,
) {
    // Update scoring matrices
    if update_scoring_matrices {
        write_all_scoring_matrices()
    }

    // Read in fasta
    let reader = from_path(fasta_file).unwrap();
    let records: Vec<FastxRecord> = reader.collect::<Result<Vec<FastxRecord>, _>>().unwrap();

    // get schemes dependent on selected numbering method
    let schemes = match scheme {
        Scheme::IMGT => vec![
            get_imgt_heavy_scheme(),
            get_imgt_kappa_scheme(),
            get_imgt_lambda_scheme(),
        ],
        Scheme::KABAT => vec![
            get_kabat_heavy_scheme(),
            get_kabat_kappa_scheme(),
            get_kabat_lambda_scheme(),
        ],
    };

    // only include schemes for selected chains
    let schemes: Vec<NumberingScheme> = schemes
        .into_iter()
        .filter(|scheme| chains.contains(&scheme.chain_type))
        .collect();

    // get pre-filter schemes
    let terminal_schemes = get_terminal_schemes(&schemes);

    let mut output_str = "".to_string();
    output_str.push_str(
        "Name\tSequence\tNumbering\tScore\tChain\tcdr1\tcdr2\tcdr3\tfmwk1\tfmwk2\tfmwk3\tfmwk4\n",
    );
    // run annotation for all sequences
    for r in records {
        let converted_sequence = r.sequence.into_bytes();

        // Select models to run using pre-scan
        let (pre_scan_output, highest_score) = run_pre_scan(&converted_sequence, &terminal_schemes);
        let pre_filter_chains = select_chains_from_pre_scan(&pre_scan_output, highest_score);
        let filtered_schemes: Vec<&NumberingScheme> = schemes
            .iter()
            .filter(|scheme| pre_filter_chains.contains(&scheme.chain_type))
            .collect();
        let output_results: Vec<Result<NumberingOutput, &str>> =
            find_all_chains(&converted_sequence, filtered_schemes);

        for output_result in output_results {
            match output_result {
                Ok(output) => {
                    //create string output
                    output_str.push_str(&r._name);
                    output_str.push('\t');
                    output_str.push_str(
                        std::str::from_utf8(output.sequence)
                            .expect("Non-UTF8 character in sequence"),
                    );
                    output_str.push('\t');
                    output_str.push_str(&output.numbering.join(","));
                    output_str.push('\t');
                    output_str.push_str(&format!("{}", output.identity));
                    output_str.push('\t');
                    output_str.push_str(match output.scheme.chain_type {
                        Chain::IGH => "H",
                        Chain::IGK => "K",
                        Chain::IGL => "L",
                        Chain::TRA => "A",
                        Chain::TRB => "B",
                        Chain::TRD => "D",
                        Chain::TRG => "G",
                    });
                    output_str.push('\t');

                    // TODO add regions
                    output_str.push_str(output.get_query_region(&output.scheme.cdr1));
                    output_str.push('\t');
                    output_str.push_str(output.get_query_region(&output.scheme.cdr2));
                    output_str.push('\t');
                    output_str.push_str(output.get_query_region(&output.scheme.cdr3));
                    output_str.push('\t');
                    output_str.push_str(output.get_query_region(&output.scheme.fr1));
                    output_str.push('\t');
                    output_str.push_str(output.get_query_region(&output.scheme.fr2));
                    output_str.push('\t');
                    output_str.push_str(output.get_query_region(&output.scheme.fr3));
                    output_str.push('\t');
                    output_str.push_str(output.get_query_region(&output.scheme.fr4));

                    output_str.push('\n');
                }
                Err(e) => {
                    output_str.push_str(&r._name);
                    output_str.push('\t');
                    output_str.push_str(
                        std::str::from_utf8(&converted_sequence)
                            .expect("Non-UTF8 character in sequence"),
                    );
                    output_str.push('\t');
                    output_str.push_str(&format!("Failed numbering {e}"));
                    output_str.push('\t');
                    output_str.push_str(&format!("{}", 0));
                    output_str.push('\t');
                    output_str.push('X');
                    for _ in 0..8 {
                        output_str.push('\t')
                    } // no regions
                    output_str.push('\n');
                }
            }
        }
    }
    fs::write(output_file, output_str).expect("Something went wrong writing to data file")
}

/// Attempts to find all antibody chains in a sequence
fn find_all_chains<'a>(
    query_sequence: &'a [u8],
    numbering_schemes: Vec<&'a NumberingScheme>,
) -> Vec<Result<NumberingOutput<'a>, &'static str>> {
    let mut chains_found: Vec<Result<NumberingOutput, &str>> = Vec::new();

    let full_query_length: u32 = query_sequence.len() as u32;
    let mut sequence_list: Vec<(&[u8], u32, u32)> = Vec::new();
    sequence_list.push((query_sequence, 0u32, full_query_length - 1));

    while let Some(item) = sequence_list.pop() {
        let (current_sequence, current_start, current_end) = item;

        let numbering_result: Result<NumberingOutput, &str> =
            find_highest_identity_chain(current_sequence, &numbering_schemes);

        match numbering_result {
            Ok(mut best_chain) => {
                if best_chain.identity > MINIMAL_CHAIN_IDENTITY {
                    // split the remaining sequence, add sequences that are long enough to the list
                    let front_sequence: &[u8] = &best_chain.sequence[0..best_chain.start as usize];
                    let end_sequence: &[u8] = &best_chain.sequence[(best_chain.end as usize + 1)..];

                    if front_sequence.len() > MINIMAL_CHAIN_LENGTH as usize {
                        sequence_list.push((front_sequence, current_start, (best_chain.start - 1)))
                    }
                    if end_sequence.len() > MINIMAL_CHAIN_LENGTH as usize {
                        sequence_list.push((end_sequence, (best_chain.end + 1), current_end))
                    }

                    // set sequence to full original sequence
                    best_chain.sequence = query_sequence;
                    best_chain.start += current_start;
                    best_chain.end += current_start;

                    //add gaps to front and end to match with length of original sequence
                    let mut start_addition = vec![String::from("-"); current_start as usize];
                    let end_addition =
                        vec![String::from("-"); (full_query_length - current_end - 1) as usize];
                    start_addition.extend(best_chain.numbering);
                    start_addition.extend(end_addition);

                    best_chain.numbering = start_addition;
                    chains_found.push(Ok(best_chain));
                }
            }
            Err(e) => chains_found.push(Err::<NumberingOutput, &str>(e)),
        }
    }

    chains_found // TODO convert to result?
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schemes::{get_imgt_heavy_scheme, get_imgt_lambda_scheme, get_kabat_kappa_scheme};
    use crate::types::Chain;

    //#[test]
    // fn number_fasta_file() {
    //     //r"C:\Antibody_Numbering\fastas\abpdseq_non_redundant.fasta"
    //     number_sequences_and_write_output(
    //         r"C:\Antibody_Numbering\fastas\abpdseq_non_redundant.fasta",
    //         //r"C:\Users/Siemen/immunum-rs/immunum/fixtures/test.fasta",
    //         Scheme::IMGT,
    //         &[Chain::IGH, Chain::IGK, Chain::IGL],
    //         r"C:\Users\Siemen\immunum-rs\immunum\fixtures\rust_output_imgt_regions.txt",
    //         true,
    //     );
    // }

    #[test]
    fn single_sequence_find_all() {
        let seq = "VLTQSPGTLSLSPGETAIISCRTSQYGSLAWYQQRPGQAPRLVIYSGSTRAAGIPDRFSGSRWGPDYNLTISNLESGDFGVYYCQQYEFFGQGTKVQVDIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLRSPVTKSFNRGEC".as_bytes();
        let mut schemes: Vec<&NumberingScheme> = Vec::new();
        let lambda_scheme = get_imgt_lambda_scheme();
        let kappa_scheme = get_kabat_kappa_scheme();
        let heavy_scheme = get_imgt_heavy_scheme();

        schemes.push(&lambda_scheme);
        schemes.push(&kappa_scheme);
        schemes.push(&heavy_scheme);

        let output = find_all_chains(seq, schemes);
        for o in output {
            println!("{:?}", o);
        }
    }

    #[test]
    fn test_correct_chain_identification() {
        let heavy_chain: &[u8] = "QVQLVQSGAVIKTPGSSVKISCRASGYNFRDYSIHWVRLIPDKGFEWIGWIKPLWGAVSYARQL\
        QGRVSMTRQLSQDPDDPDWGVAYMEFSGLTPADTAEYFCVRRGSCDYCGDFPWQYWCQGTVVVVSSASTKGPSVFPLAPSSGGTAALGCLV\
        KDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK"
            .as_bytes();
        let lambda_chain: &[u8] = "SALTQPPSASGSLGQSVTISCTGTSSDVGGYNYVSWYQQHAGKAPKVIIYEVNKRPSGVPDRF\
        SGSKSGNTASLTVSGLQAEDEADYYCSSYEGSDNFVFGTGTKVTVLGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWK\
        ADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
            .as_bytes();
        let kappa_chain: &[u8] = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTG\
        SGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSE\
        RQNGVLNSATDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
            .as_bytes();
        let mut schemes: Vec<&NumberingScheme> = Vec::new();
        let lambda_scheme = get_imgt_lambda_scheme();
        let kappa_scheme = get_kabat_kappa_scheme();
        let heavy_scheme = get_imgt_heavy_scheme();

        schemes.push(&lambda_scheme);
        schemes.push(&kappa_scheme);
        schemes.push(&heavy_scheme);

        assert_eq!(
            find_highest_identity_chain(heavy_chain, &schemes)
                .expect("")
                .scheme
                .chain_type,
            Chain::IGH
        );

        assert_eq!(
            find_highest_identity_chain(kappa_chain, &schemes)
                .expect("")
                .scheme
                .chain_type,
            Chain::IGK
        );

        assert_eq!(
            find_highest_identity_chain(lambda_chain, &schemes)
                .expect("")
                .scheme
                .chain_type,
            Chain::IGL
        );
    }
}
