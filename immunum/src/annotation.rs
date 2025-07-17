use crate::consensus_scoring::write_all_scoring_matrices;
use crate::fastx::{from_path, FastxRecord};
use crate::numbering_scheme_type::{NumberingOutput, NumberingScheme};
use crate::schemes::{get_imgt_heavy_scheme, get_imgt_lambda_scheme};
use crate::types::{Chain, Scheme};
use std::fs;

pub(crate) fn find_highest_identity_chain<'a>(
    query_sequence: &'a [u8],
    numbering_schemes: &'a Vec<NumberingScheme>,
) -> NumberingOutput<'a> {
    let mut results: Vec<NumberingOutput> = Vec::with_capacity(numbering_schemes.len());

    for scheme in numbering_schemes {
        results.push(scheme.number_sequence(query_sequence));
    }
    //Determine highest scoring model
    results
        .into_iter()
        .max_by(|a, b| {
            a.identity
                .partial_cmp(&b.identity)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .expect("No numbering scheme selected")
}

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
            get_imgt_heavy_scheme(),
            get_imgt_lambda_scheme(),
        ],
        Scheme::KABAT => vec![
            get_imgt_heavy_scheme(),
            get_imgt_heavy_scheme(),
            get_imgt_lambda_scheme(),
        ],
    };

    // filter only chains selected
    let schemes: Vec<NumberingScheme> = schemes
        .into_iter()
        .filter(|scheme| chains.contains(&scheme.chain_type))
        .collect();

    let mut output_str = "".to_string();
    // run annotation for all sequences
    for r in records {
        println!("{}", r._name);
        let converted_sequence = r.sequence.into_bytes();
        let output = find_highest_identity_chain(&converted_sequence, &schemes);

        //create string output
        output_str.push_str(&r._name);
        output_str.push('\t');
        output_str.push_str(
            std::str::from_utf8(output.sequence).expect("Non-UTF8 character in sequence"),
        );
        output_str.push('\t');
        output_str.push_str(&output.numbering.join(","));
        output_str.push('\t');
        output_str.push_str(&format!("{}", output.identity));
        output_str.push('\n');
        // TODO add regions

        // TODO remove this temporary print statement using unused variables
        println!(
            "Found {0}{1} from {2} to {3}",
            output.scheme.name, output.scheme.description, output.start, output.end
        );
        // fill in value
    }
    fs::write(output_file, output_str).expect("Should be able to write to `/foo/tmp`")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schemes::{get_imgt_heavy_scheme, get_imgt_lambda_scheme, get_kabat_kappa_scheme};
    use crate::types::Chain;

    // #[test]
    // fn number_fasta_file() {
    //     //r"C:\Antibody_Numbering\fastas\abpdseq_non_redundant.fasta"
    //     number_sequences_and_write_output(
    //         r"C:\Antibody_Numbering\fastas\abpdseq_non_redundant.fasta",
    //         Scheme::IMGT,
    //         &[Chain::IGH, Chain::IGK, Chain::IGL],
    //         r"C:\Users\Siemen\immunum-rs\immunum\fixtures\rust_output.txt",
    //         false,
    //     );
    // }

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
        let schemes = vec![
            get_imgt_lambda_scheme(),
            get_kabat_kappa_scheme(),
            get_imgt_heavy_scheme(),
        ];
        assert_eq!(
            find_highest_identity_chain(heavy_chain, &schemes)
                .scheme
                .chain_type,
            Chain::IGH
        );
        assert_eq!(
            find_highest_identity_chain(kappa_chain, &schemes)
                .scheme
                .chain_type,
            Chain::IGK
        );
        assert_eq!(
            find_highest_identity_chain(lambda_chain, &schemes)
                .scheme
                .chain_type,
            Chain::IGL
        );
        println!(
            "{:?}",
            find_highest_identity_chain(heavy_chain, &schemes).numbering
        );
        println!(
            "{:?}",
            find_highest_identity_chain(lambda_chain, &schemes).numbering
        );
        println!(
            "{:?}",
            find_highest_identity_chain(kappa_chain, &schemes).numbering
        );
    }
}
