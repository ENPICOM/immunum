use crate::fastx::from_path;
use crate::numbering::number_sequence;
use crate::schemes::{get_imgt_heavy_scheme, get_imgt_lambda_scheme};
use crate::types::{NumberingOutput, NumberingScheme, Scheme};
use std::fs;

fn find_highest_identity_chain<'a>(
    query_sequence: &'a String,
    numbering_schemes: &'a Vec<NumberingScheme>,
) -> NumberingOutput<'a> {
    let mut results: Vec<NumberingOutput> = Vec::with_capacity(numbering_schemes.len());

    for scheme in numbering_schemes {
        results.push(scheme.number_sequence(&query_sequence));
    }
    //Determine highest scoring model
    results
        .iter()
        .max_by(|a, b| {
            a.identity
                .partial_cmp(&b.identity)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap()
        .clone()
}

pub fn number_sequences_and_write_output(fasta_file: &str, scheme: Scheme, output_file: &str) {
    // Read in fasta
    let reader = from_path(fasta_file).unwrap();
    let records: Result<Vec<_>, _> = reader.collect();
    let records = records.unwrap();

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

    let mut output_str = "".to_string();
    // run annotation for all sequences
    for r in records {
        println!("{}", r._name);
        let output = find_highest_identity_chain(&r.sequence, &schemes);

        // create string output
        // output_str.push_str(&*r._name);
        // output_str.push_str(&"\t");
        // output_str.push_str(&output.sequence);
        // output_str.push_str(&"\t");
        // output_str.push_str(&output.numbering.join(","));
        // output_str.push_str(&"\t");
        // output_str.push_str(&format!("{}", output.identity));
        // output_str.push_str(&"\n");
        // TODO add regions

        // fill in value
    }
    //fs::write(&output_file, output_str).expect("Should be able to write to `/foo/tmp`")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schemes::{get_imgt_heavy_scheme, get_imgt_lambda_scheme, get_kabat_kappa_scheme};
    use crate::types::Chain;

    #[test]
    fn number_fasta_file() {
        // r"C:\Antibody_Numbering\fastas\abpdseq_non_redundant.fasta",
        number_sequences_and_write_output(
            r"C:\Users\Siemen\immunum-rs\immunum\fixtures\test_bigger.fasta",
            Scheme::IMGT,
            r"C:\Users\Siemen\immunum-rs\immunum\fixtures\rust_output.txt",
        )
    }

    #[test]
    fn test_correct_chain_identification() {
        let heavy_chain: String =
            "QVQLVQSGAVIKTPGSSVKISCRASGYNFRDYSIHWVRLIPDKGFEWIGWIKPLWGAVSYARQL\
        QGRVSMTRQLSQDPDDPDWGVAYMEFSGLTPADTAEYFCVRRGSCDYCGDFPWQYWCQGTVVVVSSASTKGPSVFPLAPSSGGTAALGCLV\
        KDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK"
                .to_string();
        let lambda_chain: String =
            "SALTQPPSASGSLGQSVTISCTGTSSDVGGYNYVSWYQQHAGKAPKVIIYEVNKRPSGVPDRF\
        SGSKSGNTASLTVSGLQAEDEADYYCSSYEGSDNFVFGTGTKVTVLGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWK\
        ADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
                .to_string();
        let kappa_chain: String =
            "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTG\
        SGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSE\
        RQNGVLNSATDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
                .to_string();
        let schemes = vec![
            get_imgt_lambda_scheme(),
            get_kabat_kappa_scheme(),
            get_imgt_heavy_scheme(),
        ];
        assert_eq!(
            find_highest_identity_chain(&heavy_chain, &schemes)
                .scheme
                .chain_type,
            Chain::IGH
        );
        assert_eq!(
            find_highest_identity_chain(&kappa_chain, &schemes)
                .scheme
                .chain_type,
            Chain::IGK
        );
        assert_eq!(
            find_highest_identity_chain(&lambda_chain, &schemes)
                .scheme
                .chain_type,
            Chain::IGL
        );
        println!(
            "{:?}",
            find_highest_identity_chain(&heavy_chain, &schemes).numbering
        );
        println!(
            "{:?}",
            find_highest_identity_chain(&lambda_chain, &schemes).numbering
        );
        println!(
            "{:?}",
            find_highest_identity_chain(&kappa_chain, &schemes).numbering
        );
    }
}
