use crate::annotation::{number_sequences_and_write_output};
use crate::types::{Chain, Scheme};

// TODO To be changed, think about entry point for user
pub fn number_sequence(fasta_file: &str, scheme: &Scheme, chains: &[Chain]) -> String{
    // TODO Temporary code, prevents any unused variable warnings by clippy
    number_sequences_and_write_output(fasta_file, scheme.clone(), chains, "numbering_output.txt", true);
    "Numbered sequences and stored in numbering_output.txt".to_string()
}
