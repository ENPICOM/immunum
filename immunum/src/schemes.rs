use crate::types::{NumberingScheme, Scheme, Chain, RegionRange};
use crate::consensus_scoring::read_consensus_file;


pub fn get_imgt_heavy_scheme() -> NumberingScheme{
    NumberingScheme{
        name: "IMGT Heavy".to_string(),
        description:"IMGT numbering scheme for heavy chains".to_string(),
        scheme_type:Scheme::IMGT,
        chain_type:Chain::IGH,
        conserved_positions:vec![23, 41, 104, 118, 119, 121],
        insertion_positions:vec![],
        gap_positions:vec![10, 73],
        consensus_amino_acids:read_consensus_file(r"C:\Anti_Num\numbering\consensus\IMGT_CONSENSUS_H.txt"),
        fr1:RegionRange{start:1, end:27},
        cdr1:RegionRange{start:27, end:39},
        fr2:RegionRange{start:39, end:56},
        cdr2:RegionRange{start:56, end:66},
        fr3:RegionRange{start:66, end:105},
        cdr3:RegionRange{start:105, end:118},
        fr4:RegionRange{start:118, end:129},
    }
}

#[cfg(test)]
mod tests {
    use wasm_bindgen::__rt::assert_not_null;
    use crate::constants::{CDR1_INSERTION_POSITION_IMGT, GAP_PEN_CDR, GAP_PEN_FR, GAP_PEN_OP, GAP_PEN_OTHER};
    use super::*;

    #[test]
    fn scheme_creation() {
        let scheme = get_imgt_heavy_scheme();
        assert_eq!(scheme.gap_positions, vec![10, 73]);
        assert_eq!(scheme.restricted_sites().len(), 88);
        assert_eq!(scheme.consensus_amino_acids.len(), 128);
        assert_eq!(scheme.consensus_amino_acids[&1], vec!['Q', 'E', 'D']);
        assert_eq!(scheme.gap_penalty(25), (GAP_PEN_FR, GAP_PEN_FR));
        assert_eq!(scheme.gap_penalty(200), (GAP_PEN_OTHER, GAP_PEN_OTHER));
        assert_eq!(scheme.gap_penalty(CDR1_INSERTION_POSITION_IMGT), (GAP_PEN_CDR, GAP_PEN_CDR));
        assert_eq!(scheme.gap_penalty(10), (GAP_PEN_FR, GAP_PEN_OP))
        //TODO add more tests
    }
}