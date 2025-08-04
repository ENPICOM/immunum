use crate::numbering_scheme_type::NumberingOutput;
use crate::types::{Chain, Scheme};
use std::collections::HashMap;

/// Result of numbering a sequence, containing all relevant information
pub struct AnnotationResult {
    /// The original sequence that was numbered
    pub sequence: String,
    /// The numbering for each position in the sequence
    pub numbers: Vec<String>,
    /// The numbering scheme used
    pub scheme: Scheme,
    /// The chain type that was matched
    pub chain: Chain,
    /// The identity/confidence score (0.0 to 1.0)
    pub identity: f64,
    /// Regions mapped to their start and end positions in the sequence
    pub regions: HashMap<String, (usize, usize)>,
}

impl AnnotationResult {
    /// Create an AnnotationResult from a NumberingOutput
    pub fn from_numbering_output(
        output: NumberingOutput,
        sequence: String,
    ) -> Result<Self, String> {
        // Extract chain type from the scheme
        let chain = output.scheme.chain_type;
        
        // Convert scheme type (we'll need to add a method to get scheme type from NumberingScheme)
        let scheme = match output.scheme.chain_type {
            Chain::IGH | Chain::IGK | Chain::IGL => {
                // We need to infer the scheme type from the numbering scheme
                // For now, we'll use IMGT as default - TODO: add scheme type to NumberingScheme
                Scheme::IMGT
            }
            _ => Scheme::IMGT,
        };
        
        // Extract regions from the numbering scheme
        let mut regions = HashMap::new();
        
        // Extract CDR and FR regions based on the scheme's region definitions
        let cdr1_seq = output.get_query_region(&output.scheme.cdr1);
        let cdr2_seq = output.get_query_region(&output.scheme.cdr2);
        let cdr3_seq = output.get_query_region(&output.scheme.cdr3);
        let fr1_seq = output.get_query_region(&output.scheme.fr1);
        let fr2_seq = output.get_query_region(&output.scheme.fr2);
        let fr3_seq = output.get_query_region(&output.scheme.fr3);
        let fr4_seq = output.get_query_region(&output.scheme.fr4);
        
        // Find positions in the original sequence
        // This is a simplified approach - in practice, we'd need more sophisticated region mapping
        if !cdr1_seq.is_empty() {
            if let Some(pos) = sequence.find(cdr1_seq) {
                regions.insert("cdr1".to_string(), (pos, pos + cdr1_seq.len()));
            }
        }
        if !cdr2_seq.is_empty() {
            if let Some(pos) = sequence.find(cdr2_seq) {
                regions.insert("cdr2".to_string(), (pos, pos + cdr2_seq.len()));
            }
        }
        if !cdr3_seq.is_empty() {
            if let Some(pos) = sequence.find(cdr3_seq) {
                regions.insert("cdr3".to_string(), (pos, pos + cdr3_seq.len()));
            }
        }
        if !fr1_seq.is_empty() {
            if let Some(pos) = sequence.find(fr1_seq) {
                regions.insert("fr1".to_string(), (pos, pos + fr1_seq.len()));
            }
        }
        if !fr2_seq.is_empty() {
            if let Some(pos) = sequence.find(fr2_seq) {
                regions.insert("fr2".to_string(), (pos, pos + fr2_seq.len()));
            }
        }
        if !fr3_seq.is_empty() {
            if let Some(pos) = sequence.find(fr3_seq) {
                regions.insert("fr3".to_string(), (pos, pos + fr3_seq.len()));
            }
        }
        if !fr4_seq.is_empty() {
            if let Some(pos) = sequence.find(fr4_seq) {
                regions.insert("fr4".to_string(), (pos, pos + fr4_seq.len()));
            }
        }
        
        Ok(AnnotationResult {
            sequence,
            numbers: output.numbering,
            scheme,
            chain,
            identity: output.identity,
            regions,
        })
    }
    
    /// Get the sequence for a specific region
    pub fn get_region_sequence(&self, region_name: &str) -> Option<String> {
        if let Some((start, end)) = self.regions.get(region_name) {
            if *end <= self.sequence.len() {
                Some(self.sequence[*start..*end].to_string())
            } else {
                None
            }
        } else {
            None
        }
    }
    
    /// Get all CDR sequences as a map
    pub fn get_cdr_sequences(&self) -> HashMap<String, String> {
        let mut cdrs = HashMap::new();
        
        for region in ["cdr1", "cdr2", "cdr3"] {
            if let Some(seq) = self.get_region_sequence(region) {
                cdrs.insert(region.to_string(), seq);
            }
        }
        
        cdrs
    }
    
    /// Get all framework sequences as a map
    pub fn get_framework_sequences(&self) -> HashMap<String, String> {
        let mut frameworks = HashMap::new();
        
        for region in ["fr1", "fr2", "fr3", "fr4"] {
            if let Some(seq) = self.get_region_sequence(region) {
                frameworks.insert(region.to_string(), seq);
            }
        }
        
        frameworks
    }
    
    /// Check if the annotation result meets a minimum identity threshold
    pub fn is_high_confidence(&self, threshold: f64) -> bool {
        self.identity >= threshold
    }
    
    /// Get a summary string of the annotation
    pub fn summary(&self) -> String {
        format!(
            "Chain: {:?}, Scheme: {:?}, Identity: {:.1}%, Regions: {}",
            self.chain,
            self.scheme,
            self.identity * 100.0,
            self.regions.len()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Chain, Scheme};
    
    #[test]
    fn test_annotation_result_creation() {
        let mut regions = HashMap::new();
        regions.insert("cdr1".to_string(), (10, 20));
        regions.insert("cdr2".to_string(), (40, 50));
        
        let result = AnnotationResult {
            sequence: "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
            numbers: vec!["1".to_string(), "2".to_string(), "3".to_string()],
            scheme: Scheme::IMGT,
            chain: Chain::IGH,
            identity: 0.95,
            regions,
        };
        
        assert_eq!(result.identity, 0.95);
        assert_eq!(result.chain, Chain::IGH);
        assert!(result.is_high_confidence(0.9));
        assert!(!result.is_high_confidence(0.99));
    }
    
    #[test]
    fn test_region_sequence_extraction() {
        let mut regions = HashMap::new();
        regions.insert("cdr1".to_string(), (0, 5));
        
        let result = AnnotationResult {
            sequence: "ATCGATCG".to_string(),
            numbers: vec!["1".to_string(), "2".to_string()],
            scheme: Scheme::IMGT,
            chain: Chain::IGH,
            identity: 0.95,
            regions,
        };
        
        assert_eq!(result.get_region_sequence("cdr1"), Some("ATCGA".to_string()));
        assert_eq!(result.get_region_sequence("cdr2"), None);
    }
}