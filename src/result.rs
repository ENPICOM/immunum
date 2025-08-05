use crate::types::{Chain, RegionRange, Scheme};
use clap::ValueEnum;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

/// Output format options for AnnotationResult
#[derive(Debug, Clone, PartialEq, ValueEnum)]
pub enum OutputFormat {
    /// Tab-separated values with full details (legacy format)
    /// Format: Name\tSequence\tNumbering\tScore\tChain\tcdr1\tcdr2\tcdr3\tfr1\tfr2\tfr3\tfr4\tStart\tEnd
    Tsv,
    /// Simple format with identity and numbering
    /// Format: Identity: X.XX%, Numbering: [...]
    Simple,
    /// JSON format for structured data
    Json,
    /// CSV format for spreadsheet compatibility
    Csv,
}

/// Result of numbering a sequence, containing all relevant information
#[derive(Debug)]
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
    /// Start position of the numbered region in the original sequence
    pub start: u32,
    /// End position of the numbered region in the original sequence
    pub end: u32,
    /// Region definitions for CDR/FR extraction
    pub cdr1: RegionRange,
    pub cdr2: RegionRange,
    pub cdr3: RegionRange,
    pub fr1: RegionRange,
    pub fr2: RegionRange,
    pub fr3: RegionRange,
    pub fr4: RegionRange,
}

impl AnnotationResult {
    /// Get slice of sequence corresponding to a framework or cdr region
    pub fn get_query_region(&self, region: &RegionRange) -> &str {
        let mut start_position = region.start;
        let mut end_position = region.end - 1;

        // get final position
        let mut final_position: usize = 0;
        for p in 0..self.numbers.len() {
            if self.numbers[p] != "-" {
                final_position = p
            }
        }

        let mut start_sequence: Option<usize> = None;
        let mut end_sequence: Option<usize> = None;

        // find start position in sequence
        while start_position <= end_position {
            let pos_option: Option<usize> = self
                .numbers
                .iter()
                .position(|num| num == &start_position.to_string());
            match pos_option {
                Some(i) => {
                    start_sequence = Some(i);
                    break;
                }
                None => {
                    start_position += 1;
                }
            }
        }
        // check if start is found, if not, region is not present
        let start_index = match start_sequence {
            Some(i) => i,
            None => return "",
        };

        // find end position
        while end_position < u32::MAX {
            // Replace with reasonable upper bound
            let pos_option: Option<usize> = self
                .numbers
                .iter()
                .position(|num| num == &end_position.to_string());
            match pos_option {
                Some(i) => {
                    end_sequence = Some(i);
                    break;
                }
                None => {
                    end_position += 1;
                }
            }
        }
        let end_index = end_sequence.unwrap_or(final_position);

        &self.sequence[start_index..=end_index]
    }

    /// Populate the regions HashMap with calculated positions
    pub fn populate_regions(&mut self) {
        let regions_to_calc = [
            ("cdr1", &self.cdr1),
            ("cdr2", &self.cdr2),
            ("cdr3", &self.cdr3),
            ("fr1", &self.fr1),
            ("fr2", &self.fr2),
            ("fr3", &self.fr3),
            ("fr4", &self.fr4),
        ];

        for (name, region) in regions_to_calc {
            let region_seq = self.get_query_region(region);
            if !region_seq.is_empty() {
                if let Some(pos) = self.sequence.find(region_seq) {
                    self.regions
                        .insert(name.to_string(), (pos, pos + region_seq.len()));
                }
            }
        }
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

    /// Convert the result to a string in the specified format
    pub fn to_string(&self, format: OutputFormat) -> String {
        match format {
            OutputFormat::Tsv => self.to_tsv(),
            OutputFormat::Simple => self.to_simple(),
            OutputFormat::Json => self.to_json(),
            OutputFormat::Csv => self.to_csv(),
        }
    }

    /// Write the result to a file in the specified format
    pub fn to_file<P: AsRef<Path>>(&self, path: P, format: OutputFormat) -> Result<(), String> {
        let content = self.to_string(format);
        fs::write(path, content).map_err(|e| format!("Failed to write file: {}", e))
    }

    /// Get the legacy TSV output string (for backward compatibility)
    pub fn get_output_string(&self) -> String {
        self.to_tsv()
    }

    /// Convert to TSV format (tab-separated values)
    fn to_tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.sequence,
            self.numbers.join(","),
            self.identity,
            self.chain_code(),
            self.get_query_region(&self.cdr1),
            self.get_query_region(&self.cdr2),
            self.get_query_region(&self.cdr3),
            self.get_query_region(&self.fr1),
            self.get_query_region(&self.fr2),
            self.get_query_region(&self.fr3),
            self.get_query_region(&self.fr4),
            self.start,
            self.end
        )
    }

    /// Convert to simple format (identity and numbering only)
    fn to_simple(&self) -> String {
        format!(
            "Identity: {:.2}%, Numbering: {:?}",
            self.identity * 100.0,
            self.numbers
        )
    }

    /// Convert to JSON format
    fn to_json(&self) -> String {
        format!(
            r#"{{"sequence": "{}", "numbers": {:?}, "scheme": "{:?}", "chain": "{:?}", "identity": {:.3}, "regions": {{"cdr1": "{}", "cdr2": "{}", "cdr3": "{}", "fr1": "{}", "fr2": "{}", "fr3": "{}", "fr4": "{}"}}, "start": {}, "end": {}}}"#,
            self.sequence,
            self.numbers,
            self.scheme,
            self.chain,
            self.identity,
            self.get_query_region(&self.cdr1),
            self.get_query_region(&self.cdr2),
            self.get_query_region(&self.cdr3),
            self.get_query_region(&self.fr1),
            self.get_query_region(&self.fr2),
            self.get_query_region(&self.fr3),
            self.get_query_region(&self.fr4),
            self.start,
            self.end
        )
    }

    /// Convert to CSV format
    fn to_csv(&self) -> String {
        format!(
            "\"{}\",\"{}\",{:.3},\"{:?}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",{},{}\n",
            self.sequence,
            self.numbers.join(","),
            self.identity,
            self.chain,
            self.get_query_region(&self.cdr1),
            self.get_query_region(&self.cdr2),
            self.get_query_region(&self.cdr3),
            self.get_query_region(&self.fr1),
            self.get_query_region(&self.fr2),
            self.get_query_region(&self.fr3),
            self.get_query_region(&self.fr4),
            self.start,
            self.end
        )
    }

    /// Get single-character chain code
    fn chain_code(&self) -> &'static str {
        match self.chain {
            Chain::IGH => "H",
            Chain::IGK => "K",
            Chain::IGL => "L",
            Chain::TRA => "A",
            Chain::TRB => "B",
            Chain::TRD => "D",
            Chain::TRG => "G",
        }
    }
}

impl OutputFormat {
    /// Get the appropriate file extension for this format
    pub fn file_extension(&self) -> &'static str {
        match self {
            OutputFormat::Tsv => "tsv",
            OutputFormat::Simple => "txt",
            OutputFormat::Json => "json",
            OutputFormat::Csv => "csv",
        }
    }

    /// Get the header line for tabular formats
    pub fn header(&self) -> Option<String> {
        match self {
            OutputFormat::Tsv => Some("Sequence\tNumbering\tIdentity\tChain\tCDR1\tCDR2\tCDR3\tFR1\tFR2\tFR3\tFR4\tStart\tEnd".to_string()),
            OutputFormat::Csv => Some("Sequence,Numbering,Identity,Chain,CDR1,CDR2,CDR3,FR1,FR2,FR3,FR4,Start,End".to_string()),
            OutputFormat::Simple | OutputFormat::Json => None,
        }
    }
}

/// Utility functions for working with multiple AnnotationResults
impl AnnotationResult {
    /// Write multiple results to a file with optional header
    pub fn write_batch_to_file<P: AsRef<Path>>(
        results: &[AnnotationResult],
        path: P,
        format: OutputFormat,
        include_header: bool,
    ) -> Result<(), String> {
        let mut content = String::new();

        // Add header if requested and format supports it
        if include_header {
            if let Some(header) = format.header() {
                content.push_str(&header);
                content.push('\n');
            }
        }

        // Add each result
        for result in results {
            content.push_str(&result.to_string(format.clone()));
        }

        fs::write(path, content).map_err(|e| format!("Failed to write file: {}", e))
    }

    /// Convert multiple results to a single string with optional header
    pub fn batch_to_string(
        results: &[AnnotationResult],
        format: OutputFormat,
        include_header: bool,
    ) -> String {
        let mut content = String::new();

        // Add header if requested and format supports it
        if include_header {
            if let Some(header) = format.header() {
                content.push_str(&header);
                content.push('\n');
            }
        }

        // Add each result
        for result in results {
            content.push_str(&result.to_string(format.clone()));
        }

        content
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
            start: 0,
            end: 39,
            cdr1: RegionRange { start: 27, end: 38 },
            cdr2: RegionRange { start: 56, end: 65 },
            cdr3: RegionRange {
                start: 105,
                end: 117,
            },
            fr1: RegionRange { start: 1, end: 26 },
            fr2: RegionRange { start: 39, end: 55 },
            fr3: RegionRange {
                start: 66,
                end: 104,
            },
            fr4: RegionRange {
                start: 118,
                end: 128,
            },
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
            start: 0,
            end: 7,
            cdr1: RegionRange { start: 27, end: 38 },
            cdr2: RegionRange { start: 56, end: 65 },
            cdr3: RegionRange {
                start: 105,
                end: 117,
            },
            fr1: RegionRange { start: 1, end: 26 },
            fr2: RegionRange { start: 39, end: 55 },
            fr3: RegionRange {
                start: 66,
                end: 104,
            },
            fr4: RegionRange {
                start: 118,
                end: 128,
            },
        };

        assert_eq!(
            result.get_region_sequence("cdr1"),
            Some("ATCGA".to_string())
        );
        assert_eq!(result.get_region_sequence("cdr2"), None);
    }

    #[test]
    fn test_output_formats() {
        let result = create_test_result();

        // Test TSV format
        let tsv = result.to_string(OutputFormat::Tsv);
        assert!(tsv.contains("\t"));
        assert!(tsv.contains("ATCGATCG"));

        // Test Simple format
        let simple = result.to_string(OutputFormat::Simple);
        assert!(simple.contains("Identity: 95.0%"));
        assert!(simple.contains("Numbering:"));

        // Test JSON format
        let json = result.to_string(OutputFormat::Json);
        assert!(json.contains("\"sequence\""));
        assert!(json.contains("\"identity\""));

        // Test CSV format
        let csv = result.to_string(OutputFormat::Csv);
        assert!(csv.contains(","));
        assert!(csv.contains("\"ATCGATCG\""));
    }

    #[test]
    fn test_output_format_extensions() {
        assert_eq!(OutputFormat::Tsv.file_extension(), "tsv");
        assert_eq!(OutputFormat::Json.file_extension(), "json");
        assert_eq!(OutputFormat::Csv.file_extension(), "csv");
        assert_eq!(OutputFormat::Simple.file_extension(), "txt");
    }

    #[test]
    fn test_headers() {
        assert!(OutputFormat::Tsv.header().is_some());
        assert!(OutputFormat::Csv.header().is_some());
        assert!(OutputFormat::Simple.header().is_none());
        assert!(OutputFormat::Json.header().is_none());
    }

    fn create_test_result() -> AnnotationResult {
        let mut regions = HashMap::new();
        regions.insert("cdr1".to_string(), (10, 20));

        AnnotationResult {
            sequence: "ATCGATCG".to_string(),
            numbers: vec!["1".to_string(), "2".to_string(), "3".to_string()],
            scheme: Scheme::IMGT,
            chain: Chain::IGH,
            identity: 0.95,
            regions,
            start: 0,
            end: 7,
            cdr1: RegionRange { start: 27, end: 38 },
            cdr2: RegionRange { start: 56, end: 65 },
            cdr3: RegionRange {
                start: 105,
                end: 117,
            },
            fr1: RegionRange { start: 1, end: 26 },
            fr2: RegionRange { start: 39, end: 55 },
            fr3: RegionRange {
                start: 66,
                end: 104,
            },
            fr4: RegionRange {
                start: 118,
                end: 128,
            },
        }
    }
}
