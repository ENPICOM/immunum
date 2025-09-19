use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SequenceError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Invalid FASTA/FASTQ format: {0}")]
    InvalidFormat(String),
}

/// Represents a sequence record from a FASTA or FASTQ file
#[derive(Debug, Clone)]
pub struct SequenceRecord {
    pub name: String,
    pub sequence: String,
}

impl SequenceRecord {
    pub fn new(name: String, sequence: String) -> Self {
        Self {
            name,
            sequence,
        }
    }
}

/// Iterator over sequence records
pub struct SequenceReader<R: BufRead> {
    reader: R,
    is_fastq: bool,
    line_buffer: String,
    finished: bool,
}

impl<R: BufRead> SequenceReader<R> {
    pub fn new(mut reader: R) -> Result<Self, SequenceError> {
        let mut line_buffer = String::new();
        reader.read_line(&mut line_buffer)?;

        let is_fastq = if line_buffer.starts_with('@') {
            true
        } else if line_buffer.starts_with('>') {
            false
        } else {
            return Err(SequenceError::InvalidFormat(
                "File must start with '>' (FASTA) or '@' (FASTQ)".to_string(),
            ));
        };

        Ok(Self {
            reader,
            is_fastq,
            line_buffer,
            finished: false,
        })
    }

    fn prepare_next_record(&mut self) -> Result<(), SequenceError> {
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            self.finished = true;
        }
        Ok(())
    }

    fn read_fasta_sequence(&mut self) -> Result<String, SequenceError> {
        let mut sequence = String::new();

        // Read sequence lines until next header or EOF
        loop {
            self.line_buffer.clear();
            if self.reader.read_line(&mut self.line_buffer)? == 0 {
                // End of file
                self.finished = true;
                break;
            }

            if self.line_buffer.starts_with('>') {
                // Found next record header, keep it in buffer for next iteration
                break;
            }

            sequence.push_str(self.line_buffer.trim());
        }

        Ok(sequence)
    }

    fn read_fastq_sequence(&mut self) -> Result<String, SequenceError> {
        // Read sequence line
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            return Err(SequenceError::InvalidFormat(
                "Unexpected end of file while reading sequence".to_string(),
            ));
        }
        let sequence = self.line_buffer.trim().to_string();

        // Read and skip separator line (no validation)
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            return Err(SequenceError::InvalidFormat(
                "Unexpected end of file while reading separator".to_string(),
            ));
        }

        // Read and skip quality line (no validation)
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            return Err(SequenceError::InvalidFormat(
                "Unexpected end of file while reading quality".to_string(),
            ));
        }

        Ok(sequence)
    }

    fn read_record(&mut self) -> Result<Option<SequenceRecord>, SequenceError> {
        if self.finished || self.line_buffer.is_empty() {
            return Ok(None);
        }

        // Extract name without header validation (trust the format detection)
        let name = self.line_buffer[1..].trim().to_string();

        // Read sequence using format-specific method
        let sequence = if self.is_fastq {
            self.read_fastq_sequence()?
        } else {
            self.read_fasta_sequence()?
        };

        // Prepare for next record (only needed for FASTQ since FASTA handles it internally)
        if self.is_fastq {
            self.prepare_next_record()?;
        }

        Ok(Some(SequenceRecord::new(name, sequence)))
    }

}

impl<R: BufRead> Iterator for SequenceReader<R> {
    type Item = Result<SequenceRecord, SequenceError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Determines if a file is gzipped based on its extension
fn is_gzipped(path: &Path) -> bool {
    path.to_string_lossy().to_lowercase().ends_with(".gz")
}

/// Creates a SequenceReader from a file path, automatically handling gzipped files
pub fn from_path<P: AsRef<Path>>(path: P) -> Result<SequenceReader<Box<dyn BufRead>>, SequenceError> {
    let path = path.as_ref();
    let file = File::open(path)?;

    let reader: Box<dyn BufRead> = if is_gzipped(path) {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    SequenceReader::new(reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn from_reader<R: BufRead + 'static>(
        reader: R,
    ) -> Result<SequenceReader<Box<dyn BufRead>>, SequenceError> {
        SequenceReader::new(Box::new(reader))
    }

    #[test]
    fn test_fasta_parsing() {
        let fasta_data = ">seq1 description\nATCG\nTGCA\n>seq2\nGGCC\n";
        let cursor = Cursor::new(fasta_data);
        let reader = from_reader(BufReader::new(cursor)).unwrap();

        let records: Result<Vec<_>, _> = reader.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "seq1 description");
        assert_eq!(records[0].sequence, "ATCGTGCA");
        assert_eq!(records[1].name, "seq2");
        assert_eq!(records[1].sequence, "GGCC");
    }

    #[test]
    fn test_fastq_parsing() {
        let fastq_data = "@seq1\nATCG\n+\nIIII\n@seq2\nGGCC\n+\n!!!!\n";
        let cursor = Cursor::new(fastq_data);
        let reader = from_reader(BufReader::new(cursor)).unwrap();

        let records: Result<Vec<_>, _> = reader.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "seq1");
        assert_eq!(records[0].sequence, "ATCG");
        assert_eq!(records[1].name, "seq2");
        assert_eq!(records[1].sequence, "GGCC");
    }

    #[test]
    fn test_fasta_from_file() {
        let reader = from_path("fixtures/test.fasta").unwrap();
        let records: Result<Vec<_>, _> = reader.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 4);

        // Test first record
        assert_eq!(
            records[0].name,
            "heavy_chain_1 Human IgG heavy chain variable region"
        );
        assert!(records[0]
            .sequence
            .starts_with("QVQLVQSGAEVKKPGASVKVSCKASGYTFTS"));

        // Test second record
        assert_eq!(
            records[1].name,
            "light_chain_1 Human kappa light chain variable region"
        );
        assert!(records[1]
            .sequence
            .starts_with("DIQMTQSPSSLSASVGDRVTITCRASQSISS"));

        // Test third record
        assert_eq!(
            records[2].name,
            "heavy_chain_2 Murine IgG heavy chain variable region"
        );
        assert!(records[2]
            .sequence
            .starts_with("EVQLLESGGGLVQPGGSLRLSCAASGFTFSS"));

        // Test fourth record
        assert_eq!(
            records[3].name,
            "light_chain_2 Human lambda light chain variable region"
        );
        assert!(records[3]
            .sequence
            .starts_with("QSALTQPASVSGSPGQSITISCTGTSSDVGG"));
    }

    #[test]
    fn test_fastq_from_file() {
        let reader = from_path("fixtures/test.fastq").unwrap();
        let records: Result<Vec<_>, _> = reader.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 3);

        // Test first record
        assert_eq!(
            records[0].name,
            "heavy_chain_1 Human IgG heavy chain variable region"
        );
        assert!(records[0]
            .sequence
            .starts_with("QVQLVQSGAEVKKPGASVKVSCKASGYTFTS"));

        // Test second record
        assert_eq!(
            records[1].name,
            "light_chain_1 Human kappa light chain variable region"
        );
        assert!(records[1]
            .sequence
            .starts_with("DIQMTQSPSSLSASVGDRVTITCRASQSISS"));

        // Test third record
        assert_eq!(
            records[2].name,
            "heavy_chain_2 Murine IgG heavy chain variable region"
        );
        assert!(records[2]
            .sequence
            .starts_with("EVQLLESGGGLVQPGGSLRLSCAASGFTFSS"));
    }

    #[test]
    fn test_fasta_gz_from_file() {
        let reader = from_path("fixtures/test.fasta.gz").unwrap();
        let records: Result<Vec<_>, _> = reader.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 4);

        // Test that gzipped content matches uncompressed content
        assert_eq!(
            records[0].name,
            "heavy_chain_1 Human IgG heavy chain variable region"
        );
        assert!(records[0]
            .sequence
            .starts_with("QVQLVQSGAEVKKPGASVKVSCKASGYTFTS"));

        assert_eq!(
            records[1].name,
            "light_chain_1 Human kappa light chain variable region"
        );
        assert!(records[1]
            .sequence
            .starts_with("DIQMTQSPSSLSASVGDRVTITCRASQSISS"));
    }

    #[test]
    fn test_fastq_gz_from_file() {
        let reader = from_path("fixtures/test.fastq.gz").unwrap();
        let records: Result<Vec<_>, _> = reader.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 3);

        // Test that gzipped content matches uncompressed content
        assert_eq!(
            records[0].name,
            "heavy_chain_1 Human IgG heavy chain variable region"
        );
        assert!(records[0]
            .sequence
            .starts_with("QVQLVQSGAEVKKPGASVKVSCKASGYTFTS"));

        assert_eq!(
            records[1].name,
            "light_chain_1 Human kappa light chain variable region"
        );
        assert!(records[1]
            .sequence
            .starts_with("DIQMTQSPSSLSASVGDRVTITCRASQSISS"));
    }

    #[test]
    fn test_nonexistent_file() {
        let result = from_path("fixtures/nonexistent.fasta");
        assert!(result.is_err());
    }

    #[test]
    fn test_empty_file() {
        use std::fs::File;
        use std::io::Write;

        // Create a temporary empty file
        let temp_path = "fixtures/temp_empty";
        File::create(temp_path).unwrap().write_all(b"").unwrap();

        let result = from_path(temp_path);
        assert!(result.is_err());

        // Clean up
        std::fs::remove_file(temp_path).unwrap();
    }

    #[test]
    fn test_sequence_stream_from_sequence() {
        let sequence = "ACDEFGHIKLMNPQRSTVWY";
        let stream = SequenceStream::from_sequence(sequence);

        let records: Result<Vec<_>, _> = stream.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "INPUT SEQUENCE");
        assert_eq!(records[0].sequence, sequence);
    }

    #[test]
    fn test_sequence_stream_new_with_sequence() {
        let sequence = "ACDEFGHIKLMNPQRSTVWY";
        let stream = SequenceStream::new(sequence).unwrap();

        let records: Result<Vec<_>, _> = stream.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, sequence);
    }

    #[test]
    fn test_sequence_stream_new_with_invalid_input() {
        let invalid_input = "not_a_file_and_not_a_sequence_123";
        let result = SequenceStream::new(invalid_input);

        assert!(result.is_err());
        match result {
            Err(SequenceError::InvalidFormat(msg)) => {
                assert!(msg.contains("neither a valid file path nor a valid sequence"));
            }
            _ => panic!("Expected InvalidFormat error"),
        }
    }

    #[test]
    fn test_is_valid_sequence() {
        assert!(is_valid_sequence("ACDEFGHIKLMNPQRSTVWY"));
        assert!(is_valid_sequence("ACG"));
        assert!(is_valid_sequence("")); // TODO - this should probably be invalid?
        assert!(!is_valid_sequence("ACGX"));
        assert!(!is_valid_sequence("123"));
        assert!(!is_valid_sequence("acg")); // lowercase not allowed
    }
}

/// Validates if a string looks like a biological sequence
fn is_valid_sequence(input: &str) -> bool {
    // Check if input string contains valid protein chars
    input.chars().all(|c| "ACDEFGHIKLMNPQRSTVWY".contains(c))
}

/// A wrapper around sequence iterators that can be created from files or direct sequence input
pub struct SequenceStream {
    inner: Box<dyn Iterator<Item = Result<SequenceRecord, SequenceError>>>,
}

impl SequenceStream {
    /// Creates a new SequenceStream by auto-detecting if the input is a file path or direct sequence
    pub fn new(input: &str) -> Result<Self, SequenceError> {
        let input_path = Path::new(input);

        if input_path.exists() {
            Self::from_file(input)
        } else if is_valid_sequence(input) {
            Ok(Self::from_sequence(input))
        } else {
            Err(SequenceError::InvalidFormat(
                "Input is neither a valid file path nor a valid sequence".to_string(),
            ))
        }
    }

    /// Creates a SequenceStream from a file path
    fn from_file(path: &str) -> Result<Self, SequenceError> {
        let reader = from_path(path)?;
        Ok(Self {
            inner: Box::new(reader),
        })
    }

    /// Creates a SequenceStream from a direct sequence string
    fn from_sequence(sequence: &str) -> Self {
        let record = SequenceRecord::new(
            "INPUT SEQUENCE".to_string(),
            sequence.to_string(),
        );
        Self {
            inner: Box::new(std::iter::once(Ok(record))),
        }
    }
}

impl Iterator for SequenceStream {
    type Item = Result<SequenceRecord, SequenceError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }
}
