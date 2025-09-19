use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FastxError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Invalid FASTA/FASTQ format: {0}")]
    InvalidFormat(String),
}

/// Represents a sequence record from a FASTA or FASTQ file
#[derive(Debug, Clone)]
pub struct FastxRecord {
    pub name: String,
    pub sequence: String,
}

impl FastxRecord {
    pub fn new(name: String, sequence: String) -> Self {
        Self {
            name,
            sequence,
        }
    }
}

/// Iterator over FASTX records
pub struct FastxReader<R: BufRead> {
    reader: R,
    is_fastq: bool,
    line_buffer: String,
    finished: bool,
}

impl<R: BufRead> FastxReader<R> {
    pub fn new(mut reader: R) -> Result<Self, FastxError> {
        let mut line_buffer = String::new();
        reader.read_line(&mut line_buffer)?;

        let is_fastq = if line_buffer.starts_with('@') {
            true
        } else if line_buffer.starts_with('>') {
            false
        } else {
            return Err(FastxError::InvalidFormat(
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

    fn prepare_next_record(&mut self) -> Result<(), FastxError> {
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            self.finished = true;
        }
        Ok(())
    }

    fn read_fasta_sequence(&mut self) -> Result<String, FastxError> {
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

    fn read_fastq_sequence(&mut self) -> Result<String, FastxError> {
        // Read sequence line
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            return Err(FastxError::InvalidFormat(
                "Unexpected end of file while reading sequence".to_string(),
            ));
        }
        let sequence = self.line_buffer.trim().to_string();

        // Read and skip separator line (no validation)
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            return Err(FastxError::InvalidFormat(
                "Unexpected end of file while reading separator".to_string(),
            ));
        }

        // Read and skip quality line (no validation)
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            return Err(FastxError::InvalidFormat(
                "Unexpected end of file while reading quality".to_string(),
            ));
        }

        Ok(sequence)
    }

    fn read_record(&mut self) -> Result<Option<FastxRecord>, FastxError> {
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

        Ok(Some(FastxRecord::new(name, sequence)))
    }

}

impl<R: BufRead> Iterator for FastxReader<R> {
    type Item = Result<FastxRecord, FastxError>;

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

/// Creates a FastxReader from a file path, automatically handling gzipped files
pub fn from_path<P: AsRef<Path>>(path: P) -> Result<FastxReader<Box<dyn BufRead>>, FastxError> {
    let path = path.as_ref();
    let file = File::open(path)?;

    let reader: Box<dyn BufRead> = if is_gzipped(path) {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    FastxReader::new(reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn from_reader<R: BufRead + 'static>(
        reader: R,
    ) -> Result<FastxReader<Box<dyn BufRead>>, FastxError> {
        FastxReader::new(Box::new(reader))
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
}
