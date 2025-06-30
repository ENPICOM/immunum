use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FastxError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Invalid FASTA/FASTQ format: {0}")]
    InvalidFormat(String),
    #[error("Unsupported file extension: {0}")]
    UnsupportedExtension(String),
}

/// Represents a sequence record from a FASTA or FASTQ file
#[derive(Debug, Clone)]
pub struct FastxRecord {
    pub name: String,
    pub sequence: String,
    pub quality: Option<String>, // Only present for FASTQ files
}

impl FastxRecord {
    pub fn new(name: String, sequence: String, quality: Option<String>) -> Self {
        Self {
            name,
            sequence,
            quality,
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

    fn read_fasta_record(&mut self) -> Result<Option<FastxRecord>, FastxError> {
        if self.finished || self.line_buffer.is_empty() {
            return Ok(None);
        }

        // Parse header line (should start with '>')
        if !self.line_buffer.starts_with('>') {
            return Err(FastxError::InvalidFormat(
                "FASTA header must start with '>'".to_string(),
            ));
        }

        let name = self.line_buffer[1..].trim().to_string();
        let mut sequence = String::new();

        // Read sequence lines until next header or EOF
        self.line_buffer.clear();
        while self.reader.read_line(&mut self.line_buffer)? > 0 {
            if self.line_buffer.starts_with('>') {
                // Found next record, don't consume this line
                break;
            }
            sequence.push_str(self.line_buffer.trim());
            self.line_buffer.clear();
        }

        if self.line_buffer.is_empty() {
            self.finished = true;
        }

        Ok(Some(FastxRecord::new(name, sequence, None)))
    }

    fn read_fastq_record(&mut self) -> Result<Option<FastxRecord>, FastxError> {
        if self.finished || self.line_buffer.is_empty() {
            return Ok(None);
        }

        // Parse header line (should start with '@')
        if !self.line_buffer.starts_with('@') {
            return Err(FastxError::InvalidFormat(
                "FASTQ header must start with '@'".to_string(),
            ));
        }

        let name = self.line_buffer[1..].trim().to_string();

        // Read sequence line
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            return Err(FastxError::InvalidFormat(
                "Unexpected end of file while reading sequence".to_string(),
            ));
        }
        let sequence = self.line_buffer.trim().to_string();

        // Read separator line (should start with '+')
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            return Err(FastxError::InvalidFormat(
                "Unexpected end of file while reading separator".to_string(),
            ));
        }
        if !self.line_buffer.starts_with('+') {
            return Err(FastxError::InvalidFormat(
                "FASTQ separator must start with '+'".to_string(),
            ));
        }

        // Read quality line
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            return Err(FastxError::InvalidFormat(
                "Unexpected end of file while reading quality".to_string(),
            ));
        }
        let quality = self.line_buffer.trim().to_string();

        // Prepare for next record
        self.line_buffer.clear();
        if self.reader.read_line(&mut self.line_buffer)? == 0 {
            self.finished = true;
        }

        Ok(Some(FastxRecord::new(name, sequence, Some(quality))))
    }
}

impl<R: BufRead> Iterator for FastxReader<R> {
    type Item = Result<FastxRecord, FastxError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let result = if self.is_fastq {
            self.read_fastq_record()
        } else {
            self.read_fasta_record()
        };

        match result {
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
        assert_eq!(records[0].quality, None);
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
        assert_eq!(records[0].quality, Some("IIII".to_string()));
        assert_eq!(records[1].name, "seq2");
        assert_eq!(records[1].sequence, "GGCC");
        assert_eq!(records[1].quality, Some("!!!!".to_string()));
    }
}
