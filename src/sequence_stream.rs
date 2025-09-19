use crate::fastx::{self, FastxError, FastxRecord};
use std::path::Path;

/// Validates if a string looks like a biological sequence
fn is_valid_sequence(input: &str) -> bool {
    // Check if input string contains valid protein chars
    input.chars().all(|c| "ACDEFGHIKLMNPQRSTVWY".contains(c))
}

/// A wrapper around sequence iterators that can be created from files or direct sequence input
pub struct SequenceStream {
    inner: Box<dyn Iterator<Item = Result<FastxRecord, FastxError>>>,
}

impl SequenceStream {
    /// Creates a new SequenceStream by auto-detecting if the input is a file path or direct sequence
    pub fn new(input: &str) -> Result<Self, FastxError> {
        let input_path = Path::new(input);

        if input_path.exists() {
            Self::from_file(input)
        } else if is_valid_sequence(input) {
            Ok(Self::from_sequence(input))
        } else {
            Err(FastxError::InvalidFormat(
                "Input is neither a valid file path nor a valid sequence".to_string(),
            ))
        }
    }

    /// Creates a SequenceStream from a file path
    fn from_file(path: &str) -> Result<Self, FastxError> {
        let reader = fastx::from_path(path)?;
        Ok(Self {
            inner: Box::new(reader),
        })
    }

    /// Creates a SequenceStream from a direct sequence string
    fn from_sequence(sequence: &str) -> Self {
        let record = FastxRecord::new(
            "INPUT SEQUENCE".to_string(),
            sequence.to_string(),
            None, // No quality scores for direct sequence input
        );
        Self {
            inner: Box::new(std::iter::once(Ok(record))),
        }
    }
}

impl Iterator for SequenceStream {
    type Item = Result<FastxRecord, FastxError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_sequence() {
        let sequence = "ACDEFGHIKLMNPQRSTVWY";
        let stream = SequenceStream::from_sequence(sequence);

        let records: Result<Vec<_>, _> = stream.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "INPUT SEQUENCE");
        assert_eq!(records[0].sequence, sequence);
        assert_eq!(records[0]._quality, None);
    }

    #[test]
    fn test_new_with_sequence() {
        let sequence = "ACDEFGHIKLMNPQRSTVWY";
        let stream = SequenceStream::new(sequence).unwrap();

        let records: Result<Vec<_>, _> = stream.collect();
        let records = records.unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, sequence);
    }

    #[test]
    fn test_new_with_invalid_input() {
        let invalid_input = "not_a_file_and_not_a_sequence_123";
        let result = SequenceStream::new(invalid_input);

        assert!(result.is_err());
        match result {
            Err(FastxError::InvalidFormat(msg)) => {
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
