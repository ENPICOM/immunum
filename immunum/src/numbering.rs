use crate::types::{Chain, Scheme};

/// Mock function for numbering sequences
/// This is a placeholder that returns mock numbering results
pub fn number_sequence(sequence: &str, scheme: &Scheme, chains: &[Chain]) -> String {
    let mut result = String::new();

    result.push_str(&format!("Sequence: {}\n", sequence));
    result.push_str(&format!("Scheme: {:?}\n", scheme));
    result.push_str(&format!("Chains: {:?}\n", chains));
    result.push('\n');

    // Mock numbering output
    result.push_str("Mock Numbering Results:\n");
    result.push_str("Position | Residue | IMGT Position\n");
    result.push_str("---------|---------|---------------\n");

    for (i, residue) in sequence.chars().enumerate() {
        let mock_position = match scheme {
            Scheme::IMGT => format!("{}", i + 1),
            Scheme::KABAT => format!("{}K", i + 1),
        };
        result.push_str(&format!(
            "{:8} | {:7} | {}\n",
            i + 1,
            residue,
            mock_position
        ));
    }

    result.push('\n');
    result.push_str("Chain-specific annotations:\n");
    for chain in chains {
        result.push_str(&format!(
            "- {:?}: Mock annotation for this chain type\n",
            chain
        ));
    }

    result
}
