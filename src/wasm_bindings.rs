#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

#[cfg(feature = "wasm")]
use crate::annotator::Annotator as RustAnnotator;
#[cfg(feature = "wasm")]
use crate::sequence::SequenceRecord;
#[cfg(feature = "wasm")]
use crate::types::{Chain, Scheme};

#[cfg(feature = "wasm")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

/// WASM wrapper for Annotator - the main entry point
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = "Annotator")]
pub struct Annotator {
    inner: RustAnnotator,
}

#[cfg(feature = "wasm")]
#[wasm_bindgen]
impl Annotator {
    #[wasm_bindgen(constructor)]
    pub fn new(
        scheme: Scheme,
        chains: Vec<Chain>,
        disable_prefiltering: Option<bool>,
        min_confidence: Option<f64>,
    ) -> Result<Annotator, JsValue> {
        use crate::constants::MINIMAL_CHAIN_IDENTITY;

        let disable_prefiltering = disable_prefiltering.unwrap_or(false);
        let min_confidence = Some(min_confidence.unwrap_or(MINIMAL_CHAIN_IDENTITY));
        match RustAnnotator::new(scheme, chains, disable_prefiltering, min_confidence) {
            Ok(annotator) => Ok(Annotator { inner: annotator }),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Number a single sequence - returns JSON string result
    #[wasm_bindgen(js_name = "numberSequence")]
    pub fn number_sequence(&self, sequence: &str, max_chains: Option<usize>) -> String {
        let max_chains = max_chains.unwrap_or(2);

        let sequence_record = SequenceRecord {
            name: b"sequence".to_vec(),
            sequence: sequence.as_bytes().to_vec(),
        };

        match self
            .inner
            .number_sequence(&sequence_record, Some(max_chains))
        {
            Ok(chains) => {
                let chain_results: Vec<String> = chains
                    .into_iter()
                    .map(|chain| {
                        format!(
                            r#"{{"numbers":[{}],"identity":{},"chain":"{}","scheme":"{}","start":{},"end":{}}}"#,
                            chain.numbers.iter().map(|n| format!(r#""{}""#, n)).collect::<Vec<_>>().join(","),
                            chain.identity,
                            chain.chain.to_short(),
                            match chain.scheme {
                                Scheme::IMGT => "IMGT",
                                Scheme::KABAT => "KABAT",
                            },
                            chain.start,
                            chain.end
                        )
                    })
                    .collect();
                format!("[{}]", chain_results.join(","))
            }
            Err(e) => format!(r#"{{"error": "{}"}}"#, e),
        }
    }

    /// Number multiple sequences - takes JSON string input and returns JSON string output
    #[wasm_bindgen(js_name = "numberSequences")]
    pub fn number_sequences(&self, sequences_json: &str, max_chains: Option<usize>) -> String {
        let max_chains = max_chains.unwrap_or(2);

        // Parse input JSON - expect array of strings or array of [name, sequence] pairs
        let sequences: Result<Vec<serde_json::Value>, _> = serde_json::from_str(sequences_json);

        let sequences = match sequences {
            Ok(seqs) => seqs,
            Err(e) => return format!(r#"{{"error": "Invalid JSON input: {}"}}"#, e),
        };

        let results: Vec<_> = sequences
            .into_iter()
            .enumerate()
            .map(|(i, seq_value)| {
                // Handle both string and [name, sequence] array formats
                let (name, sequence) = if let Some(seq_str) = seq_value.as_str() {
                    (format!("sequence_{}", i), seq_str.to_string())
                } else if let Some(arr) = seq_value.as_array() {
                    if arr.len() == 2 {
                        if let (Some(name), Some(seq)) = (arr[0].as_str(), arr[1].as_str()) {
                            (name.to_string(), seq.to_string())
                        } else {
                            return r#"{"error": "Array elements must be strings"}"#.to_string();
                        }
                    } else {
                        return r#"{"error": "Array must have exactly 2 elements [name, sequence]"}"#.to_string();
                    }
                } else {
                    return r#"{"error": "Each sequence must be a string or [name, sequence] array"}"#.to_string();
                };

                let sequence_record = SequenceRecord {
                    name: name.into_bytes(),
                    sequence: sequence.into_bytes(),
                };

                match self
                    .inner
                    .number_sequence(&sequence_record, Some(max_chains))
                {
                    Ok(chains) => {
                        let chain_results: Vec<String> = chains
                            .into_iter()
                            .map(|chain| {
                                format!(
                                    r#"{{"numbers":[{}],"identity":{},"chain":"{}","scheme":"{}","start":{},"end":{}}}"#,
                                    chain.numbers.iter().map(|n| format!(r#""{}""#, n)).collect::<Vec<_>>().join(","),
                                    chain.identity,
                                    chain.chain.to_short(),
                                    match chain.scheme {
                                        Scheme::IMGT => "IMGT",
                                        Scheme::KABAT => "KABAT",
                                    },
                                    chain.start,
                                    chain.end
                                )
                            })
                            .collect();
                        format!("[{}]", chain_results.join(","))
                    }
                    Err(e) => format!(r#"{{"error": "{}"}}"#, e),
                }
            })
            .collect();

        format!("[{}]", results.join(","))
    }
}

/// Create a new Annotator with default parameters
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = "createAnnotator")]
pub fn create_annotator(
    scheme: Option<Scheme>,
    chains: Option<Vec<Chain>>,
    disable_prefiltering: Option<bool>,
    min_confidence: Option<f64>,
) -> Result<Annotator, JsValue> {
    let scheme = scheme.unwrap_or(Scheme::IMGT);
    let chains = chains.unwrap_or_else(|| vec![Chain::IGH, Chain::IGK, Chain::IGL]);

    Annotator::new(scheme, chains, disable_prefiltering, min_confidence)
}
