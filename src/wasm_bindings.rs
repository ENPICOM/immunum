#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

use crate::annotator::Annotator as RustAnnotator;
use crate::constants::{get_scoring_params, ScoringParams};
use crate::result::AnnotationResult;
use crate::types::{Chain, Scheme};

#[cfg(feature = "wasm")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

/// WASM wrapper for AnnotationResult
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = "AnnotationResult")]
pub struct WasmAnnotationResult {
    inner: AnnotationResult,
}

#[cfg(feature = "wasm")]
#[wasm_bindgen]
impl WasmAnnotationResult {
    #[wasm_bindgen(getter)]
    pub fn sequence(&self) -> String {
        self.inner.sequence_string()
    }

    #[wasm_bindgen(getter)]
    pub fn numbers(&self) -> Vec<String> {
        self.inner.numbers.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn scheme(&self) -> Scheme {
        self.inner.scheme
    }

    #[wasm_bindgen(getter)]
    pub fn chain(&self) -> Chain {
        self.inner.chain
    }

    #[wasm_bindgen(getter)]
    pub fn identity(&self) -> f64 {
        self.inner.identity
    }

    #[wasm_bindgen(getter)]
    pub fn start(&self) -> u32 {
        self.inner.start
    }

    #[wasm_bindgen(getter)]
    pub fn end(&self) -> u32 {
        self.inner.end
    }

    #[wasm_bindgen(js_name = getRegionSequence)]
    pub fn get_region_sequence(&self, region_name: &str) -> Option<String> {
        self.inner.get_region_sequence(region_name)
    }
}

/// WASM wrapper for Annotator - the main entry point
#[cfg(feature = "wasm")]
#[wasm_bindgen]
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
        scoring_params: Option<ScoringParams>,
        use_prefiltering: Option<bool>,
    ) -> Result<Annotator, JsValue> {
        match RustAnnotator::new(scheme, chains, scoring_params, use_prefiltering) {
            Ok(annotator) => Ok(Annotator { inner: annotator }),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    #[wasm_bindgen(js_name = numberSequence)]
    pub fn number_sequence(&self, sequence: &str, sequence_id: String) -> Result<WasmAnnotationResult, JsValue> {
        match self.inner.number_sequence(sequence, sequence_id) {
            Ok(result) => Ok(WasmAnnotationResult { inner: result }),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    #[wasm_bindgen(js_name = numberSequences)]
    pub fn number_sequences(&self, sequences: Vec<String>) -> Vec<JsValue> {
        let results = self.inner.number_sequences(&sequences, false);
        results
            .into_iter()
            .map(|result| match result {
                Ok(annotation_result) => JsValue::from(WasmAnnotationResult {
                    inner: annotation_result,
                }),
                Err(e) => JsValue::from_str(&format!("Error: {}", e)),
            })
            .collect()
    }
}

/// Get default scoring parameters
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = "defaultScoringParams")]
pub fn default_scoring_params() -> ScoringParams {
    get_scoring_params()
}
