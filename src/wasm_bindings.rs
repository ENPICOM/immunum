#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

use crate::annotator::Annotator;
use crate::constants::{get_scoring_params, ScoringParams};
use crate::result::AnnotationResult;
use crate::types::{Chain, Scheme};

#[cfg(feature = "wasm")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

/// WASM wrapper for ScoringParams
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = "ScoringParams")]
pub struct WasmScoringParams {
    inner: ScoringParams,
}

#[cfg(feature = "wasm")]
#[wasm_bindgen]
impl WasmScoringParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        gap_pen_cp: Option<f64>,
        gap_pen_fr: Option<f64>,
        gap_pen_ip: Option<f64>,
        gap_pen_op: Option<f64>,
        gap_pen_cdr: Option<f64>,
        gap_pen_other: Option<f64>,
        cdr_increase: Option<f64>,
        pen_leap_insertion_point_imgt: Option<f64>,
        pen_leap_insertion_point_kabat: Option<f64>,
    ) -> WasmScoringParams {
        let default_params = ScoringParams::default();
        WasmScoringParams {
            inner: ScoringParams {
                gap_pen_cp: gap_pen_cp.unwrap_or(default_params.gap_pen_cp),
                gap_pen_fr: gap_pen_fr.unwrap_or(default_params.gap_pen_fr),
                gap_pen_ip: gap_pen_ip.unwrap_or(default_params.gap_pen_ip),
                gap_pen_op: gap_pen_op.unwrap_or(default_params.gap_pen_op),
                gap_pen_cdr: gap_pen_cdr.unwrap_or(default_params.gap_pen_cdr),
                gap_pen_other: gap_pen_other.unwrap_or(default_params.gap_pen_other),
                cdr_increase: cdr_increase.unwrap_or(default_params.cdr_increase),
                pen_leap_insertion_point_imgt: pen_leap_insertion_point_imgt
                    .unwrap_or(default_params.pen_leap_insertion_point_imgt),
                pen_leap_insertion_point_kabat: pen_leap_insertion_point_kabat
                    .unwrap_or(default_params.pen_leap_insertion_point_kabat),
            },
        }
    }

    #[wasm_bindgen(getter)]
    pub fn gap_pen_cp(&self) -> f64 {
        self.inner.gap_pen_cp
    }

    #[wasm_bindgen(setter)]
    pub fn set_gap_pen_cp(&mut self, value: f64) {
        self.inner.gap_pen_cp = value;
    }

    #[wasm_bindgen(getter)]
    pub fn gap_pen_op(&self) -> f64 {
        self.inner.gap_pen_op
    }

    #[wasm_bindgen(setter)]
    pub fn set_gap_pen_op(&mut self, value: f64) {
        self.inner.gap_pen_op = value;
    }
}

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
#[wasm_bindgen(js_name = "Annotator")]
pub struct WasmAnnotator {
    inner: Annotator,
}

#[cfg(feature = "wasm")]
#[wasm_bindgen]
impl WasmAnnotator {
    #[wasm_bindgen(constructor)]
    pub fn new(
        scheme: Scheme,
        chains: Vec<Chain>,
        scoring_params: Option<WasmScoringParams>,
        use_prefiltering: Option<bool>,
    ) -> Result<WasmAnnotator, JsValue> {
        let scoring_params_inner = scoring_params.map(|p| p.inner);

        match Annotator::new(scheme, chains, scoring_params_inner, use_prefiltering) {
            Ok(annotator) => Ok(WasmAnnotator { inner: annotator }),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    #[wasm_bindgen(js_name = numberSequence)]
    pub fn number_sequence(&self, sequence: &str) -> Result<WasmAnnotationResult, JsValue> {
        match self.inner.number_sequence(sequence) {
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
pub fn default_scoring_params() -> WasmScoringParams {
    WasmScoringParams {
        inner: get_scoring_params(),
    }
}
