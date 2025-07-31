#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

use crate::constants::ScoringParams;
use crate::numbering;
use crate::numbering_scheme_type::NumberingScheme;
use crate::schemes::get_scheme;
use crate::types::{Chain, Scheme};

/// WASM enum for numbering schemes
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = "Scheme")]
pub enum WasmScheme {
    /// IMGT numbering scheme
    IMGT,
    /// Kabat numbering scheme
    KABAT,
}

impl From<WasmScheme> for Scheme {
    fn from(wasm_scheme: WasmScheme) -> Self {
        match wasm_scheme {
            WasmScheme::IMGT => Scheme::IMGT,
            WasmScheme::KABAT => Scheme::KABAT,
        }
    }
}

/// WASM enum for chain types
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = "Chain")]
pub enum WasmChain {
    /// Immunoglobulin Heavy chain
    IGH,
    /// Immunoglobulin Kappa chain
    IGK,
    /// Immunoglobulin Lambda chain
    IGL,
    /// T-cell receptor Alpha chain
    TRA,
    /// T-cell receptor Beta chain
    TRB,
    /// T-cell receptor Gamma chain
    TRG,
    /// T-cell receptor Delta chain
    TRD,
}

impl From<WasmChain> for Chain {
    fn from(wasm_chain: WasmChain) -> Self {
        match wasm_chain {
            WasmChain::IGH => Chain::IGH,
            WasmChain::IGK => Chain::IGK,
            WasmChain::IGL => Chain::IGL,
            WasmChain::TRA => Chain::TRA,
            WasmChain::TRB => Chain::TRB,
            WasmChain::TRG => Chain::TRG,
            WasmChain::TRD => Chain::TRD,
        }
    }
}

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
    pub fn new() -> WasmScoringParams {
        WasmScoringParams {
            inner: ScoringParams::default(),
        }
    }

    #[wasm_bindgen(js_name = withCustomParams)]
    pub fn with_custom_params(
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
                pen_leap_insertion_point_imgt: pen_leap_insertion_point_imgt.unwrap_or(default_params.pen_leap_insertion_point_imgt),
                pen_leap_insertion_point_kabat: pen_leap_insertion_point_kabat.unwrap_or(default_params.pen_leap_insertion_point_kabat),
            },
        }
    }
}

/// WASM wrapper for NumberingScheme
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = "NumberingScheme")]
pub struct WasmNumberingScheme {
    inner: NumberingScheme,
}

#[cfg(feature = "wasm")]
#[wasm_bindgen]
impl WasmNumberingScheme {
    #[wasm_bindgen(js_name = getScheme)]
    pub fn get_scheme_wasm(scheme: WasmScheme, chain: WasmChain, params: Option<WasmScoringParams>) -> Result<WasmNumberingScheme, JsValue> {
        let rust_scheme: Scheme = scheme.into();
        let rust_chain: Chain = chain.into();
        let scoring_params = params.map(|p| p.inner);
        
        let numbering_scheme = get_scheme(rust_scheme, rust_chain, scoring_params);
        Ok(WasmNumberingScheme { inner: numbering_scheme })
    }

    #[wasm_bindgen(js_name = numberSequence)]
    pub fn number_sequence(&self, sequence: &str) -> Result<String, JsValue> {
        let output = self.inner.number_sequence(sequence.as_bytes());
        Ok(format!("Identity: {:.2}%, Numbering: {:?}", output.identity * 100.0, output.numbering))
    }
}

/// Batch processing function for multiple sequences
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = numberSequencesBatch)]
pub fn number_sequences_batch(
    sequences: Vec<String>,
    scheme: WasmScheme,
    chains: Vec<WasmChain>,
) -> Result<Vec<String>, JsValue> {
    let rust_scheme: Scheme = scheme.into();
    let rust_chains: Vec<Chain> = chains.into_iter().map(|c| c.into()).collect();

    let results: Vec<String> = sequences
        .iter()
        .map(|seq| numbering::number_sequence(seq, &rust_scheme, &rust_chains))
        .collect();

    Ok(results)
}
