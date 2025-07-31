#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

use crate::constants::ScoringParams;
use crate::numbering;
use crate::numbering_scheme_type::NumberingScheme;
use crate::schemes::{
    get_imgt_heavy_scheme, get_imgt_heavy_scheme_with_params,
    get_imgt_kappa_scheme, get_imgt_kappa_scheme_with_params,
    get_imgt_lambda_scheme, get_imgt_lambda_scheme_with_params,
    get_kabat_heavy_scheme, get_kabat_heavy_scheme_with_params,
    get_kabat_kappa_scheme, get_kabat_kappa_scheme_with_params,
    get_kabat_lambda_scheme, get_kabat_lambda_scheme_with_params,
};
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

impl From<Scheme> for WasmScheme {
    fn from(scheme: Scheme) -> Self {
        match scheme {
            Scheme::IMGT => WasmScheme::IMGT,
            Scheme::KABAT => WasmScheme::KABAT,
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

impl From<Chain> for WasmChain {
    fn from(chain: Chain) -> Self {
        match chain {
            Chain::IGH => WasmChain::IGH,
            Chain::IGK => WasmChain::IGK,
            Chain::IGL => WasmChain::IGL,
            Chain::TRA => WasmChain::TRA,
            Chain::TRB => WasmChain::TRB,
            Chain::TRG => WasmChain::TRG,
            Chain::TRD => WasmChain::TRD,
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

    #[wasm_bindgen(getter, js_name = gapPenCp)]
    pub fn gap_pen_cp(&self) -> f64 {
        self.inner.gap_pen_cp
    }

    #[wasm_bindgen(setter, js_name = gapPenCp)]
    pub fn set_gap_pen_cp(&mut self, value: f64) {
        self.inner.gap_pen_cp = value;
    }

    #[wasm_bindgen(getter, js_name = gapPenFr)]
    pub fn gap_pen_fr(&self) -> f64 {
        self.inner.gap_pen_fr
    }

    #[wasm_bindgen(setter, js_name = gapPenFr)]
    pub fn set_gap_pen_fr(&mut self, value: f64) {
        self.inner.gap_pen_fr = value;
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
    #[wasm_bindgen(js_name = imgtHeavy)]
    pub fn imgt_heavy(params: Option<WasmScoringParams>) -> Result<WasmNumberingScheme, JsValue> {
        let scheme = match params {
            Some(wasm_params) => get_imgt_heavy_scheme_with_params(Some(wasm_params.inner)),
            None => get_imgt_heavy_scheme(),
        };
        Ok(WasmNumberingScheme { inner: scheme })
    }

    #[wasm_bindgen(js_name = imgtKappa)]
    pub fn imgt_kappa(params: Option<WasmScoringParams>) -> Result<WasmNumberingScheme, JsValue> {
        let scheme = match params {
            Some(wasm_params) => get_imgt_kappa_scheme_with_params(Some(wasm_params.inner)),
            None => get_imgt_kappa_scheme(),
        };
        Ok(WasmNumberingScheme { inner: scheme })
    }

    #[wasm_bindgen(js_name = imgtLambda)]
    pub fn imgt_lambda(params: Option<WasmScoringParams>) -> Result<WasmNumberingScheme, JsValue> {
        let scheme = match params {
            Some(wasm_params) => get_imgt_lambda_scheme_with_params(Some(wasm_params.inner)),
            None => get_imgt_lambda_scheme(),
        };
        Ok(WasmNumberingScheme { inner: scheme })
    }

    #[wasm_bindgen(js_name = kabatHeavy)]
    pub fn kabat_heavy(params: Option<WasmScoringParams>) -> Result<WasmNumberingScheme, JsValue> {
        let scheme = match params {
            Some(wasm_params) => get_kabat_heavy_scheme_with_params(Some(wasm_params.inner)),
            None => get_kabat_heavy_scheme(),
        };
        Ok(WasmNumberingScheme { inner: scheme })
    }

    #[wasm_bindgen(js_name = kabatKappa)]
    pub fn kabat_kappa(params: Option<WasmScoringParams>) -> Result<WasmNumberingScheme, JsValue> {
        let scheme = match params {
            Some(wasm_params) => get_kabat_kappa_scheme_with_params(Some(wasm_params.inner)),
            None => get_kabat_kappa_scheme(),
        };
        Ok(WasmNumberingScheme { inner: scheme })
    }

    #[wasm_bindgen(js_name = kabatLambda)]
    pub fn kabat_lambda(params: Option<WasmScoringParams>) -> Result<WasmNumberingScheme, JsValue> {
        let scheme = match params {
            Some(wasm_params) => get_kabat_lambda_scheme_with_params(Some(wasm_params.inner)),
            None => get_kabat_lambda_scheme(),
        };
        Ok(WasmNumberingScheme { inner: scheme })
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
