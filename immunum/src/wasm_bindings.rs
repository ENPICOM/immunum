#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

use crate::numbering;
use crate::types::Chain;
use crate::types::Scheme;

#[cfg(feature = "wasm")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

/// Number a sequence using the specified scheme and chains
#[cfg(feature = "wasm")]
#[wasm_bindgen(js_name = numberSequence)]
pub fn number_sequence(
    sequence: &str,
    scheme: &str,
    chains: Vec<String>,
) -> Result<String, JsValue> {
    let valid_scheme = Scheme::parse_from_string(scheme).map_err(|e| JsValue::from_str(&e))?;
    let valid_chains = Chain::parse_vec_from_strings(chains).map_err(|e| JsValue::from_str(&e))?;

    Ok(numbering::number_sequence(
        sequence,
        &valid_scheme,
        &valid_chains,
    ))
}
