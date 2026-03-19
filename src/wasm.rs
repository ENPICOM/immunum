use js_sys::{Object, Reflect};
use std::str::FromStr;
use wasm_bindgen::prelude::*;

use crate::annotator::Annotator;
use crate::numbering::segment;
use crate::types::{Chain, Scheme};

#[wasm_bindgen]
impl Annotator {
    #[wasm_bindgen(constructor, js_name = "new")]
    pub fn wasm_new(
        chains: Vec<String>,
        scheme: String,
        min_confidence: Option<f32>,
    ) -> Result<Annotator, JsValue> {
        let parsed_chains = chains
            .iter()
            .map(|chain| {
                Chain::from_str(chain)
                    .map_err(|_| JsValue::from_str(&format!("Invalid chain: {}", chain)))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let parsed_scheme = Scheme::from_str(&scheme)
            .map_err(|_| JsValue::from_str(&format!("Invalid scheme: {}", scheme)))?;
        Annotator::new(&parsed_chains, parsed_scheme, min_confidence)
            .map_err(|e| JsValue::from_str(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "number")]
    pub fn wasm_number(&self, sequence: &str) -> Result<JsValue, JsValue> {
        let result = self
            .number(sequence)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;

        let aligned_seq = &sequence[result.query_start..=result.query_end];
        let numbering = Object::new();
        for (pos, ch) in result.positions.iter().zip(aligned_seq.chars()) {
            Reflect::set(
                &numbering,
                &JsValue::from_str(&pos.to_string()),
                &JsValue::from_str(&ch.to_string()),
            )
            .unwrap();
        }

        let dict = Object::new();
        Reflect::set(&dict, &"chain".into(), &result.chain.to_string().into()).unwrap();
        Reflect::set(&dict, &"scheme".into(), &result.scheme.to_string().into()).unwrap();
        Reflect::set(&dict, &"confidence".into(), &result.confidence.into()).unwrap();
        Reflect::set(&dict, &"numbering".into(), &numbering.into()).unwrap();
        Ok(dict.into())
    }

    #[wasm_bindgen(js_name = "segment")]
    pub fn wasm_segment(&self, sequence: &str) -> Result<JsValue, JsValue> {
        let result = self
            .number(sequence)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;

        let aligned_seq = &sequence[result.query_start..=result.query_end];
        let segments = segment(&result.positions, aligned_seq, result.scheme);

        let dict = Object::new();
        for (region, seq) in &segments {
            Reflect::set(&dict, &JsValue::from_str(region), &JsValue::from_str(seq)).unwrap();
        }
        Ok(dict.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wasm_bindgen_test::wasm_bindgen_test;

    const IGH_SEQ: &str =
        "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

    #[wasm_bindgen_test]
    fn test_number_igh() {
        let ann = Annotator::wasm_new(
            vec!["IGH".to_string(), "IGK".to_string(), "IGL".to_string()],
            "IMGT".to_string(),
            None,
        )
        .unwrap();

        let result = ann.wasm_number(IGH_SEQ).unwrap();
        let chain = Reflect::get(&result, &"chain".into()).unwrap();
        assert_eq!(chain.as_string().unwrap(), "H");
        let confidence = Reflect::get(&result, &"confidence".into()).unwrap();
        assert!(confidence.as_f64().unwrap() > 0.5);
    }

    #[wasm_bindgen_test]
    fn test_segment_igh() {
        let ann = Annotator::wasm_new(
            vec!["IGH".to_string(), "IGK".to_string(), "IGL".to_string()],
            "IMGT".to_string(),
            None,
        )
        .unwrap();

        let result = ann.wasm_segment(IGH_SEQ).unwrap();
        let fr1 = Reflect::get(&result, &"fr1".into()).unwrap();
        assert!(!fr1.as_string().unwrap().is_empty());
    }

    #[wasm_bindgen_test]
    fn test_invalid_chain_errors() {
        let err = Annotator::wasm_new(vec!["INVALID".to_string()], "IMGT".to_string(), None);
        assert!(err.is_err());
    }

    #[wasm_bindgen_test]
    fn test_invalid_scheme_errors() {
        let err = Annotator::wasm_new(vec!["IGH".to_string()], "NOTASCHEME".to_string(), None);
        assert!(err.is_err());
    }
}
