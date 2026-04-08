use js_sys::{Object, Reflect};
use std::str::FromStr;
use wasm_bindgen::prelude::*;

use crate::annotator::Annotator;
use crate::numbering::segment;
use crate::types::{Chain, Scheme};

#[wasm_bindgen(typescript_custom_section)]
const TS_TYPES: &str = r#"
/** Numbered residues keyed by IMGT/Kabat position string (e.g. `"112A"`). */
export type Numbering = Record<string, string>;

/** Result returned by {@link Annotator.number}. On failure, chain/scheme/confidence/numbering are null and error contains the reason. */
export interface NumberingResult {
    /** Detected chain type: `"H"`, `"K"`, `"L"`, `"A"`, `"B"`, `"G"`, or `"D"`. Null on failure. */
    chain: string | null;
    /** Numbering scheme used: `"imgt"` or `"kabat"`. Null on failure. */
    scheme: string | null;
    /** Alignment confidence score between 0 and 1. Null on failure. */
    confidence: number | null;
    /** Position-to-residue mapping for the aligned region. Null on failure. */
    numbering: Numbering | null;
    /** Error message if numbering failed, null on success. */
    error: string | null;
}

/** FR/CDR segments returned by {@link Annotator.segment}. On failure, all region fields are absent and error contains the reason. */
export interface SegmentationResult {
    fr1?: string;
    cdr1?: string;
    fr2?: string;
    cdr2?: string;
    fr3?: string;
    cdr3?: string;
    fr4?: string;
    /** Residues before FR1 (non-canonical N-terminal extension). */
    prefix?: string;
    /** Residues after FR4 (non-canonical C-terminal extension). */
    postfix?: string;
    /** Error message if segmentation failed, null on success. */
    error: string | null;
}

/**
 * Annotates antibody and T-cell receptor sequences with IMGT or Kabat position numbers.
 *
 * @param chains - Chain types to consider during auto-detection. Each entry is a
 *   case-insensitive string. Accepted values:
 *   - Antibody heavy chain: `"IGH"` / `"H"` / `"heavy"`
 *   - Antibody kappa chain: `"IGK"` / `"K"` / `"kappa"`
 *   - Antibody lambda chain: `"IGL"` / `"L"` / `"lambda"`
 *   - TCR alpha chain:       `"TRA"` / `"A"` / `"alpha"`
 *   - TCR beta chain:        `"TRB"` / `"B"` / `"beta"`
 *   - TCR gamma chain:       `"TRG"` / `"G"` / `"gamma"`
 *   - TCR delta chain:       `"TRD"` / `"D"` / `"delta"`
 *
 *   Pass all chains you want to consider; the annotator scores each and picks the
 *   best-matching one. To consider every supported chain pass all seven values.
 *
 * @param scheme - Numbering scheme to use for output positions. Accepted values
 *   (case-insensitive):
 *   - `"IMGT"` / `"i"` — IMGT numbering (recommended; used internally)
 *   - `"Kabat"` / `"k"` — Kabat numbering (derived from IMGT)
 *
 * @param min_confidence - Optional minimum alignment confidence threshold in the
 *   range `[0, 1]`. Sequences scoring below this value are rejected with an error.
 *   Defaults to `0.5` when `null` or omitted.
 */
export class Annotator {
    free(): void;
    [Symbol.dispose](): void;
    constructor(chains: string[], scheme: string, min_confidence?: number | null);
    number(sequence: string): NumberingResult;
    segment(sequence: string): SegmentationResult;

}
"#;

#[wasm_bindgen]
impl Annotator {
    #[wasm_bindgen(constructor, js_name = "new", skip_typescript)]
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

    #[wasm_bindgen(js_name = "number", skip_typescript)]
    pub fn wasm_number(&self, sequence: &str) -> JsValue {
        let dict = Object::new();
        match self.number(sequence) {
            Ok(result) => {
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
                Reflect::set(&dict, &"chain".into(), &result.chain.to_string().into()).unwrap();
                Reflect::set(&dict, &"scheme".into(), &result.scheme.to_string().into()).unwrap();
                Reflect::set(&dict, &"confidence".into(), &result.confidence.into()).unwrap();
                Reflect::set(&dict, &"numbering".into(), &numbering.into()).unwrap();
                Reflect::set(&dict, &"error".into(), &JsValue::NULL).unwrap();
            }
            Err(e) => {
                Reflect::set(&dict, &"chain".into(), &JsValue::NULL).unwrap();
                Reflect::set(&dict, &"scheme".into(), &JsValue::NULL).unwrap();
                Reflect::set(&dict, &"confidence".into(), &JsValue::NULL).unwrap();
                Reflect::set(&dict, &"numbering".into(), &JsValue::NULL).unwrap();
                Reflect::set(&dict, &"error".into(), &JsValue::from_str(&e.to_string())).unwrap();
            }
        }
        dict.into()
    }

    #[wasm_bindgen(js_name = "segment", skip_typescript)]
    pub fn wasm_segment(&self, sequence: &str) -> JsValue {
        let dict = Object::new();
        match self.number(sequence) {
            Ok(result) => {
                let aligned_seq = &sequence[result.query_start..=result.query_end];
                let segments = segment(&result.positions, aligned_seq, result.scheme);
                for (region, seq) in &segments {
                    Reflect::set(&dict, &JsValue::from_str(region), &JsValue::from_str(seq))
                        .unwrap();
                }
                Reflect::set(&dict, &"error".into(), &JsValue::NULL).unwrap();
            }
            Err(e) => {
                Reflect::set(&dict, &"error".into(), &JsValue::from_str(&e.to_string())).unwrap();
            }
        }
        dict.into()
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

        let result = ann.wasm_number(IGH_SEQ);
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

        let result = ann.wasm_segment(IGH_SEQ);
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
