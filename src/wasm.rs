use js_sys::{Object, Reflect};
use wasm_bindgen::prelude::*;

use crate::annotator::{Annotator, NumberingResult, SegmentResult};
use crate::types::{Chain, Scheme};

// What `Annotator.number` returns for `sequence`: the numbering, or the error with every other field
// null
fn numbering_object(sequence: &str, result: crate::Result<NumberingResult>) -> JsValue {
    let dict = Object::new();
    match result {
        Ok(result) => {
            let numbering = js_sys::Map::new();
            for (pos, ch) in result
                .residues(sequence)
                .expect("`result` numbered `sequence`")
            {
                numbering.set(
                    &JsValue::from_str(&pos.to_string()),
                    &JsValue::from_str(&ch.to_string()),
                );
            }
            Reflect::set(&dict, &"chain".into(), &result.chain.to_string().into()).unwrap();
            Reflect::set(&dict, &"scheme".into(), &result.scheme.to_string().into()).unwrap();
            Reflect::set(&dict, &"confidence".into(), &result.confidence.into()).unwrap();
            Reflect::set(&dict, &"numbering".into(), &numbering.into()).unwrap();
            Reflect::set(
                &dict,
                &"queryStart".into(),
                &(result.query_start as u32).into(),
            )
            .unwrap();
            Reflect::set(&dict, &"queryEnd".into(), &(result.query_end as u32).into()).unwrap();
            Reflect::set(&dict, &"error".into(), &JsValue::NULL).unwrap();
        }
        Err(e) => {
            Reflect::set(&dict, &"chain".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"scheme".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"confidence".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"numbering".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"queryStart".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"queryEnd".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"error".into(), &JsValue::from_str(&e.to_string())).unwrap();
        }
    }
    dict.into()
}

// What `Annotator.segment` returns: the segments, or the error with no segments
fn segment_object(result: crate::Result<SegmentResult>) -> JsValue {
    let dict = Object::new();
    match result {
        Ok(s) => {
            for (name, residues) in s.regions() {
                Reflect::set(
                    &dict,
                    &JsValue::from_str(name),
                    &JsValue::from_str(residues),
                )
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

// One object per domain, or a single error object when the sequence couldn't be searched
fn per_domain<T>(
    results: crate::Result<Vec<T>>,
    object: impl Fn(crate::Result<T>) -> JsValue,
) -> JsValue {
    let objects = js_sys::Array::new();
    match results {
        Ok(results) => {
            for result in results {
                objects.push(&object(Ok(result)));
            }
        }
        Err(e) => {
            objects.push(&object(Err(e)));
        }
    }
    objects.into()
}

#[wasm_bindgen(typescript_custom_section)]
const TS_TYPES: &str = r#"
/**
 * Ordered position → residue map, iterated in IMGT-correct order.
 * Keys are position strings (e.g. `"112A"`) and values are single-character residues.
 */
export type Numbering = Map<string, string>;

/** Result returned by {@link Annotator.number}. On failure, chain/scheme/confidence/numbering/queryStart/queryEnd are null and error contains the reason. */
export interface NumberingResult {
    /** Detected chain type: `"H"`, `"K"`, `"L"`, `"A"`, `"B"`, `"G"`, or `"D"`. Null on failure. */
    chain: string | null;
    /** Numbering scheme used: `"IMGT"`, `"Kabat"`, `"Chothia"`, `"Martin"` or `"Aho"`. Null on failure. */
    scheme: string | null;
    /** Alignment confidence score between 0 and 1. Null on failure. */
    confidence: number | null;
    /** Position-to-residue mapping for the aligned region. Null on failure. */
    numbering: Numbering | null;
    /** 0-indexed start of the aligned region in the input sequence (inclusive). Null on failure. */
    queryStart: number | null;
    /** 0-indexed end of the aligned region in the input sequence (inclusive). Null on failure. */
    queryEnd: number | null;
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
 * Annotates antibody and T-cell receptor sequences with scheme-specific position numbers.
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
 *   A group of chains is accepted too: `"ig"` (IGH, IGK, IGL), `"tcr"` (TRA, TRB, TRG,
 *   TRD) or `"all"`. Pass all chains you want to consider; the annotator scores each and
 *   picks the best-matching one.
 *
 * @param scheme - Numbering scheme to use for output positions. Accepted values
 *   (case-insensitive):
 *   - `"IMGT"` / `"i"` — IMGT numbering (recommended; used internally)
 *   - `"Kabat"` / `"k"` — Kabat numbering (derived from IMGT)
 *   - `"Chothia"` / `"c"` — Chothia numbering (derived from IMGT)
 *   - `"Martin"` / `"m"` — Martin / extended Chothia numbering (derived from IMGT)
 *   - `"Aho"` / `"a"` — AHo numbering (derived from IMGT)
 *
 *   Only IMGT supports TCR chains; the other schemes are restricted to antibody
 *   chains (IGH, IGK, IGL). {@link schemeSupportsChain} checks a pair up front.
 *
 * @param minConfidence - Optional minimum alignment confidence threshold in the
 *   range `[0, 1]`. Sequences scoring below this value are rejected with an error.
 *   Defaults to `0.5` when `null` or omitted.
 */
export class Annotator {
    free(): void;
    [Symbol.dispose](): void;
    constructor(chains: string[], scheme: string, minConfidence?: number | null);
    number(sequence: string): NumberingResult;
    /**
     * Number every variable domain in a sequence, such as both domains of an scFv. One
     * {@link NumberingResult} per domain, ordered by position, each what `number` returns for that
     * domain; empty when no domain aligns with enough confidence. When the sequence itself is
     * invalid, a single result with `error` set.
     *
     * A domain that lacks its first IMGT positions (a light chain starting at position 2, say) and
     * directly follows other residues, such as a linker, can have the residue just before it
     * numbered as its first position: IMGT position 1 is so variable that the sequence alone can't
     * tell a linker residue from the domain's own first residue.
     */
    numberDomains(sequence: string): NumberingResult[];
    segment(sequence: string): SegmentationResult;
    /**
     * Split every variable domain in a sequence into FR/CDR regions. One
     * {@link SegmentationResult} per domain, ordered by position; empty when no domain aligns with
     * enough confidence. When the sequence itself is invalid, a single result with `error` set.
     *
     * Every residue lands in exactly one domain's regions: a domain's `prefix` holds the residues
     * since the previous domain (or the start of the sequence), and only the last domain has the
     * residues after it as its `postfix`. All domains' regions in order rebuild the sequence.
     */
    segmentDomains(sequence: string): SegmentationResult[];

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
        Annotator::from_names(chains.iter().map(String::as_str), &scheme, min_confidence)
            .map_err(|e| JsValue::from_str(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "number", skip_typescript)]
    pub fn wasm_number(&self, sequence: &str) -> JsValue {
        numbering_object(sequence, self.number(sequence))
    }

    #[wasm_bindgen(js_name = "numberDomains", skip_typescript)]
    pub fn wasm_number_domains(&self, sequence: &str) -> JsValue {
        per_domain(self.number_domains(sequence), |result| {
            numbering_object(sequence, result)
        })
    }

    #[wasm_bindgen(js_name = "segment", skip_typescript)]
    pub fn wasm_segment(&self, sequence: &str) -> JsValue {
        segment_object(self.segment(sequence))
    }

    #[wasm_bindgen(js_name = "segmentDomains", skip_typescript)]
    pub fn wasm_segment_domains(&self, sequence: &str) -> JsValue {
        per_domain(self.segment_domains(sequence), segment_object)
    }
}

/// Whether `scheme` numbers `chain`: IMGT numbers every chain, the other schemes antibody chains
/// (IGH, IGK, IGL) only. Takes a scheme and a single chain by the names `Annotator` accepts, and
/// throws on an unknown one.
#[wasm_bindgen(js_name = "schemeSupportsChain")]
pub fn scheme_supports_chain(scheme: &str, chain: &str) -> Result<bool, JsValue> {
    let to_js = |e: crate::Error| JsValue::from_str(&e.to_string());
    let scheme: Scheme = scheme.parse().map_err(to_js)?;
    let chain: Chain = chain.parse().map_err(to_js)?;
    Ok(scheme.supports(chain))
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
