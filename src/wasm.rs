use js_sys::{Object, Reflect};
use wasm_bindgen::prelude::*;

use crate::annotator::{per_domain, Annotator, NumberingResult, SegmentResult};
use crate::numbering::{region_spans, SEGMENT_NAMES};
use crate::types::scheme_supports_chain;
use crate::SequenceError;

// What every function throws: an `Error` with the message and the error's `kind`
fn to_js(e: crate::Error) -> JsValue {
    let error = js_sys::Error::new(&e.to_string());
    Reflect::set(&error, &"kind".into(), &e.kind().into()).unwrap();
    error.into()
}

// Sets `error` and `errorKind` on a returned object, both null on success
fn set_error(dict: &Object, error: Option<&SequenceError>) {
    let (message, kind) = match error {
        Some(e) => (
            JsValue::from_str(&e.to_string()),
            JsValue::from_str(e.kind()),
        ),
        None => (JsValue::NULL, JsValue::NULL),
    };
    Reflect::set(dict, &"error".into(), &message).unwrap();
    Reflect::set(dict, &"errorKind".into(), &kind).unwrap();
}

// What `Annotator.number` returns for `sequence`: the numbering, or the error with every other field
// null. The fields are camelCase on purpose, as JavaScript and TypeScript code expects:
// `queryStart`, `queryEnd` and `errorKind` here are `query_start`, `query_end` and `error_kind` in
// every other interface.
fn numbering_object(sequence: &str, result: Result<NumberingResult, SequenceError>) -> JsValue {
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
            set_error(&dict, None);
        }
        Err(e) => {
            Reflect::set(&dict, &"chain".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"scheme".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"confidence".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"numbering".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"queryStart".into(), &JsValue::NULL).unwrap();
            Reflect::set(&dict, &"queryEnd".into(), &JsValue::NULL).unwrap();
            set_error(&dict, Some(&e));
        }
    }
    dict.into()
}

// What `Annotator.segment` returns: the segments, or the error with every segment null
fn segment_object(result: Result<SegmentResult, SequenceError>) -> JsValue {
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
            set_error(&dict, None);
        }
        Err(e) => {
            for name in SEGMENT_NAMES {
                Reflect::set(&dict, &JsValue::from_str(name), &JsValue::NULL).unwrap();
            }
            set_error(&dict, Some(&e));
        }
    }
    dict.into()
}

// One object per domain, as `per_domain` returns them
fn domain_objects<T>(
    results: Result<Vec<T>, SequenceError>,
    object: impl Fn(Result<T, SequenceError>) -> JsValue,
) -> JsValue {
    per_domain(results)
        .into_iter()
        .map(object)
        .collect::<js_sys::Array>()
        .into()
}

#[wasm_bindgen(typescript_custom_section)]
const TS_TYPES: &str = r#"
/**
 * Ordered position → residue map, iterated in IMGT-correct order.
 * Keys are position strings (e.g. `"112A"`) and values are single-character residues.
 */
export type Numbering = Map<string, string>;

/**
 * Error thrown when immunum is set up or called wrongly, such as with an unknown chain name.
 * `message` says what went wrong; `kind` is a stable code for it: `"invalid_chain"`,
 * `"invalid_scheme"`, `"unsupported_chain"` or `"invalid_min_confidence"`.
 */
export interface ImmunumError extends Error {
    kind: string;
}

/**
 * Result returned by {@link Annotator.number}. On failure, chain/scheme/confidence/numbering/queryStart/queryEnd are null and error and errorKind contain the reason.
 *
 * The fields are camelCase on purpose, as JavaScript and TypeScript code expects: `queryStart`,
 * `queryEnd` and `errorKind` are `query_start`, `query_end` and `error_kind` in Python, Polars and
 * the CLI.
 */
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
    /** What went wrong if numbering failed: `"invalid_sequence"` or `"low_confidence"`. Null on success. */
    errorKind: string | null;
}

/** FR/CDR segments returned by {@link Annotator.segment}. On failure, every region field is null and error and errorKind contain the reason. */
export interface SegmentationResult {
    fr1: string | null;
    cdr1: string | null;
    fr2: string | null;
    cdr2: string | null;
    fr3: string | null;
    cdr3: string | null;
    fr4: string | null;
    /** Residues before FR1 (non-canonical N-terminal extension). */
    prefix: string | null;
    /** Residues after FR4 (non-canonical C-terminal extension). */
    postfix: string | null;
    /** Error message if segmentation failed, null on success. */
    error: string | null;
    /** What went wrong if segmentation failed: `"invalid_sequence"` or `"low_confidence"`. Null on success. */
    errorKind: string | null;
}

/** Inclusive `[start, end]` position bounds of each FR/CDR region, N- to C-terminal. */
export interface RegionSpans {
    fr1: [number, number];
    cdr1: [number, number];
    fr2: [number, number];
    cdr2: [number, number];
    fr3: [number, number];
    cdr3: [number, number];
    fr4: [number, number];
}

/**
 * The FR/CDR region boundaries a scheme uses for a chain, both named as for {@link Annotator}.
 * IMGT and AHo number every chain alike; Kabat, Chothia and Martin place their CDRs differently
 * on heavy and light chains.
 *
 * @throws {ImmunumError} `invalid_scheme` or `invalid_chain` for an unknown scheme or chain, and
 *   `unsupported_chain` for a chain the scheme doesn't number (only IMGT numbers TCR chains).
 */
export function regionsFor(scheme: string, chain: string): RegionSpans;

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
 *
 * @throws {ImmunumError} `invalid_chain` for an unknown chain name or no chains,
 *   `invalid_scheme` for an unknown scheme, `unsupported_chain` for a chain the scheme
 *   doesn't number and `invalid_min_confidence` for a `minConfidence` outside `[0, 1]`.
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
            .map_err(to_js)
    }

    #[wasm_bindgen(js_name = "number", skip_typescript)]
    pub fn wasm_number(&self, sequence: &str) -> JsValue {
        numbering_object(sequence, self.number(sequence))
    }

    #[wasm_bindgen(js_name = "numberDomains", skip_typescript)]
    pub fn wasm_number_domains(&self, sequence: &str) -> JsValue {
        domain_objects(self.number_domains(sequence), |result| {
            numbering_object(sequence, result)
        })
    }

    #[wasm_bindgen(js_name = "segment", skip_typescript)]
    pub fn wasm_segment(&self, sequence: &str) -> JsValue {
        segment_object(self.segment(sequence))
    }

    #[wasm_bindgen(js_name = "segmentDomains", skip_typescript)]
    pub fn wasm_segment_domains(&self, sequence: &str) -> JsValue {
        domain_objects(self.segment_domains(sequence), segment_object)
    }
}

/// Whether `scheme` numbers `chain`: IMGT numbers every chain, the other schemes antibody chains
/// (IGH, IGK, IGL) only. Takes a scheme and a single chain by the names `Annotator` accepts.
///
/// @throws {ImmunumError} `invalid_scheme` or `invalid_chain` for an unknown scheme or chain.
#[wasm_bindgen(js_name = "schemeSupportsChain")]
pub fn wasm_scheme_supports_chain(scheme: &str, chain: &str) -> Result<bool, JsValue> {
    scheme_supports_chain(scheme, chain).map_err(to_js)
}

#[wasm_bindgen(js_name = "regionsFor", skip_typescript)]
pub fn wasm_regions_for(scheme: &str, chain: &str) -> Result<JsValue, JsValue> {
    let spans = Object::new();
    for (region, (start, end)) in region_spans(scheme, chain).map_err(to_js)? {
        let span = js_sys::Array::of2(&start.into(), &end.into());
        Reflect::set(&spans, &JsValue::from_str(region), &span).unwrap();
    }
    Ok(spans.into())
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
