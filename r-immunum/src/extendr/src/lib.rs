use extendr_api::prelude::*;
use rayon::prelude::*;
use std::str::FromStr;

use immunum::numbering::segment as segment_positions;
use immunum::{Annotator as InnerAnnotator, Chain, Scheme};

// extendr 0.7 converts Result::Err via throw_r_error → Rf_error (longjmp).
// On Windows+GNU this longjmp across Rust frames is UB and segfaults.
// Annotator::new must return Result<Self> (extendr requirement); R-side
// validates inputs so the Err path is unreachable in normal operation.
// number()/segment() return a plain List: on error the list has an
// `error` field carrying the Rust error message and all other fields NULL.
// On success, `error` is absent. This mirrors the Python wrapper's
// error-returning (non-throwing) approach added in upstream v1.1.0.

#[extendr]
fn immunum_version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

struct Annotator {
    inner: InnerAnnotator,
}

#[extendr]
impl Annotator {
    // extendr requires Result<Self> for struct constructors.
    // R-side validates chains/scheme/min_confidence before calling,
    // so the Err branch is unreachable in normal operation.
    fn new(chains: Strings, scheme: String, min_confidence: Nullable<f64>) -> Result<Self> {
        let parsed_chains: Vec<Chain> = chains
            .iter()
            .map(|c| {
                Chain::from_str(c.as_str())
                    .map_err(|_| Error::Other(format!("invalid chain {:?}", c.as_str())))
            })
            .collect::<Result<Vec<_>>>()?;

        let parsed_scheme = Scheme::from_str(&scheme)
            .map_err(|_| Error::Other(format!("invalid scheme {:?}", scheme)))?;

        let conf = match min_confidence {
            Nullable::NotNull(v) => Some(v as f32),
            Nullable::Null => None,
        };

        let inner = InnerAnnotator::new(&parsed_chains, parsed_scheme, conf)
            .map_err(|e| Error::Other(format!("{}", e)))?;

        Ok(Self { inner })
    }

    fn number(&self, sequence: &str) -> List {
        let result = match self.inner.number(sequence) {
            Ok(r) => r,
            Err(e) => return list!(error = e.to_string()),
        };

        let aligned = &sequence[result.query_start..=result.query_end];
        let positions: Vec<String> = result.positions.iter().map(|p| p.to_string()).collect();
        let residues: Vec<String> = aligned.chars().map(|c| c.to_string()).collect();

        list!(
            chain = result.chain.to_string(),
            scheme = result.scheme.to_string(),
            confidence = result.confidence as f64,
            positions = positions,
            residues = residues,
            query_start = (result.query_start + 1) as i32,
            query_end = (result.query_end + 1) as i32,
        )
    }

    fn segment(&self, sequence: &str) -> List {
        let result = match self.inner.segment(sequence) {
            Ok(r) => r,
            Err(e) => return list!(error = e.to_string()),
        };

        list!(
            prefix = result.prefix,
            fr1 = result.fr1,
            cdr1 = result.cdr1,
            fr2 = result.fr2,
            cdr2 = result.cdr2,
            fr3 = result.fr3,
            cdr3 = result.cdr3,
            fr4 = result.fr4,
            postfix = result.postfix,
        )
    }
}

// ---------------------------------------------------------------------------
// Batch primitives
// ---------------------------------------------------------------------------

fn parse_chains(chains: &Strings) -> Result<Vec<Chain>> {
    chains
        .iter()
        .map(|c| {
            Chain::from_str(c.as_str())
                .map_err(|_| Error::Other(format!("invalid chain {:?}", c.as_str())))
        })
        .collect()
}

fn parse_scheme(scheme: &str) -> Result<Scheme> {
    Scheme::from_str(scheme)
        .map_err(|_| Error::Other(format!("invalid scheme {:?}", scheme)))
}

fn build_inner_annotator(
    chains: Strings,
    scheme: String,
    min_confidence: Nullable<f64>,
) -> Result<InnerAnnotator> {
    let parsed_chains = parse_chains(&chains)?;
    let parsed_scheme = parse_scheme(&scheme)?;
    let conf = match min_confidence {
        Nullable::NotNull(v) => Some(v as f32),
        Nullable::Null => None,
    };
    InnerAnnotator::new(&parsed_chains, parsed_scheme, conf)
        .map_err(|e| Error::Other(format!("{}", e)))
}

fn materialize_inputs(seqs: Strings) -> Vec<Option<String>> {
    seqs.iter()
        .map(|s| {
            if s.is_na() {
                None
            } else {
                Some(s.as_str().to_string())
            }
        })
        .collect()
}

type NumRow = Option<std::result::Result<(String, String, f64, Vec<String>, Vec<String>), String>>;

fn process_number(ann: &InnerAnnotator, seq: &str) -> NumRow {
    let result = match ann.number(seq) {
        Ok(r) => r,
        Err(e) => return Some(Err(e.to_string())),
    };
    let (positions, residues): (Vec<String>, Vec<String>) = result
        .positions
        .iter()
        .zip(seq.chars())
        .map(|(pos, ch)| (pos.to_string(), ch.to_string()))
        .unzip();
    Some(Ok((
        result.chain.to_string(),
        result.scheme.to_string(),
        result.confidence as f64,
        positions,
        residues,
    )))
}

type SegRow = Option<std::result::Result<[String; 9], String>>;

fn process_segment(ann: &InnerAnnotator, seq: &str) -> SegRow {
    let result = match ann.number(seq) {
        Ok(r) => r,
        Err(e) => return Some(Err(e.to_string())),
    };
    let mut s = segment_positions(&result.positions, seq, result.scheme);
    let mut take = |k: &str| s.remove(k).unwrap_or_default();
    Some(Ok([
        take("prefix"),
        take("fr1"),
        take("cdr1"),
        take("fr2"),
        take("cdr2"),
        take("fr3"),
        take("cdr3"),
        take("fr4"),
        take("postfix"),
    ]))
}

fn collect_number_columns(rows: Vec<NumRow>) -> List {
    let n = rows.len();
    let mut chain: Vec<Option<String>> = Vec::with_capacity(n);
    let mut scheme: Vec<Option<String>> = Vec::with_capacity(n);
    let mut conf: Vec<Option<f64>> = Vec::with_capacity(n);
    let mut positions: Vec<Robj> = Vec::with_capacity(n);
    let mut residues: Vec<Robj> = Vec::with_capacity(n);
    let mut error: Vec<Option<String>> = Vec::with_capacity(n);
    for row in rows {
        match row {
            None => {
                chain.push(None);
                scheme.push(None);
                conf.push(None);
                positions.push(Robj::from(()));
                residues.push(Robj::from(()));
                error.push(None);
            }
            Some(Err(e)) => {
                chain.push(None);
                scheme.push(None);
                conf.push(None);
                positions.push(Robj::from(()));
                residues.push(Robj::from(()));
                error.push(Some(e));
            }
            Some(Ok((c, sc, cf, pos, res))) => {
                chain.push(Some(c));
                scheme.push(Some(sc));
                conf.push(Some(cf));
                positions.push(Robj::from(pos));
                residues.push(Robj::from(res));
                error.push(None);
            }
        }
    }
    list!(
        chain = chain,
        scheme = scheme,
        confidence = conf,
        positions = List::from_values(positions),
        residues = List::from_values(residues),
        error = error,
    )
}

fn collect_segment_columns(rows: Vec<SegRow>) -> List {
    let n = rows.len();
    let mut prefix: Vec<Option<String>> = Vec::with_capacity(n);
    let mut fr1: Vec<Option<String>> = Vec::with_capacity(n);
    let mut cdr1: Vec<Option<String>> = Vec::with_capacity(n);
    let mut fr2: Vec<Option<String>> = Vec::with_capacity(n);
    let mut cdr2: Vec<Option<String>> = Vec::with_capacity(n);
    let mut fr3: Vec<Option<String>> = Vec::with_capacity(n);
    let mut cdr3: Vec<Option<String>> = Vec::with_capacity(n);
    let mut fr4: Vec<Option<String>> = Vec::with_capacity(n);
    let mut postfix: Vec<Option<String>> = Vec::with_capacity(n);
    let mut error: Vec<Option<String>> = Vec::with_capacity(n);
    for row in rows {
        match row {
            None => {
                prefix.push(None);
                fr1.push(None);
                cdr1.push(None);
                fr2.push(None);
                cdr2.push(None);
                fr3.push(None);
                cdr3.push(None);
                fr4.push(None);
                postfix.push(None);
                error.push(None);
            }
            Some(Err(e)) => {
                prefix.push(None);
                fr1.push(None);
                cdr1.push(None);
                fr2.push(None);
                cdr2.push(None);
                fr3.push(None);
                cdr3.push(None);
                fr4.push(None);
                postfix.push(None);
                error.push(Some(e));
            }
            Some(Ok([p, f1, c1, f2, c2, f3, c3, f4, pf])) => {
                prefix.push(Some(p));
                fr1.push(Some(f1));
                cdr1.push(Some(c1));
                fr2.push(Some(f2));
                cdr2.push(Some(c2));
                fr3.push(Some(f3));
                cdr3.push(Some(c3));
                fr4.push(Some(f4));
                postfix.push(Some(pf));
                error.push(None);
            }
        }
    }
    list!(
        prefix = prefix,
        fr1 = fr1,
        cdr1 = cdr1,
        fr2 = fr2,
        cdr2 = cdr2,
        fr3 = fr3,
        cdr3 = cdr3,
        fr4 = fr4,
        postfix = postfix,
        error = error,
    )
}

// R-side validates chains/scheme before calling, so build_inner_annotator's
// Err path is unreachable in normal operation.
#[extendr]
fn numbering_batch(
    seqs: Strings,
    chains: Strings,
    scheme: String,
    min_confidence: Nullable<f64>,
) -> Result<List> {
    let annotator = build_inner_annotator(chains, scheme, min_confidence)?;
    let inputs = materialize_inputs(seqs);
    let rows: Vec<NumRow> = inputs
        .par_iter()
        .map_with(annotator, |ann, opt_seq| {
            opt_seq.as_deref().and_then(|s| process_number(ann, s))
        })
        .collect();
    Ok(collect_number_columns(rows))
}

#[extendr]
fn numbering_batch_with(seqs: Strings, annotator: &Annotator) -> List {
    let cloned = annotator.inner.clone();
    let inputs = materialize_inputs(seqs);
    let rows: Vec<NumRow> = inputs
        .par_iter()
        .map_with(cloned, |ann, opt_seq| {
            opt_seq.as_deref().and_then(|s| process_number(ann, s))
        })
        .collect();
    collect_number_columns(rows)
}

#[extendr]
fn segmentation_batch(
    seqs: Strings,
    chains: Strings,
    scheme: String,
    min_confidence: Nullable<f64>,
) -> Result<List> {
    let annotator = build_inner_annotator(chains, scheme, min_confidence)?;
    let inputs = materialize_inputs(seqs);
    let rows: Vec<SegRow> = inputs
        .par_iter()
        .map_with(annotator, |ann, opt_seq| {
            opt_seq.as_deref().and_then(|s| process_segment(ann, s))
        })
        .collect();
    Ok(collect_segment_columns(rows))
}

#[extendr]
fn segmentation_batch_with(seqs: Strings, annotator: &Annotator) -> List {
    let cloned = annotator.inner.clone();
    let inputs = materialize_inputs(seqs);
    let rows: Vec<SegRow> = inputs
        .par_iter()
        .map_with(cloned, |ann, opt_seq| {
            opt_seq.as_deref().and_then(|s| process_segment(ann, s))
        })
        .collect();
    collect_segment_columns(rows)
}

extendr_module! {
    mod immunum;
    fn immunum_version;
    fn numbering_batch;
    fn numbering_batch_with;
    fn segmentation_batch;
    fn segmentation_batch_with;
    impl Annotator;
}
