use crate::{annotator::Annotator, numbering::segment, Chain, Scheme};
use polars::chunked_array::builder::AnonymousListBuilder;
use polars::prelude::*;
use polars_core::utils::rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use polars_core::POOL;
use pyo3_polars::derive::polars_expr;
use pyo3_polars::PolarsAllocator;
use serde::{Deserialize, Deserializer, Serialize};

#[global_allocator]
static ALLOC: PolarsAllocator = PolarsAllocator::new();

fn deserialize_annotator_from_bytes<'de, D: Deserializer<'de>>(
    d: D,
) -> Result<Annotator, D::Error> {
    struct BytesVisitor;
    impl<'de> serde::de::Visitor<'de> for BytesVisitor {
        type Value = Vec<u8>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "byte array")
        }
        fn visit_bytes<E: serde::de::Error>(self, v: &[u8]) -> Result<Vec<u8>, E> {
            Ok(v.to_vec())
        }
        fn visit_byte_buf<E: serde::de::Error>(self, v: Vec<u8>) -> Result<Vec<u8>, E> {
            Ok(v)
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut seq: A) -> Result<Vec<u8>, A::Error> {
            let mut bytes = Vec::new();
            while let Some(b) = seq.next_element::<u8>()? {
                bytes.push(b);
            }
            Ok(bytes)
        }
    }
    let bytes = d.deserialize_bytes(BytesVisitor)?;
    postcard::from_bytes(&bytes).map_err(serde::de::Error::custom)
}

#[derive(Serialize, Deserialize)]
struct NumberKwargs {
    #[serde(deserialize_with = "deserialize_annotator_from_bytes")]
    annotator: Annotator,
}

#[derive(Serialize, Deserialize)]
struct NumberFuncKwargs {
    chains: Vec<Chain>,
    scheme: Scheme,
}

// ── Numbering ────────────────────────────────────────────────────────────────

fn numbering_class_struct_output(_input_fields: &[Field]) -> PolarsResult<Field> {
    let inner_fields = vec![
        Field::new("position".into(), DataType::String),
        Field::new("residue".into(), DataType::String),
    ];
    let fields = vec![
        Field::new("chain".into(), DataType::String),
        Field::new("scheme".into(), DataType::String),
        Field::new("confidence".into(), DataType::Float32),
        Field::new(
            "numbering".into(),
            DataType::List(Box::new(DataType::Struct(inner_fields))),
        ),
    ];
    Ok(Field::new("numbering".into(), DataType::Struct(fields)))
}

fn numbering_struct_output(_input_fields: &[Field]) -> PolarsResult<Field> {
    let fields = vec![
        Field::new("chain".into(), DataType::String),
        Field::new("scheme".into(), DataType::String),
        Field::new(
            "positions".into(),
            DataType::List(Box::new(DataType::String)),
        ),
        Field::new(
            "residues".into(),
            DataType::List(Box::new(DataType::String)),
        ),
    ];
    Ok(Field::new("numbering".into(), DataType::Struct(fields)))
}

#[polars_expr(output_type_func=numbering_class_struct_output)]
fn numbering_class_struct_expr(inputs: &[Series], kwargs: NumberKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let len = ca.len();
    let name = ca.name().clone();

    // Build per-row structs inside the parallel closure so Series allocation
    // is parallelized. Series is Send+Sync in Polars, so this is safe.
    type RowResult = Option<(String, String, f32, Series)>;
    let values: Vec<Option<&str>> = ca.into_iter().collect();
    let results: Vec<RowResult> = POOL.install(|| {
        values
            .par_iter()
            .map_with(kwargs.annotator, |ann, opt_v| {
                let value = (*opt_v)?;
                let result = ann.number(value).ok()?;
                let (positions, residues): (Vec<String>, Vec<String>) = result
                    .positions
                    .iter()
                    .zip(value.chars())
                    .map(|(pos, ch)| (pos.to_string(), ch.to_string()))
                    .unzip();
                let n = positions.len();
                let pos_series = Series::new("position".into(), positions);
                let res_series = Series::new("residue".into(), residues);
                let row_struct =
                    StructChunked::from_series("".into(), n, [pos_series, res_series].iter())
                        .ok()?
                        .into_series();
                Some((
                    result.chain.to_string(),
                    result.scheme.to_string(),
                    result.confidence,
                    row_struct,
                ))
            })
            .collect()
    });

    let mut chain_vec: Vec<Option<String>> = Vec::with_capacity(len);
    let mut scheme_vec: Vec<Option<String>> = Vec::with_capacity(len);
    let mut confidence_vec: Vec<Option<f32>> = Vec::with_capacity(len);
    let mut row_structs: Vec<Option<Series>> = Vec::with_capacity(len);

    for row in results {
        match row {
            None => {
                chain_vec.push(None);
                scheme_vec.push(None);
                confidence_vec.push(None);
                row_structs.push(None);
            }
            Some((chain, scheme, confidence, row_struct)) => {
                chain_vec.push(Some(chain));
                scheme_vec.push(Some(scheme));
                confidence_vec.push(Some(confidence));
                row_structs.push(Some(row_struct));
            }
        }
    }

    let mut numbering_builder = AnonymousListBuilder::new("numbering".into(), len, None);
    for opt_s in &row_structs {
        match opt_s {
            None => numbering_builder.append_null(),
            Some(s) => numbering_builder.append_series(s)?,
        }
    }

    let chain_series = Series::new("chain".into(), chain_vec);
    let scheme_series = Series::new("scheme".into(), scheme_vec);
    let confidence_series = Series::new("confidence".into(), confidence_vec);
    let numbering_series = numbering_builder.finish().into_series();

    let fields = [
        chain_series,
        scheme_series,
        confidence_series,
        numbering_series,
    ];
    StructChunked::from_series(name, len, fields.iter()).map(|ca| ca.into_series())
}

#[polars_expr(output_type_func=numbering_struct_output)]
fn numbering_struct_expr(inputs: &[Series], kwargs: NumberFuncKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let len = ca.len();
    let name = ca.name().clone();
    let annotator: Annotator = match Annotator::new(kwargs.chains.as_slice(), kwargs.scheme, None) {
        Ok(a) => a,
        Err(e) => polars_bail!(InvalidOperation: "{}", e),
    };

    type ResultType = Vec<Option<(String, String, Series, Series)>>;
    let values: Vec<Option<&str>> = ca.into_iter().collect();
    let results: ResultType = POOL.install(|| {
        values
            .par_iter()
            .map_with(annotator, |ann, opt_v| {
                let value = (*opt_v)?;
                let result = ann.number(value).ok()?;
                let (positions, residues): (Vec<String>, Vec<String>) = result
                    .positions
                    .iter()
                    .zip(value.chars())
                    .map(|(pos, ch)| (pos.to_string(), ch.to_string()))
                    .unzip();
                Some((
                    result.chain.to_string(),
                    result.scheme.to_string(),
                    Series::new("".into(), positions),
                    Series::new("".into(), residues),
                ))
            })
            .collect()
    });

    let mut chain_builder = StringChunkedBuilder::new("chain".into(), len);
    let mut scheme_builder = StringChunkedBuilder::new("scheme".into(), len);
    let mut positions_builder = ListStringChunkedBuilder::new("positions".into(), len, len);
    let mut residues_builder = ListStringChunkedBuilder::new("residues".into(), len, len);

    for row in results {
        match row {
            None => {
                chain_builder.append_null();
                scheme_builder.append_null();
                positions_builder.append_null();
                residues_builder.append_null();
            }
            Some((chain, scheme, positions, residues)) => {
                chain_builder.append_value(&chain);
                scheme_builder.append_value(&scheme);
                positions_builder.append_series(&positions)?;
                residues_builder.append_series(&residues)?;
            }
        }
    }

    let fields = [
        chain_builder.finish().into_series(),
        scheme_builder.finish().into_series(),
        positions_builder.finish().into_series(),
        residues_builder.finish().into_series(),
    ];
    StructChunked::from_series(name, len, fields.iter()).map(|ca| ca.into_series())
}

// ── Segmentation ─────────────────────────────────────────────────────────────

fn segmentation_struct_output(_input_fields: &[Field]) -> PolarsResult<Field> {
    let fields = vec![
        Field::new("prefix".into(), DataType::String),
        Field::new("fr1".into(), DataType::String),
        Field::new("cdr1".into(), DataType::String),
        Field::new("fr2".into(), DataType::String),
        Field::new("cdr2".into(), DataType::String),
        Field::new("fr3".into(), DataType::String),
        Field::new("cdr3".into(), DataType::String),
        Field::new("fr4".into(), DataType::String),
        Field::new("postfix".into(), DataType::String),
    ];
    Ok(Field::new("segmentation".into(), DataType::Struct(fields)))
}

#[polars_expr(output_type_func=segmentation_struct_output)]
fn segmentation_class_struct_expr(inputs: &[Series], kwargs: NumberKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let len = ca.len();
    let name = ca.name().clone();

    type SegResult = Option<[String; 9]>;
    let values: Vec<Option<&str>> = ca.into_iter().collect();
    let results: Vec<SegResult> = POOL.install(|| {
        values
            .par_iter()
            .map_with(kwargs.annotator, |ann, opt_v| {
                let value = (*opt_v)?;
                let result = ann.number(value).ok()?;
                let s = segment(&result.positions, value, result.scheme);
                let get = |k: &str| s.get(k).map(|v| v.as_str()).unwrap_or("").to_string();
                Some([
                    get("prefix"),
                    get("fr1"),
                    get("cdr1"),
                    get("fr2"),
                    get("cdr2"),
                    get("fr3"),
                    get("cdr3"),
                    get("fr4"),
                    get("postfix"),
                ])
            })
            .collect()
    });

    let mut prefix_b = StringChunkedBuilder::new("prefix".into(), len);
    let mut fr1_b = StringChunkedBuilder::new("fr1".into(), len);
    let mut cdr1_b = StringChunkedBuilder::new("cdr1".into(), len);
    let mut fr2_b = StringChunkedBuilder::new("fr2".into(), len);
    let mut cdr2_b = StringChunkedBuilder::new("cdr2".into(), len);
    let mut fr3_b = StringChunkedBuilder::new("fr3".into(), len);
    let mut cdr3_b = StringChunkedBuilder::new("cdr3".into(), len);
    let mut fr4_b = StringChunkedBuilder::new("fr4".into(), len);
    let mut postfix_b = StringChunkedBuilder::new("postfix".into(), len);

    for row in results {
        match row {
            None => {
                prefix_b.append_null();
                fr1_b.append_null();
                cdr1_b.append_null();
                fr2_b.append_null();
                cdr2_b.append_null();
                fr3_b.append_null();
                cdr3_b.append_null();
                fr4_b.append_null();
                postfix_b.append_null();
            }
            Some([prefix, fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4, postfix]) => {
                prefix_b.append_value(&prefix);
                fr1_b.append_value(&fr1);
                cdr1_b.append_value(&cdr1);
                fr2_b.append_value(&fr2);
                cdr2_b.append_value(&cdr2);
                fr3_b.append_value(&fr3);
                cdr3_b.append_value(&cdr3);
                fr4_b.append_value(&fr4);
                postfix_b.append_value(&postfix);
            }
        }
    }

    let fields = [
        prefix_b.finish().into_series(),
        fr1_b.finish().into_series(),
        cdr1_b.finish().into_series(),
        fr2_b.finish().into_series(),
        cdr2_b.finish().into_series(),
        fr3_b.finish().into_series(),
        cdr3_b.finish().into_series(),
        fr4_b.finish().into_series(),
        postfix_b.finish().into_series(),
    ];
    StructChunked::from_series(name, len, fields.iter()).map(|ca| ca.into_series())
}

#[polars_expr(output_type_func=segmentation_struct_output)]
fn segmentation_struct_expr(inputs: &[Series], kwargs: NumberFuncKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let len = ca.len();
    let name = ca.name().clone();
    let annotator: Annotator = match Annotator::new(kwargs.chains.as_slice(), kwargs.scheme) {
        Ok(a) => a,
        Err(e) => polars_bail!(InvalidOperation: "{}", e),
    };

    type SegResult = Option<[String; 9]>;
    let values: Vec<Option<&str>> = ca.into_iter().collect();
    let results: Vec<SegResult> = POOL.install(|| {
        values
            .par_iter()
            .map_with(annotator, |ann, opt_v| {
                let value = (*opt_v)?;
                let result = ann.number(value).ok()?;
                let s = segment(&result.positions, value, result.scheme);
                let get = |k: &str| s.get(k).map(|v| v.as_str()).unwrap_or("").to_string();
                Some([
                    get("prefix"),
                    get("fr1"),
                    get("cdr1"),
                    get("fr2"),
                    get("cdr2"),
                    get("fr3"),
                    get("cdr3"),
                    get("fr4"),
                    get("postfix"),
                ])
            })
            .collect()
    });

    let mut prefix_b = StringChunkedBuilder::new("prefix".into(), len);
    let mut fr1_b = StringChunkedBuilder::new("fr1".into(), len);
    let mut cdr1_b = StringChunkedBuilder::new("cdr1".into(), len);
    let mut fr2_b = StringChunkedBuilder::new("fr2".into(), len);
    let mut cdr2_b = StringChunkedBuilder::new("cdr2".into(), len);
    let mut fr3_b = StringChunkedBuilder::new("fr3".into(), len);
    let mut cdr3_b = StringChunkedBuilder::new("cdr3".into(), len);
    let mut fr4_b = StringChunkedBuilder::new("fr4".into(), len);
    let mut postfix_b = StringChunkedBuilder::new("postfix".into(), len);

    for row in results {
        match row {
            None => {
                prefix_b.append_null();
                fr1_b.append_null();
                cdr1_b.append_null();
                fr2_b.append_null();
                cdr2_b.append_null();
                fr3_b.append_null();
                cdr3_b.append_null();
                fr4_b.append_null();
                postfix_b.append_null();
            }
            Some([prefix, fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4, postfix]) => {
                prefix_b.append_value(&prefix);
                fr1_b.append_value(&fr1);
                cdr1_b.append_value(&cdr1);
                fr2_b.append_value(&fr2);
                cdr2_b.append_value(&cdr2);
                fr3_b.append_value(&fr3);
                cdr3_b.append_value(&cdr3);
                fr4_b.append_value(&fr4);
                postfix_b.append_value(&postfix);
            }
        }
    }

    let fields = [
        prefix_b.finish().into_series(),
        fr1_b.finish().into_series(),
        cdr1_b.finish().into_series(),
        fr2_b.finish().into_series(),
        cdr2_b.finish().into_series(),
        fr3_b.finish().into_series(),
        cdr3_b.finish().into_series(),
        fr4_b.finish().into_series(),
        postfix_b.finish().into_series(),
    ];
    StructChunked::from_series(name, len, fields.iter()).map(|ca| ca.into_series())
}
