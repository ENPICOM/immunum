use crate::{annotator::Annotator, numbering::segment, Chain, Scheme};
use polars::chunked_array::builder::AnonymousListBuilder;
use polars::prelude::*;
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

    // Collect per-row data; row_structs holds owned Series to keep them alive
    // while AnonymousListBuilder borrows them.
    let mut chain_vec: Vec<Option<String>> = Vec::with_capacity(len);
    let mut scheme_vec: Vec<Option<String>> = Vec::with_capacity(len);
    let mut confidence_vec: Vec<Option<f32>> = Vec::with_capacity(len);
    let mut row_structs: Vec<Option<Series>> = Vec::with_capacity(len);

    for opt_v in ca.into_iter() {
        match opt_v {
            None => {
                chain_vec.push(None);
                scheme_vec.push(None);
                confidence_vec.push(None);
                row_structs.push(None);
            }
            Some(value) => match kwargs.annotator.number(value) {
                Err(_) => {
                    chain_vec.push(None);
                    scheme_vec.push(None);
                    confidence_vec.push(None);
                    row_structs.push(None);
                }
                Ok(result) => {
                    chain_vec.push(Some(result.chain.to_string()));
                    scheme_vec.push(Some(result.scheme.to_string()));
                    confidence_vec.push(Some(result.confidence));
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
                        StructChunked::from_series("".into(), n, [pos_series, res_series].iter())?
                            .into_series();
                    row_structs.push(Some(row_struct));
                }
            },
        }
    }

    let inner_dtype = DataType::Struct(vec![
        Field::new("position".into(), DataType::String),
        Field::new("residue".into(), DataType::String),
    ]);
    let mut numbering_builder =
        AnonymousListBuilder::new("numbering".into(), len, Some(inner_dtype));
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
    StructChunked::from_series(ca.name().clone(), len, fields.iter()).map(|ca| ca.into_series())
}

#[polars_expr(output_type_func=numbering_struct_output)]
fn numbering_struct_expr(inputs: &[Series], kwargs: NumberFuncKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let len = ca.len();
    let mut chain_builder = StringChunkedBuilder::new("chain".into(), len);
    let mut scheme_builder = StringChunkedBuilder::new("scheme".into(), len);
    let mut positions_builder = ListStringChunkedBuilder::new("positions".into(), len, len);
    let mut residues_builder = ListStringChunkedBuilder::new("residues".into(), len, len);
    let annotator: Annotator = match Annotator::new(kwargs.chains.as_slice(), kwargs.scheme) {
        Ok(a) => a,
        Err(e) => polars_bail!(InvalidOperation: "{}", e),
    };

    ca.into_iter().try_for_each(|opt_v| -> PolarsResult<()> {
        match opt_v {
            None => {
                chain_builder.append_null();
                scheme_builder.append_null();
                positions_builder.append_null();
                residues_builder.append_null();
            }
            Some(value) => match annotator.number(value) {
                Err(_) => {
                    chain_builder.append_null();
                    scheme_builder.append_null();
                    positions_builder.append_null();
                    residues_builder.append_null();
                }
                Ok(result) => {
                    chain_builder.append_value(result.chain.to_string());
                    scheme_builder.append_value(result.scheme.to_string());
                    let (positions, residues): (Vec<String>, Vec<String>) = result
                        .positions
                        .iter()
                        .zip(value.chars())
                        .map(|(pos, ch)| (pos.to_string(), ch.to_string()))
                        .unzip();
                    positions_builder.append_series(&Series::new("".into(), positions))?;
                    residues_builder.append_series(&Series::new("".into(), residues))?;
                }
            },
        }
        Ok(())
    })?;

    let fields = [
        chain_builder.finish().into_series(),
        scheme_builder.finish().into_series(),
        positions_builder.finish().into_series(),
        residues_builder.finish().into_series(),
    ];
    StructChunked::from_series(ca.name().clone(), len, fields.iter()).map(|ca| ca.into_series())
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
    let mut prefix_b = StringChunkedBuilder::new("prefix".into(), len);
    let mut fr1_b = StringChunkedBuilder::new("fr1".into(), len);
    let mut cdr1_b = StringChunkedBuilder::new("cdr1".into(), len);
    let mut fr2_b = StringChunkedBuilder::new("fr2".into(), len);
    let mut cdr2_b = StringChunkedBuilder::new("cdr2".into(), len);
    let mut fr3_b = StringChunkedBuilder::new("fr3".into(), len);
    let mut cdr3_b = StringChunkedBuilder::new("cdr3".into(), len);
    let mut fr4_b = StringChunkedBuilder::new("fr4".into(), len);
    let mut postfix_b = StringChunkedBuilder::new("postfix".into(), len);

    ca.into_iter().try_for_each(|opt_v| -> PolarsResult<()> {
        match opt_v {
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
            Some(value) => match kwargs.annotator.number(value) {
                Err(_) => {
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
                Ok(result) => {
                    let s = segment(&result.positions, value, result.scheme);
                    prefix_b.append_value(s.get("prefix").map(|v| v.as_str()).unwrap_or(""));
                    fr1_b.append_value(s.get("fr1").map(|v| v.as_str()).unwrap_or(""));
                    cdr1_b.append_value(s.get("cdr1").map(|v| v.as_str()).unwrap_or(""));
                    fr2_b.append_value(s.get("fr2").map(|v| v.as_str()).unwrap_or(""));
                    cdr2_b.append_value(s.get("cdr2").map(|v| v.as_str()).unwrap_or(""));
                    fr3_b.append_value(s.get("fr3").map(|v| v.as_str()).unwrap_or(""));
                    cdr3_b.append_value(s.get("cdr3").map(|v| v.as_str()).unwrap_or(""));
                    fr4_b.append_value(s.get("fr4").map(|v| v.as_str()).unwrap_or(""));
                    postfix_b.append_value(s.get("postfix").map(|v| v.as_str()).unwrap_or(""));
                }
            },
        }
        Ok(())
    })?;

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
    StructChunked::from_series(ca.name().clone(), len, fields.iter()).map(|ca| ca.into_series())
}

#[polars_expr(output_type_func=segmentation_struct_output)]
fn segmentation_struct_expr(inputs: &[Series], kwargs: NumberFuncKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let len = ca.len();
    let mut prefix_b = StringChunkedBuilder::new("prefix".into(), len);
    let mut fr1_b = StringChunkedBuilder::new("fr1".into(), len);
    let mut cdr1_b = StringChunkedBuilder::new("cdr1".into(), len);
    let mut fr2_b = StringChunkedBuilder::new("fr2".into(), len);
    let mut cdr2_b = StringChunkedBuilder::new("cdr2".into(), len);
    let mut fr3_b = StringChunkedBuilder::new("fr3".into(), len);
    let mut cdr3_b = StringChunkedBuilder::new("cdr3".into(), len);
    let mut fr4_b = StringChunkedBuilder::new("fr4".into(), len);
    let mut postfix_b = StringChunkedBuilder::new("postfix".into(), len);
    let annotator: Annotator = match Annotator::new(kwargs.chains.as_slice(), kwargs.scheme) {
        Ok(a) => a,
        Err(e) => polars_bail!(InvalidOperation: "{}", e),
    };

    ca.into_iter().try_for_each(|opt_v| -> PolarsResult<()> {
        match opt_v {
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
            Some(value) => match annotator.number(value) {
                Err(_) => {
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
                Ok(result) => {
                    let s = segment(&result.positions, value, result.scheme);
                    prefix_b.append_value(s.get("prefix").map(|v| v.as_str()).unwrap_or(""));
                    fr1_b.append_value(s.get("fr1").map(|v| v.as_str()).unwrap_or(""));
                    cdr1_b.append_value(s.get("cdr1").map(|v| v.as_str()).unwrap_or(""));
                    fr2_b.append_value(s.get("fr2").map(|v| v.as_str()).unwrap_or(""));
                    cdr2_b.append_value(s.get("cdr2").map(|v| v.as_str()).unwrap_or(""));
                    fr3_b.append_value(s.get("fr3").map(|v| v.as_str()).unwrap_or(""));
                    cdr3_b.append_value(s.get("cdr3").map(|v| v.as_str()).unwrap_or(""));
                    fr4_b.append_value(s.get("fr4").map(|v| v.as_str()).unwrap_or(""));
                    postfix_b.append_value(s.get("postfix").map(|v| v.as_str()).unwrap_or(""));
                }
            },
        }
        Ok(())
    })?;

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
    StructChunked::from_series(ca.name().clone(), len, fields.iter()).map(|ca| ca.into_series())
}
