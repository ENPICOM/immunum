use crate::{annotator::Annotator, Chain, Scheme};
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

#[polars_expr(output_type_func=numbering_struct_output)]
fn numbering_class_struct_expr(inputs: &[Series], kwargs: NumberKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let len = ca.len();
    let name = ca.name().clone();

    let values: Vec<Option<&str>> = ca.into_iter().collect();
    let results: Vec<Option<(String, String, Vec<String>, Vec<String>)>> = POOL.install(|| {
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
                Some((
                    result.chain.to_string(),
                    result.scheme.to_string(),
                    positions,
                    residues,
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
                positions_builder.append_series(&Series::new("".into(), positions))?;
                residues_builder.append_series(&Series::new("".into(), residues))?;
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

#[derive(Serialize, Deserialize)]
struct NumberFuncKwargs {
    chains: Vec<Chain>,
    scheme: Scheme,
}

#[polars_expr(output_type_func=numbering_struct_output)]
fn numbering_struct_expr(inputs: &[Series], kwargs: NumberFuncKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let len = ca.len();
    let name = ca.name().clone();
    let annotator: Annotator = match Annotator::new(kwargs.chains.as_slice(), kwargs.scheme) {
        Ok(a) => a,
        Err(e) => polars_bail!(InvalidOperation: "{}", e),
    };

    let values: Vec<Option<&str>> = ca.into_iter().collect();
    let results: Vec<Option<(String, String, Vec<String>, Vec<String>)>> = POOL.install(|| {
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
                    positions,
                    residues,
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
                positions_builder.append_series(&Series::new("".into(), positions))?;
                residues_builder.append_series(&Series::new("".into(), residues))?;
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
