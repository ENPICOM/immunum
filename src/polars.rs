use crate::annotator::Annotator;
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

#[polars_expr(output_type=Int64)]
fn numbering_end_expr(inputs: &[Series], kwargs: NumberKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let mut builder = PrimitiveChunkedBuilder::<Int64Type>::new(ca.name().clone(), ca.len());
    ca.into_iter().for_each(|opt_v: Option<&str>| match opt_v {
        None => builder.append_null(),
        Some(value) => match kwargs.annotator.number(value) {
            Ok(result) => builder.append_value(result.end as i64),
            Err(_) => builder.append_null(),
        },
    });

    Ok(builder.finish().into_series())
}
