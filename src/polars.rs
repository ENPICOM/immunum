use crate::annotator::Annotator;
use crate::python::AnnotatorSerializationWrapper;
use polars::prelude::*;
use postcard::{from_bytes, to_allocvec};
use pyo3_polars::derive::polars_expr;
use pyo3_polars::PolarsAllocator;

use serde::{Deserialize, Serialize};

#[global_allocator]
static ALLOC: PolarsAllocator = PolarsAllocator::new();

#[derive(Serialize, Deserialize)]
struct NumberKwargs {
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

// // the code below has been taken from: https://github.com/MarcoGorelli/polars-plugins-tutorial/issues/75
impl serde::Serialize for Annotator {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let bytes = to_allocvec(&self).map_err(serde::ser::Error::custom)?;
        serializer.serialize_bytes(&bytes)
    }
}

impl<'de> serde::Deserialize<'de> for Annotator {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        struct FieldVisitor;

        impl<'de> serde::de::Visitor<'de> for FieldVisitor {
            type Value = Annotator;
            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str(concat!(
                    "a byte array containing bincode-serialized ",
                    stringify!(Encoder),
                    " data"
                ))
            }
            fn visit_bytes<E>(self, v: &[u8]) -> Result<Self::Value, E>
            where
                E: serde::de::Error,
            {
                let annotator: Annotator = from_bytes(v).map_err(serde::de::Error::custom)?;

                Ok(annotator)
            }

            fn visit_byte_buf<E>(self, v: Vec<u8>) -> Result<Self::Value, E>
            where
                E: serde::de::Error,
            {
                self.visit_bytes(&v)
            }
        }
        deserializer.deserialize_bytes(FieldVisitor)
    }
}
