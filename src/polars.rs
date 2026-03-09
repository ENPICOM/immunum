use crate::annotator::Annotator;
use polars::prelude::*;
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
