use crate::annotator::Annotator;
use crate::types::{Chain, Scheme};
use polars::prelude::*;
use pyo3_polars::derive::polars_expr;
use pyo3_polars::PolarsAllocator;
use serde::{Deserialize, Serialize};
use std::str::FromStr;

#[global_allocator]
static ALLOC: PolarsAllocator = PolarsAllocator::new();

#[derive(Serialize, Deserialize)]
struct NumberKwargs {
    chains: Vec<String>,
    scheme: String,
}

#[polars_expr(output_type=Int64)]
fn numbering_end_expr(inputs: &[Series], kwargs: NumberKwargs) -> PolarsResult<Series> {
    let chains: Vec<Chain> = kwargs
        .chains
        .iter()
        .map(|s| Chain::from_str(s))
        .collect::<Result<_, _>>()
        .map_err(|e| PolarsError::ComputeError(e.to_string().into()))?;
    let scheme = Scheme::from_str(&kwargs.scheme)
        .map_err(|e| PolarsError::ComputeError(e.to_string().into()))?;
    let annotator = Annotator::new(&chains, scheme)
        .map_err(|e| PolarsError::ComputeError(e.to_string().into()))?;

    let ca = inputs[0].str()?;
    let mut builder = PrimitiveChunkedBuilder::<Int64Type>::new(ca.name().clone(), ca.len());
    ca.into_iter().for_each(|opt_v: Option<&str>| match opt_v {
        None => builder.append_null(),
        Some(value) => match annotator.number(value) {
            Ok(result) => builder.append_value(result.end as i64),
            Err(_) => builder.append_null(),
        },
    });

    Ok(builder.finish().into_series())
}
