use polars::prelude::*;
use pyo3_polars::PolarsAllocator;

#[global_allocator]
static ALLOC: PolarsAllocator = PolarsAllocator::new();

#[polars_expr(output_type=String)]
fn tokenize_expr(inputs: &[Series], kwargs: TokenizeExprKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;

    let mut builder = ListPrimitiveChunkedBuilder::<UInt32Type>::new(
        ca.name().clone(),
        ca.len(),
        128,
        DataType::UInt32,
    );
    ca.into_iter().for_each(|opt_v: Option<&str>| match opt_v {
        None => builder.append_null(),
        Some(value) => builder.append_slice(
            kwargs
                .encoder
                .encode(value, kwargs.seed)
                .unwrap()
                .as_slice(),
        ),
    });
    let out = builder.finish().into_series();

    Ok(out)
}
