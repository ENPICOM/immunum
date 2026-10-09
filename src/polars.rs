use crate::annotator::{Annotator, NumberingResult, SegmentResult};
use crate::numbering::SEGMENT_NAMES;
use polars::prelude::*;
use polars_arrow::bitmap::MutableBitmap;
use polars_arrow::offset::Offsets;
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
    chains: Vec<String>,
    scheme: String,
    min_confidence: Option<f32>,
}

impl NumberFuncKwargs {
    fn annotator(&self) -> PolarsResult<Annotator> {
        Annotator::from_names(
            self.chains.iter().map(String::as_str),
            &self.scheme,
            self.min_confidence,
        )
        .map_err(|e| polars_err!(InvalidOperation: "{}", e))
    }
}

// ── Numbering ────────────────────────────────────────────────────────────────

// One numbered residue: its position and its amino acid
fn residue_dtype() -> DataType {
    DataType::Struct(vec![
        Field::new("position".into(), DataType::String),
        Field::new("residue".into(), DataType::String),
    ])
}

// The fields `Annotator.number` returns in Python and JavaScript, as Polars types
fn numbering_struct_output(_input_fields: &[Field]) -> PolarsResult<Field> {
    let fields = vec![
        Field::new("chain".into(), DataType::String),
        Field::new("scheme".into(), DataType::String),
        Field::new("confidence".into(), DataType::Float32),
        Field::new(
            "numbering".into(),
            DataType::List(Box::new(residue_dtype())),
        ),
        Field::new("query_start".into(), DataType::UInt32),
        Field::new("query_end".into(), DataType::UInt32),
        Field::new("error".into(), DataType::String),
    ];
    Ok(Field::new("numbering".into(), DataType::Struct(fields)))
}

#[polars_expr(output_type_func=numbering_struct_output)]
fn numbering_class_struct_expr(inputs: &[Series], kwargs: NumberKwargs) -> PolarsResult<Series> {
    numbering_series(inputs[0].str()?, kwargs.annotator)
}

#[polars_expr(output_type_func=numbering_struct_output)]
fn numbering_struct_expr(inputs: &[Series], kwargs: NumberFuncKwargs) -> PolarsResult<Series> {
    numbering_series(inputs[0].str()?, kwargs.annotator()?)
}

// A sequence's numbering, with its positions and residues formatted for the output columns
struct NumberedRow {
    result: NumberingResult,
    positions: Vec<String>,
    residues: Vec<String>,
}

// One struct per sequence, with the fields of `numbering_struct_output`
fn numbering_series(ca: &StringChunked, annotator: Annotator) -> PolarsResult<Series> {
    let len = ca.len();
    let values: Vec<Option<&str>> = ca.into_iter().collect();
    let rows: Vec<Option<Result<NumberedRow, String>>> = POOL.install(|| {
        values
            .par_iter()
            .map_with(annotator, |ann, opt_v| {
                let value = (*opt_v)?;
                Some(
                    ann.number(value)
                        .map(|result| {
                            let (positions, residues) = result
                                .residues(value)
                                .expect("`result` numbered `value`")
                                .map(|(pos, ch)| (pos.to_string(), ch.to_string()))
                                .unzip();
                            NumberedRow {
                                result,
                                positions,
                                residues,
                            }
                        })
                        .map_err(|e| e.to_string()),
                )
            })
            .collect()
    });

    // Every numbered residue of every row in one flat struct column, which the rows' lists offset into
    let numbered_residues: usize = rows
        .iter()
        .flatten()
        .flatten()
        .map(|row| row.positions.len())
        .sum();
    let mut positions = StringChunkedBuilder::new("position".into(), numbered_residues);
    let mut residues = StringChunkedBuilder::new("residue".into(), numbered_residues);
    for row in rows.iter().flatten().flatten() {
        row.positions.iter().for_each(|p| positions.append_value(p));
        row.residues.iter().for_each(|r| residues.append_value(r));
    }
    let flat = StructChunked::from_series(
        "".into(),
        numbered_residues,
        [
            positions.finish().into_series(),
            residues.finish().into_series(),
        ]
        .iter(),
    )?
    .into_series();

    let mut chain = StringChunkedBuilder::new("chain".into(), len);
    let mut scheme = StringChunkedBuilder::new("scheme".into(), len);
    let mut confidence = PrimitiveChunkedBuilder::<Float32Type>::new("confidence".into(), len);
    let mut offsets = Offsets::<i64>::with_capacity(len);
    let mut numbered = MutableBitmap::with_capacity(len);
    let mut query_start = PrimitiveChunkedBuilder::<UInt32Type>::new("query_start".into(), len);
    let mut query_end = PrimitiveChunkedBuilder::<UInt32Type>::new("query_end".into(), len);
    let mut error = StringChunkedBuilder::new("error".into(), len);
    for row in &rows {
        match row {
            Some(Ok(row)) => {
                chain.append_value(row.result.chain.to_string());
                scheme.append_value(row.result.scheme.to_string());
                confidence.append_value(row.result.confidence);
                offsets.try_push(row.positions.len())?;
                numbered.push(true);
                query_start.append_value(row.result.query_start as u32);
                query_end.append_value(row.result.query_end as u32);
                error.append_null();
            }
            _ => {
                chain.append_null();
                scheme.append_null();
                confidence.append_null();
                offsets.try_push(0)?;
                numbered.push(false);
                query_start.append_null();
                query_end.append_null();
                error.append_option(row.as_ref().and_then(|r| r.as_ref().err()));
            }
        }
    }

    let values = flat.rechunk().chunks()[0].clone();
    let numbering = ListChunked::with_chunk(
        "numbering".into(),
        LargeListArray::try_new(
            LargeListArray::default_datatype(values.dtype().clone()),
            offsets.into(),
            values,
            Some(numbered.into()),
        )?,
    );

    let fields = [
        chain.finish().into_series(),
        scheme.finish().into_series(),
        confidence.finish().into_series(),
        numbering.into_series(),
        query_start.finish().into_series(),
        query_end.finish().into_series(),
        error.finish().into_series(),
    ];
    StructChunked::from_series(ca.name().clone(), len, fields.iter()).map(|ca| ca.into_series())
}

// ── Segmentation ─────────────────────────────────────────────────────────────

fn segmentation_struct_output(_input_fields: &[Field]) -> PolarsResult<Field> {
    let fields = SEGMENT_NAMES
        .into_iter()
        .chain(["error"])
        .map(|name| Field::new(name.into(), DataType::String))
        .collect();
    Ok(Field::new("segmentation".into(), DataType::Struct(fields)))
}

#[polars_expr(output_type_func=segmentation_struct_output)]
fn segmentation_class_struct_expr(inputs: &[Series], kwargs: NumberKwargs) -> PolarsResult<Series> {
    segmentation_series(inputs[0].str()?, kwargs.annotator)
}

#[polars_expr(output_type_func=segmentation_struct_output)]
fn segmentation_struct_expr(inputs: &[Series], kwargs: NumberFuncKwargs) -> PolarsResult<Series> {
    segmentation_series(inputs[0].str()?, kwargs.annotator()?)
}

// One struct per sequence: a field per segment, and the error when the sequence couldn't be segmented
fn segmentation_series(ca: &StringChunked, annotator: Annotator) -> PolarsResult<Series> {
    let len = ca.len();
    let values: Vec<Option<&str>> = ca.into_iter().collect();
    let results: Vec<Option<Result<SegmentResult, String>>> = POOL.install(|| {
        values
            .par_iter()
            .map_with(annotator, |ann, opt_v| {
                let value = (*opt_v)?;
                Some(ann.segment(value).map_err(|e| e.to_string()))
            })
            .collect()
    });

    let mut segments = SEGMENT_NAMES.map(|name| StringChunkedBuilder::new(name.into(), len));
    let mut errors = StringChunkedBuilder::new("error".into(), len);
    for row in &results {
        match row {
            Some(Ok(s)) => {
                for (builder, (_, residues)) in segments.iter_mut().zip(s.regions()) {
                    builder.append_value(residues);
                }
                errors.append_null();
            }
            _ => {
                segments
                    .iter_mut()
                    .for_each(|builder| builder.append_null());
                errors.append_option(row.as_ref().and_then(|r| r.as_ref().err()));
            }
        }
    }

    let fields: Vec<Series> = segments
        .into_iter()
        .chain([errors])
        .map(|builder| builder.finish().into_series())
        .collect();
    StructChunked::from_series(ca.name().clone(), len, fields.iter()).map(|ca| ca.into_series())
}
