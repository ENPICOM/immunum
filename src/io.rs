//! Input parsing and output formatting for sequence records

use crate::annotator::{NumberingResult, SegmentResult};
use crate::error::SequenceError;
use crate::numbering::SEGMENT_NAMES;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::str::FromStr;

/// A raw input record (id + sequence)
pub struct Record {
    pub id: String,
    pub sequence: String,
}

/// A numbered record: input record paired with its numbering result
pub struct NumberedRecord {
    pub id: String,
    pub sequence: String,
    /// Which of the sequence's domains `result` numbers, 0-based, when every domain was numbered
    pub domain: Option<usize>,
    pub result: Result<NumberingResult, SequenceError>,
}

impl NumberedRecord {
    pub fn new(
        id: String,
        sequence: String,
        result: Result<NumberingResult, SequenceError>,
    ) -> Self {
        Self {
            id,
            sequence,
            domain: None,
            result,
        }
    }
    /// This record as the `domain`th domain of its sequence
    pub fn in_domain(self, domain: usize) -> Self {
        Self {
            domain: Some(domain),
            ..self
        }
    }
}

/// A segmented record: input record id paired with its FR/CDR split
pub struct SegmentedRecord {
    pub id: String,
    /// Which of the sequence's domains `result` splits, 0-based, when every domain was segmented
    pub domain: Option<usize>,
    pub result: Result<SegmentResult, SequenceError>,
}

impl SegmentedRecord {
    pub fn new(id: String, result: Result<SegmentResult, SequenceError>) -> Self {
        Self {
            id,
            domain: None,
            result,
        }
    }
    /// This record as the `domain`th domain of its sequence
    pub fn in_domain(self, domain: usize) -> Self {
        Self {
            domain: Some(domain),
            ..self
        }
    }
}

/// Output format
#[derive(Clone, Copy)]
pub enum OutputFormat {
    Tsv,
    Json,
    Jsonl,
}

impl FromStr for OutputFormat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "tsv" => Ok(Self::Tsv),
            "json" => Ok(Self::Json),
            "jsonl" => Ok(Self::Jsonl),
            _ => Err(format!(
                "unknown format '{}' (options: tsv, json, jsonl)",
                s
            )),
        }
    }
}

impl OutputFormat {
    /// Write numbered records in this format, one per sequence
    pub fn write(&self, writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
        match self {
            Self::Tsv => write_tsv(writer, records),
            Self::Json => write_json(writer, records),
            Self::Jsonl => write_jsonl(writer, records),
        }
    }

    /// Write format header (e.g. TSV column names, JSON array opening) for numbered records. With
    /// `all_domains`, records carry the index of the domain they number.
    pub fn write_header(&self, writer: &mut impl Write, all_domains: bool) -> io::Result<()> {
        match self {
            Self::Tsv => write_tsv_header(writer, all_domains),
            Self::Json => writeln!(writer, "["),
            Self::Jsonl => Ok(()),
        }
    }

    /// Write format header (e.g. TSV column names, JSON array opening) for segmented records. With
    /// `all_domains`, records carry the index of the domain they split.
    pub fn write_segment_header(
        &self,
        writer: &mut impl Write,
        all_domains: bool,
    ) -> io::Result<()> {
        match self {
            Self::Tsv => {
                write_tsv_id_header(writer, all_domains)?;
                for name in SEGMENT_NAMES {
                    write!(writer, "\t{name}")?;
                }
                writeln!(writer, "\terror\terror_kind")
            }
            Self::Json => writeln!(writer, "["),
            Self::Jsonl => Ok(()),
        }
    }

    /// Write a single numbered record. With `all_domains`, it carries the index of the domain it
    /// numbers.
    pub fn write_record(
        &self,
        writer: &mut impl Write,
        record: &NumberedRecord,
        index: usize,
        all_domains: bool,
    ) -> io::Result<()> {
        match self {
            Self::Tsv => write_tsv_record(writer, record, all_domains),
            _ => self.write_json_record(writer, &record_to_json(record, all_domains)?, index),
        }
    }

    /// Write a single segmented record, one TSV row per record. With `all_domains`, it carries the
    /// index of the domain it splits.
    pub fn write_segment_record(
        &self,
        writer: &mut impl Write,
        record: &SegmentedRecord,
        index: usize,
        all_domains: bool,
    ) -> io::Result<()> {
        match self {
            Self::Tsv => {
                write_tsv_id(writer, &record.id, record.domain, all_domains)?;
                match &record.result {
                    Ok(segments) => {
                        for (_, residues) in segments.regions() {
                            write!(writer, "\t{residues}")?;
                        }
                        writeln!(writer, "\t\t")
                    }
                    Err(e) => {
                        for _ in SEGMENT_NAMES {
                            write!(writer, "\t")?;
                        }
                        writeln!(writer, "\t{e}\t{}", e.kind())
                    }
                }
            }
            _ => self.write_json_record(writer, &segments_to_json(record, all_domains), index),
        }
    }

    // One JSON record: an element of the JSON array, or a JSONL line
    fn write_json_record(
        &self,
        writer: &mut impl Write,
        json: &serde_json::Value,
        index: usize,
    ) -> io::Result<()> {
        if matches!(self, Self::Json) {
            if index > 0 {
                writeln!(writer, ",")?;
            }
            serde_json::to_writer_pretty(&mut *writer, json).map_err(io::Error::other)
        } else {
            serde_json::to_writer(&mut *writer, json).map_err(io::Error::other)?;
            writeln!(writer)
        }
    }

    /// Write format footer (e.g. JSON array closing)
    pub fn write_footer(&self, writer: &mut impl Write) -> io::Result<()> {
        match self {
            Self::Json => writeln!(writer, "\n]"),
            _ => Ok(()),
        }
    }
}

/// Read input records: auto-detects FASTA file, stdin, or raw sequence string. An argument is a
/// raw sequence only when it isn't an existing path and consists solely of ASCII letters.
pub fn read_input(input: Option<&str>) -> Result<Vec<Record>, String> {
    match input {
        None | Some("-") => {
            let stdin = io::stdin();
            read_auto(BufReader::new(stdin.lock()))
        }
        Some(s) => {
            let path = Path::new(s);
            if path.exists() {
                let file = File::open(path).map_err(|e| format!("cannot open '{}': {}", s, e))?;
                read_auto(BufReader::new(file))
            } else if !s.is_empty() && s.bytes().all(|b| b.is_ascii_alphabetic()) {
                Ok(vec![Record {
                    id: "seq_1".to_string(),
                    sequence: s.to_string(),
                }])
            } else {
                Err(format!(
                    "cannot open '{s}': no such file (a raw sequence may only contain letters)"
                ))
            }
        }
    }
}

/// Read records, auto-detecting FASTA (starts with '>') or raw sequence lines
fn read_auto(reader: impl BufRead) -> Result<Vec<Record>, String> {
    let mut lines = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("read error: {}", e))?;
        let trimmed = line.trim().to_string();
        if !trimmed.is_empty() {
            lines.push(trimmed);
        }
    }
    if lines.is_empty() {
        return Ok(Vec::new());
    }
    if lines[0].starts_with('>') {
        read_fasta(io::Cursor::new(lines.join("\n")))
    } else {
        Ok(lines
            .into_iter()
            .enumerate()
            .map(|(i, seq)| Record {
                id: format!("seq_{}", i + 1),
                sequence: seq,
            })
            .collect())
    }
}

/// Parse FASTA records from a buffered reader
pub fn read_fasta(reader: impl BufRead) -> Result<Vec<Record>, String> {
    let mut records = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("read error: {}", e))?;
        let line = line.trim_end();
        if let Some(header) = line.strip_prefix('>') {
            if !current_id.is_empty() && !current_seq.is_empty() {
                records.push(Record {
                    id: current_id,
                    sequence: current_seq,
                });
                current_seq = String::new();
            }
            current_id = header
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
        } else if !line.is_empty() {
            current_seq.push_str(line);
        }
    }
    if !current_id.is_empty() && !current_seq.is_empty() {
        records.push(Record {
            id: current_id,
            sequence: current_seq,
        });
    }
    Ok(records)
}

/// Write records in TSV long format (one row per position)
pub fn write_tsv(writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
    write_tsv_header(writer, false)?;
    for rec in records {
        write_tsv_record(writer, rec, false)?;
    }
    Ok(())
}

fn write_tsv_header(writer: &mut impl Write, all_domains: bool) -> io::Result<()> {
    write_tsv_id_header(writer, all_domains)?;
    writeln!(
        writer,
        "\tchain\tscheme\tconfidence\tposition\tresidue\terror\terror_kind"
    )
}

// The columns naming a record: its sequence, and its domain when every domain was processed
fn write_tsv_id_header(writer: &mut impl Write, all_domains: bool) -> io::Result<()> {
    write!(writer, "sequence_id")?;
    if all_domains {
        write!(writer, "\tdomain")?;
    }
    Ok(())
}

fn write_tsv_id(
    writer: &mut impl Write,
    id: &str,
    domain: Option<usize>,
    all_domains: bool,
) -> io::Result<()> {
    write!(writer, "{id}")?;
    match (all_domains, domain) {
        (true, Some(domain)) => write!(writer, "\t{domain}"),
        (true, None) => write!(writer, "\t"),
        (false, _) => Ok(()),
    }
}

/// Write a single record in TSV format (without header)
fn write_tsv_record(
    writer: &mut impl Write,
    rec: &NumberedRecord,
    all_domains: bool,
) -> io::Result<()> {
    match &rec.result {
        Ok(result) => {
            for (pos, ch) in result.residues(&rec.sequence).map_err(io::Error::other)? {
                write_tsv_id(writer, &rec.id, rec.domain, all_domains)?;
                writeln!(
                    writer,
                    "\t{}\t{}\t{:.4}\t{}\t{}\t\t",
                    result.chain, result.scheme, result.confidence, pos, ch
                )?;
            }
        }
        Err(e) => {
            write_tsv_id(writer, &rec.id, rec.domain, all_domains)?;
            writeln!(writer, "\t\t\t\t\t\t{e}\t{}", e.kind())?;
        }
    }
    Ok(())
}

/// Write records as a JSON array
pub fn write_json(writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
    let json_records = records
        .iter()
        .map(|rec| record_to_json(rec, false))
        .collect::<io::Result<Vec<_>>>()?;
    serde_json::to_writer_pretty(&mut *writer, &json_records).map_err(io::Error::other)?;
    writeln!(writer)?;
    Ok(())
}

/// Write records as JSON lines (one object per line)
pub fn write_jsonl(writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
    for rec in records {
        let json = record_to_json(rec, false)?;
        serde_json::to_writer(&mut *writer, &json).map_err(io::Error::other)?;
        writeln!(writer)?;
    }
    Ok(())
}

// The keys naming a record: its sequence, and its domain when every domain was processed
fn json_id(
    id: &str,
    domain: Option<usize>,
    all_domains: bool,
) -> serde_json::Map<String, serde_json::Value> {
    let mut record = serde_json::Map::new();
    record.insert("sequence_id".into(), id.into());
    if all_domains {
        record.insert("domain".into(), domain.into());
    }
    record
}

fn segments_to_json(rec: &SegmentedRecord, all_domains: bool) -> serde_json::Value {
    let mut record = json_id(&rec.id, rec.domain, all_domains);
    match &rec.result {
        Ok(segments) => {
            for (name, residues) in segments.regions() {
                record.insert(name.into(), residues.into());
            }
            record.insert("error".into(), serde_json::Value::Null);
            record.insert("error_kind".into(), serde_json::Value::Null);
        }
        Err(e) => {
            for name in SEGMENT_NAMES {
                record.insert(name.into(), serde_json::Value::Null);
            }
            record.insert("error".into(), e.to_string().into());
            record.insert("error_kind".into(), e.kind().into());
        }
    }
    record.into()
}

fn record_to_json(rec: &NumberedRecord, all_domains: bool) -> io::Result<serde_json::Value> {
    let mut record = json_id(&rec.id, rec.domain, all_domains);
    let fields = match &rec.result {
        Ok(result) => {
            let numbering: serde_json::Map<String, serde_json::Value> = result
                .residues(&rec.sequence)
                .map_err(io::Error::other)?
                .map(|(pos, ch)| (pos.to_string(), serde_json::Value::String(ch.to_string())))
                .collect();
            serde_json::json!({
                "chain": result.chain.to_string(),
                "scheme": result.scheme.to_string(),
                "confidence": result.confidence,
                "numbering": numbering,
                "query_start": result.query_start,
                "query_end": result.query_end,
                "error": null,
                "error_kind": null,
            })
        }
        Err(e) => serde_json::json!({
            "chain": null,
            "scheme": null,
            "confidence": null,
            "numbering": null,
            "query_start": null,
            "query_end": null,
            "error": e.to_string(),
            "error_kind": e.kind(),
        }),
    };
    if let serde_json::Value::Object(fields) = fields {
        record.extend(fields);
    }
    Ok(record.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Chain, Position, Scheme};
    use std::io::Cursor;

    const LOW_CONFIDENCE: SequenceError = SequenceError::LowConfidence {
        confidence: 0.1,
        threshold: 0.5,
    };

    fn simple_test_result(positions: Vec<Position>) -> NumberingResult {
        let query_end = positions.len().saturating_sub(1);
        NumberingResult {
            chain: Chain::IGH,
            scheme: Scheme::IMGT,
            positions,
            cons_start: 0,
            cons_end: 0,
            confidence: 1.0,
            query_start: 0,
            query_end,
        }
    }

    #[test]
    fn test_read_fasta_single() {
        let input = b">seq1\nEVQLVES\n";
        let records = read_fasta(Cursor::new(input)).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence, "EVQLVES");
    }

    #[test]
    fn test_read_fasta_multi() {
        let input = b">seq1 some description\nEVQL\nVES\n\n>seq2\nDIQMT\n";
        let records = read_fasta(Cursor::new(input)).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence, "EVQLVES");
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].sequence, "DIQMT");
    }

    #[test]
    fn test_read_fasta_empty() {
        let input = b"";
        let records = read_fasta(Cursor::new(input)).unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn test_write_tsv() {
        let result = simple_test_result(vec![
            Position {
                number: 1,
                insertion: None,
            },
            Position {
                number: 2,
                insertion: None,
            },
        ]);
        let records = vec![NumberedRecord::new(
            "s1".to_string(),
            "EV".to_string(),
            Ok(result),
        )];
        let mut buf = Vec::new();
        write_tsv(&mut buf, &records).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(
            lines[0],
            "sequence_id\tchain\tscheme\tconfidence\tposition\tresidue\terror\terror_kind"
        );
        assert_eq!(lines[1], "s1\tH\tIMGT\t1.0000\t1\tE\t\t");
        assert_eq!(lines[2], "s1\tH\tIMGT\t1.0000\t2\tV\t\t");
    }

    #[test]
    fn test_write_tsv_error() {
        let records = vec![NumberedRecord::new(
            "bad".to_string(),
            "AAAAA".to_string(),
            Err(LOW_CONFIDENCE),
        )];
        let mut buf = Vec::new();
        write_tsv(&mut buf, &records).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(
            lines[0],
            "sequence_id\tchain\tscheme\tconfidence\tposition\tresidue\terror\terror_kind"
        );
        assert_eq!(
            lines[1],
            "bad\t\t\t\t\t\talignment confidence 0.1000 is below min_confidence 0.5000\tlow_confidence"
        );
    }

    #[test]
    fn test_write_jsonl() {
        let result = simple_test_result(vec![Position {
            number: 1,
            insertion: None,
        }]);
        let records = vec![NumberedRecord::new(
            "s1".to_string(),
            "E".to_string(),
            Ok(result),
        )];
        let mut buf = Vec::new();
        write_jsonl(&mut buf, &records).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(output.trim()).unwrap();
        assert_eq!(parsed["sequence_id"], "s1");
        assert_eq!(parsed["numbering"]["1"], "E");
        assert!(parsed["error"].is_null());
        assert!(parsed["error_kind"].is_null());
    }

    #[test]
    fn test_write_jsonl_error() {
        let records = vec![NumberedRecord::new(
            "bad".to_string(),
            "AAAAA".to_string(),
            Err(LOW_CONFIDENCE),
        )];
        let mut buf = Vec::new();
        write_jsonl(&mut buf, &records).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(output.trim()).unwrap();
        assert_eq!(parsed["sequence_id"], "bad");
        assert!(parsed["chain"].is_null());
        assert_eq!(
            parsed["error"],
            "alignment confidence 0.1000 is below min_confidence 0.5000"
        );
        assert_eq!(parsed["error_kind"], "low_confidence");
    }

    #[test]
    fn test_write_json() {
        let result = simple_test_result(vec![Position {
            number: 1,
            insertion: None,
        }]);
        let records = vec![NumberedRecord::new(
            "s1".to_string(),
            "E".to_string(),
            Ok(result),
        )];
        let mut buf = Vec::new();
        write_json(&mut buf, &records).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let parsed: Vec<serde_json::Value> = serde_json::from_str(&output).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0]["sequence_id"], "s1");
        assert!(parsed[0]["error"].is_null());
    }

    #[test]
    fn test_output_format_from_str() {
        assert!(matches!(
            "tsv".parse::<OutputFormat>().unwrap(),
            OutputFormat::Tsv
        ));
        assert!(matches!(
            "JSON".parse::<OutputFormat>().unwrap(),
            OutputFormat::Json
        ));
        assert!(matches!(
            "jsonl".parse::<OutputFormat>().unwrap(),
            OutputFormat::Jsonl
        ));
        assert!("xml".parse::<OutputFormat>().is_err());
    }
}
