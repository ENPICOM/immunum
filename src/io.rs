//! Input parsing and output formatting for sequence records

use crate::annotator::NumberingResult;
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
    pub result: NumberingResult,
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
    /// Write numbered records in this format
    pub fn write(&self, writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
        match self {
            Self::Tsv => write_tsv(writer, records),
            Self::Json => write_json(writer, records),
            Self::Jsonl => write_jsonl(writer, records),
        }
    }
}

/// Read input records: auto-detects FASTA file, stdin, or raw sequence string
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
            } else {
                Ok(vec![Record {
                    id: "seq_1".to_string(),
                    sequence: s.to_string(),
                }])
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
    writeln!(writer, "sequence_id\tchain\tscheme\tposition\tresidue")?;
    for rec in records {
        for (pos, ch) in rec.result.positions.iter().zip(rec.sequence.chars()) {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                rec.id, rec.result.chain, rec.result.scheme, pos, ch
            )?;
        }
    }
    Ok(())
}

/// Write records as a JSON array
pub fn write_json(writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
    let json_records: Vec<serde_json::Value> = records.iter().map(record_to_json).collect();
    serde_json::to_writer_pretty(&mut *writer, &json_records).map_err(io::Error::other)?;
    writeln!(writer)?;
    Ok(())
}

/// Write records as JSON lines (one object per line)
pub fn write_jsonl(writer: &mut impl Write, records: &[NumberedRecord]) -> io::Result<()> {
    for rec in records {
        let json = record_to_json(rec);
        serde_json::to_writer(&mut *writer, &json).map_err(io::Error::other)?;
        writeln!(writer)?;
    }
    Ok(())
}

fn record_to_json(rec: &NumberedRecord) -> serde_json::Value {
    let numbering: serde_json::Map<String, serde_json::Value> = rec
        .result
        .positions
        .iter()
        .zip(rec.sequence.chars())
        .map(|(pos, ch)| (pos.to_string(), serde_json::Value::String(ch.to_string())))
        .collect();

    serde_json::json!({
        "sequence_id": rec.id,
        "chain": rec.result.chain.to_string(),
        "scheme": rec.result.scheme.to_string(),
        "numbering": numbering,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Chain, Position, Scheme};
    use std::io::Cursor;

    fn simple_test_result(positions: Vec<Position>) -> NumberingResult {
        NumberingResult {
            chain: Chain::IGH,
            scheme: Scheme::IMGT,
            positions,
            start: 0,
            end: 0,
            confidence: 1.0,
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
        let records = vec![NumberedRecord {
            id: "s1".to_string(),
            sequence: "EV".to_string(),
            result,
        }];
        let mut buf = Vec::new();
        write_tsv(&mut buf, &records).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines[0], "sequence_id\tchain\tscheme\tposition\tresidue");
        assert_eq!(lines[1], "s1\tH\tIMGT\t1\tE");
        assert_eq!(lines[2], "s1\tH\tIMGT\t2\tV");
    }

    #[test]
    fn test_write_jsonl() {
        let result = simple_test_result(vec![Position {
            number: 1,
            insertion: None,
        }]);
        let records = vec![NumberedRecord {
            id: "s1".to_string(),
            sequence: "E".to_string(),
            result,
        }];
        let mut buf = Vec::new();
        write_jsonl(&mut buf, &records).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(output.trim()).unwrap();
        assert_eq!(parsed["sequence_id"], "s1");
        assert_eq!(parsed["numbering"]["1"], "E");
    }

    #[test]
    fn test_write_json() {
        let result = simple_test_result(vec![Position {
            number: 1,
            insertion: None,
        }]);
        let records = vec![NumberedRecord {
            id: "s1".to_string(),
            sequence: "E".to_string(),
            result,
        }];
        let mut buf = Vec::new();
        write_json(&mut buf, &records).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let parsed: Vec<serde_json::Value> = serde_json::from_str(&output).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0]["sequence_id"], "s1");
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
