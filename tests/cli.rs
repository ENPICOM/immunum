#![cfg(feature = "cli")]

use assert_cmd::cargo;
use predicates::prelude::*;
use std::fs;

fn immunum() -> assert_cmd::Command {
    cargo::cargo_bin_cmd!("immunum")
}

// --- Basic input modes ---

#[test]
fn raw_sequence_argument() {
    immunum()
        .args(["number", "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS"])
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "sequence_id\tchain\tscheme\tconfidence\tposition\tresidue\terror",
        ));
}

#[test]
fn fasta_file_input() {
    immunum()
        .args(["number", "fixtures/ig.fasta"])
        .assert()
        .success()
        .stdout(predicate::str::contains("4qo1_A|Heavy|A"));
}

#[test]
fn stdin_fasta_pipe() {
    let fasta = fs::read("fixtures/ig.fasta").unwrap();
    immunum()
        .args(["number"])
        .write_stdin(fasta)
        .assert()
        .success()
        .stdout(predicate::str::contains("4qo1_A|Heavy|A"));
}

#[test]
fn stdin_dash_arg() {
    let fasta = fs::read("fixtures/ig.fasta").unwrap();
    immunum()
        .args(["number", "-"])
        .write_stdin(fasta)
        .assert()
        .success()
        .stdout(predicate::str::contains("4qo1_A|Heavy|A"));
}

#[test]
fn stdin_raw_sequence() {
    immunum()
        .args(["number"])
        .write_stdin("EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS\n")
        .assert()
        .success()
        .stdout(predicate::str::contains("seq_1"));
}

// --- Output formats ---

#[test]
fn output_tsv_default() {
    immunum()
        .args(["number", "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS"])
        .assert()
        .success()
        .stdout(predicate::str::contains("\t"));
}

#[test]
fn output_json() {
    immunum()
        .args([
            "number",
            "-f",
            "json",
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS",
        ])
        .assert()
        .success()
        .stdout(predicate::str::starts_with("["));
}

#[test]
fn output_jsonl() {
    immunum()
        .args([
            "number",
            "-f",
            "jsonl",
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("\"sequence_id\""));
}

// --- Options ---

#[test]
fn scheme_kabat() {
    immunum()
        .args([
            "number",
            "-s",
            "kabat",
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Kabat"));
}

#[test]
fn chain_filter_tcr() {
    immunum()
        .args(["number", "-c", "tcr", "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS"])
        .assert()
        .success();
}

#[test]
fn chain_aliases_case_insensitive() {
    // All these should be equivalent
    for alias in ["h", "heavy", "igh", "H", "Heavy", "IGH"] {
        immunum()
            .args(["number", "-c", alias, "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS"])
            .assert()
            .success()
            .stdout(predicate::str::contains("\tH\t"));
    }
}

// --- Output to file ---

#[test]
fn output_to_file() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("output.tsv");

    immunum()
        .args([
            "number",
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS",
            out_path.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::is_empty());

    let contents = fs::read_to_string(&out_path).unwrap();
    assert!(contents.contains("sequence_id\tchain\tscheme\tconfidence\tposition\tresidue\terror"));
}

// --- TSV piping via stdin ---

#[test]
fn stdin_multiple_raw_sequences() {
    immunum()
        .args(["number"])
        .write_stdin("EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS\nDIQMTQSPSSLSASVGDRVTITC\n")
        .assert()
        .success()
        .stdout(predicate::str::contains("seq_1").and(predicate::str::contains("seq_2")));
}

// --- Error field in output ---

#[test]
fn invalid_sequence_emits_error_record_jsonl() {
    // A garbage sequence should appear as an error record, not abort the batch
    let output = immunum()
        .args(["number", "-f", "jsonl"])
        .write_stdin("AAAAAAAAAA\n")
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(stdout.trim()).expect("valid jsonl");
    assert!(
        parsed["error"].is_string(),
        "error field should be a string"
    );
    assert!(parsed["chain"].is_null(), "chain should be null on error");
    assert!(
        parsed["numbering"].is_null(),
        "numbering should be null on error"
    );
}

#[test]
fn valid_sequence_has_null_error_jsonl() {
    let output = immunum()
        .args([
            "number",
            "-f",
            "jsonl",
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS",
        ])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(stdout.trim()).expect("valid jsonl");
    assert!(parsed["error"].is_null(), "error should be null on success");
    assert!(
        parsed["chain"].is_string(),
        "chain should be set on success"
    );
}

#[test]
fn mixed_batch_always_emits_one_record_per_input() {
    // Two sequences: one valid IGH, one garbage
    let input = "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS\nAAAAAAAAAAAAAAAAA\n";
    let output = immunum()
        .args(["number", "-f", "jsonl"])
        .write_stdin(input)
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.trim().lines().collect();
    assert_eq!(lines.len(), 2, "one output per input");

    let first: serde_json::Value = serde_json::from_str(lines[0]).unwrap();
    let second: serde_json::Value = serde_json::from_str(lines[1]).unwrap();
    assert!(first["error"].is_null());
    assert!(second["error"].is_string());
}

#[test]
fn error_record_appears_in_tsv() {
    let output = immunum()
        .args(["number", "-f", "tsv"])
        .write_stdin("AAAAAAAAAA\n")
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.trim().lines().collect();
    assert_eq!(lines.len(), 2); // header + one error row
    assert!(
        lines[1].contains("seq_1"),
        "error row should have sequence id"
    );
    // last column (error) should be non-empty
    let cols: Vec<&str> = lines[1].split('\t').collect();
    assert!(
        !cols.last().unwrap().is_empty(),
        "error column should be non-empty"
    );
}

// --- Error cases ---

#[test]
fn no_input_on_terminal_shows_error() {
    // Without stdin and no argument, should fail
    // (assert_cmd doesn't set is_terminal, so stdin will be piped — send empty)
    immunum()
        .args(["number"])
        .write_stdin("")
        .assert()
        .success(); // empty input produces no output, no error
}

#[test]
fn invalid_format_shows_error() {
    immunum()
        .args(["number", "-f", "xml", "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("unknown format"));
}

#[test]
fn invalid_scheme_shows_error() {
    immunum()
        .args([
            "number",
            "-s",
            "chothia",
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS",
        ])
        .assert()
        .failure();
}

// --- JSON output is valid ---

#[test]
fn json_output_is_valid_json() {
    let output = immunum()
        .args(["number", "-f", "json", "fixtures/ig.fasta"])
        .output()
        .unwrap();

    assert!(output.status.success());
    let parsed: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("output should be valid JSON");
    assert!(parsed.is_array());
    assert_eq!(parsed.as_array().unwrap().len(), 3);
}

#[test]
fn jsonl_output_has_one_object_per_line() {
    let output = immunum()
        .args(["number", "-f", "jsonl", "fixtures/ig.fasta"])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = stdout.trim().lines().collect();
    assert_eq!(lines.len(), 3);
    for line in &lines {
        let parsed: serde_json::Value =
            serde_json::from_str(line).expect("each line should be valid JSON");
        assert!(parsed.is_object());
    }
}
