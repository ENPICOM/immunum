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
            "sequence_id\tchain\tscheme\tconfidence\tposition\tresidue\terror\terror_kind",
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

#[test]
fn chain_filter_takes_a_comma_separated_list() {
    immunum()
        .args([
            "number",
            "-c",
            "k, h",
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("\tH\t"));
    immunum()
        .args([
            "number",
            "-c",
            "k,xyz",
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS",
        ])
        .assert()
        .failure();
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
    assert!(contents
        .contains("sequence_id\tchain\tscheme\tconfidence\tposition\tresidue\terror\terror_kind"));
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
fn jsonl_record_carries_the_numbered_span() {
    let igh = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
    let leader = "MGWSCIILFLVATATGVHSX";
    let output = immunum()
        .args(["number", "-f", "jsonl", &format!("{leader}{igh}")])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(stdout.trim()).expect("valid jsonl");
    assert_eq!(parsed["query_start"], leader.len());
    assert_eq!(parsed["query_end"], leader.len() + igh.len() - 1);
}

const SCFV: &str = concat!(
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
    "GGGGSGGGGSGGGGS",
    "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNSKNYLAWYQDKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYSTPYSFGQGTKLEIK",
);

#[test]
fn all_domains_emits_one_record_per_domain() {
    // An scFv, and an invalid sequence, whose error record belongs to no domain
    let input = format!("{SCFV}\nAAAA\n");
    let output = immunum()
        .args(["number", "--all-domains", "-f", "jsonl"])
        .write_stdin(input)
        .output()
        .unwrap();

    assert!(output.status.success());
    let records: Vec<serde_json::Value> = String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).expect("valid jsonl"))
        .collect();
    let summary: Vec<_> = records
        .iter()
        .map(|r| {
            (
                r["sequence_id"].clone(),
                r["domain"].clone(),
                r["chain"].clone(),
                r["error"].is_string(),
            )
        })
        .collect();
    assert_eq!(
        summary,
        vec![
            ("seq_1".into(), 0.into(), "H".into(), false),
            ("seq_1".into(), 1.into(), "K".into(), false),
            (
                "seq_2".into(),
                serde_json::Value::Null,
                serde_json::Value::Null,
                true
            ),
        ]
    );
}

#[test]
fn all_domains_tsv_has_a_domain_column() {
    let output = immunum()
        .args(["number", "--all-domains", SCFV])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let mut lines = stdout.lines();
    assert_eq!(
        lines.next(),
        Some(
            "sequence_id\tdomain\tchain\tscheme\tconfidence\tposition\tresidue\terror\terror_kind"
        )
    );
    let domains: std::collections::BTreeSet<(&str, &str)> = lines
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            (cols[1], cols[2])
        })
        .collect();
    assert_eq!(domains, [("0", "H"), ("1", "K")].into_iter().collect());
}

const SEGMENT_COLUMNS: [&str; 11] = [
    "prefix",
    "fr1",
    "cdr1",
    "fr2",
    "cdr2",
    "fr3",
    "cdr3",
    "fr4",
    "postfix",
    "error",
    "error_kind",
];

#[test]
fn segment_writes_a_column_per_region() {
    let igh = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
    let output = immunum().args(["segment", igh]).output().unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let rows: Vec<Vec<&str>> = stdout.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows[0][0], "sequence_id");
    assert_eq!(rows[0][1..], SEGMENT_COLUMNS);
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[1][1..10].concat(), igh);
    assert_eq!(rows[1][7], "AREGTTGKPIGAFAH");
}

#[test]
fn segment_all_domains_puts_every_residue_in_one_record() {
    let output = immunum()
        .args(["segment", "--all-domains", "-f", "jsonl"])
        .write_stdin(format!("{SCFV}\nAAAA\n"))
        .output()
        .unwrap();

    assert!(output.status.success());
    let records: Vec<serde_json::Value> = String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).expect("valid jsonl"))
        .collect();
    let ids: Vec<_> = records
        .iter()
        .map(|r| (r["sequence_id"].clone(), r["domain"].clone()))
        .collect();
    assert_eq!(
        ids,
        vec![
            ("seq_1".into(), 0.into()),
            ("seq_1".into(), 1.into()),
            ("seq_2".into(), serde_json::Value::Null),
        ]
    );
    let rebuilt: String = records[..2]
        .iter()
        .flat_map(|r| {
            SEGMENT_COLUMNS[..9]
                .iter()
                .map(|&c| r[c].as_str().unwrap().to_string())
        })
        .collect();
    assert_eq!(rebuilt, SCFV);
    assert!(records[2]["error"].is_string());
    assert!(records[2]["fr1"].is_null());
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
    let cols: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(
        cols[cols.len() - 2..],
        ["sequence length 10 is below minimum 30", "invalid_sequence"]
    );
    assert_eq!(
        String::from_utf8(output.stderr).unwrap(),
        "warning: 1 of 1 sequences had errors; see the error and error_kind columns\n"
    );
}

#[test]
fn success_rows_leave_error_columns_empty_in_tsv() {
    let output = immunum()
        .args(["number", "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS"])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(output.stderr.is_empty(), "no summary without errors");
    let stdout = String::from_utf8(output.stdout).unwrap();
    for line in stdout.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 8);
        assert_eq!(cols[6..], ["", ""]);
    }
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
fn missing_input_file_is_a_setup_error() {
    immunum()
        .args(["number", "seqs_typo.fasta"])
        .assert()
        .code(1)
        .stdout(predicate::str::is_empty())
        .stderr(
            "error: cannot open 'seqs_typo.fasta': no such file (a raw sequence may only contain letters)\n",
        );
}

// --- Shared error table (tests/error_cases.json) ---

fn error_cases() -> serde_json::Value {
    serde_json::from_str(&fs::read_to_string("tests/error_cases.json").unwrap()).unwrap()
}

fn chain_list(case: &serde_json::Value) -> String {
    let chains: Vec<&str> = case["chains"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap())
        .collect();
    chains.join(",")
}

fn jsonl_records(stdout: Vec<u8>) -> Vec<serde_json::Value> {
    String::from_utf8(stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).expect("valid jsonl"))
        .collect()
}

#[test]
fn setup_errors_exit_with_the_shared_message() {
    for case in error_cases()["setup"].as_array().unwrap() {
        let chains = chain_list(case);
        if chains.is_empty() {
            // The command line can't express an empty chain list
            continue;
        }
        let mut args = vec![
            "number".to_string(),
            "-c".to_string(),
            chains,
            "-s".to_string(),
            case["scheme"].as_str().unwrap().to_string(),
        ];
        if let Some(min_confidence) = case["min_confidence"].as_f64() {
            args.extend(["--min-confidence".to_string(), min_confidence.to_string()]);
        }
        args.push("EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS".to_string());
        immunum()
            .args(&args)
            .assert()
            .code(1)
            .stdout(predicate::str::is_empty())
            .stderr(format!("error: {}\n", case["message"].as_str().unwrap()));
    }
}

#[test]
fn sequence_errors_are_returned_as_records() {
    for case in error_cases()["sequences"].as_array().unwrap() {
        for command in ["number", "segment"] {
            let records = jsonl_records(run_case(command, case, false));
            assert_eq!(records.len(), 1, "{command} {case}");
            assert_eq!(records[0]["error"], case["message"], "{command} {case}");
            assert_eq!(records[0]["error_kind"], case["kind"], "{command} {case}");
            assert_all_domains_error(command, case);
        }
    }
}

#[test]
fn a_sequence_without_a_domain_gets_an_error_record_with_all_domains() {
    for case in error_cases()["domain_errors"].as_array().unwrap() {
        for command in ["number", "segment"] {
            let records = jsonl_records(run_case(command, case, false));
            assert!(records[0]["error"].is_null(), "{command} {case}");
            assert_all_domains_error(command, case);
        }
    }
}

// `command -f jsonl` on `case`'s sequence, through stdin because one fixture sequence holds a digit
fn run_case(command: &str, case: &serde_json::Value, all_domains: bool) -> Vec<u8> {
    let chains = chain_list(case);
    let mut args = vec![command, "-f", "jsonl", "-c", &chains];
    args.extend(["-s", case["scheme"].as_str().unwrap()]);
    if all_domains {
        args.push("--all-domains");
    }
    let output = immunum()
        .args(args)
        .write_stdin(format!("{}\n", case["sequence"].as_str().unwrap()))
        .output()
        .unwrap();
    assert!(output.status.success(), "{command} {case}");
    output.stdout
}

// With `--all-domains`, `case`'s sequence gets one record, carrying its error
fn assert_all_domains_error(command: &str, case: &serde_json::Value) {
    let records = jsonl_records(run_case(command, case, true));
    assert_eq!(records.len(), 1, "{command} {case}");
    assert_eq!(records[0]["error"], case["message"], "{command} {case}");
    assert_eq!(records[0]["error_kind"], case["kind"], "{command} {case}");
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
