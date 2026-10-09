//! The core reports every case in `error_cases.json` with its kind and message; every other
//! interface's tests read the same file.

use immunum::numbering::region_spans;
use immunum::{per_domain, scheme_supports_chain, Annotator};
use serde_json::Value;

fn cases(section: &str) -> Vec<Value> {
    let all: Value = serde_json::from_str(include_str!("error_cases.json")).unwrap();
    all[section].as_array().unwrap().clone()
}

fn str_list(value: &Value) -> Vec<&str> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect()
}

fn annotator(case: &Value) -> Result<Annotator, immunum::Error> {
    Annotator::from_names(
        str_list(&case["chains"]),
        case["scheme"].as_str().unwrap(),
        case["min_confidence"].as_f64().map(|c| c as f32),
    )
}

#[test]
fn setup_errors() {
    for case in cases("setup") {
        let error = annotator(&case).err().expect("setup should fail");
        assert_eq!(error.kind(), case["kind"], "{case}");
        assert_eq!(error.to_string(), case["message"], "{case}");
    }
}

#[test]
fn lookup_errors() {
    for case in cases("lookups") {
        let (scheme, chain) = (
            case["scheme"].as_str().unwrap(),
            case["chain"].as_str().unwrap(),
        );

        let error = region_spans(scheme, chain).unwrap_err();
        assert_eq!(error.kind(), case["regions_for"]["kind"], "{case}");
        assert_eq!(error.to_string(), case["regions_for"]["message"], "{case}");

        let expected = &case["scheme_supports_chain"];
        match scheme_supports_chain(scheme, chain) {
            Ok(supported) => assert_eq!(Value::Bool(supported), *expected, "{case}"),
            Err(error) => {
                assert_eq!(error.kind(), expected["kind"], "{case}");
                assert_eq!(error.to_string(), expected["message"], "{case}");
            }
        }
    }
}

#[test]
fn sequence_errors() {
    for case in cases("sequences") {
        let annotator = annotator(&case).unwrap();
        let sequence = case["sequence"].as_str().unwrap();

        let number = annotator.number(sequence).unwrap_err();
        let segment = annotator.segment(sequence).unwrap_err();
        for error in [number, segment] {
            assert_eq!(error.kind(), case["kind"], "{case}");
            assert_eq!(error.to_string(), case["message"], "{case}");
        }

        for domains in [
            per_domain(annotator.number_domains(sequence))
                .into_iter()
                .map(|r| r.map(|_| ()))
                .collect::<Vec<_>>(),
            per_domain(annotator.segment_domains(sequence))
                .into_iter()
                .map(|r| r.map(|_| ()))
                .collect(),
        ] {
            match case["domains"].as_str().unwrap() {
                "empty" => assert!(domains.is_empty(), "{case}"),
                "error" => {
                    let [Err(error)] = &domains[..] else {
                        panic!("expected one error for {case}")
                    };
                    assert_eq!(error.kind(), case["kind"], "{case}");
                    assert_eq!(error.to_string(), case["message"], "{case}");
                }
                other => panic!("unknown domains expectation {other}"),
            }
        }
    }
}
