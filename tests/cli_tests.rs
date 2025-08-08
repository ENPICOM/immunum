use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;

#[test]
fn cli_basic_file_input() {
    let mut cmd = Command::cargo_bin("immunum-cli").expect("binary exists");
    cmd.args(["-i", "fixtures/test.fasta", "-s", "imgt", "-c", "igh"]);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("# "));
}

#[test]
fn cli_all_chains_flag() {
    let mut cmd = Command::cargo_bin("immunum-cli").expect("binary exists");
    cmd.args([
        "-i",
        "fixtures/test.fasta",
        "-s",
        "imgt",
        "-c",
        "igh",
        "--all-chains",
    ]);
    cmd.assert().success();
}

#[test]
fn cli_parallel_flag_on_file() {
    let mut cmd = Command::cargo_bin("immunum-cli").expect("binary exists");
    cmd.args([
        "-i",
        "fixtures/test.fasta",
        "-s",
        "imgt",
        "-c",
        "igh",
        "--parallel",
    ]);
    cmd.assert().success();
}

#[test]
fn cli_nonexistent_file_fails() {
    let mut cmd = Command::cargo_bin("immunum-cli").expect("binary exists");
    cmd.args(["-i", "fixtures/does_not_exist.fasta", "-s", "imgt", "-c", "igh"]);
    cmd.assert().failure();
}

#[test]
fn cli_direct_sequence_input() {
    // A realistic heavy chain fragment (letters only)
    let seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR";
    let mut cmd = Command::cargo_bin("immunum-cli").expect("binary exists");
    cmd.args(["-i", seq, "-s", "imgt", "-c", "igh"]);
    cmd.assert().success();
}

