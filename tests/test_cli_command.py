"""The `immunum` command the Python package installs, as `pip install immunum` and `uvx immunum`
run it."""

import json
import shutil
import signal
import subprocess
import sys
import sysconfig
import time
from pathlib import Path

import pytest

IMMUNUM = shutil.which("immunum", path=sysconfig.get_path("scripts"))
FASTA = Path(__file__).parent.parent / "fixtures" / "ig.fasta"
IGH_SEQ = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"


def immunum(*args: str, stdin: str | None = None) -> subprocess.CompletedProcess[str]:
    assert IMMUNUM is not None, "the immunum script isn't installed in this environment"
    return subprocess.run(
        [IMMUNUM, *args], input=stdin, capture_output=True, text=True, timeout=60
    )


def test_segments_a_sequence():
    result = immunum("segment", "-f", "jsonl", IGH_SEQ)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["cdr3"] == "AREGTTGKPIGAFAH"


def test_reads_fasta_from_stdin():
    result = immunum("number", "-f", "json", stdin=FASTA.read_text())
    assert result.returncode == 0, result.stderr
    assert [r["sequence_id"] for r in json.loads(result.stdout)] == [
        "4qo1_A|Heavy|A",
        "4olu_L|Light|L",
        "4ktd_L|Light|L",
    ]


def test_exit_status_is_the_cli_s():
    setup_error = immunum("number", "-c", "IGX", IGH_SEQ)
    assert setup_error.returncode == 1
    assert setup_error.stderr.startswith("error: unknown chain 'IGX'")

    assert immunum("number", "--no-such-flag", IGH_SEQ).returncode == 2

    help = immunum("--help")
    assert help.returncode == 0
    assert "segment" in help.stdout


@pytest.mark.skipif(sys.platform == "win32", reason="sends SIGINT")
def test_ctrl_c_stops_a_run():
    # Without input on stdin the command waits for it inside Rust, where Python's own Ctrl-C
    # handling never runs: a run that ignored SIGINT would wait until the timeout
    assert IMMUNUM is not None
    process = subprocess.Popen(
        [IMMUNUM, "number"],
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        time.sleep(1)
        process.send_signal(signal.SIGINT)
        assert process.wait(timeout=10) != 0
    finally:
        process.kill()
        process.wait()
