"""The integration-test comparator (test/common/almost_diff.py)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "common" / "almost_diff.py"


def _compare(tmp_path, ref, out):
    (tmp_path / "ref").write_text(ref)
    (tmp_path / "out").write_text(out)
    return subprocess.run(
        [sys.executable, str(_SCRIPT), str(tmp_path / "ref"), str(tmp_path / "out")],
        capture_output=True, text=True).returncode


@pytest.mark.parametrize("ref, out", [
    ("a 1\nb 2\n", "a 1\nb 2\n"),
    ("x 1.000000000000000e+00\n", "x 1.000000000000001e+00\n"),   # float drift
    ("x 0.000000000000000\n", "x -0.000000000000000\n"),           # signed zero
    ("0  0      0  1\n", "0      0  0      1\n"),                  # spacing
    ("x nan\n", "x -nan\n"),
    ("x 07\n", "x 7\n"),
], ids=["equal", "drift", "signed_zero", "spacing", "nan", "int_form"])
def test_equivalent_files_pass(tmp_path, ref, out):
    assert _compare(tmp_path, ref, out) == 0


@pytest.mark.parametrize("ref, out", [
    ("x 1.0\n", "x 1.1\n"),
    ("x 1\n", "y 1\n"),
    ("a 1\nb 2\nc 3\n", "a 1\nb 2\n"),          # output is truncated
    ("a 1\nb 2\n", "a 1\nb 2\nc 3\n"),          # output has extra lines
    ("a 1\nb 2\n", ""),                         # output is empty
    ("a 1 2\n", "a 1\n"),                       # a token is missing
    ("a 1\n", "a 1 2\n"),                       # an extra token
], ids=["value", "keyword", "truncated", "extra_lines", "empty", "missing_token", "extra_token"])
def test_different_files_fail(tmp_path, ref, out):
    assert _compare(tmp_path, ref, out) != 0


def test_missing_argument_is_an_error():
    assert subprocess.run([sys.executable, str(_SCRIPT), "only_one"],
                          capture_output=True).returncode != 0
