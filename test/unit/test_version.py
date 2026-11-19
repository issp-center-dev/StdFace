"""Unit tests for the Python package version.

The Python package is versioned independently of the C generator; its
version is written once, in ``stdface.core.version.__version__``.
"""
from __future__ import annotations

import pathlib
import sys

import pytest

import stdface
import stdface.core.version as version
from stdface.__main__ import main

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Python 3.10
    import tomli as tomllib

_PYPROJECT = pathlib.Path(__file__).resolve().parents[2] / "python" / "pyproject.toml"


class TestVersionString:
    def test_is_a_valid_pep440_version(self):
        packaging_version = pytest.importorskip("packaging.version")
        packaging_version.Version(version.__version__)

    def test_package_exposes_the_same_value(self):
        assert stdface.__version__ == version.__version__


class TestPrintVersion:
    def test_prints_the_python_label_and_version(self, capsys):
        version.print_version()
        assert capsys.readouterr().out.strip() == (
            f"StdFace(python) version {version.__version__}")

    def test_follows_the_single_constant(self, capsys, monkeypatch):
        monkeypatch.setattr(version, "__version__", "9.8.7")
        version.print_version()
        assert capsys.readouterr().out.strip() == "StdFace(python) version 9.8.7"

    def test_command_line_flag(self, capsys):
        assert main(["-v"]) == 0
        assert capsys.readouterr().out.strip() == (
            f"StdFace(python) version {version.__version__}")


class TestSingleSource:
    """pyproject.toml must not carry its own copy of the version."""

    def test_metadata_reads_the_version_from_the_module(self):
        with open(_PYPROJECT, "rb") as fp:
            data = tomllib.load(fp)
        assert "version" not in data["project"]
        assert "version" in data["project"]["dynamic"]
        assert data["tool"]["setuptools"]["dynamic"]["version"] == {
            "attr": "stdface.core.version.__version__"}
