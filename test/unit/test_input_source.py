"""Unit tests for the input_source module (D4)."""
from __future__ import annotations

import json

import pytest

from stdface.core.input_source import (
    StanFileSource,
    TOMLSource,
    JSONSource,
    DictSource,
    _normalize_keys,
    _TOML_BARE_ALIASES,
)


class TestStanFileSource:
    """Tests for StanFileSource.load()."""

    def test_parses_keyword_value(self, tmp_path):
        stan = tmp_path / "stan.in"
        stan.write_text("model = hubbard\nlattice = chain\nL = 4\n")
        data = StanFileSource(stan).load()
        assert data == {"model": "hubbard", "lattice": "chain", "l": "4"}

    def test_lowercases_keywords(self, tmp_path):
        stan = tmp_path / "stan.in"
        stan.write_text("Model = Hubbard\nLATTICE = chain\n")
        data = StanFileSource(stan).load()
        assert data["model"] == "Hubbard"  # value kept verbatim
        assert data["lattice"] == "chain"

    def test_skips_comments_and_blanks(self, tmp_path):
        stan = tmp_path / "stan.in"
        stan.write_text("// a comment\n\nmodel = spin\n\n// another\n")
        data = StanFileSource(stan).load()
        assert data == {"model": "spin"}

    def test_missing_equals_raises(self, tmp_path):
        stan = tmp_path / "stan.in"
        stan.write_text("model hubbard\n")
        with pytest.raises(ValueError):
            StanFileSource(stan).load()

    def test_duplicate_keyword_raises(self, tmp_path):
        stan = tmp_path / "stan.in"
        stan.write_text("t = 1.0\nt = 2.0\n")
        with pytest.raises(ValueError):
            StanFileSource(stan).load()

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            StanFileSource(tmp_path / "nope.in").load()


class TestTOMLSource:
    """Tests for TOMLSource.load() (with bare-key aliases)."""

    def test_parses_toml(self, tmp_path):
        toml = tmp_path / "in.toml"
        toml.write_text('model = "hubbard"\nL = 4\nt = 1.0\n')
        data = TOMLSource(toml).load()
        assert data["model"] == "hubbard"
        assert data["L"] == 4
        assert data["t"] == 1.0

    def test_bare_alias_translation(self, tmp_path):
        toml = tmp_path / "in.toml"
        toml.write_text("tp = 1.0\ntpp = 2.0\nj0p = 3.0\n")
        data = TOMLSource(toml).load()
        assert data["t'"] == 1.0
        assert data["t''"] == 2.0
        assert data["j0'"] == 3.0

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            TOMLSource(tmp_path / "nope.toml").load()


class TestJSONSource:
    """Tests for JSONSource.load()."""

    def test_parses_json(self, tmp_path):
        js = tmp_path / "in.json"
        js.write_text(json.dumps({"model": "hubbard", "L": 4}))
        data = JSONSource(js).load()
        assert data == {"model": "hubbard", "L": 4}

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            JSONSource(tmp_path / "nope.json").load()


class TestDictSource:
    """Tests for DictSource.load()."""

    def test_returns_copy(self):
        original = {"model": "hubbard", "L": 4}
        src = DictSource(original)
        data = src.load()
        assert data == original
        data["L"] = 99
        assert original["L"] == 4  # load() returns an independent copy


class TestNormalizeKeys:
    """Tests for the TOML bare-alias table and _normalize_keys."""

    def test_alias_table_examples(self):
        assert _TOML_BARE_ALIASES["tp"] == "t'"
        assert _TOML_BARE_ALIASES["tpp"] == "t''"
        assert _TOML_BARE_ALIASES["j0p"] == "j0'"
        assert _TOML_BARE_ALIASES["v2pp"] == "v2''"

    def test_normalize_case_insensitive(self):
        assert _normalize_keys({"Tp": 1.0}) == {"t'": 1.0}

    def test_non_alias_passthrough(self):
        assert _normalize_keys({"W": 4, "model": "spin"}) == {"W": 4, "model": "spin"}
