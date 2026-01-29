"""Unit tests for helper functions extracted from stdface_main.

Tests for ``_parse_input_file`` and ``_resolve_model_and_method``.
"""
from __future__ import annotations

import os

import pytest

from stdface.core.stdface_vals import StdIntList, ModelType, SolverType, MethodType, NaN_d
from stdface.core.stdface_main import _parse_input_file, _resolve_model_and_method


class TestParseInputFile:
    """Tests for _parse_input_file."""

    def test_parses_simple_keywords(self, tmp_path):
        """Test that basic keyword = value lines are parsed."""
        infile = tmp_path / "stan.in"
        infile.write_text("model = hubbard\nlattice = chain\n")
        StdI = StdIntList()
        _parse_input_file(str(infile), StdI, SolverType.HPhi)
        assert StdI.model == "hubbard"
        assert StdI.lattice == "chain"

    def test_skips_comment_lines(self, tmp_path):
        """Test that // comment lines are skipped."""
        infile = tmp_path / "stan.in"
        infile.write_text("// this is a comment\nmodel = spin\nlattice = square\n")
        StdI = StdIntList()
        _parse_input_file(str(infile), StdI, SolverType.HPhi)
        assert StdI.model == "spin"

    def test_skips_blank_lines(self, tmp_path):
        """Test that blank lines are skipped."""
        infile = tmp_path / "stan.in"
        infile.write_text("\nmodel = kondo\n\nlattice = chain\n")
        StdI = StdIntList()
        _parse_input_file(str(infile), StdI, SolverType.HPhi)
        assert StdI.model == "kondo"

    def test_exits_on_missing_file(self, tmp_path):
        """Test that missing file causes SystemExit."""
        StdI = StdIntList()
        with pytest.raises(SystemExit):
            _parse_input_file(str(tmp_path / "nonexistent.in"), StdI, SolverType.HPhi)

    def test_exits_on_missing_equals(self, tmp_path):
        """Test that a line without '=' causes SystemExit."""
        infile = tmp_path / "stan.in"
        infile.write_text("model hubbard\n")
        StdI = StdIntList()
        with pytest.raises(SystemExit):
            _parse_input_file(str(infile), StdI, SolverType.HPhi)

    def test_exits_on_unknown_keyword(self, tmp_path):
        """Test that an unrecognised keyword causes SystemExit."""
        infile = tmp_path / "stan.in"
        infile.write_text("totally_fake_keyword = 42\n")
        StdI = StdIntList()
        with pytest.raises(SystemExit):
            _parse_input_file(str(infile), StdI, SolverType.HPhi)

    def test_case_insensitive_keywords(self, tmp_path):
        """Test that keywords are case-insensitive."""
        infile = tmp_path / "stan.in"
        infile.write_text("Model = hubbard\nLATTICE = chain\n")
        StdI = StdIntList()
        _parse_input_file(str(infile), StdI, SolverType.HPhi)
        assert StdI.model == "hubbard"
        assert StdI.lattice == "chain"


class TestResolveModelAndMethod:
    """Tests for _resolve_model_and_method."""

    def test_hubbard_canonical(self):
        """Test hubbard resolves to ModelType.HUBBARD with lGC=0."""
        StdI = StdIntList()
        StdI.model = "hubbard"
        StdI.lattice = "chain"
        _resolve_model_and_method(StdI, SolverType.HPhi)
        assert StdI.model == ModelType.HUBBARD
        assert StdI.lGC == 0
        assert StdI.lBoost == 0

    def test_spingc(self):
        """Test spingc resolves to ModelType.SPIN with lGC=1."""
        StdI = StdIntList()
        StdI.model = "spingc"
        StdI.lattice = "chain"
        _resolve_model_and_method(StdI, SolverType.mVMC)
        assert StdI.model == ModelType.SPIN
        assert StdI.lGC == 1
        assert StdI.lBoost == 0

    def test_spingcboost_hphi_only(self):
        """Test spingcboost resolves only for HPhi solver."""
        StdI = StdIntList()
        StdI.model = "spingcboost"
        StdI.lattice = "chain"
        _resolve_model_and_method(StdI, SolverType.HPhi)
        assert StdI.model == ModelType.SPIN
        assert StdI.lGC == 1
        assert StdI.lBoost == 1

    def test_spingcboost_non_hphi_exits(self):
        """Test spingcboost exits for non-HPhi solver."""
        StdI = StdIntList()
        StdI.model = "spingcboost"
        StdI.lattice = "chain"
        with pytest.raises(SystemExit):
            _resolve_model_and_method(StdI, SolverType.mVMC)

    def test_unknown_model_exits(self):
        """Test unknown model name causes SystemExit."""
        StdI = StdIntList()
        StdI.model = "nosuchmodel"
        StdI.lattice = "chain"
        with pytest.raises(SystemExit):
            _resolve_model_and_method(StdI, SolverType.HPhi)

    def test_method_alias_direct(self):
        """Test that method 'direct' is normalised to 'fulldiag' for HPhi."""
        StdI = StdIntList()
        StdI.model = "hubbard"
        StdI.lattice = "chain"
        StdI.method = "direct"
        _resolve_model_and_method(StdI, SolverType.HPhi)
        assert StdI.method == MethodType.FULLDIAG

    def test_method_not_changed_for_non_hphi(self):
        """Test that method aliases are not applied for non-HPhi solvers."""
        StdI = StdIntList()
        StdI.model = "hubbard"
        StdI.lattice = "chain"
        StdI.method = "direct"
        _resolve_model_and_method(StdI, SolverType.mVMC)
        # Method should remain unchanged for non-HPhi
        assert StdI.method == "direct"

    def test_kondo_gc(self):
        """Test kondogc resolves correctly."""
        StdI = StdIntList()
        StdI.model = "kondogc"
        StdI.lattice = "chain"
        _resolve_model_and_method(StdI, SolverType.HPhi)
        assert StdI.model == ModelType.KONDO
        assert StdI.lGC == 1
