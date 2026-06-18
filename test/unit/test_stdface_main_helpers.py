"""Unit tests for helper functions extracted from stdface_main.

Tests for ``_parse_input_file`` and ``_resolve_model_and_method``.
"""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from stdface.core.stdface_vals import (
    StdIntList,
    ModelType,
    SolverType,
    MethodType,
)
from stdface.core.stdface_main import (
    BOOST_DISPATCH,
    _build_lattice_and_boost,
    _parse_input_file,
    _parse_solver_keyword_via_plugin,
    _resolve_model_and_method,
    stdface_main,
)


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
        """Test that a missing file raises FileNotFoundError."""
        StdI = StdIntList()
        with pytest.raises(FileNotFoundError):
            _parse_input_file(str(tmp_path / "nonexistent.in"), StdI, SolverType.HPhi)

    def test_exits_on_missing_equals(self, tmp_path):
        """Test that a line without '=' raises ValueError."""
        infile = tmp_path / "stan.in"
        infile.write_text("model hubbard\n")
        StdI = StdIntList()
        with pytest.raises(ValueError):
            _parse_input_file(str(infile), StdI, SolverType.HPhi)

    def test_exits_on_unknown_keyword(self, tmp_path):
        """Test that an unrecognised keyword raises ValueError."""
        infile = tmp_path / "stan.in"
        infile.write_text("totally_fake_keyword = 42\n")
        StdI = StdIntList()
        with pytest.raises(ValueError):
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
        with pytest.raises(ValueError):
            _resolve_model_and_method(StdI, SolverType.mVMC)

    def test_unknown_model_exits(self):
        """Test unknown model name causes ValueError."""
        StdI = StdIntList()
        StdI.model = "nosuchmodel"
        StdI.lattice = "chain"
        with pytest.raises(ValueError):
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

    def test_time_evolution_calls_vector_potential(self):
        """HPhi + timeevolution triggers vector_potential (Boost / gauge prep)."""
        StdI = StdIntList()
        StdI.model = "hubbard"
        StdI.lattice = "chain"
        StdI.method = MethodType.TIME_EVOLUTION
        with patch("stdface.core.stdface_main._vector_potential") as vp:
            _resolve_model_and_method(StdI, SolverType.HPhi)
            vp.assert_called_once_with(StdI)


class TestParseSolverKeywordViaPlugin:
    """Tests for _parse_solver_keyword_via_plugin plugin / legacy fallback."""

    def test_fallback_when_solver_plugin_missing(self):
        """KeyError from get_plugin falls back to parse_solver_keyword."""
        StdI = StdIntList()
        result = _parse_solver_keyword_via_plugin(
            "lanczos_max", "100", StdI, "__no_such_solver__"
        )
        assert result is False


class TestStdfaceMainIntegration:
    """Light integration test for stdface_main (parse → lattice → plugin.write)."""

    def test_minimal_hubbard_chain_hphi_writes_namelist(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        infile = tmp_path / "stan.in"
        infile.write_text(
            "model = hubbard\n"
            "lattice = chain\n"
            "L = 4\n"
            "nelec = 4\n"
            "method = lanczos\n"
        )
        stdface_main(str(infile), SolverType.HPhi)
        assert (tmp_path / "namelist.def").is_file()

    def test_cdatafilehead_set_skips_default_branch(self, tmp_path, monkeypatch):
        """When CDataFileHead is set, the default ``zvo`` branch is not taken."""
        monkeypatch.chdir(tmp_path)
        infile = tmp_path / "stan.in"
        infile.write_text(
            "model = hubbard\n"
            "lattice = chain\n"
            "L = 4\n"
            "nelec = 4\n"
            "method = lanczos\n"
            "cdatafilehead = myrun\n"
        )
        stdface_main(str(infile), SolverType.HPhi)
        nml = (tmp_path / "namelist.def").read_text()
        assert "myrun" in nml


class TestBuildLatticeAndBoost:
    """Tests for ``_build_lattice_and_boost`` error and plugin edge paths."""

    def test_unknown_lattice_exits(self):
        StdI = StdIntList()
        StdI.model = "hubbard"
        StdI.lattice = "__not_a_registered_lattice__"
        with pytest.raises(ValueError):
            _build_lattice_and_boost(StdI, SolverType.HPhi)

    def test_unknown_solver_skips_post_lattice(self):
        """``get_plugin`` KeyError is swallowed (no post_lattice hook)."""
        StdI = StdIntList()
        StdI.model = ModelType.HUBBARD
        StdI.lattice = "chain"
        lattice_plugin = MagicMock()
        with patch(
            "stdface.core.stdface_main._get_lattice",
            return_value=lattice_plugin,
        ):
            _build_lattice_and_boost(StdI, "__no_solver_plugin__")
        lattice_plugin.setup.assert_called_once_with(StdI)


class TestBoostDispatchItemsDefensive:
    """Cover ``except KeyError: pass`` inside ``BOOST_DISPATCH.items()``."""

    def test_items_skips_alias_when_get_lattice_raises(self):
        import stdface.core.stdface_main as sm

        real = sm._get_lattice

        def _flaky(name: str):
            if name == "chain":
                raise KeyError(name)
            return real(name)

        with patch.object(sm, "_get_lattice", side_effect=_flaky):
            entries = list(BOOST_DISPATCH.items())
        assert isinstance(entries, list)
