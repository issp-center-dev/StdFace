"""Unit tests for hphi_writer module.

Tests for the HPhi solver-specific output functions extracted from
``stdface_main.py``.
"""
from __future__ import annotations

import math
import os
import tempfile

import pytest

from stdface.core.stdface_vals import StdIntList, MethodType, ModelType
from stdface.solvers.hphi.writer import (
    large_value,
    print_calc_mod,
    print_excitation,
    vector_potential,
    print_pump,
    _resolve_string_param,
    _configure_spectrum_ops,
    _spectrum_szsz,
    _spectrum_spsm,
    _spectrum_density,
    _spectrum_up,
    _spectrum_down,
    _SPECTRUM_HANDLERS,
    _compute_fourier_coefficients,
    _write_excitation_file,
    _validate_ngpu_scalapack,
    _write_calcmod_file,
    _CalcModParams,
    _pulselaser_At_Et,
    _aclaser_At_Et,
    _dclaser_At_Et,
    _PUMP_TYPE_HANDLERS,
    METHOD_TO_CALC_TYPE,
    RESTART_TO_INT,
    CALC_SPEC_TO_INT,
    MODEL_GC_TO_CALC_MODEL,
    INITIAL_VEC_TYPE_TO_INT,
    EIGENVEC_IO_TO_FLAGS,
    HAM_IO_TO_FLAGS,
    OUTPUT_EX_VEC_TO_INT,
)

# Sentinel values matching what _reset_vals sets at runtime


def _make_stdi_for_large_value(**overrides) -> StdIntList:
    """Create a StdIntList with fields needed by large_value."""
    StdI = StdIntList()
    StdI.LargeValue = float("nan")  # as set by _reset_vals
    StdI.nsite = overrides.get("nsite", 4)
    StdI.trans_list = [(v, 0, 0, 0, 0) for v in overrides.get("trans", [])]
    StdI.intr_list = [(v, 0, 0, 0, 0, 0, 0, 0, 0) for v in overrides.get("intr", [])]
    StdI.Cintra_list = [(v, 0) for v in overrides.get("Cintra", [])]
    StdI.Cinter_list = [(v, 0, 1) for v in overrides.get("Cinter", [])]
    StdI.Ex_list = [(v, 0, 1) for v in overrides.get("Ex", [])]
    StdI.PairLift_list = [(v, 0, 1) for v in overrides.get("PairLift", [])]
    StdI.Hund_list = [(v, 0, 1) for v in overrides.get("Hund", [])]
    return StdI


def _make_stdi_for_calcmod(**overrides) -> StdIntList:
    """Create a StdIntList with fields needed by print_calc_mod."""
    StdI = StdIntList()
    StdI.method = overrides.get("method", "fulldiag")
    StdI.model = overrides.get("model", "hubbard")
    StdI.lGC = overrides.get("lGC", 0)
    StdI.Restart = overrides.get("Restart", None)
    StdI.InitialVecType = overrides.get("InitialVecType", None)
    StdI.EigenVecIO = overrides.get("EigenVecIO", None)
    StdI.HamIO = overrides.get("HamIO", None)
    StdI.CalcSpec = overrides.get("CalcSpec", None)
    StdI.OutputExVec = overrides.get("OutputExVec", None)
    StdI.NGPU = overrides.get("NGPU", None)
    StdI.Scalapack = overrides.get("Scalapack", None)
    return StdI


class TestLargeValue:
    """Tests for the large_value function."""

    def test_basic_computation(self):
        """Test LargeValue computation with known transfer terms."""
        StdI = _make_stdi_for_large_value(
            nsite=4,
            ntrans=4,
            trans=[1.0 + 0j, -1.0 + 0j, 0.5 + 0j, -0.5 + 0j],
        )

        large_value(StdI)

        # sum(|t|) = 1+1+0.5+0.5 = 3.0, divided by nsite=4 => 0.75
        assert StdI.LargeValue == pytest.approx(0.75)

    def test_with_interactions(self):
        """Test LargeValue includes exchange terms with 2x weight."""
        StdI = _make_stdi_for_large_value(
            nsite=2,
            NEx=1, Ex=[0.5],
            NHund=1, Hund=[0.3],
        )

        large_value(StdI)

        # 2*|Ex| + 2*|Hund| = 2*0.5 + 2*0.3 = 1.6, / nsite=2 => 0.8
        assert StdI.LargeValue == pytest.approx(0.8)

    def test_user_specified_value_preserved(self):
        """Test that user-specified LargeValue is preserved."""
        StdI = _make_stdi_for_large_value(nsite=4)
        StdI.LargeValue = 42.0  # user specified (not NaN)

        large_value(StdI)

        assert StdI.LargeValue == pytest.approx(42.0)


class TestPrintCalcMod:
    """Tests for the print_calc_mod function."""

    def test_writes_calcmod_def(self):
        """Test that calcmod.def is created with correct content."""
        StdI = _make_stdi_for_calcmod(method="fulldiag", model="hubbard")

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_calc_mod(StdI)
                assert os.path.exists("calcmod.def")
                content = open("calcmod.def").read()
                assert "CalcType   2" in content  # fulldiag = 2
                assert "CalcModel   0" in content  # hubbard, non-GC = 0
            finally:
                os.chdir(orig)

    def test_spin_gc_model(self):
        """Test SpinGC model mapping."""
        StdI = _make_stdi_for_calcmod(method="lanczos", model="spin", lGC=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_calc_mod(StdI)
                content = open("calcmod.def").read()
                assert "CalcType   0" in content  # lanczos = 0
                assert "CalcModel   4" in content  # spinGC = 4
            finally:
                os.chdir(orig)

    def test_eigenvecio_out(self):
        """Test EigenVecIO=out setting."""
        StdI = _make_stdi_for_calcmod(EigenVecIO="out")

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_calc_mod(StdI)
                content = open("calcmod.def").read()
                assert "OutputEigenVec   1" in content
                assert "InputEigenVec   0" in content
            finally:
                os.chdir(orig)

    def test_kondo_gc_model(self):
        """Test KondoGC model mapping."""
        StdI = _make_stdi_for_calcmod(method="cg", model="kondo", lGC=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_calc_mod(StdI)
                content = open("calcmod.def").read()
                assert "CalcType   3" in content  # cg = 3
                assert "CalcModel   5" in content  # kondoGC = 5
            finally:
                os.chdir(orig)


class TestVectorPotential:
    """Tests for the vector_potential function."""

    def test_quench_default(self):
        """Test that default PumpType is quench with PumpBody=2."""
        StdI = StdIntList()
        StdI.PumpType = None
        StdI.VecPot = [float("nan")] * 3
        StdI.Lanczos_max = None
        StdI.dt = float("nan")
        StdI.freq = float("nan")
        StdI.tshift = float("nan")
        StdI.tdump = float("nan")
        StdI.Uquench = float("nan")
        StdI.ExpandCoef = None

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                vector_potential(StdI)
                assert StdI.PumpBody == 2
                assert StdI.PumpType == "quench"
                assert not os.path.exists("potential.dat")
            finally:
                os.chdir(orig)

    def test_aclaser_writes_potential(self):
        """Test that aclaser pump writes potential.dat."""
        StdI = StdIntList()
        StdI.PumpType = "aclaser"
        StdI.VecPot = [1.0, 0.0, 0.0]
        StdI.Lanczos_max = 10
        StdI.dt = 0.1
        StdI.freq = 1.0
        StdI.tshift = 0.0
        StdI.tdump = float("nan")
        StdI.Uquench = float("nan")
        StdI.ExpandCoef = None

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                vector_potential(StdI)
                assert StdI.PumpBody == 1
                assert os.path.exists("potential.dat")
                lines = open("potential.dat").readlines()
                assert len(lines) == 11  # header + 10 timesteps
            finally:
                os.chdir(orig)


class TestDispatchDicts:
    """Tests for the module-level dispatch dictionaries."""

    # -- METHOD_TO_CALC_TYPE --

    def test_method_lanczos(self):
        assert METHOD_TO_CALC_TYPE["lanczos"] == 0

    def test_method_lanczos_energy(self):
        assert METHOD_TO_CALC_TYPE[MethodType.LANCZOS_ENERGY] == 0

    def test_method_tpq(self):
        assert METHOD_TO_CALC_TYPE["tpq"] == 1

    def test_method_fulldiag(self):
        assert METHOD_TO_CALC_TYPE["fulldiag"] == 2

    def test_method_cg(self):
        assert METHOD_TO_CALC_TYPE["cg"] == 3

    def test_method_time_evolution(self):
        assert METHOD_TO_CALC_TYPE[MethodType.TIME_EVOLUTION] == 4

    def test_method_ctpq(self):
        assert METHOD_TO_CALC_TYPE["ctpq"] == 5

    def test_method_unknown_returns_none(self):
        assert METHOD_TO_CALC_TYPE.get("bogus") is None

    # -- RESTART_TO_INT --

    def test_restart_none(self):
        assert RESTART_TO_INT["none"] == 0

    def test_restart_out(self):
        assert RESTART_TO_INT["restart_out"] == 1

    def test_restart_save_alias(self):
        assert RESTART_TO_INT["save"] == 1

    def test_restart_restartsave(self):
        assert RESTART_TO_INT["restartsave"] == 2

    def test_restart_alias(self):
        assert RESTART_TO_INT["restart"] == 2

    def test_restart_in(self):
        assert RESTART_TO_INT["restart_in"] == 3

    def test_restart_unknown_returns_none(self):
        assert RESTART_TO_INT.get("bogus") is None

    # -- CALC_SPEC_TO_INT --

    def test_calcspec_none(self):
        assert CALC_SPEC_TO_INT["none"] == 0

    def test_calcspec_normal(self):
        assert CALC_SPEC_TO_INT["normal"] == 1

    def test_calcspec_noiteration(self):
        assert CALC_SPEC_TO_INT["noiteration"] == 2

    def test_calcspec_restart_out(self):
        assert CALC_SPEC_TO_INT["restart_out"] == 3

    def test_calcspec_restart_in(self):
        assert CALC_SPEC_TO_INT["restart_in"] == 4

    def test_calcspec_restartsave(self):
        assert CALC_SPEC_TO_INT["restartsave"] == 5

    def test_calcspec_restart_alias(self):
        assert CALC_SPEC_TO_INT["restart"] == 5

    def test_calcspec_unknown_returns_none(self):
        assert CALC_SPEC_TO_INT.get("bogus") is None

    # -- MODEL_GC_TO_CALC_MODEL --

    def test_model_hubbard_non_gc(self):
        assert MODEL_GC_TO_CALC_MODEL[(ModelType.HUBBARD, 0)] == 0

    def test_model_hubbard_gc(self):
        assert MODEL_GC_TO_CALC_MODEL[(ModelType.HUBBARD, 1)] == 3

    def test_model_spin_non_gc(self):
        assert MODEL_GC_TO_CALC_MODEL[(ModelType.SPIN, 0)] == 1

    def test_model_spin_gc(self):
        assert MODEL_GC_TO_CALC_MODEL[(ModelType.SPIN, 1)] == 4

    def test_model_kondo_non_gc(self):
        assert MODEL_GC_TO_CALC_MODEL[(ModelType.KONDO, 0)] == 2

    def test_model_kondo_gc(self):
        assert MODEL_GC_TO_CALC_MODEL[(ModelType.KONDO, 1)] == 5

    def test_model_unknown_returns_none(self):
        assert MODEL_GC_TO_CALC_MODEL.get(("bogus", 0)) is None

    # -- INITIAL_VEC_TYPE_TO_INT --

    def test_initial_vec_c(self):
        assert INITIAL_VEC_TYPE_TO_INT["c"] == 0

    def test_initial_vec_r(self):
        assert INITIAL_VEC_TYPE_TO_INT["r"] == 1

    def test_initial_vec_unknown_returns_none(self):
        assert INITIAL_VEC_TYPE_TO_INT.get("bogus") is None

    # -- EIGENVEC_IO_TO_FLAGS --

    def test_eigenvec_none(self):
        assert EIGENVEC_IO_TO_FLAGS["none"] == (0, 0)

    def test_eigenvec_in(self):
        assert EIGENVEC_IO_TO_FLAGS["in"] == (1, 0)

    def test_eigenvec_out(self):
        assert EIGENVEC_IO_TO_FLAGS["out"] == (0, 1)

    def test_eigenvec_inout(self):
        assert EIGENVEC_IO_TO_FLAGS["inout"] == (1, 1)

    def test_eigenvec_unknown_returns_none(self):
        assert EIGENVEC_IO_TO_FLAGS.get("bogus") is None

    # -- HAM_IO_TO_FLAGS --

    def test_hamio_none(self):
        assert HAM_IO_TO_FLAGS["none"] == (0, 0)

    def test_hamio_out(self):
        assert HAM_IO_TO_FLAGS["out"] == (1, 0)

    def test_hamio_in(self):
        assert HAM_IO_TO_FLAGS["in"] == (0, 1)

    def test_hamio_unknown_returns_none(self):
        assert HAM_IO_TO_FLAGS.get("bogus") is None

    # -- OUTPUT_EX_VEC_TO_INT --

    def test_output_ex_vec_none(self):
        assert OUTPUT_EX_VEC_TO_INT["none"] == 0

    def test_output_ex_vec_out(self):
        assert OUTPUT_EX_VEC_TO_INT["out"] == 1

    def test_output_ex_vec_unknown_returns_none(self):
        assert OUTPUT_EX_VEC_TO_INT.get("bogus") is None


# -------------------------------------------------------------------
#  _resolve_string_param
# -------------------------------------------------------------------


class TestResolveStringParam:
    """Tests for the _resolve_string_param helper."""

    _SAMPLE_DISPATCH = {"alpha": 10, "beta": 20}

    def test_unset_returns_default_value(self):
        """Test that None field returns the default value."""
        StdI = StdIntList()
        StdI.Restart = None
        result = _resolve_string_param(
            StdI, "Restart", "Restart", "none", 0, RESTART_TO_INT)
        assert result == 0

    def test_unset_sets_field_to_default(self):
        """Test that None field is overwritten with default string."""
        StdI = StdIntList()
        StdI.Restart = None
        _resolve_string_param(
            StdI, "Restart", "Restart", "none", 0, RESTART_TO_INT)
        assert StdI.Restart == "none"

    def test_set_value_returns_dispatched_int(self):
        """Test that a set value returns its dispatched integer."""
        StdI = StdIntList()
        StdI.Restart = "save"
        result = _resolve_string_param(
            StdI, "Restart", "Restart", "none", 0, RESTART_TO_INT)
        assert result == 1  # "save" → 1

    def test_set_value_field_unchanged(self):
        """Test that a set value is not overwritten."""
        StdI = StdIntList()
        StdI.Restart = "save"
        _resolve_string_param(
            StdI, "Restart", "Restart", "none", 0, RESTART_TO_INT)
        assert StdI.Restart == "save"

    def test_unknown_value_raises_exit(self):
        """Test that an unknown value raises ValueError."""
        StdI = StdIntList()
        StdI.Restart = "bogus_value"
        with pytest.raises(ValueError):
            _resolve_string_param(
                StdI, "Restart", "Restart", "none", 0, RESTART_TO_INT)

    def test_tuple_dispatch(self):
        """Test that tuple return values are passed through correctly."""
        StdI = StdIntList()
        StdI.EigenVecIO = "inout"
        result = _resolve_string_param(
            StdI, "EigenVecIO", "EigenVecIO", "none",
            (0, 0), EIGENVEC_IO_TO_FLAGS)
        assert result == (1, 1)

    def test_tuple_default(self):
        """Test that tuple default values are returned for unset fields."""
        StdI = StdIntList()
        StdI.EigenVecIO = None
        result = _resolve_string_param(
            StdI, "EigenVecIO", "EigenVecIO", "none",
            (0, 0), EIGENVEC_IO_TO_FLAGS)
        assert result == (0, 0)

    def test_prints_default_message(self, caplog):
        """Test that unset fields produce a DEFAULT VALUE message."""
        import logging

        StdI = StdIntList()
        StdI.CalcSpec = None
        with caplog.at_level(logging.INFO, logger="stdface.solvers.hphi.writer"):
            _resolve_string_param(
                StdI, "CalcSpec", "CalcSpec", "none", 0, CALC_SPEC_TO_INT)
        assert "DEFAULT VALUE IS USED" in caplog.text
        assert "CalcSpec" in caplog.text

    def test_prints_set_value(self, caplog):
        """Test that set fields log their value without DEFAULT."""
        import logging

        StdI = StdIntList()
        StdI.CalcSpec = "normal"
        with caplog.at_level(logging.INFO, logger="stdface.solvers.hphi.writer"):
            _resolve_string_param(
                StdI, "CalcSpec", "CalcSpec", "none", 0, CALC_SPEC_TO_INT)
        assert "normal" in caplog.text
        assert "DEFAULT VALUE IS USED" not in caplog.text


# -------------------------------------------------------------------
#  _configure_spectrum_ops
# -------------------------------------------------------------------


class TestConfigureSpectrumOps:
    """Tests for the _configure_spectrum_ops helper."""

    def _make_arrays(self, n=2):
        """Create coef and spin arrays of size *n*."""
        coef = [0.0] * n
        spin = [[0, 0] for _ in range(n)]
        return coef, spin

    # -- szsz (non-spin model) --

    def test_szsz_hubbard_numop(self):
        """Test szsz with Hubbard model returns NumOp=2."""
        coef, spin = self._make_arrays(2)
        NumOp, body = _configure_spectrum_ops("szsz", ModelType.HUBBARD, 1, coef, spin)
        assert NumOp == 2

    def test_szsz_hubbard_body(self):
        """Test szsz returns SpectrumBody=2 (pair)."""
        coef, spin = self._make_arrays(2)
        _, body = _configure_spectrum_ops("szsz", ModelType.HUBBARD, 1, coef, spin)
        assert body == 2

    def test_szsz_hubbard_coef(self):
        """Test szsz Hubbard coefficients are +0.5 and -0.5."""
        coef, spin = self._make_arrays(2)
        _configure_spectrum_ops("szsz", ModelType.HUBBARD, 1, coef, spin)
        assert coef[0] == pytest.approx(0.5)
        assert coef[1] == pytest.approx(-0.5)

    def test_szsz_hubbard_spin(self):
        """Test szsz Hubbard spin indices are diagonal."""
        coef, spin = self._make_arrays(2)
        _configure_spectrum_ops("szsz", ModelType.HUBBARD, 1, coef, spin)
        assert spin[0] == [0, 0]
        assert spin[1] == [1, 1]

    # -- szsz (spin model, S2=1) --

    def test_szsz_spin_s2_1_numop(self):
        """Test szsz with spin model S2=1 returns NumOp=S2+1=2."""
        coef, spin = self._make_arrays(2)
        NumOp, _ = _configure_spectrum_ops("szsz", ModelType.SPIN, 1, coef, spin)
        assert NumOp == 2

    def test_szsz_spin_s2_1_coef(self):
        """Test szsz spin S2=1 coefficients are Sz values (-0.5, +0.5)."""
        coef, spin = self._make_arrays(2)
        _configure_spectrum_ops("szsz", ModelType.SPIN, 1, coef, spin)
        assert coef[0] == pytest.approx(-0.5)
        assert coef[1] == pytest.approx(0.5)

    # -- szsz (spin model, S2=2 → S=1, three operators) --

    def test_szsz_spin_s2_2_numop(self):
        """Test szsz with spin model S2=2 returns NumOp=3."""
        coef, spin = self._make_arrays(3)
        NumOp, _ = _configure_spectrum_ops("szsz", ModelType.SPIN, 2, coef, spin)
        assert NumOp == 3

    def test_szsz_spin_s2_2_coef(self):
        """Test szsz spin S2=2 coefficients are Sz = -1, 0, +1."""
        coef, spin = self._make_arrays(3)
        _configure_spectrum_ops("szsz", ModelType.SPIN, 2, coef, spin)
        assert coef[0] == pytest.approx(-1.0)
        assert coef[1] == pytest.approx(0.0)
        assert coef[2] == pytest.approx(1.0)

    # -- s+s- (non-spin or S2=1) --

    def test_spsm_hubbard_numop(self):
        """Test s+s- with Hubbard model returns NumOp=1."""
        coef, spin = self._make_arrays(2)
        NumOp, body = _configure_spectrum_ops("s+s-", ModelType.HUBBARD, 1, coef, spin)
        assert NumOp == 1
        assert body == 2

    def test_spsm_hubbard_coef_and_spin(self):
        """Test s+s- Hubbard has coef=1 and spin=[1,0]."""
        coef, spin = self._make_arrays(2)
        _configure_spectrum_ops("s+s-", ModelType.HUBBARD, 1, coef, spin)
        assert coef[0] == pytest.approx(1.0)
        assert spin[0] == [1, 0]

    # -- s+s- (spin model, S2=2 → high spin) --

    def test_spsm_spin_s2_2_numop(self):
        """Test s+s- with spin model S2=2 returns NumOp=S2=2."""
        coef, spin = self._make_arrays(3)
        NumOp, _ = _configure_spectrum_ops("s+s-", ModelType.SPIN, 2, coef, spin)
        assert NumOp == 2

    def test_spsm_spin_s2_2_coef(self):
        """Test s+s- spin S2=2 coefficients are sqrt(S(S+1)-Sz(Sz+1))."""
        coef, spin = self._make_arrays(3)
        _configure_spectrum_ops("s+s-", ModelType.SPIN, 2, coef, spin)
        # S=1; ispin=1: Sz=0 → sqrt(2-0)=sqrt(2); ispin=2: Sz=-1 → sqrt(2-0)=sqrt(2)
        assert coef[0] == pytest.approx(math.sqrt(2.0))
        assert coef[1] == pytest.approx(math.sqrt(2.0))

    # -- density --

    def test_density_returns(self):
        """Test density returns NumOp=2, SpectrumBody=2."""
        coef, spin = self._make_arrays(2)
        NumOp, body = _configure_spectrum_ops("density", ModelType.HUBBARD, 1, coef, spin)
        assert NumOp == 2
        assert body == 2

    def test_density_coef(self):
        """Test density coefficients are both 1.0."""
        coef, spin = self._make_arrays(2)
        _configure_spectrum_ops("density", ModelType.HUBBARD, 1, coef, spin)
        assert coef[0] == pytest.approx(1.0)
        assert coef[1] == pytest.approx(1.0)

    # -- up --

    def test_up_returns(self):
        """Test up returns NumOp=1, SpectrumBody=1 (single)."""
        coef, spin = self._make_arrays(2)
        NumOp, body = _configure_spectrum_ops("up", ModelType.HUBBARD, 1, coef, spin)
        assert NumOp == 1
        assert body == 1

    def test_up_spin(self):
        """Test up sets spin[0][0]=0."""
        coef, spin = self._make_arrays(2)
        _configure_spectrum_ops("up", ModelType.HUBBARD, 1, coef, spin)
        assert spin[0][0] == 0

    # -- down --

    def test_down_returns(self):
        """Test down returns NumOp=1, SpectrumBody=1 (single)."""
        coef, spin = self._make_arrays(2)
        NumOp, body = _configure_spectrum_ops("down", ModelType.HUBBARD, 1, coef, spin)
        assert NumOp == 1
        assert body == 1

    def test_down_spin(self):
        """Test down sets spin[0][0]=1."""
        coef, spin = self._make_arrays(2)
        _configure_spectrum_ops("down", ModelType.HUBBARD, 1, coef, spin)
        assert spin[0][0] == 1

    # -- unknown --

    def test_unknown_spectrum_type_exits(self):
        """Test that an unknown spectrum type raises ValueError."""
        coef, spin = self._make_arrays(2)
        with pytest.raises(ValueError):
            _configure_spectrum_ops("bogus", ModelType.HUBBARD, 1, coef, spin)


# -------------------------------------------------------------------
#  Spectrum-type handler functions
# -------------------------------------------------------------------


class TestSpectrumHandlersSzSz:
    """Tests for _spectrum_szsz handler."""

    def _make_arrays(self, n=2):
        coef = [0.0] * n
        spin = [[0, 0] for _ in range(n)]
        return coef, spin

    def test_matches_dispatcher(self):
        """Handler result matches _configure_spectrum_ops for szsz."""
        coef1, spin1 = self._make_arrays(3)
        coef2, spin2 = self._make_arrays(3)
        r1 = _spectrum_szsz(ModelType.SPIN, 2, coef1, spin1)
        r2 = _configure_spectrum_ops("szsz", ModelType.SPIN, 2, coef2, spin2)
        assert r1 == r2
        assert coef1 == coef2
        assert spin1 == spin2


class TestSpectrumHandlersSpsm:
    """Tests for _spectrum_spsm handler."""

    def _make_arrays(self, n=2):
        coef = [0.0] * n
        spin = [[0, 0] for _ in range(n)]
        return coef, spin

    def test_matches_dispatcher(self):
        """Handler result matches _configure_spectrum_ops for s+s-."""
        coef1, spin1 = self._make_arrays(3)
        coef2, spin2 = self._make_arrays(3)
        r1 = _spectrum_spsm(ModelType.SPIN, 2, coef1, spin1)
        r2 = _configure_spectrum_ops("s+s-", ModelType.SPIN, 2, coef2, spin2)
        assert r1 == r2
        assert coef1 == coef2
        assert spin1 == spin2


class TestSpectrumHandlersDensity:
    """Tests for _spectrum_density handler."""

    def _make_arrays(self, n=2):
        coef = [0.0] * n
        spin = [[0, 0] for _ in range(n)]
        return coef, spin

    def test_matches_dispatcher(self):
        """Handler result matches _configure_spectrum_ops for density."""
        coef1, spin1 = self._make_arrays(2)
        coef2, spin2 = self._make_arrays(2)
        r1 = _spectrum_density(ModelType.HUBBARD, 1, coef1, spin1)
        r2 = _configure_spectrum_ops("density", ModelType.HUBBARD, 1, coef2, spin2)
        assert r1 == r2
        assert coef1 == coef2


class TestSpectrumHandlersUpDown:
    """Tests for _spectrum_up and _spectrum_down handlers."""

    def _make_arrays(self, n=2):
        coef = [0.0] * n
        spin = [[0, 0] for _ in range(n)]
        return coef, spin

    def test_up_matches_dispatcher(self):
        """_spectrum_up matches _configure_spectrum_ops for up."""
        coef1, spin1 = self._make_arrays(2)
        coef2, spin2 = self._make_arrays(2)
        r1 = _spectrum_up(ModelType.HUBBARD, 1, coef1, spin1)
        r2 = _configure_spectrum_ops("up", ModelType.HUBBARD, 1, coef2, spin2)
        assert r1 == r2
        assert coef1 == coef2

    def test_down_matches_dispatcher(self):
        """_spectrum_down matches _configure_spectrum_ops for down."""
        coef1, spin1 = self._make_arrays(2)
        coef2, spin2 = self._make_arrays(2)
        r1 = _spectrum_down(ModelType.HUBBARD, 1, coef1, spin1)
        r2 = _configure_spectrum_ops("down", ModelType.HUBBARD, 1, coef2, spin2)
        assert r1 == r2
        assert coef1 == coef2


class TestSpectrumHandlersDispatch:
    """Tests for the _SPECTRUM_HANDLERS dispatch table."""

    def test_all_five_types_registered(self):
        """All five spectrum types are in the dispatch table."""
        assert set(_SPECTRUM_HANDLERS.keys()) == {
            "szsz", "s+s-", "density", "up", "down",
        }

    def test_szsz_points_to_handler(self):
        """szsz maps to _spectrum_szsz."""
        assert _SPECTRUM_HANDLERS["szsz"] is _spectrum_szsz

    def test_spsm_points_to_handler(self):
        """s+s- maps to _spectrum_spsm."""
        assert _SPECTRUM_HANDLERS["s+s-"] is _spectrum_spsm

    def test_density_points_to_handler(self):
        """density maps to _spectrum_density."""
        assert _SPECTRUM_HANDLERS["density"] is _spectrum_density


# -------------------------------------------------------------------
#  Pump-type handler functions
# -------------------------------------------------------------------


class TestPulseLaserAtEt:
    """Tests for _pulselaser_At_Et helper."""

    def test_at_zero_at_center(self):
        """Test A(t) at t=tshift is V*1*1 = V (cos(0)=1, gauss peak)."""
        At, Et = _pulselaser_At_Et(time=1.0, V=2.0, freq=1.0, tshift=1.0, tdump=0.5)
        assert At == pytest.approx(2.0)

    def test_at_is_zero_far_from_center(self):
        """Test A(t) decays to ~0 far from tshift."""
        At, Et = _pulselaser_At_Et(time=100.0, V=2.0, freq=1.0, tshift=0.0, tdump=0.1)
        assert abs(At) < 1e-10

    def test_et_at_center_tshift_zero(self):
        """Test E(t) at t=tshift=0: dt_shift=0, gauss=1, cos=1, sin=0.
        Et = -V * (0 * 1 - freq * 0) * 1 = 0."""
        _, Et = _pulselaser_At_Et(time=0.0, V=1.0, freq=2.0, tshift=0.0, tdump=1.0)
        assert Et == pytest.approx(0.0)

    def test_symmetry_with_negative_v(self):
        """Test that negative V flips sign of both A and E."""
        At_pos, Et_pos = _pulselaser_At_Et(0.5, 1.0, 1.0, 0.0, 0.5)
        At_neg, Et_neg = _pulselaser_At_Et(0.5, -1.0, 1.0, 0.0, 0.5)
        assert At_neg == pytest.approx(-At_pos)
        assert Et_neg == pytest.approx(-Et_pos)

    def test_matches_original_formula(self):
        """Test against the original inline expressions."""
        time, V, freq, tshift, tdump = 0.3, 1.5, 2.0, 0.1, 0.4
        dt = time - tshift
        gauss = math.exp(-0.5 * dt ** 2 / (tdump * tdump))
        expected_At = V * math.cos(freq * dt) * gauss
        expected_Et = (-V * ((tshift - time) / (tdump * tdump)
                             * math.cos(freq * dt)
                             - freq * math.sin(freq * dt)) * gauss)
        At, Et = _pulselaser_At_Et(time, V, freq, tshift, tdump)
        assert At == pytest.approx(expected_At)
        assert Et == pytest.approx(expected_Et)


class TestAcLaserAtEt:
    """Tests for _aclaser_At_Et helper."""

    def test_at_at_tshift(self):
        """Test A(tshift) = V * sin(0) = 0."""
        At, Et = _aclaser_At_Et(time=1.0, V=2.0, freq=1.0, tshift=1.0, tdump=0.0)
        assert At == pytest.approx(0.0)

    def test_et_at_tshift(self):
        """Test E(tshift) = V * cos(0) * freq = V * freq."""
        _, Et = _aclaser_At_Et(time=1.0, V=2.0, freq=3.0, tshift=1.0, tdump=0.0)
        assert Et == pytest.approx(6.0)

    def test_quarter_period(self):
        """Test A at quarter period: sin(pi/2) = 1."""
        freq = 1.0
        quarter_T = math.pi / (2.0 * freq)
        At, _ = _aclaser_At_Et(time=quarter_T, V=1.0, freq=freq, tshift=0.0, tdump=0.0)
        assert At == pytest.approx(1.0)


class TestDcLaserAtEt:
    """Tests for _dclaser_At_Et helper."""

    def test_at_is_linear(self):
        """Test A(t) = V * t."""
        At, _ = _dclaser_At_Et(time=3.0, V=2.0, freq=0.0, tshift=0.0, tdump=0.0)
        assert At == pytest.approx(6.0)

    def test_et_is_negative_v(self):
        """Test E(t) = -V (constant)."""
        _, Et = _dclaser_At_Et(time=3.0, V=2.0, freq=0.0, tshift=0.0, tdump=0.0)
        assert Et == pytest.approx(-2.0)

    def test_at_zero_at_origin(self):
        """Test A(0) = 0."""
        At, _ = _dclaser_At_Et(time=0.0, V=5.0, freq=0.0, tshift=0.0, tdump=0.0)
        assert At == pytest.approx(0.0)


class TestPumpTypeHandlers:
    """Tests for _PUMP_TYPE_HANDLERS dispatch table."""

    def test_quench_entry(self):
        """Test quench has PumpBody=2 and no handler."""
        body, fn = _PUMP_TYPE_HANDLERS["quench"]
        assert body == 2
        assert fn is None

    def test_pulselaser_entry(self):
        """Test pulselaser has PumpBody=1 and correct handler."""
        body, fn = _PUMP_TYPE_HANDLERS["pulselaser"]
        assert body == 1
        assert fn is _pulselaser_At_Et

    def test_aclaser_entry(self):
        """Test aclaser has PumpBody=1 and correct handler."""
        body, fn = _PUMP_TYPE_HANDLERS["aclaser"]
        assert body == 1
        assert fn is _aclaser_At_Et

    def test_dclaser_entry(self):
        """Test dclaser has PumpBody=1 and correct handler."""
        body, fn = _PUMP_TYPE_HANDLERS["dclaser"]
        assert body == 1
        assert fn is _dclaser_At_Et

    def test_unknown_returns_none(self):
        """Test unknown PumpType returns None from dict."""
        assert _PUMP_TYPE_HANDLERS.get("bogus") is None

    def test_all_one_body_have_handlers(self):
        """Test all PumpBody=1 entries have non-None handlers."""
        for name, (body, fn) in _PUMP_TYPE_HANDLERS.items():
            if body == 1:
                assert fn is not None, f"{name} has PumpBody=1 but no handler"

    def test_unknown_pump_type_exits(self):
        """Test that vector_potential rejects unknown PumpType."""
        StdI = StdIntList()
        StdI.PumpType = "bogus_laser"
        StdI.VecPot = [0.0] * 3
        StdI.Lanczos_max = None
        StdI.dt = float("nan")
        StdI.freq = float("nan")
        StdI.tshift = float("nan")
        StdI.tdump = float("nan")
        StdI.Uquench = float("nan")
        StdI.ExpandCoef = None
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                with pytest.raises(ValueError):
                    vector_potential(StdI)
            finally:
                os.chdir(orig)


# ---------------------------------------------------------------------------
#  Tests: _compute_fourier_coefficients
# ---------------------------------------------------------------------------


def _make_spectrum_StdI(
    model: str = "hubbard",
    nsite: int = 4,
    NCell: int = 2,
    NsiteUC: int = 2,
) -> StdIntList:
    """Create a minimal StdIntList for Fourier-coefficient tests."""
    import numpy as np
    s = StdIntList()
    s.model = model
    s.nsite = nsite
    s.NCell = NCell
    s.NsiteUC = NsiteUC
    s.pi = math.acos(-1.0)
    s.Cell = np.array([[0, 0, 0], [1, 0, 0]], dtype=int)
    s.tau = np.zeros((NsiteUC, 3))
    s.tau[0] = [0.0, 0.0, 0.0]
    s.tau[1] = [0.5, 0.0, 0.0]
    s.SpectrumQ = [0.0, 0.0, 0.0]
    return s


class TestComputeFourierCoefficients:
    """Tests for _compute_fourier_coefficients."""

    def test_zero_q_gives_all_ones(self):
        """With Q=0, all Fourier real parts should be 1, imaginary 0."""
        s = _make_spectrum_StdI()
        s.SpectrumQ = [0.0, 0.0, 0.0]
        fr, fi = _compute_fourier_coefficients(s)
        assert len(fr) == s.nsite
        assert len(fi) == s.nsite
        for i in range(s.nsite):
            assert abs(fr[i] - 1.0) < 1e-12
            assert abs(fi[i]) < 1e-12

    def test_nonzero_q(self):
        """With Q!=0, Fourier coefficients should be non-trivial."""
        s = _make_spectrum_StdI()
        s.SpectrumQ = [0.5, 0.0, 0.0]
        fr, fi = _compute_fourier_coefficients(s)
        # Cell[0]=(0,0,0), tau[0]=(0,0,0) => Cphase=0 => cos=1, sin=0
        assert abs(fr[0] - 1.0) < 1e-12
        # Cell[0]=(0,0,0), tau[1]=(0.5,0,0) => Cphase=0.25 => cos(pi/2)=0
        assert abs(fr[1]) < 1e-12

    def test_kondo_duplicates_coefficients(self):
        """For Kondo model, second half should mirror first half."""
        s = _make_spectrum_StdI(model="kondo", nsite=8, NCell=2, NsiteUC=2)
        # Expand Cell and tau for 4 itinerant sites + 4 local-spin
        import numpy as np
        s.Cell = np.array([[0, 0, 0], [1, 0, 0]], dtype=int)
        s.tau = np.zeros((2, 3))
        s.SpectrumQ = [0.25, 0.0, 0.0]
        fr, fi = _compute_fourier_coefficients(s)
        half = s.nsite // 2
        for i in range(half):
            assert abs(fr[i] - fr[i + half]) < 1e-12
            assert abs(fi[i] - fi[i + half]) < 1e-12

    def test_returns_correct_length(self):
        """Output lists should have length nsite."""
        s = _make_spectrum_StdI(nsite=6, NCell=3, NsiteUC=2)
        import numpy as np
        s.Cell = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=int)
        s.SpectrumQ = [0.0, 0.0, 0.0]
        fr, fi = _compute_fourier_coefficients(s)
        assert len(fr) == 6
        assert len(fi) == 6


# ---------------------------------------------------------------------------
#  Tests: _write_excitation_file
# ---------------------------------------------------------------------------


class TestWriteExcitationFile:
    """Tests for _write_excitation_file."""

    def test_writes_pair_def(self, tmp_path):
        """Should create pair.def for SpectrumBody=2."""
        os.chdir(tmp_path)
        s = StdIntList()
        s.model = "hubbard"
        s.nsite = 2
        s.SpectrumBody = 2
        NumOp = 2
        coef = [1.0, -1.0]
        spin = [[0, 0], [1, 1]]
        fr = [1.0, 0.5]
        fi = [0.0, 0.1]
        _write_excitation_file(s, NumOp, coef, spin, fr, fi)
        assert (tmp_path / "pair.def").exists()
        content = (tmp_path / "pair.def").read_text()
        assert "NPair 4" in content
        assert "Pair Excitation" in content

    def test_writes_single_def(self, tmp_path):
        """Should create single.def for SpectrumBody=1."""
        os.chdir(tmp_path)
        s = StdIntList()
        s.model = "hubbard"
        s.nsite = 2
        s.SpectrumBody = 1
        NumOp = 1
        coef = [1.0]
        spin = [[0, 0]]
        fr = [1.0, 0.5]
        fi = [0.0, 0.1]
        _write_excitation_file(s, NumOp, coef, spin, fr, fi)
        assert (tmp_path / "single.def").exists()
        content = (tmp_path / "single.def").read_text()
        assert "NSingle 2" in content
        assert "Single Excitation" in content

    def test_kondo_single_uses_half_sites(self, tmp_path):
        """For Kondo model, NSingle should use nsite//2."""
        os.chdir(tmp_path)
        s = StdIntList()
        s.model = "kondo"
        s.nsite = 8
        s.SpectrumBody = 1
        NumOp = 1
        coef = [1.0]
        spin = [[0, 0]]
        fr = [1.0] * 8
        fi = [0.0] * 8
        _write_excitation_file(s, NumOp, coef, spin, fr, fi)
        content = (tmp_path / "single.def").read_text()
        assert "NSingle 4" in content

    def test_pair_def_line_count(self, tmp_path):
        """pair.def should have nsite*NumOp data lines."""
        os.chdir(tmp_path)
        s = StdIntList()
        s.model = "hubbard"
        s.nsite = 3
        s.SpectrumBody = 2
        NumOp = 2
        coef = [1.0, -1.0]
        spin = [[0, 0], [1, 1]]
        fr = [1.0, 0.5, 0.3]
        fi = [0.0, 0.1, 0.2]
        _write_excitation_file(s, NumOp, coef, spin, fr, fi)
        content = (tmp_path / "pair.def").read_text()
        # 3 header lines + nsite*NumOp data lines = 3 + 6 = 9 total
        lines = [l for l in content.strip().split("\n") if l.strip()]
        data_lines = [l for l in lines if "=" not in l and "NPair" not in l]
        assert len(data_lines) == 3 * 2  # nsite * NumOp


# ===================================================================
#  _validate_ngpu_scalapack
# ===================================================================


class TestValidateNGPUScalapack:
    """Tests for the _validate_ngpu_scalapack validation helper."""

    def test_no_validation_when_unset(self):
        """Test that no error occurs when NGPU and Scalapack are unset."""
        StdI = StdIntList()
        StdI.NGPU = None
        StdI.Scalapack = None
        _validate_ngpu_scalapack(StdI)  # should not raise

    def test_valid_ngpu(self):
        """Test that valid NGPU passes validation."""
        StdI = StdIntList()
        StdI.NGPU = 2
        StdI.Scalapack = None
        _validate_ngpu_scalapack(StdI)  # should not raise

    def test_invalid_ngpu_exits(self):
        """Test that NGPU < 1 causes exit."""
        StdI = StdIntList()
        StdI.NGPU = 0
        StdI.Scalapack = None
        with pytest.raises(ValueError):
            _validate_ngpu_scalapack(StdI)

    def test_valid_scalapack(self):
        """Test that valid Scalapack (0 or 1) passes."""
        StdI = StdIntList()
        StdI.NGPU = None
        StdI.Scalapack = 1
        _validate_ngpu_scalapack(StdI)  # should not raise

    def test_invalid_scalapack_exits(self):
        """Test that Scalapack > 1 causes exit."""
        StdI = StdIntList()
        StdI.NGPU = None
        StdI.Scalapack = 2
        with pytest.raises(ValueError):
            _validate_ngpu_scalapack(StdI)

    def test_negative_scalapack_exits(self):
        """Test that Scalapack < 0 causes exit."""
        StdI = StdIntList()
        StdI.NGPU = None
        StdI.Scalapack = -1
        with pytest.raises(ValueError):
            _validate_ngpu_scalapack(StdI)


# ===================================================================
#  _write_calcmod_file
# ===================================================================


class TestWriteCalcmodFile:
    """Tests for the _write_calcmod_file output helper."""

    def test_writes_calcmod_def(self, tmp_path):
        """Test that calcmod.def is created."""
        os.chdir(tmp_path)
        StdI = StdIntList()
        StdI.NGPU = None
        StdI.Scalapack = None
        _write_calcmod_file(
            StdI, _CalcModParams(
                iCalcType=0, iCalcModel=0, iCalcEigenvec=0,
                iRestart=0, iCalcSpec=0, iInitialVecType=0,
                InputEigenVec=0, OutputEigenVec=0,
                iInputHam=0, iOutputHam=0, iOutputExVec=0,
            ),
        )
        assert (tmp_path / "calcmod.def").exists()

    def test_contains_calc_type(self, tmp_path):
        """Test that CalcType value appears in output."""
        os.chdir(tmp_path)
        StdI = StdIntList()
        StdI.NGPU = None
        StdI.Scalapack = None
        _write_calcmod_file(
            StdI, _CalcModParams(
                iCalcType=2, iCalcModel=3, iCalcEigenvec=1,
                iRestart=1, iCalcSpec=0, iInitialVecType=0,
                InputEigenVec=0, OutputEigenVec=0,
                iInputHam=0, iOutputHam=0, iOutputExVec=0,
            ),
        )
        content = (tmp_path / "calcmod.def").read_text()
        assert "CalcType   2" in content
        assert "CalcModel   3" in content
        assert "CalcEigenVec   1" in content

    def test_ngpu_line_present(self, tmp_path):
        """Test that NGPU line appears when NGPU is set."""
        os.chdir(tmp_path)
        StdI = StdIntList()
        StdI.NGPU = 4
        StdI.Scalapack = None
        _write_calcmod_file(
            StdI, _CalcModParams(
                iCalcType=2, iCalcModel=0, iCalcEigenvec=0,
                iRestart=0, iCalcSpec=0, iInitialVecType=0,
                InputEigenVec=0, OutputEigenVec=0,
                iInputHam=0, iOutputHam=0, iOutputExVec=0,
            ),
        )
        content = (tmp_path / "calcmod.def").read_text()
        assert "NGPU   4" in content

    def test_scalapack_line_present(self, tmp_path):
        """Test that Scalapack line appears when Scalapack is set."""
        os.chdir(tmp_path)
        StdI = StdIntList()
        StdI.NGPU = None
        StdI.Scalapack = 1
        _write_calcmod_file(
            StdI, _CalcModParams(
                iCalcType=0, iCalcModel=0, iCalcEigenvec=0,
                iRestart=0, iCalcSpec=0, iInitialVecType=0,
                InputEigenVec=0, OutputEigenVec=0,
                iInputHam=0, iOutputHam=0, iOutputExVec=0,
            ),
        )
        content = (tmp_path / "calcmod.def").read_text()
        assert "Scalapack   1" in content

    def test_no_ngpu_line_when_unset(self, tmp_path):
        """Test that NGPU line is absent when NGPU is unset."""
        os.chdir(tmp_path)
        StdI = StdIntList()
        StdI.NGPU = None
        StdI.Scalapack = None
        _write_calcmod_file(
            StdI, _CalcModParams(
                iCalcType=0, iCalcModel=0, iCalcEigenvec=0,
                iRestart=0, iCalcSpec=0, iInitialVecType=0,
                InputEigenVec=0, OutputEigenVec=0,
                iInputHam=0, iOutputHam=0, iOutputExVec=0,
            ),
        )
        content = (tmp_path / "calcmod.def").read_text()
        assert "NGPU" not in content
        assert "Scalapack" not in content
