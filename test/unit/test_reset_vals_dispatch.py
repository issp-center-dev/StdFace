"""Unit tests for data-driven field reset tables and functions.

Tests for the common (solver-independent) reset tables
``_COMMON_RESET_SCALARS``, ``_COMMON_RESET_ARRAYS``, the solver-specific
tables ``_SOLVER_RESET_SCALARS``, ``_SOLVER_RESET_ARRAYS``,
``_UHF_BASE_SCALARS``, ``_UHF_BASE_ARRAYS``, the ``_apply_field_resets``
function, and the refactored ``_reset_vals`` function in ``stdface_main``.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from stdface.core.stdface_main import (
    _COMMON_RESET_SCALARS,
    _COMMON_RESET_ARRAYS,
    _SOLVER_RESET_SCALARS,
    _SOLVER_RESET_ARRAYS,
    _UHF_BASE_SCALARS,
    _UHF_BASE_ARRAYS,
    _apply_field_resets,
    _reset_vals,
    NaN_i,
)
from stdface.core.stdface_vals import StdIntList, SolverType


def _make_stdi(solver: str) -> StdIntList:
    """Create a minimal StdIntList for testing reset functions."""
    from stdface.core.stdface_main import _attach_solver_config
    StdI = StdIntList()
    StdI.solver = solver
    _attach_solver_config(StdI)  # C3: solver-specific fields live on the config
    StdI.pi = math.acos(-1.0)
    return StdI


# -------------------------------------------------------------------
#  Table-structure tests
# -------------------------------------------------------------------


class TestResetTableStructure:
    """Tests for the reset-field data tables."""

    def test_scalars_has_four_solvers(self):
        """Test that scalar table has exactly 4 solver entries."""
        assert len(_SOLVER_RESET_SCALARS) == 4

    def test_arrays_has_four_solvers(self):
        """Test that array table has exactly 4 solver entries."""
        assert len(_SOLVER_RESET_ARRAYS) == 4

    def test_all_solvers_present_in_scalars(self):
        """Test that HPhi, mVMC, UHF, HWAVE are all in scalars table."""
        for solver in (SolverType.HPhi, SolverType.mVMC,
                       SolverType.UHF, SolverType.HWAVE):
            assert solver in _SOLVER_RESET_SCALARS

    def test_all_solvers_present_in_arrays(self):
        """Test that HPhi, mVMC, UHF, HWAVE are all in arrays table."""
        for solver in (SolverType.HPhi, SolverType.mVMC,
                       SolverType.UHF, SolverType.HWAVE):
            assert solver in _SOLVER_RESET_ARRAYS

    def test_scalar_entries_are_name_value_tuples(self):
        """Test that every scalar entry is a (str, value) tuple."""
        for solver, entries in _SOLVER_RESET_SCALARS.items():
            for entry in entries:
                assert isinstance(entry, tuple), f"{solver}: {entry}"
                assert len(entry) == 2, f"{solver}: {entry}"
                assert isinstance(entry[0], str), f"{solver}: {entry}"

    def test_array_entries_are_name_value_tuples(self):
        """Test that every array entry is a (str, value) tuple."""
        for solver, entries in _SOLVER_RESET_ARRAYS.items():
            for entry in entries:
                assert isinstance(entry, tuple), f"{solver}: {entry}"
                assert len(entry) == 2, f"{solver}: {entry}"
                assert isinstance(entry[0], str), f"{solver}: {entry}"


# -------------------------------------------------------------------
#  UHF/HWAVE base sharing tests
# -------------------------------------------------------------------


class TestUHFBaseSharing:
    """Tests for UHF/HWAVE shared base tables."""

    def test_uhf_uses_base_scalars(self):
        """Test that UHF solver scalars is exactly _UHF_BASE_SCALARS."""
        assert _SOLVER_RESET_SCALARS[SolverType.UHF] is _UHF_BASE_SCALARS

    def test_uhf_uses_base_arrays(self):
        """Test that UHF solver arrays is exactly _UHF_BASE_ARRAYS."""
        assert _SOLVER_RESET_ARRAYS[SolverType.UHF] is _UHF_BASE_ARRAYS

    def test_hwave_scalars_extend_uhf(self):
        """Test that HWAVE scalars start with the UHF base entries."""
        hwave = _SOLVER_RESET_SCALARS[SolverType.HWAVE]
        uhf = _UHF_BASE_SCALARS
        assert hwave[:len(uhf)] == uhf

    def test_hwave_has_extra_scalars(self):
        """Test that HWAVE adds export_all and lattice_gp beyond UHF."""
        hwave = _SOLVER_RESET_SCALARS[SolverType.HWAVE]
        uhf = _UHF_BASE_SCALARS
        extra_names = [name for name, _ in hwave[len(uhf):]]
        assert "export_all" in extra_names
        assert "lattice_gp" in extra_names

    def test_hwave_arrays_same_as_uhf(self):
        """Test that HWAVE arrays is exactly _UHF_BASE_ARRAYS."""
        assert _SOLVER_RESET_ARRAYS[SolverType.HWAVE] is _UHF_BASE_ARRAYS


# -------------------------------------------------------------------
#  _apply_field_resets — HPhi
# -------------------------------------------------------------------


class TestApplyFieldResetsHPhi:
    """Tests for _apply_field_resets with HPhi solver."""

    def test_sets_nan_d_fields(self):
        """Test that float scalar fields are set to None (A2)."""
        StdI = _make_stdi("HPhi")
        _apply_field_resets(StdI, SolverType.HPhi)
        assert StdI.LargeValue is None
        assert StdI.OmegaMax is None
        assert StdI.OmegaMin is None
        assert StdI.dt is None

    def test_sets_nan_i_fields(self):
        """Test that integer fields are set to None."""
        StdI = _make_stdi("HPhi")
        _apply_field_resets(StdI, SolverType.HPhi)
        assert StdI.Nomega is None
        assert StdI.Lanczos_max is None
        assert StdI.exct is None
        assert StdI.NGPU is None

    def test_sets_flg_temp(self):
        """Test that FlgTemp is set to 1 (not NaN)."""
        StdI = _make_stdi("HPhi")
        _apply_field_resets(StdI, SolverType.HPhi)
        assert StdI.FlgTemp == 1

    def test_fills_spectrum_q_array(self):
        """Test that SpectrumQ array is filled with NaN."""
        StdI = _make_stdi("HPhi")
        _apply_field_resets(StdI, SolverType.HPhi)
        assert all(math.isnan(v) for v in StdI.SpectrumQ)

    def test_fills_vecpot_array(self):
        """Test that VecPot array is filled with NaN."""
        StdI = _make_stdi("HPhi")
        _apply_field_resets(StdI, SolverType.HPhi)
        assert all(math.isnan(v) for v in StdI.VecPot)


# -------------------------------------------------------------------
#  _apply_field_resets — mVMC
# -------------------------------------------------------------------


class TestApplyFieldResetsMVMC:
    """Tests for _apply_field_resets with mVMC solver."""

    def test_sets_nan_i_fields(self):
        """Test that integer fields are set to None."""
        StdI = _make_stdi("mVMC")
        _apply_field_resets(StdI, SolverType.mVMC)
        assert StdI.NVMCCalMode is None
        assert StdI.NVMCSample is None
        assert StdI.ComplexType is None
        assert StdI.RndSeed is None

    def test_sets_nan_d_fields(self):
        """Test that float scalar fields are set to None (A2)."""
        StdI = _make_stdi("mVMC")
        _apply_field_resets(StdI, SolverType.mVMC)
        assert StdI.DSROptRedCut is None
        assert StdI.DSROptStaDel is None

    def test_sets_boxsub_scalars(self):
        """Test that Hsub/Lsub/Wsub are set to None."""
        StdI = _make_stdi("mVMC")
        _apply_field_resets(StdI, SolverType.mVMC)
        assert StdI.Hsub is None
        assert StdI.Lsub is None
        assert StdI.Wsub is None

    def test_fills_boxsub_array(self):
        """Test that boxsub array is filled with NaN_i."""
        StdI = _make_stdi("mVMC")
        _apply_field_resets(StdI, SolverType.mVMC)
        assert (StdI.boxsub == NaN_i).all()


# -------------------------------------------------------------------
#  _apply_field_resets — UHF
# -------------------------------------------------------------------


class TestApplyFieldResetsUHF:
    """Tests for _apply_field_resets with UHF solver."""

    def test_sets_nan_i_fields(self):
        """Test that integer fields are set to None."""
        StdI = _make_stdi("UHF")
        _apply_field_resets(StdI, SolverType.UHF)
        assert StdI.NMPTrans is None
        assert StdI.RndSeed is None
        assert StdI.Iteration_max is None

    def test_sets_nan_d_fields(self):
        """Test that float scalar fields are set to None (A2)."""
        StdI = _make_stdi("UHF")
        _apply_field_resets(StdI, SolverType.UHF)
        assert StdI.mix is None

    def test_fills_boxsub_array(self):
        """Test that boxsub array is filled with NaN_i."""
        StdI = _make_stdi("UHF")
        _apply_field_resets(StdI, SolverType.UHF)
        assert (StdI.boxsub == NaN_i).all()

    def test_sets_boxsub_scalars(self):
        """Test that Hsub/Lsub/Wsub are set to None."""
        StdI = _make_stdi("UHF")
        _apply_field_resets(StdI, SolverType.UHF)
        assert StdI.Hsub is None
        assert StdI.Lsub is None
        assert StdI.Wsub is None


# -------------------------------------------------------------------
#  _apply_field_resets — HWAVE
# -------------------------------------------------------------------


class TestApplyFieldResetsHWAVE:
    """Tests for _apply_field_resets with HWAVE solver."""

    def test_sets_nan_i_fields(self):
        """Test that integer fields are set to None."""
        StdI = _make_stdi("HWAVE")
        _apply_field_resets(StdI, SolverType.HWAVE)
        assert StdI.NMPTrans is None
        assert StdI.RndSeed is None
        assert StdI.export_all is None
        assert StdI.lattice_gp is None

    def test_sets_nan_d_fields(self):
        """Test that float scalar fields are set to None (A2)."""
        StdI = _make_stdi("HWAVE")
        _apply_field_resets(StdI, SolverType.HWAVE)
        assert StdI.mix is None

    def test_fills_boxsub_array(self):
        """Test that boxsub array is filled with NaN_i."""
        StdI = _make_stdi("HWAVE")
        _apply_field_resets(StdI, SolverType.HWAVE)
        assert (StdI.boxsub == NaN_i).all()

    def test_has_extra_fields_vs_uhf(self):
        """Test that HWAVE resets export_all/lattice_gp beyond UHF."""
        StdI = _make_stdi("HWAVE")
        _apply_field_resets(StdI, SolverType.HWAVE)
        assert StdI.export_all is None
        assert StdI.lattice_gp is None


# -------------------------------------------------------------------
#  _apply_field_resets — unknown solver
# -------------------------------------------------------------------


class TestApplyFieldResetsUnknown:
    """Tests for _apply_field_resets with an unknown solver."""

    def test_unknown_solver_no_error(self):
        """Test that unknown solver type does not raise."""
        StdI = _make_stdi("unknown")
        _apply_field_resets(StdI, "unknown")  # should be a no-op

    def test_unknown_solver_leaves_fields_unchanged(self):
        """Test that unknown solver does not modify any fields."""
        StdI = _make_stdi("unknown")
        original_h = StdI.h
        _apply_field_resets(StdI, "unknown")
        assert StdI.h == original_h


# -------------------------------------------------------------------
#  Common reset table — structure tests
# -------------------------------------------------------------------


class TestCommonResetTableStructure:
    """Tests for _COMMON_RESET_SCALARS and _COMMON_RESET_ARRAYS tables."""

    def test_scalars_is_nonempty_list(self):
        """Test that common scalars table is a non-empty list."""
        assert isinstance(_COMMON_RESET_SCALARS, list)
        assert len(_COMMON_RESET_SCALARS) > 0

    def test_arrays_is_nonempty_list(self):
        """Test that common arrays table is a non-empty list."""
        assert isinstance(_COMMON_RESET_ARRAYS, list)
        assert len(_COMMON_RESET_ARRAYS) > 0

    def test_scalar_entries_are_name_value_tuples(self):
        """Test that every common scalar entry is a (str, value) tuple."""
        for entry in _COMMON_RESET_SCALARS:
            assert isinstance(entry, tuple), entry
            assert len(entry) == 2, entry
            assert isinstance(entry[0], str), entry

    def test_array_entries_are_name_value_tuples(self):
        """Test that every common array entry is a (str, value) tuple."""
        for entry in _COMMON_RESET_ARRAYS:
            assert isinstance(entry, tuple), entry
            assert len(entry) == 2, entry
            assert isinstance(entry[0], str), entry

    def test_scalar_fields_exist_on_stdi(self):
        """Test that every scalar field name exists as a StdIntList attribute."""
        StdI = StdIntList()
        for name, _ in _COMMON_RESET_SCALARS:
            assert hasattr(StdI, name), f"StdIntList missing field: {name}"

    def test_array_fields_exist_on_stdi(self):
        """Test that every array field name exists as a StdIntList attribute."""
        StdI = StdIntList()
        for name, _ in _COMMON_RESET_ARRAYS:
            assert hasattr(StdI, name), f"StdIntList missing field: {name}"

    def test_no_duplicate_scalar_names(self):
        """Test that scalar table has no duplicate field names."""
        names = [name for name, _ in _COMMON_RESET_SCALARS]
        assert len(names) == len(set(names))

    def test_no_duplicate_array_names(self):
        """Test that array table has no duplicate field names."""
        names = [name for name, _ in _COMMON_RESET_ARRAYS]
        assert len(names) == len(set(names))

    def test_scalar_values_are_sentinels(self):
        """Scalar reset values are None (int scalars), NaN_d, or NaN_c."""
        for name, value in _COMMON_RESET_SCALARS:
            if value is None:
                continue  # A2: integer scalars are unset as None
            elif isinstance(value, complex):
                assert math.isnan(value.real), f"{name}: {value} not NaN_c"
            elif isinstance(value, float):
                assert math.isnan(value), f"{name}: {value} not NaN_d"
            else:
                pytest.fail(f"{name}: unexpected type {type(value)}")

    def test_array_values_are_sentinels(self):
        """Array fill values are NaN_d (float matrices) or NaN_i (int matrices)."""
        for name, value in _COMMON_RESET_ARRAYS:
            if isinstance(value, float):
                assert math.isnan(value), f"{name}: {value} not NaN_d"
            elif isinstance(value, int):
                assert value == NaN_i, f"{name}: {value} not NaN_i"
            else:
                pytest.fail(f"{name}: unexpected type {type(value)}")


# -------------------------------------------------------------------
#  Common reset table — content checks
# -------------------------------------------------------------------


class TestCommonResetTableContent:
    """Tests for specific entries in the common reset tables."""

    def test_hopping_parameters_present(self):
        """Test that all 12 hopping parameters are in scalars table."""
        names = {name for name, _ in _COMMON_RESET_SCALARS}
        for base in ("t", "tp", "tpp", "t0", "t0p", "t0pp",
                     "t1", "t1p", "t1pp", "t2", "t2p", "t2pp"):
            assert base in names, f"missing hopping: {base}"

    def test_coulomb_parameters_present(self):
        """Test that all 13 Coulomb parameters are in scalars table."""
        names = {name for name, _ in _COMMON_RESET_SCALARS}
        for base in ("U", "V", "Vp", "Vpp", "V0", "V0p", "V0pp",
                     "V1", "V1p", "V1pp", "V2", "V2p", "V2pp"):
            assert base in names, f"missing Coulomb: {base}"

    def test_spin_coupling_scalars_present(self):
        """Test that all 12 J*All coupling scalars are in table."""
        names = {name for name, _ in _COMMON_RESET_SCALARS}
        for base in ("JAll", "JpAll", "JppAll", "J0All", "J0pAll", "J0ppAll",
                     "J1All", "J1pAll", "J1ppAll", "J2All", "J2pAll", "J2ppAll"):
            assert base in names, f"missing coupling: {base}"

    def test_spin_coupling_matrices_present(self):
        """Test that all 12 J coupling 3x3 matrices are in arrays table."""
        names = {name for name, _ in _COMMON_RESET_ARRAYS}
        for base in ("J", "Jp", "Jpp", "J0", "J0p", "J0pp",
                     "J1", "J1p", "J1pp", "J2", "J2p", "J2pp"):
            assert base in names, f"missing matrix: {base}"

    def test_wannier90_cutoff_scalars_present(self):
        """Test that Wannier90 cutoff scalars are in table."""
        names = {name for name, _ in _COMMON_RESET_SCALARS}
        for field in ("cutoff_t", "cutoff_u", "cutoff_j",
                      "cutoff_length_t", "cutoff_length_U", "cutoff_length_J",
                      "lambda_", "lambda_U", "lambda_J", "alpha"):
            assert field in names, f"missing cutoff scalar: {field}"

    def test_wannier90_cutoff_arrays_present(self):
        """Test that Wannier90 cutoff arrays are in table."""
        names = {name for name, _ in _COMMON_RESET_ARRAYS}
        for field in ("cutoff_tR", "cutoff_UR", "cutoff_JR",
                      "cutoff_tVec", "cutoff_UVec", "cutoff_JVec"):
            assert field in names, f"missing cutoff array: {field}"


# -------------------------------------------------------------------
#  _reset_vals — end-to-end tests
# -------------------------------------------------------------------


class TestResetVals:
    """Tests for _reset_vals using the data-driven tables."""

    def test_common_scalars_reset_to_none(self):
        """A2: all common scalar fields (int/float/complex) reset to None."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        _reset_vals(StdI)
        for name, value in _COMMON_RESET_SCALARS:
            assert value is None, f"{name} table value should be None, got {value}"
            actual = getattr(StdI, name)
            assert actual is None, f"{name} = {actual}, expected None"

    def test_common_arrays_filled(self):
        """Test that common array fields are filled with sentinel values."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        _reset_vals(StdI)
        for name, value in _COMMON_RESET_ARRAYS:
            arr = getattr(StdI, name)
            if isinstance(value, float):
                assert np.all(np.isnan(arr)), f"{name} not all NaN"
            else:
                assert np.all(arr == value), f"{name} not all {value}"

    def test_d_matrix_special_case(self):
        """Test that D matrix is zeros except D[2,2] = NaN."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        _reset_vals(StdI)
        for i in range(3):
            for j in range(3):
                if i == 2 and j == 2:
                    assert math.isnan(StdI.D[i, j])
                else:
                    assert StdI.D[i, j] == 0.0

    def test_lboost_is_zero(self):
        """Test that lBoost is set to 0 (not NaN)."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        _reset_vals(StdI)
        assert StdI.lBoost == 0

    def test_solver_fields_also_set(self):
        """Test that solver-specific fields are also set via _apply_field_resets."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        _reset_vals(StdI)
        assert StdI.LargeValue is None
        assert StdI.FlgTemp == 1
