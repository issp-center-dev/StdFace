"""Unit tests for input_params module.

Tests for the input parameter resolution helpers: ``input_spin_nn``,
``input_spin``, ``input_coulomb_v``, and ``input_hopp``.
"""
from __future__ import annotations


import numpy as np
import pytest

from stdface.lattice.input_params import (
    input_spin_nn,
    input_spin,
    input_coulomb_v,
    input_hopp,
    _has_set_elements,
    _first_set_index,
    _check_scalar_vs_matrix,
    _check_matrix_vs_matrix,
    _resolve_spin_matrix,
)


NaN_d = float("nan")
NaN_c = complex(NaN_d, 0.0)


def _nan_3x3() -> np.ndarray:
    """Return a 3x3 NaN matrix."""
    return np.full((3, 3), NaN_d)


# ===================================================================
#  input_spin_nn
# ===================================================================


class TestInputSpinNN:
    """Tests for input_spin_nn."""

    def test_all_nan_gives_zeros(self):
        """Test that all-NaN inputs produce a zero J0 matrix."""
        J = _nan_3x3()
        J0 = _nan_3x3()
        input_spin_nn(J, NaN_d, J0, NaN_d, "J0")
        assert np.allclose(J0, 0.0)

    def test_j0all_sets_diagonal(self):
        """Test that J0All sets diagonal elements of J0."""
        J = _nan_3x3()
        J0 = _nan_3x3()
        input_spin_nn(J, NaN_d, J0, 1.5, "J0")
        for i in range(3):
            assert J0[i, i] == 1.5
        # Off-diagonal should be 0
        assert J0[0, 1] == 0.0
        assert J0[1, 0] == 0.0

    def test_jall_sets_diagonal(self):
        """Test that JAll sets diagonal elements when J0All is NaN."""
        J = _nan_3x3()
        J0 = _nan_3x3()
        input_spin_nn(J, 2.0, J0, NaN_d, "J0")
        for i in range(3):
            assert J0[i, i] == 2.0

    def test_j_matrix_copies_to_j0(self):
        """Test that J matrix elements are copied to J0."""
        J = _nan_3x3()
        J[0, 0] = 3.0
        J[1, 2] = 0.5
        J0 = _nan_3x3()
        input_spin_nn(J, NaN_d, J0, NaN_d, "J0")
        assert J0[0, 0] == 3.0
        assert J0[1, 2] == 0.5

    def test_j0_explicit_takes_priority(self):
        """Test that explicit J0 elements are kept over J and JAll."""
        J = _nan_3x3()
        J0 = _nan_3x3()
        J0[0, 0] = 5.0
        input_spin_nn(J, NaN_d, J0, NaN_d, "J0")
        assert J0[0, 0] == 5.0

    def test_jall_j0all_conflict_raises(self):
        """Test that both JAll and J0All being set raises ValueError."""
        J = _nan_3x3()
        J0 = _nan_3x3()
        with pytest.raises(ValueError):
            input_spin_nn(J, 1.0, J0, 2.0, "J0")

    def test_jall_j_element_conflict_raises(self):
        """Test that JAll set with J element set raises ValueError."""
        J = _nan_3x3()
        J[0, 0] = 1.0
        J0 = _nan_3x3()
        with pytest.raises(ValueError):
            input_spin_nn(J, 2.0, J0, NaN_d, "J0")

    def test_j0_element_j_element_conflict_raises(self):
        """Test that J0 element and J element both set raises ValueError."""
        J = _nan_3x3()
        J[0, 0] = 1.0
        J0 = _nan_3x3()
        J0[1, 1] = 2.0
        with pytest.raises(ValueError):
            input_spin_nn(J, NaN_d, J0, NaN_d, "J0")


# ===================================================================
#  input_spin
# ===================================================================


class TestInputSpin:
    """Tests for input_spin."""

    def test_all_nan_gives_zeros(self):
        """Test that all-NaN inputs produce a zero matrix."""
        Jp = _nan_3x3()
        input_spin(Jp, NaN_d, "J'")
        assert np.allclose(Jp, 0.0)

    def test_jpall_sets_diagonal(self):
        """Test that JpAll sets diagonal elements."""
        Jp = _nan_3x3()
        input_spin(Jp, 1.0, "J'")
        for i in range(3):
            assert Jp[i, i] == 1.0
        assert Jp[0, 1] == 0.0

    def test_explicit_elements_kept(self):
        """Test that explicitly set elements are preserved."""
        Jp = _nan_3x3()
        Jp[0, 1] = 0.3
        input_spin(Jp, NaN_d, "J'")
        assert Jp[0, 1] == 0.3

    def test_jpall_jp_element_conflict_raises(self):
        """Test that JpAll and Jp element conflict raises ValueError."""
        Jp = _nan_3x3()
        Jp[0, 0] = 1.0
        with pytest.raises(ValueError):
            input_spin(Jp, 2.0, "J'")


# ===================================================================
#  input_coulomb_v
# ===================================================================


class TestInputCoulombV:
    """Tests for input_coulomb_v."""

    def test_both_nan_gives_zero(self):
        """Test that both NaN returns 0."""
        result = input_coulomb_v(NaN_d, NaN_d, "V0")
        assert result == 0.0

    def test_v0_explicit(self):
        """Test that explicit V0 is returned."""
        result = input_coulomb_v(NaN_d, 1.5, "V0")
        assert result == 1.5

    def test_v_fallback(self):
        """Test that V is used when V0 is NaN."""
        result = input_coulomb_v(2.0, NaN_d, "V0")
        assert result == 2.0

    def test_conflict_raises(self):
        """Test that both V and V0 set raises ValueError."""
        with pytest.raises(ValueError):
            input_coulomb_v(1.0, 2.0, "V0")


# ===================================================================
#  input_hopp
# ===================================================================


class TestInputHopp:
    """Tests for input_hopp."""

    def test_both_nan_gives_zero(self):
        """Test that both NaN returns 0+0j."""
        result = input_hopp(NaN_c, NaN_c, "t0")
        assert result == 0.0 + 0j

    def test_t0_explicit(self):
        """Test that explicit t0 is returned."""
        result = input_hopp(NaN_c, complex(1.0, 0.5), "t0")
        assert result == complex(1.0, 0.5)

    def test_t_fallback(self):
        """Test that t is used when t0 is NaN."""
        result = input_hopp(complex(2.0, 0.0), NaN_c, "t0")
        assert result == complex(2.0, 0.0)

    def test_conflict_raises(self):
        """Test that both t and t0 set raises ValueError."""
        with pytest.raises(ValueError):
            input_hopp(complex(1.0, 0.0), complex(2.0, 0.0), "t0")


# ===================================================================
#  Backward compatibility
# ===================================================================


class TestHasSetElements:
    """Tests for the _has_set_elements vectorised helper."""

    def test_all_nan_returns_false(self):
        """Test that an all-NaN matrix returns False."""
        assert _has_set_elements(_nan_3x3()) is False

    def test_one_set_returns_true(self):
        """Test that a matrix with one non-NaN element returns True."""
        m = _nan_3x3()
        m[1, 2] = 0.5
        assert _has_set_elements(m) is True

    def test_all_set_returns_true(self):
        """Test that a fully set matrix returns True."""
        assert _has_set_elements(np.zeros((3, 3))) is True

    def test_zero_is_set(self):
        """Test that 0.0 counts as a set element."""
        m = _nan_3x3()
        m[0, 0] = 0.0
        assert _has_set_elements(m) is True


# ===================================================================
#  _first_set_index
# ===================================================================


class TestFirstSetIndex:
    """Tests for the _first_set_index vectorised helper."""

    def test_all_nan_returns_none(self):
        """Test that an all-NaN matrix returns None."""
        assert _first_set_index(_nan_3x3()) is None

    def test_single_element(self):
        """Test that the index of a single set element is returned."""
        m = _nan_3x3()
        m[2, 1] = 1.0
        assert _first_set_index(m) == (2, 1)

    def test_first_of_multiple(self):
        """Test that the first (row-major) set element index is returned."""
        m = _nan_3x3()
        m[1, 0] = 1.0
        m[2, 2] = 2.0
        assert _first_set_index(m) == (1, 0)

    def test_top_left_corner(self):
        """Test that (0,0) is returned when it is set."""
        m = _nan_3x3()
        m[0, 0] = 0.0
        assert _first_set_index(m) == (0, 0)


# ===================================================================
#  _check_scalar_vs_matrix
# ===================================================================


class TestCheckScalarVsMatrix:
    """Tests for the _check_scalar_vs_matrix conflict detector."""

    def test_nan_scalar_no_error(self):
        """Test that NaN scalar never triggers conflict."""
        m = np.zeros((3, 3))
        _check_scalar_vs_matrix(NaN_d, m, "J", "J")  # no error

    def test_scalar_vs_nan_matrix_no_error(self):
        """Test that set scalar with all-NaN matrix is fine."""
        _check_scalar_vs_matrix(1.0, _nan_3x3(), "J", "J")  # no error

    def test_scalar_vs_set_matrix_raises(self):
        """Test that set scalar with set matrix element aborts."""
        m = _nan_3x3()
        m[0, 0] = 1.0
        with pytest.raises(ValueError):
            _check_scalar_vs_matrix(2.0, m, "J", "J")


# ===================================================================
#  _check_matrix_vs_matrix
# ===================================================================


class TestCheckMatrixVsMatrix:
    """Tests for the _check_matrix_vs_matrix conflict detector."""

    def test_both_nan_no_error(self):
        """Test that two all-NaN matrices don't conflict."""
        _check_matrix_vs_matrix(_nan_3x3(), _nan_3x3(), "J0", "J")

    def test_one_nan_no_error(self):
        """Test that one set matrix with one NaN matrix is fine."""
        m = _nan_3x3()
        m[0, 0] = 1.0
        _check_matrix_vs_matrix(m, _nan_3x3(), "J0", "J")

    def test_both_set_raises(self):
        """Test that two matrices with set elements aborts."""
        a = _nan_3x3(); a[0, 0] = 1.0
        b = _nan_3x3(); b[1, 1] = 2.0
        with pytest.raises(ValueError):
            _check_matrix_vs_matrix(a, b, "J0", "J")


# ===================================================================
#  _resolve_spin_matrix
# ===================================================================


class TestResolveSpinMatrix:
    """Tests for the _resolve_spin_matrix value cascade."""

    def test_all_nan_gives_zeros(self):
        """Test that all-NaN inputs produce zeros."""
        J0 = _nan_3x3()
        _resolve_spin_matrix(J0, "J0")
        assert np.allclose(J0, 0.0)

    def test_j0all_sets_diagonal(self):
        """Test that J0All sets diagonal elements."""
        J0 = _nan_3x3()
        _resolve_spin_matrix(J0, "J0", J0All=1.5)
        for i in range(3):
            assert J0[i, i] == 1.5
        assert J0[0, 1] == 0.0

    def test_jall_sets_diagonal(self):
        """Test that JAll sets diagonal when J0All is NaN."""
        J0 = _nan_3x3()
        _resolve_spin_matrix(J0, "J0", JAll=2.0)
        for i in range(3):
            assert J0[i, i] == 2.0

    def test_j0all_takes_priority_over_jall(self):
        """Test that J0All is used before JAll for diagonal."""
        J0 = _nan_3x3()
        _resolve_spin_matrix(J0, "J0", J0All=1.0, JAll=2.0)
        for i in range(3):
            assert J0[i, i] == 1.0

    def test_j_fallback_copies(self):
        """Test that J matrix elements are copied as fallback."""
        J = _nan_3x3()
        J[0, 1] = 0.7
        J0 = _nan_3x3()
        _resolve_spin_matrix(J0, "J0", J=J)
        assert J0[0, 1] == 0.7

    def test_explicit_j0_preserved(self):
        """Test that already-set J0 elements are not overwritten."""
        J0 = _nan_3x3()
        J0[0, 0] = 5.0
        _resolve_spin_matrix(J0, "J0", J0All=1.0)
        assert J0[0, 0] == 5.0
