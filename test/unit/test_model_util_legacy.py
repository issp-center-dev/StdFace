"""Legacy mixed tests over param_check / input_params / interaction_builder / site_util.

These cases predate the module split (they came with the former
``stdface_model_util`` shim, removed in FR-2) and complement the
per-module test files.

Tests for the Python translation of StdFace_ModelUtil.c/.h.
"""
from __future__ import annotations


import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.core.param_check import (
    print_val_c, print_val_d, print_val_dd, print_val_i,
    not_used_d, not_used_i, not_used_j, required_val_i,
)
from stdface.lattice.input_params import (
    input_spin_nn, input_spin, input_coulomb_v, input_hopp,
)
from stdface.lattice.interaction_builder import (
    trans, intr, hopping, hubbard_local, mag_field, general_j, coulomb,
    malloc_interactions,
)
from stdface.lattice.site_util import _fold_site


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def _make_allocated(ntrans: int = 100) -> StdIntList:
    """Return an StdIntList with interaction arrays pre-allocated."""
    s = StdIntList()
    malloc_interactions(s, ntrans)
    return s


# ---------------------------------------------------------------------------
#  trans
# ---------------------------------------------------------------------------


class TestTrans:
    """Tests for the trans function."""

    def test_adds_entry(self):
        """trans should add an entry to trans/transindx arrays."""
        s = _make_allocated()
        trans(s, 1.0 + 0.5j, 0, 0, 1, 1)
        assert len(s.trans_list) == 1
        assert s.trans_list[0][0] == pytest.approx(1.0 + 0.5j)
        assert list(s.trans_list[0][1:]) == [0, 0, 1, 1]

    def test_skip_tiny(self):
        """trans should skip entries with |trans0| < 1e-12."""
        s = _make_allocated()
        trans(s, 1e-13, 0, 0, 1, 1)
        assert len(s.trans_list) == 0

    def test_multiple_entries(self):
        """trans should increment ntrans correctly."""
        s = _make_allocated()
        trans(s, 1.0, 0, 0, 1, 0)
        trans(s, 2.0, 1, 1, 0, 1)
        assert len(s.trans_list) == 2
        assert s.trans_list[0][0] == 1.0
        assert s.trans_list[1][0] == 2.0


# ---------------------------------------------------------------------------
#  hopping
# ---------------------------------------------------------------------------


class TestHopping:
    """Tests for the hopping function."""

    def test_normal_hopping(self):
        """Hopping should add 4 transfer entries (2 spins × 2 directions)."""
        s = _make_allocated()
        dR = np.zeros(3)
        hopping(s, 1.0 + 0j, 0, 1, dR)
        assert len(s.trans_list) == 4

    def test_zero_hopping(self):
        """Zero hopping should add no entries."""
        s = _make_allocated()
        dR = np.zeros(3)
        hopping(s, 0.0, 0, 1, dR)
        assert len(s.trans_list) == 0


# ---------------------------------------------------------------------------
#  hubbard_local
# ---------------------------------------------------------------------------


class TestHubbardLocal:
    """Tests for hubbard_local."""

    def test_basic_terms(self):
        """hubbard_local should add transfer and intra-Coulomb terms."""
        s = _make_allocated()
        hubbard_local(s, mu0=1.0, h0=0.0, Gamma0=0.0,
                          Gamma0_y=0.0, U0=4.0, isite=0)
        # mu contributes 2 transfers (spin up, spin down)
        assert len(s.trans_list) == 2
        assert len(s.Cintra_list) == 1
        assert s.Cintra_list[0][0] == 4.0
        assert s.Cintra_list[0][1] == 0

    def test_with_magnetic_field(self):
        """hubbard_local with h and Gamma should add more transfers."""
        s = _make_allocated()
        hubbard_local(s, mu0=0.0, h0=1.0, Gamma0=0.5,
                          Gamma0_y=0.3, U0=0.0, isite=2)
        # h: 2 transfers, Gamma: 2, Gamma_y: 2 = 6 total
        assert len(s.trans_list) == 6
        assert len(s.Cintra_list) == 1


# ---------------------------------------------------------------------------
#  mag_field
# ---------------------------------------------------------------------------


class TestMagField:
    """Tests for mag_field."""

    def test_spin_half(self):
        """S=1/2 (S2=1) should produce longitudinal and transverse terms."""
        s = _make_allocated()
        mag_field(s, S2=1, h=1.0, Gamma=0.5, Gamma_y=0.0, isite=0)
        # S2=1: ispin=0 (Sz=0.5 -> 1 longitudinal), ispin=1 (Sz=-0.5 -> 1 longitudinal + 2 transverse)
        assert len(s.trans_list) == 4

    def test_zero_field(self):
        """Zero field should produce no transfers."""
        s = _make_allocated()
        mag_field(s, S2=1, h=0.0, Gamma=0.0, Gamma_y=0.0, isite=0)
        assert len(s.trans_list) == 0

    def test_spin_one(self):
        """S=1 (S2=2) should produce more terms."""
        s = _make_allocated()
        mag_field(s, S2=2, h=1.0, Gamma=0.5, Gamma_y=0.0, isite=0)
        assert len(s.trans_list) > 0


# ---------------------------------------------------------------------------
#  intr
# ---------------------------------------------------------------------------


class TestIntr:
    """Tests for the intr function."""

    def test_adds_interaction(self):
        """intr should add an interaction entry."""
        s = _make_allocated()
        intr(s, 1.0 + 0j, 0, 0, 0, 1, 1, 0, 1, 1)
        assert len(s.intr_list) == 1
        assert s.intr_list[0][0] == 1.0
        assert list(s.intr_list[0][1:]) == [0, 0, 0, 1, 1, 0, 1, 1]

    def test_skip_tiny(self):
        """intr should skip entries with |intr0| < 1e-12."""
        s = _make_allocated()
        intr(s, 1e-13, 0, 0, 0, 1, 1, 0, 1, 1)
        assert len(s.intr_list) == 0


# ---------------------------------------------------------------------------
#  coulomb
# ---------------------------------------------------------------------------


class TestCoulomb:
    """Tests for the coulomb function."""

    def test_adds_coulomb(self):
        """coulomb should add a Coulomb interaction."""
        s = _make_allocated()
        coulomb(s, V=2.5, isite=0, jsite=1)
        assert len(s.Cinter_list) == 1
        assert s.Cinter_list[0][0] == 2.5
        assert s.Cinter_list[0][1] == 0
        assert s.Cinter_list[0][2] == 1


# ---------------------------------------------------------------------------
#  general_j
# ---------------------------------------------------------------------------


class TestGeneralJ:
    """Tests for general_j."""

    def test_spin_half_diagonal(self):
        """S=1/2 with diagonal J should set Hund, Cinter, Ex, PairLift."""
        s = _make_allocated()
        J = np.zeros((3, 3))
        J[0, 0] = 1.0
        J[1, 1] = 1.0
        J[2, 2] = 1.0
        general_j(s, J, Si2=1, Sj2=1, isite=0, jsite=1)
        assert len(s.Hund_list) == 1
        assert len(s.Cinter_list) == 1
        assert len(s.Ex_list) == 1
        assert len(s.PairLift_list) == 1

    def test_spin_half_off_diagonal(self):
        """S=1/2 with off-diagonal J should use InterAll."""
        s = _make_allocated()
        J = np.zeros((3, 3))
        J[0, 0] = 1.0
        J[0, 1] = 0.5  # off-diagonal
        J[1, 1] = 1.0
        J[2, 2] = 1.0
        general_j(s, J, Si2=1, Sj2=1, isite=0, jsite=1)
        assert len(s.Hund_list) == 1
        assert len(s.Cinter_list) == 1
        # Off-diagonal means ExGeneral stays 1, so nintr is used
        assert len(s.intr_list) > 0


# ---------------------------------------------------------------------------
#  PrintVal functions
# ---------------------------------------------------------------------------


class TestPrintValD:
    """Tests for print_val_d."""

    def test_nan_uses_default(self, caplog):
        """NaN input should return default and log DEFAULT tag."""
        import logging

        with caplog.at_level(logging.INFO, logger="stdface.core.param_check"):
            result = print_val_d("mu", float("nan"), 0.5)
        assert result == 0.5
        assert "DEFAULT VALUE IS USED" in caplog.text

    def test_specified_value(self, caplog):
        """Non-NaN input should be returned as-is."""
        import logging

        with caplog.at_level(logging.INFO, logger="stdface.core.param_check"):
            result = print_val_d("mu", 1.5, 0.5)
        assert result == 1.5
        assert "DEFAULT" not in caplog.text


class TestPrintValDd:
    """Tests for print_val_dd."""

    def test_nan_with_primary_default(self):
        """When val=NaN and val0 is specified, use val0."""
        result = print_val_dd("V", float("nan"), 2.0, 0.0)
        assert result == 2.0

    def test_nan_with_secondary_default(self):
        """When val=NaN and val0=NaN, use val1."""
        result = print_val_dd("V", float("nan"), float("nan"), 3.0)
        assert result == 3.0


class TestPrintValC:
    """Tests for print_val_c."""

    def test_nan_uses_default(self):
        """NaN real part should trigger default."""
        result = print_val_c("t", complex(float("nan"), 0), 1.0 + 0.5j)
        assert result == 1.0 + 0.5j

    def test_specified_value(self):
        """Non-NaN should be returned as-is."""
        result = print_val_c("t", 2.0 + 1.0j, 0.0 + 0j)
        assert result == 2.0 + 1.0j


class TestPrintValI:
    """Tests for print_val_i."""

    def test_sentinel_uses_default(self):
        """Sentinel value should trigger default."""
        result = print_val_i("L", 2147483647, 4)
        assert result == 4

    def test_specified_value(self):
        """Non-sentinel should be returned as-is."""
        result = print_val_i("L", 8, 4)
        assert result == 8


# ---------------------------------------------------------------------------
#  NotUsed functions
# ---------------------------------------------------------------------------


class TestNotUsed:
    """Tests for not_used_d, not_used_j, not_used_i."""

    def test_not_used_d_nan_ok(self):
        """NaN value should not trigger exit."""
        not_used_d("x", float("nan"))

    def test_not_used_d_specified_exits(self):
        """Specified value should raise ValueError."""
        with pytest.raises(ValueError):
            not_used_d("x", 1.0)

    def test_not_used_d_complex_nan_ok(self):
        """Complex NaN real part should not trigger exit via not_used_d."""
        not_used_d("t", complex(float("nan"), 0))

    def test_not_used_d_complex_specified_exits(self):
        """Complex specified value should raise ValueError via not_used_d."""
        with pytest.raises(ValueError):
            not_used_d("t", 1.0 + 0j)

    def test_not_used_i_sentinel_ok(self):
        """Sentinel value should not trigger exit."""
        not_used_i("W", 2147483647)

    def test_not_used_i_specified_exits(self):
        """Specified value should raise ValueError."""
        with pytest.raises(ValueError):
            not_used_i("W", 5)

    def test_not_used_j_all_nan_ok(self):
        """All NaN values should not trigger exit."""
        J = np.full((3, 3), float("nan"))
        not_used_j("J", float("nan"), J)

    def test_not_used_j_specified_exits(self):
        """Specified scalar should raise ValueError."""
        J = np.full((3, 3), float("nan"))
        with pytest.raises(ValueError):
            not_used_j("J", 1.0, J)


# ---------------------------------------------------------------------------
#  required_val_i
# ---------------------------------------------------------------------------


class TestRequiredValI:
    """Tests for required_val_i."""

    def test_missing_exits(self):
        """Sentinel value should raise ValueError."""
        with pytest.raises(ValueError):
            required_val_i("nsite", 2147483647)

    def test_specified_ok(self, caplog):
        """Non-sentinel should just log the value."""
        import logging

        with caplog.at_level(logging.INFO, logger="stdface.core.param_check"):
            required_val_i("nsite", 16)
        assert "16" in caplog.text


# ---------------------------------------------------------------------------
#  _fold_site
# ---------------------------------------------------------------------------


class TestFoldSite:
    """Tests for _fold_site."""

    def test_identity_fold(self):
        """A site inside the cell should map to itself."""
        s = StdIntList()
        s.NCell = 4
        s.box = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]], dtype=int)
        # rbox = cofactor matrix of box
        s.rbox = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 4]], dtype=int)
        nBox, fold = _fold_site(s, [0, 0, 0])
        assert nBox == [0, 0, 0]
        assert fold == [0, 0, 0]

    def test_periodic_fold(self):
        """A site outside should be folded back."""
        s = StdIntList()
        s.NCell = 4
        s.box = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]], dtype=int)
        s.rbox = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 4]], dtype=int)
        nBox, fold = _fold_site(s, [2, 0, 0])
        assert nBox[0] == 1


# ---------------------------------------------------------------------------
#  malloc_interactions
# ---------------------------------------------------------------------------


class TestMallocInteractions:
    """Tests for malloc_interactions."""

    def test_arrays_allocated(self):
        """All interaction arrays should be properly allocated."""
        s = StdIntList()
        malloc_interactions(s, ntransMax=50)
        assert len(s.trans_list) == 0
        assert len(s.intr_list) == 0
        assert len(s.Cintra_list) == 0
        assert len(s.Cinter_list) == 0
        assert len(s.Hund_list) == 0
        assert len(s.Ex_list) == 0
        assert len(s.PairLift_list) == 0
        assert len(s.PairHopp_list) == 0

    def test_pump_arrays_hphi(self):
        """HPhi time-evolution mode should allocate pump arrays."""
        s = StdIntList()
        s.solver = "HPhi"
        s.method = "timeevolution"
        s.PumpBody = 1
        s.Lanczos_max = 5
        malloc_interactions(s, ntransMax=20)
        assert s.npump.shape == (5,)
        assert s.pumpindx.shape == (5, 20, 4)
        assert s.pump.shape == (5, 20)


# ---------------------------------------------------------------------------
#  input_spin_nn / input_spin
# ---------------------------------------------------------------------------


class TestInputSpinNN:
    """Tests for input_spin_nn."""

    def test_isotropic_fills_diagonal(self):
        """Isotropic JAll should fill diagonal of J0."""
        J = np.full((3, 3), float("nan"))
        J0 = np.full((3, 3), float("nan"))
        input_spin_nn(J, JAll=1.5, J0=J0, J0All=float("nan"), J0name="J0")
        assert J0[0, 0] == pytest.approx(1.5)
        assert J0[1, 1] == pytest.approx(1.5)
        assert J0[2, 2] == pytest.approx(1.5)
        assert J0[0, 1] == pytest.approx(0.0)

    def test_conflict_exits(self):
        """Conflicting JAll and J0All should exit."""
        J = np.full((3, 3), float("nan"))
        J0 = np.full((3, 3), float("nan"))
        with pytest.raises(ValueError):
            input_spin_nn(J, JAll=1.0, J0=J0, J0All=2.0, J0name="J0")


class TestInputSpin:
    """Tests for input_spin."""

    def test_isotropic_fills_diagonal(self):
        """Isotropic JpAll should fill diagonal of Jp."""
        Jp = np.full((3, 3), float("nan"))
        input_spin(Jp, JpAll=0.5, Jpname="J'")
        assert Jp[0, 0] == pytest.approx(0.5)
        assert Jp[1, 1] == pytest.approx(0.5)
        assert Jp[2, 2] == pytest.approx(0.5)
        assert Jp[0, 1] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
#  input_coulomb_v / input_hopp
# ---------------------------------------------------------------------------


class TestInputCoulombV:
    """Tests for input_coulomb_v."""

    def test_specified_v0(self):
        """When V0 is specified, it should be returned."""
        result = input_coulomb_v(float("nan"), 2.5, "V1")
        assert result == 2.5

    def test_inherit_from_v(self):
        """When V0=NaN but V is specified, inherit V."""
        result = input_coulomb_v(1.0, float("nan"), "V1")
        assert result == 1.0

    def test_both_nan(self):
        """Both NaN should default to 0."""
        result = input_coulomb_v(float("nan"), float("nan"), "V1")
        assert result == 0.0

    def test_conflict_exits(self):
        """Both specified should exit."""
        with pytest.raises(ValueError):
            input_coulomb_v(1.0, 2.0, "V1")


class TestInputHopp:
    """Tests for input_hopp."""

    def test_specified_t0(self):
        """When t0 is specified, it should be returned."""
        result = input_hopp(complex(float("nan"), 0), 0.5 + 0.1j, "t1")
        assert result == 0.5 + 0.1j

    def test_inherit_from_t(self):
        """When t0=NaN but t is specified, inherit t."""
        result = input_hopp(1.0 + 0j, complex(float("nan"), 0), "t1")
        assert result == 1.0 + 0j

    def test_both_nan(self):
        """Both NaN should default to 0."""
        result = input_hopp(complex(float("nan"), 0),
                                complex(float("nan"), 0), "t1")
        assert result == 0.0 + 0j


# lattice.xsf / geometry.dat coverage now lives in test_geometry_output
# (build_xsf / build_geometry); the stdface_model_util re-exports of
# print_xsf / print_geometry were removed.
