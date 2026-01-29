"""Unit tests for stdface_model_util module.

Tests for the Python translation of StdFace_ModelUtil.c/.h.
"""
from __future__ import annotations

import math
import io

import numpy as np
import pytest

from stdface_vals import StdIntList
import stdface_model_util as smu


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def _make_allocated(ntrans: int = 100, nintr: int = 100) -> StdIntList:
    """Return an StdIntList with interaction arrays pre-allocated."""
    s = StdIntList()
    smu.malloc_interactions(s, ntrans, nintr)
    return s


# ---------------------------------------------------------------------------
#  exit_program
# ---------------------------------------------------------------------------


class TestExitProgram:
    """Tests for exit_program."""

    def test_exit_raises(self):
        """exit_program should call sys.exit."""
        with pytest.raises(SystemExit) as exc_info:
            smu.exit_program(42)
        assert exc_info.value.code == 42

    def test_exit_negative(self):
        """exit_program with negative code."""
        with pytest.raises(SystemExit) as exc_info:
            smu.exit_program(-1)
        assert exc_info.value.code == -1


# ---------------------------------------------------------------------------
#  trans
# ---------------------------------------------------------------------------


class TestTrans:
    """Tests for the trans function."""

    def test_adds_entry(self):
        """trans should add an entry to trans/transindx arrays."""
        s = _make_allocated()
        smu.trans(s, 1.0 + 0.5j, 0, 0, 1, 1)
        assert s.ntrans == 1
        assert s.trans[0] == pytest.approx(1.0 + 0.5j)
        assert list(s.transindx[0]) == [0, 0, 1, 1]

    def test_skip_tiny(self):
        """trans should skip entries with |trans0| < 1e-12."""
        s = _make_allocated()
        smu.trans(s, 1e-13, 0, 0, 1, 1)
        assert s.ntrans == 0

    def test_multiple_entries(self):
        """trans should increment ntrans correctly."""
        s = _make_allocated()
        smu.trans(s, 1.0, 0, 0, 1, 0)
        smu.trans(s, 2.0, 1, 1, 0, 1)
        assert s.ntrans == 2
        assert s.trans[0] == 1.0
        assert s.trans[1] == 2.0


# ---------------------------------------------------------------------------
#  hopping
# ---------------------------------------------------------------------------


class TestHopping:
    """Tests for the hopping function."""

    def test_normal_hopping(self):
        """Hopping should add 4 transfer entries (2 spins × 2 directions)."""
        s = _make_allocated()
        dR = np.zeros(3)
        smu.hopping(s, 1.0 + 0j, 0, 1, dR)
        assert s.ntrans == 4

    def test_zero_hopping(self):
        """Zero hopping should add no entries."""
        s = _make_allocated()
        dR = np.zeros(3)
        smu.hopping(s, 0.0, 0, 1, dR)
        assert s.ntrans == 0


# ---------------------------------------------------------------------------
#  hubbard_local
# ---------------------------------------------------------------------------


class TestHubbardLocal:
    """Tests for hubbard_local."""

    def test_basic_terms(self):
        """hubbard_local should add transfer and intra-Coulomb terms."""
        s = _make_allocated()
        smu.hubbard_local(s, mu0=1.0, h0=0.0, Gamma0=0.0,
                          Gamma0_y=0.0, U0=4.0, isite=0)
        # mu contributes 2 transfers (spin up, spin down)
        assert s.ntrans == 2
        assert s.NCintra == 1
        assert s.Cintra[0] == 4.0
        assert s.CintraIndx[0, 0] == 0

    def test_with_magnetic_field(self):
        """hubbard_local with h and Gamma should add more transfers."""
        s = _make_allocated()
        smu.hubbard_local(s, mu0=0.0, h0=1.0, Gamma0=0.5,
                          Gamma0_y=0.3, U0=0.0, isite=2)
        # h: 2 transfers, Gamma: 2, Gamma_y: 2 = 6 total
        assert s.ntrans == 6
        assert s.NCintra == 1


# ---------------------------------------------------------------------------
#  mag_field
# ---------------------------------------------------------------------------


class TestMagField:
    """Tests for mag_field."""

    def test_spin_half(self):
        """S=1/2 (S2=1) should produce longitudinal and transverse terms."""
        s = _make_allocated()
        smu.mag_field(s, S2=1, h=1.0, Gamma=0.5, Gamma_y=0.0, isite=0)
        # S2=1: ispin=0 (Sz=0.5 -> 1 longitudinal), ispin=1 (Sz=-0.5 -> 1 longitudinal + 2 transverse)
        assert s.ntrans == 4

    def test_zero_field(self):
        """Zero field should produce no transfers."""
        s = _make_allocated()
        smu.mag_field(s, S2=1, h=0.0, Gamma=0.0, Gamma_y=0.0, isite=0)
        assert s.ntrans == 0

    def test_spin_one(self):
        """S=1 (S2=2) should produce more terms."""
        s = _make_allocated()
        smu.mag_field(s, S2=2, h=1.0, Gamma=0.5, Gamma_y=0.0, isite=0)
        assert s.ntrans > 0


# ---------------------------------------------------------------------------
#  intr
# ---------------------------------------------------------------------------


class TestIntr:
    """Tests for the intr function."""

    def test_adds_interaction(self):
        """intr should add an interaction entry."""
        s = _make_allocated()
        smu.intr(s, 1.0 + 0j, 0, 0, 0, 1, 1, 0, 1, 1)
        assert s.nintr == 1
        assert s.intr[0] == 1.0
        assert list(s.intrindx[0]) == [0, 0, 0, 1, 1, 0, 1, 1]

    def test_skip_tiny(self):
        """intr should skip entries with |intr0| < 1e-12."""
        s = _make_allocated()
        smu.intr(s, 1e-13, 0, 0, 0, 1, 1, 0, 1, 1)
        assert s.nintr == 0


# ---------------------------------------------------------------------------
#  coulomb
# ---------------------------------------------------------------------------


class TestCoulomb:
    """Tests for the coulomb function."""

    def test_adds_coulomb(self):
        """coulomb should add a Coulomb interaction."""
        s = _make_allocated()
        smu.coulomb(s, V=2.5, isite=0, jsite=1)
        assert s.NCinter == 1
        assert s.Cinter[0] == 2.5
        assert s.CinterIndx[0, 0] == 0
        assert s.CinterIndx[0, 1] == 1


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
        smu.general_j(s, J, Si2=1, Sj2=1, isite=0, jsite=1)
        assert s.NHund == 1
        assert s.NCinter == 1
        assert s.NEx == 1
        assert s.NPairLift == 1

    def test_spin_half_off_diagonal(self):
        """S=1/2 with off-diagonal J should use InterAll."""
        s = _make_allocated()
        J = np.zeros((3, 3))
        J[0, 0] = 1.0
        J[0, 1] = 0.5  # off-diagonal
        J[1, 1] = 1.0
        J[2, 2] = 1.0
        smu.general_j(s, J, Si2=1, Sj2=1, isite=0, jsite=1)
        assert s.NHund == 1
        assert s.NCinter == 1
        # Off-diagonal means ExGeneral stays 1, so nintr is used
        assert s.nintr > 0


# ---------------------------------------------------------------------------
#  PrintVal functions
# ---------------------------------------------------------------------------


class TestPrintValD:
    """Tests for print_val_d."""

    def test_nan_uses_default(self, capsys):
        """NaN input should return default and print DEFAULT tag."""
        result = smu.print_val_d("mu", float("nan"), 0.5)
        assert result == 0.5
        captured = capsys.readouterr()
        assert "DEFAULT VALUE IS USED" in captured.out

    def test_specified_value(self, capsys):
        """Non-NaN input should be returned as-is."""
        result = smu.print_val_d("mu", 1.5, 0.5)
        assert result == 1.5
        captured = capsys.readouterr()
        assert "DEFAULT" not in captured.out


class TestPrintValDd:
    """Tests for print_val_dd."""

    def test_nan_with_primary_default(self, capsys):
        """When val=NaN and val0 is specified, use val0."""
        result = smu.print_val_dd("V", float("nan"), 2.0, 0.0)
        assert result == 2.0

    def test_nan_with_secondary_default(self, capsys):
        """When val=NaN and val0=NaN, use val1."""
        result = smu.print_val_dd("V", float("nan"), float("nan"), 3.0)
        assert result == 3.0


class TestPrintValC:
    """Tests for print_val_c."""

    def test_nan_uses_default(self, capsys):
        """NaN real part should trigger default."""
        result = smu.print_val_c("t", complex(float("nan"), 0), 1.0 + 0.5j)
        assert result == 1.0 + 0.5j

    def test_specified_value(self, capsys):
        """Non-NaN should be returned as-is."""
        result = smu.print_val_c("t", 2.0 + 1.0j, 0.0 + 0j)
        assert result == 2.0 + 1.0j


class TestPrintValI:
    """Tests for print_val_i."""

    def test_sentinel_uses_default(self, capsys):
        """Sentinel value should trigger default."""
        result = smu.print_val_i("L", 2147483647, 4)
        assert result == 4

    def test_specified_value(self, capsys):
        """Non-sentinel should be returned as-is."""
        result = smu.print_val_i("L", 8, 4)
        assert result == 8


# ---------------------------------------------------------------------------
#  NotUsed functions
# ---------------------------------------------------------------------------


class TestNotUsed:
    """Tests for not_used_d, not_used_j, not_used_i."""

    def test_not_used_d_nan_ok(self):
        """NaN value should not trigger exit."""
        smu.not_used_d("x", float("nan"))

    def test_not_used_d_specified_exits(self):
        """Specified value should trigger exit."""
        with pytest.raises(SystemExit):
            smu.not_used_d("x", 1.0)

    def test_not_used_d_complex_nan_ok(self):
        """Complex NaN real part should not trigger exit via not_used_d."""
        smu.not_used_d("t", complex(float("nan"), 0))

    def test_not_used_d_complex_specified_exits(self):
        """Complex specified value should trigger exit via not_used_d."""
        with pytest.raises(SystemExit):
            smu.not_used_d("t", 1.0 + 0j)

    def test_not_used_i_sentinel_ok(self):
        """Sentinel value should not trigger exit."""
        smu.not_used_i("W", 2147483647)

    def test_not_used_i_specified_exits(self):
        """Specified value should trigger exit."""
        with pytest.raises(SystemExit):
            smu.not_used_i("W", 5)

    def test_not_used_j_all_nan_ok(self):
        """All NaN values should not trigger exit."""
        J = np.full((3, 3), float("nan"))
        smu.not_used_j("J", float("nan"), J)

    def test_not_used_j_specified_exits(self):
        """Specified scalar should trigger exit."""
        J = np.full((3, 3), float("nan"))
        with pytest.raises(SystemExit):
            smu.not_used_j("J", 1.0, J)


# ---------------------------------------------------------------------------
#  required_val_i
# ---------------------------------------------------------------------------


class TestRequiredValI:
    """Tests for required_val_i."""

    def test_missing_exits(self):
        """Sentinel value should trigger exit."""
        with pytest.raises(SystemExit):
            smu.required_val_i("nsite", 2147483647)

    def test_specified_ok(self, capsys):
        """Non-sentinel should just print the value."""
        smu.required_val_i("nsite", 16)
        captured = capsys.readouterr()
        assert "16" in captured.out


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
        nBox, fold = smu._fold_site(s, [0, 0, 0])
        assert nBox == [0, 0, 0]
        assert fold == [0, 0, 0]

    def test_periodic_fold(self):
        """A site outside should be folded back."""
        s = StdIntList()
        s.NCell = 4
        s.box = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]], dtype=int)
        s.rbox = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 4]], dtype=int)
        nBox, fold = smu._fold_site(s, [2, 0, 0])
        assert nBox[0] == 1


# ---------------------------------------------------------------------------
#  malloc_interactions
# ---------------------------------------------------------------------------


class TestMallocInteractions:
    """Tests for malloc_interactions."""

    def test_arrays_allocated(self):
        """All interaction arrays should be properly allocated."""
        s = StdIntList()
        smu.malloc_interactions(s, ntransMax=50, nintrMax=30)
        assert s.transindx.shape == (50, 4)
        assert s.trans.shape == (50,)
        assert s.ntrans == 0
        assert s.intrindx.shape == (30, 8)
        assert s.intr.shape == (30,)
        assert s.nintr == 0
        assert s.CintraIndx.shape == (30, 1)
        assert s.Cintra.shape == (30,)
        assert s.NCintra == 0
        assert s.CinterIndx.shape == (30, 2)
        assert s.Cinter.shape == (30,)
        assert s.NCinter == 0
        assert s.HundIndx.shape == (30, 2)
        assert s.Hund.shape == (30,)
        assert s.NHund == 0
        assert s.ExIndx.shape == (30, 2)
        assert s.Ex.shape == (30,)
        assert s.NEx == 0
        assert s.PLIndx.shape == (30, 2)
        assert s.PairLift.shape == (30,)
        assert s.NPairLift == 0
        assert s.PHIndx.shape == (30, 2)
        assert s.PairHopp.shape == (30,)
        assert s.NPairHopp == 0

    def test_pump_arrays_hphi(self):
        """HPhi time-evolution mode should allocate pump arrays."""
        s = StdIntList()
        s.solver = "HPhi"
        s.method = "timeevolution"
        s.PumpBody = 1
        s.Lanczos_max = 5
        smu.malloc_interactions(s, ntransMax=20, nintrMax=10)
        assert s.npump.shape == (5,)
        assert s.pumpindx.shape == (5, 20, 4)
        assert s.pump.shape == (5, 20)


# ---------------------------------------------------------------------------
#  input_spin_nn / input_spin
# ---------------------------------------------------------------------------


class TestInputSpinNN:
    """Tests for input_spin_nn."""

    def test_isotropic_fills_diagonal(self, capsys):
        """Isotropic JAll should fill diagonal of J0."""
        J = np.full((3, 3), float("nan"))
        J0 = np.full((3, 3), float("nan"))
        smu.input_spin_nn(J, JAll=1.5, J0=J0, J0All=float("nan"), J0name="J0")
        assert J0[0, 0] == pytest.approx(1.5)
        assert J0[1, 1] == pytest.approx(1.5)
        assert J0[2, 2] == pytest.approx(1.5)
        assert J0[0, 1] == pytest.approx(0.0)

    def test_conflict_exits(self):
        """Conflicting JAll and J0All should exit."""
        J = np.full((3, 3), float("nan"))
        J0 = np.full((3, 3), float("nan"))
        with pytest.raises(SystemExit):
            smu.input_spin_nn(J, JAll=1.0, J0=J0, J0All=2.0, J0name="J0")


class TestInputSpin:
    """Tests for input_spin."""

    def test_isotropic_fills_diagonal(self, capsys):
        """Isotropic JpAll should fill diagonal of Jp."""
        Jp = np.full((3, 3), float("nan"))
        smu.input_spin(Jp, JpAll=0.5, Jpname="J'")
        assert Jp[0, 0] == pytest.approx(0.5)
        assert Jp[1, 1] == pytest.approx(0.5)
        assert Jp[2, 2] == pytest.approx(0.5)
        assert Jp[0, 1] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
#  input_coulomb_v / input_hopp
# ---------------------------------------------------------------------------


class TestInputCoulombV:
    """Tests for input_coulomb_v."""

    def test_specified_v0(self, capsys):
        """When V0 is specified, it should be returned."""
        result = smu.input_coulomb_v(float("nan"), 2.5, "V1")
        assert result == 2.5

    def test_inherit_from_v(self, capsys):
        """When V0=NaN but V is specified, inherit V."""
        result = smu.input_coulomb_v(1.0, float("nan"), "V1")
        assert result == 1.0

    def test_both_nan(self):
        """Both NaN should default to 0."""
        result = smu.input_coulomb_v(float("nan"), float("nan"), "V1")
        assert result == 0.0

    def test_conflict_exits(self):
        """Both specified should exit."""
        with pytest.raises(SystemExit):
            smu.input_coulomb_v(1.0, 2.0, "V1")


class TestInputHopp:
    """Tests for input_hopp."""

    def test_specified_t0(self, capsys):
        """When t0 is specified, it should be returned."""
        result = smu.input_hopp(complex(float("nan"), 0), 0.5 + 0.1j, "t1")
        assert result == 0.5 + 0.1j

    def test_inherit_from_t(self, capsys):
        """When t0=NaN but t is specified, inherit t."""
        result = smu.input_hopp(1.0 + 0j, complex(float("nan"), 0), "t1")
        assert result == 1.0 + 0j

    def test_both_nan(self):
        """Both NaN should default to 0."""
        result = smu.input_hopp(complex(float("nan"), 0),
                                complex(float("nan"), 0), "t1")
        assert result == 0.0 + 0j


# ---------------------------------------------------------------------------
#  print_geometry (smoke test)
# ---------------------------------------------------------------------------


class TestPrintGeometry:
    """Tests for print_geometry."""

    def test_hwave_uhfk_suppressed(self, tmp_path, monkeypatch):
        """HWAVE uhfk mode should not create geometry.dat."""
        monkeypatch.chdir(tmp_path)
        s = StdIntList()
        s.solver = "HWAVE"
        s.calcmode = "uhfk"
        smu.print_geometry(s)
        assert not (tmp_path / "geometry.dat").exists()


# ---------------------------------------------------------------------------
#  print_xsf (smoke test)
# ---------------------------------------------------------------------------


class TestPrintXSF:
    """Tests for print_xsf."""

    def test_creates_xsf(self, tmp_path, monkeypatch):
        """print_xsf should create lattice.xsf."""
        monkeypatch.chdir(tmp_path)
        s = StdIntList()
        s.lattice = "chain"
        s.NCell = 2
        s.NsiteUC = 1
        s.direct = np.eye(3)
        s.box = np.array([[2, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
        s.Cell = np.array([[0, 0, 0], [1, 0, 0]], dtype=int)
        s.tau = np.zeros((1, 3))
        s.length = np.array([1.0, 1.0, 1.0])
        smu.print_xsf(s)
        assert (tmp_path / "lattice.xsf").exists()
        content = (tmp_path / "lattice.xsf").read_text()
        assert "CRYSTAL" in content
        assert "PRIMVEC" in content
        assert "PRIMCOORD" in content
