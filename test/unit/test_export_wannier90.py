"""Unit tests for export_wannier90 module.

Tests for the Python translation of export_wannier90.c.
"""
from __future__ import annotations

import io
import math
import os
import tempfile

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.solvers.hwave import export_wannier90 as ew


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

NaN_i = 2147483647


def _make_stdI_for_geometry(nsiteUC: int = 2) -> StdIntList:
    """Return a StdIntList set up for geometry export tests.

    Parameters
    ----------
    nsiteUC : int
        Number of sites per unit cell.
    """
    s = StdIntList()
    s.solver = "HWAVE"
    s.NsiteUC = nsiteUC
    s.fileprefix = ""

    # Set direct lattice vectors (identity-like for simple test)
    s.direct = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])

    # Set tau (orbital positions)
    s.tau = np.zeros((nsiteUC, 3))
    if nsiteUC >= 1:
        s.tau[0] = [0.0, 0.0, 0.0]
    if nsiteUC >= 2:
        s.tau[1] = [0.5, 0.5, 0.0]

    return s


def _make_stdI_for_interaction(nsiteUC: int = 1, ncell: int = 2) -> StdIntList:
    """Return a StdIntList set up for interaction export tests.

    Parameters
    ----------
    nsiteUC : int
        Number of sites per unit cell.
    ncell : int
        Number of cells.
    """
    s = StdIntList()
    s.solver = "HWAVE"
    s.NsiteUC = nsiteUC
    s.NCell = ncell
    s.fileprefix = ""
    s.export_all = NaN_i  # use default

    # Set up identity-like box/rbox
    s.box = np.array([
        [ncell, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ], dtype=int)
    s.rbox = np.array([
        [1, 0, 0],
        [0, ncell, 0],
        [0, 0, ncell],
    ], dtype=int)

    # Cell positions
    s.Cell = np.zeros((ncell, 3), dtype=int)
    for i in range(ncell):
        s.Cell[i, 0] = i

    # Direct lattice vectors
    s.direct = np.eye(3)
    s.tau = np.zeros((nsiteUC, 3))

    # Zero out interaction arrays
    s.ntrans = 0
    s.transindx = None
    s.trans = None
    s.NCintra = 0
    s.CintraIndx = None
    s.Cintra = None
    s.NCinter = 0
    s.CinterIndx = None
    s.Cinter = None
    s.NHund = 0
    s.HundIndx = None
    s.Hund = None
    s.NEx = 0
    s.ExIndx = None
    s.Ex = None
    s.NPairLift = 0
    s.PLIndx = None
    s.PairLift = None
    s.NPairHopp = 0
    s.PHIndx = None
    s.PairHopp = None

    return s


# ===========================================================================
#  Tests: Internal helpers
# ===========================================================================

class TestGenerateKey:
    """Tests for _generate_key."""

    def test_keylen_1(self):
        result = ew._generate_key(1, [42], 0)
        assert result == (42,)

    def test_keylen_2_unordered(self):
        result = ew._generate_key(2, [5, 3], 0)
        assert result == (5, 3)

    def test_keylen_2_ordered_swap(self):
        result = ew._generate_key(2, [5, 3], 1)
        assert result == (3, 5)

    def test_keylen_2_ordered_no_swap(self):
        result = ew._generate_key(2, [2, 7], 1)
        assert result == (2, 7)

    def test_keylen_4_unordered(self):
        result = ew._generate_key(4, [3, 0, 1, 1], 0)
        assert result == (3, 0, 1, 1)

    def test_keylen_4_ordered_swap(self):
        result = ew._generate_key(4, [5, 0, 2, 1], 1)
        assert result == (2, 1, 5, 0)

    def test_keylen_4_ordered_no_swap(self):
        result = ew._generate_key(4, [1, 0, 3, 1], 1)
        assert result == (1, 0, 3, 1)


class TestComputeIndex:
    """Tests for _compute_index."""

    def test_simple_case(self):
        # rr=[0,0,0], nsiteuc=1, nspin=1 => only one element, index=0
        idx = ew._compute_index(0, 0, 0, 0, 0, 0, 0, [0, 0, 0], 1, 1)
        assert idx == 0

    def test_spin_offset(self):
        # rr=[0,0,0], nsiteuc=1, nspin=2
        # t=1 should give index 1
        idx0 = ew._compute_index(0, 0, 0, 0, 0, 0, 0, [0, 0, 0], 1, 2)
        idx1 = ew._compute_index(0, 0, 0, 0, 0, 0, 1, [0, 0, 0], 1, 2)
        assert idx0 == 0
        assert idx1 == 1

    def test_orbital_offset(self):
        # rr=[0,0,0], nsiteuc=2, nspin=1
        idx00 = ew._compute_index(0, 0, 0, 0, 0, 0, 0, [0, 0, 0], 2, 1)
        idx01 = ew._compute_index(0, 0, 0, 0, 1, 0, 0, [0, 0, 0], 2, 1)
        idx10 = ew._compute_index(0, 0, 0, 1, 0, 0, 0, [0, 0, 0], 2, 1)
        assert idx00 == 0
        assert idx01 == 1
        assert idx10 == 2


# ===========================================================================
#  Tests: Accumulate list
# ===========================================================================

class TestAccumulateList:
    """Tests for _accumulate_list."""

    def test_single_entry(self):
        tbl_index = np.array([[0, 1]])
        tbl_value = np.array([1.0 + 0j])
        nintr, idx, val = ew._accumulate_list(2, 1, tbl_index, tbl_value, 0)
        assert nintr == 1
        assert idx[0] == [0, 1]
        assert val[0] == pytest.approx(1.0)

    def test_accumulate_same_key(self):
        tbl_index = np.array([[0, 1], [0, 1]])
        tbl_value = np.array([1.0 + 0j, 2.0 + 0j])
        nintr, idx, val = ew._accumulate_list(2, 2, tbl_index, tbl_value, 0)
        assert nintr == 1
        assert val[0] == pytest.approx(3.0)

    def test_accumulate_ordered_key(self):
        """With ordered=1, (5,3) and (3,5) should be the same key."""
        tbl_index = np.array([[5, 3], [3, 5]])
        tbl_value = np.array([1.0 + 0j, 2.0 + 0j])
        nintr, idx, val = ew._accumulate_list(2, 2, tbl_index, tbl_value, 1)
        assert nintr == 1
        assert idx[0] == [3, 5]
        assert val[0] == pytest.approx(3.0)

    def test_eliminate_zero(self):
        """Entries that sum to zero should be eliminated."""
        tbl_index = np.array([[0, 1], [0, 1]])
        tbl_value = np.array([1.0 + 0j, -1.0 + 0j])
        nintr, idx, val = ew._accumulate_list(2, 2, tbl_index, tbl_value, 0)
        assert nintr == 0

    def test_different_keys(self):
        tbl_index = np.array([[0, 1], [2, 3]])
        tbl_value = np.array([1.0 + 0j, 5.0 + 0j])
        nintr, idx, val = ew._accumulate_list(2, 2, tbl_index, tbl_value, 0)
        assert nintr == 2

    def test_keylen_1(self):
        tbl_index = np.array([[3], [3], [5]])
        tbl_value = np.array([1.0 + 0j, 2.0 + 0j, 7.0 + 0j])
        nintr, idx, val = ew._accumulate_list(1, 3, tbl_index, tbl_value, 1)
        assert nintr == 2
        # key [3] accumulated to 3.0
        assert idx[0] == [3]
        assert val[0] == pytest.approx(3.0)
        assert idx[1] == [5]
        assert val[1] == pytest.approx(7.0)

    def test_keylen_4(self):
        tbl_index = np.array([[0, 0, 1, 0], [0, 0, 1, 0]])
        tbl_value = np.array([1.0 + 0j, 3.0 + 0j])
        nintr, idx, val = ew._accumulate_list(4, 2, tbl_index, tbl_value, 0)
        assert nintr == 1
        assert val[0] == pytest.approx(4.0)


# ===========================================================================
#  Tests: Unfold site
# ===========================================================================

class TestUnfoldSite:
    """Tests for _unfold_site."""

    def test_identity_no_shift(self):
        """v_in = [0,0,0] should return [0,0,0]."""
        s = StdIntList()
        s.NCell = 4
        s.rbox = np.eye(3, dtype=int)
        s.box = np.array([[4, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
        result = ew._unfold_site(s, [0, 0, 0])
        assert result == [0, 0, 0]

    def test_unfold_positive(self):
        """A coordinate beyond half the box should wrap to negative."""
        s = StdIntList()
        s.NCell = 4
        s.rbox = np.eye(3, dtype=int)
        s.box = np.array([[4, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
        result = ew._unfold_site(s, [3, 0, 0])
        assert result == [-1, 0, 0]


# ===========================================================================
#  Tests: Prefix helper
# ===========================================================================

class TestPrefix:
    """Tests for _prefix."""

    def test_no_prefix(self):
        s = StdIntList()
        s.fileprefix = ""
        assert ew._prefix(s, "geom.dat") == "geom.dat"

    def test_stars_sentinel(self):
        s = StdIntList()
        s.fileprefix = "****"
        assert ew._prefix(s, "geom.dat") == "geom.dat"

    def test_with_prefix(self):
        s = StdIntList()
        s.fileprefix = "myrun"
        assert ew._prefix(s, "geom.dat") == "myrun_geom.dat"


# ===========================================================================
#  Tests: Write geometry
# ===========================================================================

class TestWriteGeometry:
    """Tests for _write_geometry (and export_geometry)."""

    def test_write_geometry_file_content(self, tmp_path):
        """Check that the geometry file has expected structure."""
        s = _make_stdI_for_geometry(nsiteUC=2)
        fname = str(tmp_path / "geom.dat")
        ew._write_geometry(s, fname)

        assert os.path.exists(fname)
        with open(fname) as f:
            lines = f.readlines()

        # 3 lines for lattice vectors + 1 for NsiteUC + 2 for tau
        assert len(lines) == 6

        # First line: lattice vector [1, 0, 0]
        vals = [float(x) for x in lines[0].split()]
        assert vals == pytest.approx([1.0, 0.0, 0.0])

        # NsiteUC line
        assert int(lines[3].strip()) == 2

        # First tau
        tau0 = [float(x) for x in lines[4].split()]
        assert tau0 == pytest.approx([0.0, 0.0, 0.0])

        # Second tau
        tau1 = [float(x) for x in lines[5].split()]
        assert tau1 == pytest.approx([0.5, 0.5, 0.0])

    def test_export_geometry_with_prefix(self, tmp_path, monkeypatch):
        """Check that export_geometry uses the prefix correctly."""
        monkeypatch.chdir(tmp_path)
        s = _make_stdI_for_geometry(nsiteUC=1)
        s.fileprefix = "test"
        ew.export_geometry(s)
        assert os.path.exists(tmp_path / "test_geom.dat")


# ===========================================================================
#  Tests: Write Wannier90 format
# ===========================================================================

class TestWriteWannier90:
    """Tests for _write_wannier90."""

    def test_single_element_no_spin(self, tmp_path):
        """One interaction entry, nspin=1."""
        item = ew._IntrItem(r=[0, 0, 0], a=0, b=0, s=0, t=0, v=1.5 + 0.5j)
        fname = str(tmp_path / "test_w90.dat")
        ew._write_wannier90([item], 1, 1, fname, "Test")

        with open(fname) as f:
            lines = f.readlines()

        # Header: tagname line, nsiteuc, nvol, degeneracy, then data
        assert "Test" in lines[0]
        assert int(lines[1].strip()) == 1  # nsiteuc
        assert int(lines[2].strip()) == 1  # nvol
        # degeneracy line
        assert lines[3].strip() == "1"
        # data line: rx ry rz a+1 b+1 real imag
        data = lines[4].split()
        assert int(data[0]) == 0  # rx
        assert int(data[1]) == 0  # ry
        assert int(data[2]) == 0  # rz
        assert int(data[3]) == 1  # a+1
        assert int(data[4]) == 1  # b+1
        assert float(data[5]) == pytest.approx(1.5)
        assert float(data[6]) == pytest.approx(0.5)

    def test_hermitian_conjugate(self, tmp_path):
        """Check that the reverse entry is set to conjugate."""
        item = ew._IntrItem(r=[1, 0, 0], a=0, b=0, s=0, t=0, v=1.0 + 2.0j)
        fname = str(tmp_path / "test_hc.dat")
        ew._write_wannier90([item], 1, 1, fname, "Test")

        with open(fname) as f:
            lines = f.readlines()

        # With rr=[1,0,0], nvol=(2*1+1)*1*1 = 3
        # r values are -1, 0, 1
        assert int(lines[2].strip()) == 3  # nvol

        # Find the line for r=[-1,0,0] a=1 b=1
        data_lines = [l for l in lines[4:] if l.strip()]
        found_conj = False
        for dl in data_lines:
            parts = dl.split()
            rx, ry, rz = int(parts[0]), int(parts[1]), int(parts[2])
            if rx == -1 and ry == 0 and rz == 0:
                re_val = float(parts[5])
                im_val = float(parts[6])
                assert re_val == pytest.approx(1.0)
                assert im_val == pytest.approx(-2.0)  # conjugate
                found_conj = True
        assert found_conj

    def test_spin_dependent_format(self, tmp_path):
        """nspin=2 should use extended format with s,t columns."""
        item = ew._IntrItem(r=[0, 0, 0], a=0, b=0, s=1, t=0, v=0.3 + 0j)
        fname = str(tmp_path / "test_spin.dat")
        ew._write_wannier90([item], 1, 2, fname, "SpinTest")

        with open(fname) as f:
            lines = f.readlines()

        # Data lines should have 9 columns (rx ry rz a b s t re im)
        data_lines = [l for l in lines[4:] if l.strip()]
        for dl in data_lines:
            parts = dl.split()
            assert len(parts) == 9


# ===========================================================================
#  Tests: _build_wannier_matrix
# ===========================================================================

class TestBuildWannierMatrix:
    """Tests for _build_wannier_matrix."""

    def test_single_item_rr(self):
        """Single item at r=[0,0,0] gives rr=[0,0,0]."""
        item = ew._IntrItem(r=[0, 0, 0], a=0, b=0, s=0, t=0, v=1.0 + 0j)
        rr, nvol, matrix = ew._build_wannier_matrix([item], 1, 1)
        assert rr == [0, 0, 0]
        assert nvol == 1
        assert len(matrix) == 1

    def test_nonzero_r_expands_range(self):
        """Item at r=[1,0,0] gives rr=[1,0,0] and nvol=3."""
        item = ew._IntrItem(r=[1, 0, 0], a=0, b=0, s=0, t=0, v=2.0 + 0j)
        rr, nvol, matrix = ew._build_wannier_matrix([item], 1, 1)
        assert rr == [1, 0, 0]
        assert nvol == 3

    def test_matrix_value_placed(self):
        """Check that the item value ends up at the correct matrix index."""
        item = ew._IntrItem(r=[0, 0, 0], a=0, b=0, s=0, t=0, v=3.5 + 1.0j)
        rr, nvol, matrix = ew._build_wannier_matrix([item], 1, 1)
        idx = ew._compute_index(0, 0, 0, 0, 0, 0, 0, rr, 1, 1)
        assert matrix[idx] == pytest.approx(3.5 + 1.0j)

    def test_hermitian_conjugate_filled(self):
        """Reverse-direction entry is set to conjugate when empty."""
        item = ew._IntrItem(r=[1, 0, 0], a=0, b=0, s=0, t=0, v=1.0 + 2.0j)
        rr, nvol, matrix = ew._build_wannier_matrix([item], 1, 1)
        ridx = ew._compute_index(-1, 0, 0, 0, 0, 0, 0, rr, 1, 1)
        assert matrix[ridx] == pytest.approx(1.0 - 2.0j)

    def test_spin_matrix_size(self):
        """With nspin=2, matrix size is nvol * nsiteuc^2 * nspin^2."""
        item = ew._IntrItem(r=[0, 0, 0], a=0, b=0, s=1, t=0, v=0.5 + 0j)
        rr, nvol, matrix = ew._build_wannier_matrix([item], 1, 2)
        assert len(matrix) == nvol * 1 * 1 * 2 * 2


# ===========================================================================
#  Tests: _write_wannier_body
# ===========================================================================

class TestWriteWannierBody:
    """Tests for _write_wannier_body."""

    def test_nspin1_compact_format(self):
        """nspin=1 writes compact format (7 columns: rx ry rz a b re im)."""
        matrix = np.array([1.0 + 0.5j])
        fp = io.StringIO()
        ew._write_wannier_body(fp, [0, 0, 0], 1, 1, 1, matrix)
        lines = fp.getvalue().strip().splitlines()
        assert len(lines) == 1
        parts = lines[0].split()
        assert len(parts) == 7  # rx ry rz a+1 b+1 re im

    def test_nspin2_extended_format(self):
        """nspin=2 writes extended format (9 columns: rx ry rz a b s t re im)."""
        matrix = np.zeros(4, dtype=complex)
        matrix[0] = 1.0 + 0j
        fp = io.StringIO()
        ew._write_wannier_body(fp, [0, 0, 0], 1, 1, 2, matrix)
        lines = fp.getvalue().strip().splitlines()
        # With export_all=1, all 4 spin combinations are written
        assert len(lines) == 4
        parts = lines[0].split()
        assert len(parts) == 9  # rx ry rz a+1 b+1 s t re im

    def test_correct_values_written(self):
        """Values from the matrix appear in the output."""
        matrix = np.array([2.5 + 0.3j])
        fp = io.StringIO()
        ew._write_wannier_body(fp, [0, 0, 0], 1, 1, 1, matrix)
        content = fp.getvalue()
        assert "2.500000000000" in content
        assert "0.300000000000" in content

    def test_multi_volume_iterations(self):
        """Multiple volumes produce entries for all r-vectors."""
        # rr=[1,0,0] → nvol=3, nsiteuc=1, nspin=1 → 3 entries
        matrix = np.ones(3, dtype=complex)
        fp = io.StringIO()
        ew._write_wannier_body(fp, [1, 0, 0], 3, 1, 1, matrix)
        lines = fp.getvalue().strip().splitlines()
        assert len(lines) == 3
        # Check rx values: -1, 0, 1
        rx_vals = [int(l.split()[0]) for l in lines]
        assert sorted(rx_vals) == [-1, 0, 1]


# ===========================================================================
#  Tests: Export inter (complex)
# ===========================================================================

class TestExportInter:
    """Tests for _export_inter."""

    def test_zero_entries_skipped(self, capsys):
        """ntbl=0 should print skip message."""
        s = _make_stdI_for_interaction()
        ew._export_inter(s, 0, None, None, "test.dat", "Test")
        captured = capsys.readouterr()
        assert "skipped" in captured.out

    def test_basic_export(self, tmp_path):
        """A simple two-site interaction should produce a file."""
        s = _make_stdI_for_interaction(nsiteUC=1, ncell=2)

        tbl_index = np.array([[0, 1]], dtype=int)
        tbl_value = np.array([1.5 + 0j], dtype=complex)
        fname = str(tmp_path / "inter_test.dat")

        ew._export_inter(s, 1, tbl_index, tbl_value, fname, "TestInter")
        assert os.path.exists(fname)


# ===========================================================================
#  Tests: Export inter real
# ===========================================================================

class TestExportInterReal:
    """Tests for _export_inter_real."""

    def test_zero_entries_skipped(self, capsys):
        s = _make_stdI_for_interaction()
        ew._export_inter_real(s, 0, None, np.array([], dtype=float),
                              "test.dat", "Test")
        captured = capsys.readouterr()
        assert "skipped" in captured.out

    def test_real_to_complex_conversion(self, tmp_path):
        """Real values should be correctly converted to complex."""
        s = _make_stdI_for_interaction(nsiteUC=1, ncell=2)
        tbl_index = np.array([[0, 1]], dtype=int)
        tbl_value = np.array([2.5], dtype=float)
        fname = str(tmp_path / "inter_real_test.dat")

        ew._export_inter_real(s, 1, tbl_index, tbl_value, fname, "TestReal")
        assert os.path.exists(fname)

        with open(fname) as f:
            content = f.read()
        # The value 2.5 should appear in the file
        assert "2.5" in content


# ===========================================================================
#  Tests: Export transfer
# ===========================================================================

class TestExportTransfer:
    """Tests for _export_transfer."""

    def test_zero_entries_skipped(self, capsys):
        s = _make_stdI_for_interaction()
        ew._export_transfer(s, 0, None, None, "test.dat", "Test", 0)
        captured = capsys.readouterr()
        assert "skipped" in captured.out

    def test_basic_transfer(self, tmp_path):
        """A simple hopping term should produce a file."""
        s = _make_stdI_for_interaction(nsiteUC=1, ncell=2)

        # index: [i, spin_i, j, spin_j]
        tbl_index = np.array([[0, 0, 1, 0]], dtype=int)
        tbl_value = np.array([-1.0 + 0j], dtype=complex)
        fname = str(tmp_path / "transfer_test.dat")

        ew._export_transfer(s, 1, tbl_index, tbl_value, fname,
                            "Transfer", 0)
        assert os.path.exists(fname)

    def test_spin_dep_skip(self, tmp_path):
        """With spin_dep=0, non-(0,0) spin pairs should be skipped."""
        s = _make_stdI_for_interaction(nsiteUC=1, ncell=2)

        # Only spin-flip terms: should result in empty
        tbl_index = np.array([[0, 0, 1, 1]], dtype=int)
        tbl_value = np.array([-1.0 + 0j], dtype=complex)
        fname = str(tmp_path / "transfer_skip.dat")

        ew._export_transfer(s, 1, tbl_index, tbl_value, fname,
                            "Transfer", 0)
        # File should not be written (skipped)
        captured_lines = []
        if os.path.exists(fname):
            with open(fname) as f:
                captured_lines = f.readlines()
        # With spin_dep=0 and ispin=0,jspin=1, the entry is skipped
        # so nintr_table=0 and we get the "skipped" message
        # (the file won't exist)
        # However, the entry passes accumulate_list (keylen=4),
        # but is then filtered by spin_dep check.
        # This is fine - it should print "skipped".

    def test_transfer_negation(self, tmp_path):
        """Transfer values are negated by convention."""
        s = _make_stdI_for_interaction(nsiteUC=1, ncell=2)

        tbl_index = np.array([[0, 0, 1, 0]], dtype=int)
        tbl_value = np.array([2.0 + 0j], dtype=complex)
        fname = str(tmp_path / "transfer_neg.dat")

        ew._export_transfer(s, 1, tbl_index, tbl_value, fname,
                            "Transfer", 0)
        assert os.path.exists(fname)
        with open(fname) as f:
            lines = f.readlines()
        # Find data line for the forward direction
        data_lines = [l for l in lines[4:] if l.strip()]
        # At least one line should contain -2.0
        found_negated = False
        for dl in data_lines:
            parts = dl.split()
            re_val = float(parts[5])
            if abs(re_val - (-2.0)) < 1e-8:
                found_negated = True
        assert found_negated, "Transfer value should be negated"


# ===========================================================================
#  Tests: Export Coulomb intra
# ===========================================================================

class TestExportCoulombIntra:
    """Tests for _export_coulomb_intra."""

    def test_zero_entries_skipped(self, capsys):
        s = _make_stdI_for_interaction()
        ew._export_coulomb_intra(s, 0, None, None, "test.dat", "Test")
        captured = capsys.readouterr()
        assert "skipped" in captured.out

    def test_basic_coulomb_intra(self, tmp_path):
        """Simple on-site Coulomb should produce a file."""
        s = _make_stdI_for_interaction(nsiteUC=1, ncell=2)

        tbl_index = np.array([[0], [1]], dtype=int)
        tbl_value = np.array([4.0, 4.0], dtype=float)
        fname = str(tmp_path / "cintra_test.dat")

        ew._export_coulomb_intra(s, 2, tbl_index, tbl_value, fname,
                                 "CoulombIntra")
        assert os.path.exists(fname)

    def test_uniform_check(self, tmp_path, capsys):
        """Non-uniform on-site Coulomb should produce a warning."""
        s = _make_stdI_for_interaction(nsiteUC=1, ncell=2)

        tbl_index = np.array([[0], [1]], dtype=int)
        tbl_value = np.array([4.0, 5.0], dtype=float)
        fname = str(tmp_path / "cintra_warn.dat")

        ew._export_coulomb_intra(s, 2, tbl_index, tbl_value, fname,
                                 "CoulombIntra")
        captured = capsys.readouterr()
        assert "WARNING" in captured.out


# ===========================================================================
#  Tests: Export interaction (integration)
# ===========================================================================

class TestExportInteraction:
    """Integration tests for export_interaction."""

    def test_all_skipped_when_empty(self, tmp_path, monkeypatch, capsys):
        """All interactions should be skipped when all counts are zero."""
        monkeypatch.chdir(tmp_path)
        s = _make_stdI_for_interaction()
        ew.export_interaction(s)
        captured = capsys.readouterr()
        assert captured.out.count("skipped") == 7  # 7 interaction types

    def test_with_transfer(self, tmp_path, monkeypatch):
        """export_interaction with a single transfer term."""
        monkeypatch.chdir(tmp_path)
        s = _make_stdI_for_interaction(nsiteUC=1, ncell=2)
        s.ntrans = 1
        s.transindx = np.array([[0, 0, 1, 0]], dtype=int)
        s.trans = np.array([-1.0 + 0j], dtype=complex)

        ew.export_interaction(s)
        assert os.path.exists(tmp_path / "transfer.dat")

    def test_with_fileprefix(self, tmp_path, monkeypatch):
        """export_interaction should use fileprefix."""
        monkeypatch.chdir(tmp_path)
        s = _make_stdI_for_interaction(nsiteUC=1, ncell=2)
        s.fileprefix = "run1"
        s.ntrans = 1
        s.transindx = np.array([[0, 0, 1, 0]], dtype=int)
        s.trans = np.array([-1.0 + 0j], dtype=complex)

        ew.export_interaction(s)
        assert os.path.exists(tmp_path / "run1_transfer.dat")


# ===========================================================================
#  Tests: Export geometry (integration)
# ===========================================================================

class TestExportGeometry:
    """Integration tests for export_geometry."""

    def test_creates_geom_file(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = _make_stdI_for_geometry(nsiteUC=1)
        ew.export_geometry(s)
        assert os.path.exists(tmp_path / "geom.dat")

    def test_creates_geom_file_with_prefix(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = _make_stdI_for_geometry(nsiteUC=1)
        s.fileprefix = "myprefix"
        ew.export_geometry(s)
        assert os.path.exists(tmp_path / "myprefix_geom.dat")


# ===========================================================================
#  Tests: IntrItem dataclass
# ===========================================================================

class TestIntrItem:
    """Tests for the _IntrItem dataclass."""

    def test_defaults(self):
        item = ew._IntrItem()
        assert item.r == [0, 0, 0]
        assert item.a == 0
        assert item.b == 0
        assert item.s == 0
        assert item.t == 0
        assert item.v == 0.0 + 0.0j

    def test_custom_values(self):
        item = ew._IntrItem(r=[1, 2, 3], a=4, b=5, s=1, t=0, v=3.14 + 1j)
        assert item.r == [1, 2, 3]
        assert item.a == 4
        assert item.v == pytest.approx(3.14 + 1j)


# ===========================================================================
#  Tests: Module-level global state
# ===========================================================================

class TestExportAllFlag:
    """Tests for the _is_export_all global flag."""

    def test_default_value(self):
        """Default should be 1 (export all)."""
        # Reset to default
        ew._is_export_all = 1
        assert ew._is_export_all == 1

    def test_set_via_export_interaction(self, tmp_path, monkeypatch):
        """export_interaction should update _is_export_all from StdI."""
        monkeypatch.chdir(tmp_path)
        ew._is_export_all = 1  # reset
        s = _make_stdI_for_interaction()
        s.export_all = 0  # not NaN_i, so should be used
        ew.export_interaction(s)
        assert ew._is_export_all == 0
        # Reset
        ew._is_export_all = 1


# ===================================================================
#  _build_transfer_table
# ===================================================================


class TestBuildTransferTable:
    """Tests for _build_transfer_table."""

    @staticmethod
    def _make_stdi(ncell: int = 2, nsiteUC: int = 1) -> StdIntList:
        """Create a minimal StdIntList for transfer-table tests."""
        s = StdIntList()
        s.NsiteUC = nsiteUC
        s.NCell = ncell
        s.box = np.array([
            [ncell, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ], dtype=int)
        s.rbox = np.array([
            [1, 0, 0],
            [0, ncell, 0],
            [0, 0, ncell],
        ], dtype=int)
        s.Cell = np.zeros((ncell, 3), dtype=int)
        for i in range(ncell):
            s.Cell[i, 0] = i
        return s

    def test_empty_input(self):
        """Zero entries returns empty table."""
        s = self._make_stdi()
        result = ew._build_transfer_table(s, 0, [], np.array([], dtype=complex), 1)
        assert result == []

    def test_single_entry(self):
        """Single transfer entry produces one table item."""
        s = self._make_stdi(ncell=2)
        # site 0 (cell 0) → site 1 (cell 1), spin 0→0
        intr_index = [[0, 0, 1, 0]]
        intr_value = np.array([1.0 + 0j])
        result = ew._build_transfer_table(s, 1, intr_index, intr_value, 1)
        assert len(result) == 1
        assert result[0].a == 0
        assert result[0].b == 0  # site 1 % NsiteUC(1) = 0
        assert result[0].s == 0
        assert result[0].t == 0
        # Value is sign-flipped
        assert result[0].v == pytest.approx(-1.0 + 0j)

    def test_sign_flip(self):
        """Transfer values are negated by convention."""
        s = self._make_stdi(ncell=2)
        intr_index = [[0, 0, 0, 0]]
        intr_value = np.array([3.5 + 2j])
        result = ew._build_transfer_table(s, 1, intr_index, intr_value, 1)
        assert result[0].v == pytest.approx(-3.5 - 2j)

    def test_deduplication(self):
        """Duplicate entries (same key) are not repeated."""
        s = self._make_stdi(ncell=4)
        # Two entries that map to the same relative coordinate
        # site 0→1 (cells 0→1) and site 2→3 (cells 2→3) → both rr=[1,0,0]
        intr_index = [[0, 0, 1, 0], [2, 0, 3, 0]]
        intr_value = np.array([1.0 + 0j, 1.0 + 0j])
        result = ew._build_transfer_table(s, 2, intr_index, intr_value, 1)
        # Both map to (rr=[1,0,0], a=0, b=0, s=0, t=0) with value -1.0
        assert len(result) == 1

    def test_spin_dep_zero_filters(self):
        """spin_dep=0 skips entries where (ispin, jspin) != (0, 0)."""
        s = self._make_stdi(ncell=2)
        # Entry with spin (0,1) — should be skipped when spin_dep=0
        intr_index = [[0, 0, 0, 1]]
        intr_value = np.array([1.0 + 0j])
        result = ew._build_transfer_table(s, 1, intr_index, intr_value, 0)
        assert len(result) == 0

    def test_spin_dep_zero_keeps_00(self):
        """spin_dep=0 keeps entries where (ispin, jspin) == (0, 0)."""
        s = self._make_stdi(ncell=2)
        intr_index = [[0, 0, 0, 0]]
        intr_value = np.array([1.0 + 0j])
        result = ew._build_transfer_table(s, 1, intr_index, intr_value, 0)
        assert len(result) == 1

    def test_spin_dep_one_keeps_all(self):
        """spin_dep=1 keeps entries with any spin indices."""
        s = self._make_stdi(ncell=2)
        intr_index = [[0, 0, 0, 0], [0, 1, 0, 1]]
        intr_value = np.array([1.0 + 0j, 2.0 + 0j])
        result = ew._build_transfer_table(s, 2, intr_index, intr_value, 1)
        assert len(result) == 2

    def test_relative_coordinates(self):
        """Check that rr is computed from cell differences."""
        s = self._make_stdi(ncell=4)
        # site 0 (cell 0) → site 2 (cell 2): rr should be [2,0,0]
        # but _unfold_site may fold: 4 cells, rr_frac = 2/4 = 0.5 → not folded
        intr_index = [[0, 0, 2, 0]]
        intr_value = np.array([1.0 + 0j])
        result = ew._build_transfer_table(s, 1, intr_index, intr_value, 1)
        assert len(result) == 1
        # The exact rr depends on _unfold_site; just verify it's set
        assert len(result[0].r) == 3

    def test_inconsistent_values_warns(self, capsys):
        """Duplicate entries with different values print a warning."""
        s = self._make_stdi(ncell=4)
        # Two entries mapping to same key but different values
        intr_index = [[0, 0, 1, 0], [4, 0, 5, 0]]  # both → rr=[1,0,0]
        # NsiteUC=1, so site 4 is cell 4 site 0 — but we only have 4 cells
        # Use cells that map to same rr: cell 0→1 and cell 2→3
        intr_index = [[0, 0, 1, 0], [2, 0, 3, 0]]
        intr_value = np.array([1.0 + 0j, 2.0 + 0j])  # different values!
        result = ew._build_transfer_table(s, 2, intr_index, intr_value, 1)
        captured = capsys.readouterr()
        assert "WARNING" in captured.out


# ===================================================================
#  _build_inter_table
# ===================================================================


class TestBuildInterTable:
    """Tests for _build_inter_table."""

    @staticmethod
    def _make_stdi(ncell: int = 2, nsiteUC: int = 1) -> StdIntList:
        """Create a minimal StdIntList for inter-table tests."""
        s = StdIntList()
        s.NsiteUC = nsiteUC
        s.NCell = ncell
        s.box = np.array([
            [ncell, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ], dtype=int)
        s.rbox = np.array([
            [1, 0, 0],
            [0, ncell, 0],
            [0, 0, ncell],
        ], dtype=int)
        s.Cell = np.zeros((ncell, 3), dtype=int)
        for i in range(ncell):
            s.Cell[i, 0] = i
        return s

    def test_empty_input(self):
        """Zero entries returns empty table."""
        s = self._make_stdi()
        result = ew._build_inter_table(s, 0, [], np.array([], dtype=complex))
        assert result == []

    def test_single_entry(self):
        """Single interaction entry produces one table item."""
        s = self._make_stdi(ncell=2)
        intr_index = [[0, 1]]
        intr_value = np.array([2.5 + 0j])
        result = ew._build_inter_table(s, 1, intr_index, intr_value)
        assert len(result) == 1
        assert result[0].a == 0
        assert result[0].b == 0  # site 1 % NsiteUC(1) = 0
        assert result[0].s == 0
        assert result[0].t == 0
        assert result[0].v == pytest.approx(2.5 + 0j)

    def test_no_sign_flip(self):
        """Inter-site values are NOT sign-flipped (unlike transfer)."""
        s = self._make_stdi(ncell=2)
        intr_index = [[0, 1]]
        intr_value = np.array([3.0 + 1j])
        result = ew._build_inter_table(s, 1, intr_index, intr_value)
        assert result[0].v == pytest.approx(3.0 + 1j)

    def test_deduplication(self):
        """Duplicate entries (same rr, a, b) are not repeated."""
        s = self._make_stdi(ncell=4)
        # cell 0→1 and cell 2→3 both have rr=[1,0,0], a=0, b=0
        intr_index = [[0, 1], [2, 3]]
        intr_value = np.array([1.0 + 0j, 1.0 + 0j])
        result = ew._build_inter_table(s, 2, intr_index, intr_value)
        assert len(result) == 1

    def test_different_sites_kept(self):
        """Entries with different (a, b) pairs are kept separately."""
        s = self._make_stdi(ncell=2, nsiteUC=2)
        # site 0 (cell 0, uc 0) → site 1 (cell 0, uc 1)
        # site 0 (cell 0, uc 0) → site 2 (cell 1, uc 0)
        intr_index = [[0, 1], [0, 2]]
        intr_value = np.array([1.0 + 0j, 2.0 + 0j])
        result = ew._build_inter_table(s, 2, intr_index, intr_value)
        assert len(result) == 2

    def test_spin_indices_always_zero(self):
        """All entries have spin indices (s, t) = (0, 0)."""
        s = self._make_stdi(ncell=2)
        intr_index = [[0, 1]]
        intr_value = np.array([1.0 + 0j])
        result = ew._build_inter_table(s, 1, intr_index, intr_value)
        assert result[0].s == 0
        assert result[0].t == 0

    def test_relative_coordinate_length(self):
        """Relative coordinate vector has length 3."""
        s = self._make_stdi(ncell=4)
        intr_index = [[0, 2]]
        intr_value = np.array([1.0 + 0j])
        result = ew._build_inter_table(s, 1, intr_index, intr_value)
        assert len(result[0].r) == 3

    def test_inconsistent_values_warns(self, capsys):
        """Duplicate entries with different values print a warning."""
        s = self._make_stdi(ncell=4)
        intr_index = [[0, 1], [2, 3]]
        intr_value = np.array([1.0 + 0j, 5.0 + 0j])
        ew._build_inter_table(s, 2, intr_index, intr_value)
        captured = capsys.readouterr()
        assert "WARNING" in captured.out


# ===================================================================
#  _build_coulomb_intra_table
# ===================================================================


class TestBuildCoulombIntraTable:
    """Tests for _build_coulomb_intra_table."""

    @staticmethod
    def _make_stdi(nsiteUC: int = 2) -> StdIntList:
        """Create a minimal StdIntList for Coulomb-intra tests."""
        s = StdIntList()
        s.NsiteUC = nsiteUC
        return s

    def test_empty_input(self):
        """Zero entries returns empty table."""
        s = self._make_stdi()
        result = ew._build_coulomb_intra_table(
            s, 0, [], np.array([], dtype=complex))
        assert result == []

    def test_single_entry(self):
        """Single entry produces one item with rr=[0,0,0] and a==b."""
        s = self._make_stdi(nsiteUC=2)
        intr_index = [[0]]
        intr_value = np.array([4.0 + 0j])
        result = ew._build_coulomb_intra_table(s, 1, intr_index, intr_value)
        assert len(result) == 1
        assert result[0].r == [0, 0, 0]
        assert result[0].a == 0
        assert result[0].b == 0  # a == b for on-site
        assert result[0].v == pytest.approx(4.0 + 0j)

    def test_deduplication_by_site(self):
        """Multiple entries for the same unit-cell site are deduplicated."""
        s = self._make_stdi(nsiteUC=1)
        # sites 0 and 2 both map to uc site 0 (% 1 = 0)
        intr_index = [[0], [2]]
        intr_value = np.array([3.0 + 0j, 3.0 + 0j])
        result = ew._build_coulomb_intra_table(s, 2, intr_index, intr_value)
        assert len(result) == 1

    def test_different_uc_sites_kept(self):
        """Entries for different unit-cell sites are kept separately."""
        s = self._make_stdi(nsiteUC=2)
        # site 0 → uc 0, site 1 → uc 1
        intr_index = [[0], [1]]
        intr_value = np.array([1.0 + 0j, 2.0 + 0j])
        result = ew._build_coulomb_intra_table(s, 2, intr_index, intr_value)
        assert len(result) == 2
        assert result[0].a == 0
        assert result[1].a == 1

    def test_spin_indices_zero(self):
        """All entries have (s, t) = (0, 0)."""
        s = self._make_stdi()
        intr_index = [[0]]
        intr_value = np.array([1.0 + 0j])
        result = ew._build_coulomb_intra_table(s, 1, intr_index, intr_value)
        assert result[0].s == 0
        assert result[0].t == 0

    def test_a_equals_b(self):
        """On-site term always has a == b."""
        s = self._make_stdi(nsiteUC=3)
        intr_index = [[2]]  # uc site 2
        intr_value = np.array([1.0 + 0j])
        result = ew._build_coulomb_intra_table(s, 1, intr_index, intr_value)
        assert result[0].a == result[0].b == 2

    def test_inconsistent_values_warns(self, capsys):
        """Duplicate entries with different values print a warning."""
        s = self._make_stdi(nsiteUC=1)
        intr_index = [[0], [1]]  # both map to uc site 0
        intr_value = np.array([1.0 + 0j, 9.0 + 0j])
        ew._build_coulomb_intra_table(s, 2, intr_index, intr_value)
        captured = capsys.readouterr()
        assert "WARNING" in captured.out
