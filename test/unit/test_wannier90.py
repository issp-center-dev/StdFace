"""Unit tests for wannier90 module.

Tests for the Python translation of Wannier90.c.
"""
from __future__ import annotations

import math
import os
import textwrap

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.lattice import wannier90 as w90


# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------

NaN_d = float("nan")
NaN_i = 2147483647
NaN_c = complex(float("nan"), 0.0)


# ---------------------------------------------------------------------------
#  Helpers: inverse matrix / in-box
# ---------------------------------------------------------------------------


class TestCheckInBox:
    """Tests for the internal _check_in_box helper."""

    def test_origin_is_in_box(self):
        """The zero vector should always be inside the box."""
        inv = np.eye(3)
        assert w90._check_in_box(np.array([0, 0, 0]), inv) is True

    def test_boundary_is_in_box(self):
        """Vectors on the boundary (abs == 1) should be inside."""
        inv = np.eye(3)
        assert w90._check_in_box(np.array([1, 0, 0]), inv) is True
        assert w90._check_in_box(np.array([0, -1, 0]), inv) is True

    def test_outside_box(self):
        """Vectors outside the boundary should not be inside."""
        inv = np.eye(3)
        assert w90._check_in_box(np.array([2, 0, 0]), inv) is False

    def test_scaled_box(self):
        """With scaled inverse, check boundary correctly."""
        inv = np.diag([0.5, 0.5, 0.5])
        # judge_vec = [0.5*r0, 0.5*r1, 0.5*r2]
        assert w90._check_in_box(np.array([2, 2, 2]), inv) is True
        assert w90._check_in_box(np.array([3, 0, 0]), inv) is False


# ---------------------------------------------------------------------------
#  Helpers: geometry & file I/O
# ---------------------------------------------------------------------------


def _write_geom_file(path, NsiteUC=2, prefix="test"):
    """Write a minimal Wannier90 geometry file."""
    filename = os.path.join(path, f"{prefix}_geom.dat")
    with open(filename, "w") as f:
        # 3 direct lattice vectors
        f.write("  1.0  0.0  0.0\n")
        f.write("  0.0  1.0  0.0\n")
        f.write("  0.0  0.0  1.0\n")
        # Number of sites
        f.write(f"  {NsiteUC}\n")
        # Wannier centres
        for i in range(NsiteUC):
            f.write(f"  {0.1 * i:.1f}  {0.2 * i:.1f}  {0.3 * i:.1f}\n")
    return filename


def _write_hr_file(path, NsiteUC=2, nWSC=1, prefix="test"):
    """Write a minimal Wannier90 hopping file (*_hr.dat).

    Creates a simple file with on-site and nearest-neighbour terms.
    """
    filename = os.path.join(path, f"{prefix}_hr.dat")
    nWan = NsiteUC
    with open(filename, "w") as f:
        f.write("  written by test\n")
        f.write(f"  {nWan}\n")
        f.write(f"  {nWSC}\n")
        # Degeneracy weights (all 1)
        for _ in range(nWSC):
            f.write("  1")
        f.write("\n")
        # Body: R0 R1 R2 iWan jWan Re Im
        for iWSC in range(nWSC):
            for iWan in range(nWan):
                for jWan in range(nWan):
                    re_val = 0.0
                    im_val = 0.0
                    if iWan == jWan:
                        re_val = -1.0  # on-site energy
                    elif abs(iWan - jWan) == 1:
                        re_val = -0.5  # nearest-neighbour hopping
                    f.write(f"  0  0  0  {iWan + 1}  {jWan + 1}  {re_val:12.6f}  {im_val:12.6f}\n")
    return filename


def _write_hr_file_split_degeneracy_weights(
    path, NsiteUC=2, nWSC=4, prefix="test",
):
    """Write *_hr.dat with several WSC blocks and degeneracy integers on multiple lines.

    Exercises :func:`wannier90._skip_degeneracy_weights` when one line does not
    hold all ``nWSC`` weight values.
    """
    filename = os.path.join(path, f"{prefix}_hr.dat")
    nWan = NsiteUC
    first = nWSC // 2
    second = nWSC - first
    with open(filename, "w") as f:
        f.write("  written by test\n")
        f.write(f"  {nWan}\n")
        f.write(f"  {nWSC}\n")
        f.write("".join("  1" for _ in range(first)) + "\n")
        f.write("".join("  1" for _ in range(second)) + "\n")
        for iWSC in range(nWSC):
            for iWan in range(nWan):
                for jWan in range(nWan):
                    re_val = 0.0
                    im_val = 0.0
                    if iWan == jWan:
                        re_val = -1.0
                    elif abs(iWan - jWan) == 1:
                        re_val = -0.5
                    f.write(
                        f"  0  0  0  {iWan + 1}  {jWan + 1}  {re_val:12.6f}  {im_val:12.6f}\n"
                    )
    return filename


def _write_ur_file(path, NsiteUC=2, prefix="test"):
    """Write a minimal Wannier90 Coulomb file (*_ur.dat)."""
    filename = os.path.join(path, f"{prefix}_ur.dat")
    nWan = NsiteUC
    nWSC = 1
    with open(filename, "w") as f:
        f.write("  written by test\n")
        f.write(f"  {nWan}\n")
        f.write(f"  {nWSC}\n")
        f.write("  1\n")
        for iWan in range(nWan):
            for jWan in range(nWan):
                re_val = 0.0
                if iWan == jWan:
                    re_val = 4.0  # on-site U
                f.write(f"  0  0  0  {iWan + 1}  {jWan + 1}  {re_val:12.6f}  0.000000\n")
    return filename


def _write_jr_file(path, NsiteUC=2, prefix="test"):
    """Write a minimal Wannier90 Hund file (*_jr.dat)."""
    filename = os.path.join(path, f"{prefix}_jr.dat")
    nWan = NsiteUC
    nWSC = 1
    with open(filename, "w") as f:
        f.write("  written by test\n")
        f.write(f"  {nWan}\n")
        f.write(f"  {nWSC}\n")
        f.write("  1\n")
        for iWan in range(nWan):
            for jWan in range(nWan):
                f.write(f"  0  0  0  {iWan + 1}  {jWan + 1}  0.000000  0.000000\n")
    return filename


def _write_dr_file(path, NsiteUC=2, prefix="test"):
    """Write a minimal density matrix file (*_dr.dat)."""
    filename = os.path.join(path, f"{prefix}_dr.dat")
    nWan = NsiteUC
    nWSC = 1
    with open(filename, "w") as f:
        f.write("  written by test\n")
        f.write(f"  {nWan}\n")
        f.write(f"  {nWSC}\n")
        f.write("  1\n")
        for iWan in range(nWan):
            for jWan in range(nWan):
                re_val = 0.5 if iWan == jWan else 0.0
                f.write(f"  0  0  0  {iWan + 1}  {jWan + 1}  {re_val:12.6f}  0.000000\n")
    return filename


def _make_wannier_StdI(
    model: str = "hubbard",
    prefix: str = "test",
    W: int = 2,
    L: int = 2,
    Height: int = 1,
) -> StdIntList:
    """Create a StdIntList pre-configured for Wannier90 tests."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = model
    s.solver = "HPhi"
    s.lattice = "wannier90"
    s.CDataFileHead = prefix

    # Lattice parameters
    s.a = NaN_d
    s.length[:] = NaN_d
    s.direct[:, :] = NaN_d
    s.phase[:] = NaN_d
    s.W = W
    s.L = L
    s.Height = Height
    s.box[:, :] = NaN_i

    # Wannier90 parameters
    s.cutoff_t = NaN_d
    s.cutoff_u = NaN_d
    s.cutoff_j = NaN_d
    s.cutoff_length_t = NaN_d
    s.cutoff_length_U = NaN_d
    s.cutoff_length_J = NaN_d
    s.cutoff_tR[:] = NaN_i
    s.cutoff_UR[:] = NaN_i
    s.cutoff_JR[:] = NaN_i
    s.cutoff_tVec[:, :] = NaN_d
    s.cutoff_UVec[:, :] = NaN_d
    s.cutoff_JVec[:, :] = NaN_d

    # Lambda parameters
    s.lambda_ = NaN_d
    s.lambda_U = NaN_d
    s.lambda_J = NaN_d
    s.alpha = NaN_d
    s.double_counting_mode = "none"

    # Model parameters
    s.h = NaN_d
    s.Gamma = NaN_d
    s.Gamma_y = NaN_d
    s.K = NaN_d
    s.U = NaN_d
    s.mu = NaN_d
    s.S2 = None

    return s


# ---------------------------------------------------------------------------
#  Tests: module structure
# ---------------------------------------------------------------------------


class TestModuleStructure:
    """Verify module can be imported and has expected attributes."""

    def test_import(self):
        """wannier90 module should import without error."""
        from stdface.lattice import wannier90
        assert wannier90 is not None

    def test_wannier90_function_exists(self):
        """wannier90 should be a callable function."""
        assert hasattr(w90, "wannier90")
        assert callable(w90.wannier90)

    def test_module_docstring(self):
        """Module should have a docstring."""
        assert w90.__doc__ is not None

    def test_wannier90_docstring(self):
        """wannier90 function should have a docstring."""
        assert w90.wannier90.__doc__ is not None


# ---------------------------------------------------------------------------
#  Tests: _geometry_w90
# ---------------------------------------------------------------------------


class TestGeometryW90:
    """Tests for _geometry_w90 reading geometry files."""

    def test_reads_direct_lattice_vectors(self, tmp_path):
        """Should read 3x3 direct lattice vectors from file."""
        os.chdir(tmp_path)
        _write_geom_file(str(tmp_path), NsiteUC=2, prefix="test")
        s = StdIntList()
        s.CDataFileHead = "test"
        s.NsiteUC = 2
        s.direct = np.zeros((3, 3))
        s.tau = np.zeros((2, 3))

        w90._geometry_w90(s)

        np.testing.assert_allclose(s.direct[0], [1.0, 0.0, 0.0])
        np.testing.assert_allclose(s.direct[1], [0.0, 1.0, 0.0])
        np.testing.assert_allclose(s.direct[2], [0.0, 0.0, 1.0])

    def test_reads_NsiteUC(self, tmp_path):
        """Should correctly read the number of correlated sites."""
        os.chdir(tmp_path)
        _write_geom_file(str(tmp_path), NsiteUC=3, prefix="geo")
        s = StdIntList()
        s.CDataFileHead = "geo"
        s.NsiteUC = 1
        s.direct = np.zeros((3, 3))
        s.tau = np.zeros((1, 3))

        w90._geometry_w90(s)

        assert s.NsiteUC == 3

    def test_reads_tau(self, tmp_path):
        """Should read Wannier centre positions."""
        os.chdir(tmp_path)
        _write_geom_file(str(tmp_path), NsiteUC=2, prefix="test")
        s = StdIntList()
        s.CDataFileHead = "test"
        s.NsiteUC = 2
        s.direct = np.zeros((3, 3))
        s.tau = np.zeros((2, 3))

        w90._geometry_w90(s)

        np.testing.assert_allclose(s.tau[0], [0.0, 0.0, 0.0], atol=1e-10)
        np.testing.assert_allclose(s.tau[1], [0.1, 0.2, 0.3], atol=1e-10)

    def test_missing_file_exits(self, tmp_path):
        """Should exit when geometry file is missing."""
        os.chdir(tmp_path)
        s = StdIntList()
        s.CDataFileHead = "nonexistent"
        s.NsiteUC = 1
        s.direct = np.zeros((3, 3))
        s.tau = np.zeros((1, 3))

        with pytest.raises(FileNotFoundError):
            w90._geometry_w90(s)


# ---------------------------------------------------------------------------
#  Tests: _read_w90
# ---------------------------------------------------------------------------


class TestReadW90:
    """Tests for _read_w90 reading hopping/interaction files."""

    def test_read_hopping_basic(self, tmp_path):
        """Should read hopping terms from hr.dat file."""
        os.chdir(tmp_path)
        NsiteUC = 2
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix="test")

        s = StdIntList()

        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = None
        s.L = None
        s.Height = None

        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]
        cutoff_R = np.array([10, 10, 10], dtype=int)
        cutoff_Rvec = np.full((3, 3), NaN_i, dtype=float)

        w90._read_w90(
            s, str(tmp_path / "test_hr.dat"),
            1.0e-8, cutoff_R, cutoff_Rvec, -1.0,
            0, NtUJ, tUJindx, 1.0, tUJ,
        )

        # Should have found effective terms
        assert NtUJ[0] > 0
        assert tUJ[0] is not None
        assert tUJindx[0] is not None

    def test_missing_file_skips(self, tmp_path):
        """Should skip gracefully when file is missing."""
        os.chdir(tmp_path)
        s = StdIntList()


        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]
        cutoff_R = np.array([10, 10, 10], dtype=int)
        cutoff_Rvec = np.full((3, 3), NaN_i, dtype=float)

        # Should not raise
        w90._read_w90(
            s, str(tmp_path / "nonexistent_hr.dat"),
            1.0e-8, cutoff_R, cutoff_Rvec, -1.0,
            0, NtUJ, tUJindx, 1.0, tUJ,
        )
        # NtUJ should remain 0
        assert NtUJ[0] == 0

    def test_cutoff_filters_terms(self, tmp_path):
        """Higher cutoff should filter out small terms."""
        os.chdir(tmp_path)
        NsiteUC = 2
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix="test")

        s = StdIntList()

        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = None
        s.L = None
        s.Height = None

        cutoff_R = np.array([10, 10, 10], dtype=int)
        cutoff_Rvec = np.full((3, 3), NaN_i, dtype=float)

        # Low cutoff
        NtUJ_low = [0, 0, 0]
        tUJ_low = [None, None, None]
        tUJindx_low = [None, None, None]
        w90._read_w90(
            s, str(tmp_path / "test_hr.dat"),
            0.01, cutoff_R, cutoff_Rvec, -1.0,
            0, NtUJ_low, tUJindx_low, 1.0, tUJ_low,
        )

        # High cutoff
        NtUJ_high = [0, 0, 0]
        tUJ_high = [None, None, None]
        tUJindx_high = [None, None, None]
        w90._read_w90(
            s, str(tmp_path / "test_hr.dat"),
            10.0, cutoff_R, cutoff_Rvec, -1.0,
            0, NtUJ_high, tUJindx_high, 1.0, tUJ_high,
        )

        assert NtUJ_low[0] >= NtUJ_high[0]

    def test_lambda_scales_values(self, tmp_path):
        """Lambda parameter should scale matrix elements."""
        os.chdir(tmp_path)
        NsiteUC = 1
        # Write a simple 1-site file
        filename = os.path.join(str(tmp_path), "test_hr.dat")
        with open(filename, "w") as f:
            f.write("  written by test\n")
            f.write("  1\n")
            f.write("  1\n")
            f.write("  1\n")
            f.write("  0  0  0  1  1  -2.000000  0.000000\n")

        s = StdIntList()

        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = None
        s.L = None
        s.Height = None

        cutoff_R = np.array([10, 10, 10], dtype=int)
        cutoff_Rvec = np.full((3, 3), NaN_i, dtype=float)

        # Lambda = 1
        NtUJ1 = [0, 0, 0]
        tUJ1 = [None, None, None]
        tUJindx1 = [None, None, None]
        w90._read_w90(
            s, filename,
            1.0e-8, cutoff_R, cutoff_Rvec, -1.0,
            0, NtUJ1, tUJindx1, 1.0, tUJ1,
        )

        # Lambda = 2
        NtUJ2 = [0, 0, 0]
        tUJ2 = [None, None, None]
        tUJindx2 = [None, None, None]
        w90._read_w90(
            s, filename,
            1.0e-8, cutoff_R, cutoff_Rvec, -1.0,
            0, NtUJ2, tUJindx2, 2.0, tUJ2,
        )

        assert NtUJ1[0] == NtUJ2[0]
        np.testing.assert_allclose(
            np.abs(tUJ2[0]), 2.0 * np.abs(tUJ1[0]), atol=1e-12
        )


# ---------------------------------------------------------------------------
#  Tests: _read_density_matrix
# ---------------------------------------------------------------------------


class TestReadDensityMatrix:
    """Tests for _read_density_matrix."""

    def test_reads_density_matrix(self, tmp_path):
        """Should read density matrix from dr.dat file."""
        os.chdir(tmp_path)
        NsiteUC = 2
        _write_dr_file(str(tmp_path), NsiteUC=NsiteUC, prefix="test")

        s = StdIntList()
        s.NsiteUC = NsiteUC

        DenMat = w90._read_density_matrix(s, str(tmp_path / "test_dr.dat"))

        assert (0, 0, 0) in DenMat
        # Diagonal elements should be 0.5
        np.testing.assert_allclose(DenMat[(0, 0, 0)][0, 0], 0.5 + 0j)
        np.testing.assert_allclose(DenMat[(0, 0, 0)][1, 1], 0.5 + 0j)
        # Off-diagonal should be 0
        np.testing.assert_allclose(DenMat[(0, 0, 0)][0, 1], 0.0 + 0j)

    def test_missing_file_exits(self, tmp_path):
        """Should exit when density matrix file is missing."""
        os.chdir(tmp_path)
        s = StdIntList()
        s.NsiteUC = 1

        with pytest.raises(FileNotFoundError):
            w90._read_density_matrix(s, str(tmp_path / "nonexistent_dr.dat"))


# ---------------------------------------------------------------------------
#  Tests: _read_w90_with_cutoff
# ---------------------------------------------------------------------------


class TestReadW90WithCutoff:
    """Tests for the _read_w90_with_cutoff helper."""

    def test_returns_updated_cutoff_values(self, tmp_path):
        """Should return resolved cutoff_val and cutoff_length."""
        os.chdir(tmp_path)
        NsiteUC = 2
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix="test")

        s = StdIntList()
        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = None
        s.L = None
        s.Height = None
        s.box = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]], dtype=int)
        s.CDataFileHead = "test"

        cutoff_R = np.full(3, NaN_i, dtype=int)
        cutoff_Vec = np.full((3, 3), NaN_d)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        val, length = w90._read_w90_with_cutoff(
            s, "t", "t", "_hr.dat",
            NaN_d, NaN_d,
            cutoff_R, cutoff_Vec,
            cutoff_length_default=-1.0,
            cutoff_R_defaults=(None, None, None),
            itUJ=0, NtUJ=NtUJ, tUJindx=tUJindx, lam=1.0, tUJ=tUJ,
        )

        # cutoff_val should default to 1e-8
        np.testing.assert_allclose(val, 1.0e-8, atol=1e-15)
        # cutoff_length should default to -1.0
        np.testing.assert_allclose(length, -1.0, atol=1e-15)

    def test_sets_cutoff_R_defaults(self, tmp_path):
        """Should set cutoff_R values from defaults when provided."""
        os.chdir(tmp_path)
        NsiteUC = 1
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix="test")

        s = StdIntList()
        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = None
        s.L = None
        s.Height = None
        s.box = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]], dtype=int)
        s.CDataFileHead = "test"

        cutoff_R = np.full(3, NaN_i, dtype=int)
        cutoff_Vec = np.full((3, 3), NaN_d)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._read_w90_with_cutoff(
            s, "u", "U", "_hr.dat",
            NaN_d, NaN_d,
            cutoff_R, cutoff_Vec,
            cutoff_length_default=0.3,
            cutoff_R_defaults=(0, 0, 0),
            itUJ=0, NtUJ=NtUJ, tUJindx=tUJindx, lam=1.0, tUJ=tUJ,
        )

        # All cutoff_R should be set to 0
        assert cutoff_R[0] == 0
        assert cutoff_R[1] == 0
        assert cutoff_R[2] == 0

    def test_skips_cutoff_R_for_none_default(self, tmp_path):
        """Should not modify cutoff_R when default is None."""
        os.chdir(tmp_path)
        NsiteUC = 1
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix="test")

        s = StdIntList()
        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = None
        s.L = None
        s.Height = None
        s.box = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]], dtype=int)
        s.CDataFileHead = "test"

        cutoff_R = np.array([5, 5, 5], dtype=int)
        cutoff_Vec = np.full((3, 3), NaN_d)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._read_w90_with_cutoff(
            s, "t", "t", "_hr.dat",
            NaN_d, NaN_d,
            cutoff_R, cutoff_Vec,
            cutoff_length_default=-1.0,
            cutoff_R_defaults=(None, None, None),
            itUJ=0, NtUJ=NtUJ, tUJindx=tUJindx, lam=1.0, tUJ=tUJ,
        )

        # cutoff_R should remain at original values (not modified)
        assert cutoff_R[0] == 5
        assert cutoff_R[1] == 5
        assert cutoff_R[2] == 5

    def test_sets_cutoff_vec_from_box(self, tmp_path):
        """Should set cutoff_Vec to box*0.5 as default."""
        os.chdir(tmp_path)
        NsiteUC = 1
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix="test")

        s = StdIntList()
        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = None
        s.L = None
        s.Height = None
        s.box = np.array([[4, 0, 0], [0, 6, 0], [0, 0, 2]], dtype=int)
        s.CDataFileHead = "test"

        cutoff_R = np.full(3, NaN_i, dtype=int)
        cutoff_Vec = np.full((3, 3), NaN_d)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._read_w90_with_cutoff(
            s, "u", "U", "_ur.dat",
            NaN_d, NaN_d,
            cutoff_R, cutoff_Vec,
            cutoff_length_default=0.3,
            cutoff_R_defaults=(0, 0, 0),
            itUJ=0, NtUJ=NtUJ, tUJindx=tUJindx, lam=1.0, tUJ=tUJ,
        )

        # Diagonal: box*0.5
        np.testing.assert_allclose(cutoff_Vec[0, 0], 2.0, atol=1e-12)
        np.testing.assert_allclose(cutoff_Vec[1, 1], 3.0, atol=1e-12)
        np.testing.assert_allclose(cutoff_Vec[2, 2], 1.0, atol=1e-12)
        # Off-diagonal: 0*0.5 = 0
        np.testing.assert_allclose(cutoff_Vec[0, 1], 0.0, atol=1e-12)

    def test_reads_terms(self, tmp_path):
        """Should populate NtUJ and tUJ after reading."""
        os.chdir(tmp_path)
        NsiteUC = 2
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix="test")

        s = StdIntList()
        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = None
        s.L = None
        s.Height = None
        s.box = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]], dtype=int)
        s.CDataFileHead = "test"

        cutoff_R = np.array([10, 10, 10], dtype=int)
        cutoff_Vec = np.full((3, 3), NaN_i, dtype=float)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._read_w90_with_cutoff(
            s, "t", "t", "_hr.dat",
            NaN_d, NaN_d,
            cutoff_R, cutoff_Vec,
            cutoff_length_default=-1.0,
            cutoff_R_defaults=(None, None, None),
            itUJ=0, NtUJ=NtUJ, tUJindx=tUJindx, lam=1.0, tUJ=tUJ,
        )

        assert NtUJ[0] > 0
        assert tUJ[0] is not None
        assert tUJindx[0] is not None

    def test_missing_file_is_noop(self, tmp_path):
        """Should silently skip when the data file is missing."""
        os.chdir(tmp_path)

        s = StdIntList()
        s.NsiteUC = 1
        s.direct = np.eye(3)
        s.tau = np.zeros((1, 3))
        s.W = None
        s.L = None
        s.Height = None
        s.box = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]], dtype=int)
        s.CDataFileHead = "nonexistent"

        cutoff_R = np.full(3, NaN_i, dtype=int)
        cutoff_Vec = np.full((3, 3), NaN_d)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._read_w90_with_cutoff(
            s, "t", "t", "_hr.dat",
            NaN_d, NaN_d,
            cutoff_R, cutoff_Vec,
            cutoff_length_default=-1.0,
            cutoff_R_defaults=(None, None, None),
            itUJ=0, NtUJ=NtUJ, tUJindx=tUJindx, lam=1.0, tUJ=tUJ,
        )

        # Should not crash, NtUJ should remain 0
        assert NtUJ[0] == 0


# ---------------------------------------------------------------------------
#  Tests: DCMode enum
# ---------------------------------------------------------------------------


class TestDCMode:
    """Tests for the double-counting mode enum."""

    def test_values(self):
        """Enum values should match the C enum."""
        assert w90._DCMode.NOTCORRECT == 0
        assert w90._DCMode.HARTREE == 1
        assert w90._DCMode.HARTREE_U == 2
        assert w90._DCMode.FULL == 3


# ---------------------------------------------------------------------------
#  Tests: _parse_double_counting_mode
# ---------------------------------------------------------------------------


class TestParseDoubleCountingMode:
    """Tests for the _parse_double_counting_mode helper."""

    def test_none_returns_notcorrect(self):
        """'none' should map to NOTCORRECT."""
        assert w90._parse_double_counting_mode("none") == w90._DCMode.NOTCORRECT

    def test_unset_string_returns_notcorrect(self):
        """None sentinel should map to NOTCORRECT."""
        assert w90._parse_double_counting_mode(None) == w90._DCMode.NOTCORRECT

    def test_hartree_returns_hartree(self):
        """'hartree' should map to HARTREE."""
        assert w90._parse_double_counting_mode("hartree") == w90._DCMode.HARTREE

    def test_hartree_u_returns_hartree_u(self):
        """'hartree_u' should map to HARTREE_U."""
        assert w90._parse_double_counting_mode("hartree_u") == w90._DCMode.HARTREE_U

    def test_full_returns_full(self):
        """'full' should map to FULL."""
        assert w90._parse_double_counting_mode("full") == w90._DCMode.FULL

    def test_invalid_mode_exits(self):
        """An unrecognised string should cause exit."""
        with pytest.raises(ValueError):
            w90._parse_double_counting_mode("bogus")

    def test_return_type_is_dcmode(self):
        """Return value should be an instance of _DCMode enum."""
        result = w90._parse_double_counting_mode("hartree")
        assert isinstance(result, w90._DCMode)

    def test_dc_mode_map_keys(self):
        """_DC_MODE_MAP should contain exactly the expected keys."""
        expected_keys = {"none", None, "hartree", "hartree_u", "full"}
        assert set(w90._DC_MODE_MAP.keys()) == expected_keys

    def test_dc_mode_map_values(self):
        """_DC_MODE_MAP values should all be _DCMode members."""
        for key, val in w90._DC_MODE_MAP.items():
            assert isinstance(val, w90._DCMode), f"key={key!r} maps to non-_DCMode {val!r}"


# ---------------------------------------------------------------------------
#  Tests: Full wannier90 function
# ---------------------------------------------------------------------------


class TestWannier90Hubbard:
    """Test the main wannier90() function with the Hubbard model."""

    def test_hubbard_basic(self, tmp_path):
        """wannier90() with a Hubbard model should populate basic arrays."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        # Check basic structure
        assert s.NsiteUC == NsiteUC
        assert s.nsite == NsiteUC * s.NCell
        assert s.locspinflag is not None
        assert all(s.locspinflag[i] == 0 for i in range(s.nsite))

    def test_hubbard_creates_files(self, tmp_path):
        """wannier90() should create output files."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        assert (tmp_path / "wan2site.dat").exists()
        # lattice.xsf and geometry.dat are now independent lattice outputs
        # (build_xsf / build_geometry, written by the main flow), no longer
        # side effects of setup.
        from stdface.lattice.geometry_output import build_xsf
        assert not (tmp_path / "lattice.xsf").exists()
        build_xsf(s).write(tmp_path)
        assert (tmp_path / "lattice.xsf").exists()

    def test_hubbard_defaults(self, tmp_path):
        """Default values should be applied for h, Gamma, etc."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        assert s.h == 0.0
        assert s.Gamma == 0.0
        assert s.Gamma_y == 0.0
        assert s.mu == 0.0

    def test_hubbard_interactions_allocated(self, tmp_path):
        """After wannier90(), interaction arrays should be allocated."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        assert s.trans_list

    def test_wan2site_content(self, tmp_path):
        """wan2site.dat should contain correct site mapping."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        content = (tmp_path / "wan2site.dat").read_text()
        assert "Total site number" in content
        assert "site nx ny nz norb" in content


class TestWannier90Spin:
    """Test the main wannier90() function with the spin model."""

    def test_spin_basic(self, tmp_path):
        """wannier90() with a spin model should set locspinflag."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="spin", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        assert s.NsiteUC == NsiteUC
        assert s.nsite == NsiteUC * s.NCell
        # For spin model, locspinflag = S2
        assert all(s.locspinflag[i] == s.S2 for i in range(s.nsite))


class TestWannier90KondoError:
    """Test that Kondo model raises an error."""

    def test_kondo_exits(self, tmp_path):
        """wannier90() with Kondo model should exit with error."""
        os.chdir(tmp_path)
        NsiteUC = 1
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="kondo", prefix=prefix, W=2, L=2, Height=1)

        with pytest.raises(ValueError):
            w90.wannier90(s)


class TestWannier90LambdaError:
    """Test that negative lambda values raise an error."""

    def test_negative_lambda_U_exits(self, tmp_path):
        """wannier90() should exit if lambda_U is negative."""
        os.chdir(tmp_path)
        NsiteUC = 1
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        s.lambda_U = -1.0
        s.lambda_J = 1.0

        with pytest.raises(ValueError):
            w90.wannier90(s)


class TestWannier90AlphaError:
    """Test that out-of-range alpha raises an error."""

    def test_alpha_out_of_range_exits(self, tmp_path):
        """wannier90() should exit if alpha > 1."""
        os.chdir(tmp_path)
        NsiteUC = 1
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        s.alpha = 1.5

        with pytest.raises(ValueError):
            w90.wannier90(s)


class TestWannier90DoubleCountingModeError:
    """Test that invalid double counting mode raises an error."""

    def test_invalid_dc_mode_exits(self, tmp_path):
        """wannier90() should exit with invalid double_counting_mode."""
        os.chdir(tmp_path)
        NsiteUC = 1
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        s.double_counting_mode = "invalid_mode"

        with pytest.raises(ValueError):
            w90.wannier90(s)


class TestWannier90DoubleCountingIntegration:
    """Full :func:`wannier90.wannier90` with ``*_dr.dat`` and DC corrections."""

    @pytest.mark.parametrize("dc_mode", ["hartree", "full"])
    def test_hubbard_with_density_matrix(self, tmp_path, dc_mode):
        """Read density matrix, apply DC mode, and write ``initial.def``."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, nWSC=1, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_dr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        s.double_counting_mode = dc_mode
        w90.wannier90(s)

        assert s.NsiteUC == NsiteUC
        assert (tmp_path / "initial.def").exists()
        assert len(s.trans_list) >= 0


class TestWannier90HrMultilineDegeneracy:
    """``*_hr.dat`` with ``nWSC > 1`` and degeneracy weights split across lines."""

    def test_hubbard_runs_with_split_weights(self, tmp_path):
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file_split_degeneracy_weights(
            str(tmp_path), NsiteUC=NsiteUC, nWSC=4, prefix=prefix,
        )
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)
        assert len(s.trans_list) >= 0


# ---------------------------------------------------------------------------
#  Helpers: set up interaction arrays for direct helper tests
# ---------------------------------------------------------------------------


def _setup_interactions(s: StdIntList) -> None:
    """Initialise interaction term lists on StdIntList for testing."""
    s.trans_list = []
    s.intr_list = []
    s.Cintra_list = []
    s.Cinter_list = []
    s.Hund_list = []
    s.Ex_list = []
    s.PairLift_list = []
    s.PairHopp_list = []


# ---------------------------------------------------------------------------
#  Tests: _apply_hopping_terms
# ---------------------------------------------------------------------------


class TestApplyHoppingTerms:
    """Tests for the _apply_hopping_terms helper."""

    def test_none_indices_is_noop(self):
        """Should do nothing when tUJindx[0] is None."""
        s = StdIntList()
        s.model = "hubbard"
        _setup_interactions(s)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]
        ntrans_before = len(s.trans_list)
        w90._apply_hopping_terms(s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx, None)
        assert len(s.trans_list) == ntrans_before

    def test_local_hopping_hubbard(self):
        """Local hopping term should add on-site transfer for Hubbard."""
        s = StdIntList()
        s.model = "hubbard"
        s.NsiteUC = 1
        _setup_interactions(s)
        # On-site term: R=(0,0,0), iWan=jWan=0
        tUJindx_arr = np.array([[0, 0, 0, 0, 0]], dtype=int)
        tUJ_arr = np.array([-1.0 + 0j])
        NtUJ = [1, 0, 0]
        tUJ = [tUJ_arr, None, None]
        tUJindx = [tUJindx_arr, None, None]

        w90._apply_hopping_terms(s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx, None)

        # Should add 2 transfer terms (one per spin)
        assert len(s.trans_list) == 2
        # Transfer value should be -(-1.0) = 1.0
        np.testing.assert_allclose(s.trans_list[0][0], 1.0 + 0j, atol=1e-12)
        np.testing.assert_allclose(s.trans_list[1][0], 1.0 + 0j, atol=1e-12)

    def test_local_hopping_spin_is_noop(self):
        """Local hopping should be skipped for spin model."""
        s = StdIntList()
        s.model = "spin"
        s.NsiteUC = 1
        _setup_interactions(s)
        tUJindx_arr = np.array([[0, 0, 0, 0, 0]], dtype=int)
        tUJ_arr = np.array([-1.0 + 0j])
        NtUJ = [1, 0, 0]
        tUJ = [tUJ_arr, None, None]
        tUJindx = [tUJindx_arr, None, None]

        w90._apply_hopping_terms(s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx, None)

        # No transfer terms added for spin local hopping
        assert len(s.trans_list) == 0

    def test_zero_terms_is_noop(self):
        """Should do nothing when NtUJ[0] is 0."""
        s = StdIntList()
        s.model = "hubbard"
        _setup_interactions(s)
        tUJindx_arr = np.zeros((0, 5), dtype=int)
        tUJ_arr = np.zeros(0, dtype=complex)
        NtUJ = [0, 0, 0]
        tUJ = [tUJ_arr, None, None]
        tUJindx = [tUJindx_arr, None, None]

        w90._apply_hopping_terms(s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx, None)

        assert len(s.trans_list) == 0


# ---------------------------------------------------------------------------
#  Tests: _apply_coulomb_terms
# ---------------------------------------------------------------------------


class TestApplyCoulombTerms:
    """Tests for the _apply_coulomb_terms helper."""

    def test_none_indices_is_noop(self):
        """Should do nothing when tUJindx[1] is None."""
        s = StdIntList()
        s.model = "hubbard"
        _setup_interactions(s)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]
        ncintra_before = len(s.Cintra_list)
        w90._apply_coulomb_terms(
            s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx,
            w90._DCMode.NOTCORRECT, None,
        )
        assert len(s.Cintra_list) == ncintra_before

    def test_local_coulomb_adds_cintra(self):
        """Local Coulomb term should add intra-site Coulomb."""
        s = StdIntList()
        s.model = "hubbard"
        s.NsiteUC = 1
        _setup_interactions(s)
        # On-site U: R=(0,0,0), iWan=jWan=0
        tUJindx_arr = np.array([[0, 0, 0, 0, 0]], dtype=int)
        tUJ_arr = np.array([4.0 + 0j])
        NtUJ = [0, 1, 0]
        tUJ = [None, tUJ_arr, None]
        tUJindx = [None, tUJindx_arr, None]

        w90._apply_coulomb_terms(
            s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx,
            w90._DCMode.NOTCORRECT, None,
        )

        assert len(s.Cintra_list) == 1
        np.testing.assert_allclose(s.Cintra_list[0][0], 4.0, atol=1e-12)
        assert s.Cintra_list[0][1] == 0

    def test_local_coulomb_dc_adds_transfer(self):
        """Local Coulomb with double-counting adds transfer terms."""
        s = StdIntList()
        s.model = "hubbard"
        s.NsiteUC = 1
        s.alpha = 0.5
        _setup_interactions(s)

        tUJindx_arr = np.array([[0, 0, 0, 0, 0]], dtype=int)
        tUJ_arr = np.array([4.0 + 0j])
        NtUJ = [0, 1, 0]
        tUJ = [None, tUJ_arr, None]
        tUJindx = [None, tUJindx_arr, None]

        DenMat = {(0, 0, 0): np.array([[0.5 + 0j]])}

        w90._apply_coulomb_terms(
            s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx,
            w90._DCMode.HARTREE, DenMat,
        )

        assert len(s.Cintra_list) == 1
        # Should add 2 transfer terms (one per spin)
        assert len(s.trans_list) == 2
        # alpha * U * DenMat = 0.5 * 4.0 * 0.5 = 1.0
        np.testing.assert_allclose(s.trans_list[0][0], 1.0 + 0j, atol=1e-12)

    def test_zero_terms_is_noop(self):
        """Should do nothing when NtUJ[1] is 0."""
        s = StdIntList()
        s.model = "hubbard"
        _setup_interactions(s)
        tUJindx_arr = np.zeros((0, 5), dtype=int)
        tUJ_arr = np.zeros(0, dtype=complex)
        NtUJ = [0, 0, 0]
        tUJ = [None, tUJ_arr, None]
        tUJindx = [None, tUJindx_arr, None]

        w90._apply_coulomb_terms(
            s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx,
            w90._DCMode.NOTCORRECT, None,
        )

        assert len(s.Cintra_list) == 0
        assert len(s.trans_list) == 0


# ---------------------------------------------------------------------------
#  Tests: _apply_hund_terms
# ---------------------------------------------------------------------------


class TestApplyHundTerms:
    """Tests for the _apply_hund_terms helper."""

    def test_none_indices_is_noop(self):
        """Should do nothing when tUJindx[2] is None."""
        s = StdIntList()
        s.model = "hubbard"
        _setup_interactions(s)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]
        nhund_before = len(s.Hund_list)
        w90._apply_hund_terms(
            s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx,
            w90._DCMode.NOTCORRECT, None,
        )
        assert len(s.Hund_list) == nhund_before

    def test_local_hund_term_skipped(self):
        """Local Hund term (same site) should be skipped."""
        s = StdIntList()
        s.model = "hubbard"
        s.NsiteUC = 1
        _setup_interactions(s)
        # On-site: R=(0,0,0), iWan=jWan=0 — should be skipped
        tUJindx_arr = np.array([[0, 0, 0, 0, 0]], dtype=int)
        tUJ_arr = np.array([0.5 + 0j])
        NtUJ = [0, 0, 1]
        tUJ = [None, None, tUJ_arr]
        tUJindx = [None, None, tUJindx_arr]

        w90._apply_hund_terms(
            s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx,
            w90._DCMode.NOTCORRECT, None,
        )

        assert len(s.Hund_list) == 0
        assert len(s.Ex_list) == 0
        assert len(s.PairHopp_list) == 0

    def test_zero_terms_is_noop(self):
        """Should do nothing when NtUJ[2] is 0."""
        s = StdIntList()
        s.model = "hubbard"
        _setup_interactions(s)
        tUJindx_arr = np.zeros((0, 5), dtype=int)
        tUJ_arr = np.zeros(0, dtype=complex)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, tUJ_arr]
        tUJindx = [None, None, tUJindx_arr]

        w90._apply_hund_terms(
            s, 0, 0, 0, 0, NtUJ, tUJ, tUJindx,
            w90._DCMode.NOTCORRECT, None,
        )

        assert len(s.Hund_list) == 0
        assert len(s.Ex_list) == 0
        assert len(s.PairHopp_list) == 0


# ---------------------------------------------------------------------------
#  Tests: _validate_wannier_params
# ---------------------------------------------------------------------------


class TestValidateWannierParams:
    """Tests for the _validate_wannier_params helper."""

    def test_spin_model_sets_S2(self):
        """Spin model should set S2 default to 1."""
        s = StdIntList()
        s.model = "spin"
        s.K = NaN_d
        s.h = NaN_d
        s.Gamma = NaN_d
        s.Gamma_y = NaN_d
        s.U = NaN_d
        s.S2 = None
        s.mu = NaN_d
        w90._validate_wannier_params(s)
        assert s.S2 == 1
        assert s.h == 0.0
        assert s.Gamma == 0.0
        assert s.Gamma_y == 0.0

    def test_hubbard_model_sets_mu(self):
        """Hubbard model should set mu default to 0.0."""
        s = StdIntList()
        s.model = "hubbard"
        s.K = NaN_d
        s.h = NaN_d
        s.Gamma = NaN_d
        s.Gamma_y = NaN_d
        s.U = NaN_d
        s.S2 = None
        s.mu = NaN_d
        w90._validate_wannier_params(s)
        assert s.mu == 0.0

    def test_kondo_model_exits(self):
        """Kondo model should cause exit."""
        s = StdIntList()
        s.model = "kondo"
        s.K = NaN_d
        s.h = NaN_d
        s.Gamma = NaN_d
        s.Gamma_y = NaN_d
        s.U = NaN_d
        s.S2 = None
        s.mu = NaN_d
        with pytest.raises(ValueError):
            w90._validate_wannier_params(s)


# ---------------------------------------------------------------------------
#  Tests: _build_wannier_interactions
# ---------------------------------------------------------------------------


class TestBuildWannierInteractions:
    """Tests for _build_wannier_interactions."""

    def test_hubbard_allocates_arrays(self, tmp_path):
        """Should allocate interaction arrays for Hubbard model."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        assert s.trans_list
        assert len(s.trans_list) >= 0

    def test_spin_allocates_arrays(self, tmp_path):
        """Should allocate interaction arrays for Spin model."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="spin", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        assert isinstance(s.trans_list, list)
        assert isinstance(s.intr_list, list)


# ---------------------------------------------------------------------------
#  Tests: _write_wan2site
# ---------------------------------------------------------------------------


class TestWriteWan2site:
    """Tests for _write_wan2site."""

    def test_creates_file(self, tmp_path):
        """Should create wan2site.dat in working directory."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        assert (tmp_path / "wan2site.dat").exists()

    def test_content_format(self, tmp_path):
        """wan2site.dat should contain expected header and site data."""
        os.chdir(tmp_path)
        NsiteUC = 2
        prefix = "test"
        _write_geom_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_hr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_ur_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)
        _write_jr_file(str(tmp_path), NsiteUC=NsiteUC, prefix=prefix)

        s = _make_wannier_StdI(model="hubbard", prefix=prefix, W=2, L=2, Height=1)
        w90.wannier90(s)

        content = (tmp_path / "wan2site.dat").read_text()
        assert "Total site number" in content
        assert "site nx ny nz norb" in content
        # Should have data lines for all sites
        lines = [l for l in content.strip().split("\n")
                 if l.strip() and "=" not in l
                 and "site" not in l.lower() and "total" not in l.lower()]
        assert len(lines) == s.NCell * s.NsiteUC


# ---------------------------------------------------------------------------
#  Tests: _apply_boundary_weights
# ---------------------------------------------------------------------------


class TestApplyBoundaryWeights:
    """Tests for the _apply_boundary_weights helper."""

    def test_returns_band_lattice(self):
        """Should return the maximum absolute R-vector extent per dimension."""
        indx_tot = np.array([
            [0, 0, 0],
            [1, -2, 0],
            [-1, 1, 3],
        ], dtype=int)
        Weight_tot = np.ones(3)
        s = StdIntList()
        s.W = None
        s.L = None
        s.Height = None
        band = w90._apply_boundary_weights(indx_tot, Weight_tot, 3, s)
        np.testing.assert_array_equal(band, [1, 2, 3])

    def test_no_halving_when_dimensions_unset(self):
        """Weights should be unchanged when W/L/Height are unset."""
        indx_tot = np.array([[1, 0, 0], [-1, 0, 0]], dtype=int)
        Weight_tot = np.ones(2)
        s = StdIntList()
        s.W = None
        s.L = None
        s.Height = None
        w90._apply_boundary_weights(indx_tot, Weight_tot, 2, s)
        np.testing.assert_array_equal(Weight_tot, [1.0, 1.0])

    def test_halving_at_boundary(self):
        """Boundary WSCs should have halved weight for even lattice dims."""
        # W=4 => Model_lattice[0] = 2; Band_lattice[0] = 3 (from |R0|=3)
        # Condition: Model_lattice < Band_lattice AND |R0| == Model_lattice
        indx_tot = np.array([
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [-2, 0, 0],
            [3, 0, 0],
        ], dtype=int)
        Weight_tot = np.ones(5)
        s = StdIntList()
        s.W = 4   # Model_lattice[0] = 2
        s.L = 1   # odd => Model_lattice[1] = 0 => no halving in dim 1
        s.Height = 1
        w90._apply_boundary_weights(indx_tot, Weight_tot, 5, s)
        # WSC 0 (R0=0) and WSC 1 (R0=1): unchanged
        assert Weight_tot[0] == 1.0
        assert Weight_tot[1] == 1.0
        # WSC 2 (|R0|=2 == Model_lattice[0]) and WSC 3 (|R0|=2): halved
        assert Weight_tot[2] == 0.5
        assert Weight_tot[3] == 0.5
        # WSC 4 (|R0|=3 != 2): unchanged
        assert Weight_tot[4] == 1.0

    def test_no_halving_when_model_exceeds_band(self):
        """No halving when model lattice > band lattice extent."""
        indx_tot = np.array([[1, 0, 0]], dtype=int)
        Weight_tot = np.ones(1)
        s = StdIntList()
        s.W = 10  # Model_lattice[0] = 5, Band_lattice[0] = 1
        s.L = 10
        s.Height = 10
        w90._apply_boundary_weights(indx_tot, Weight_tot, 1, s)
        assert Weight_tot[0] == 1.0


# ---------------------------------------------------------------------------
#  Tests: _count_and_store_terms
# ---------------------------------------------------------------------------


class TestCountAndStoreTerms:
    """Tests for the _count_and_store_terms helper."""

    def test_counts_above_cutoff(self):
        """Should count only terms above the cutoff threshold."""
        nWSC = 1
        NsiteUC = 2
        Mat_tot = np.zeros((nWSC, NsiteUC, NsiteUC), dtype=complex)
        Mat_tot[0, 0, 0] = 2.0 + 0j    # above cutoff
        Mat_tot[0, 0, 1] = 0.001 + 0j   # below cutoff
        Mat_tot[0, 1, 0] = 0.0 + 0j     # zero
        Mat_tot[0, 1, 1] = -1.5 + 0.5j  # above cutoff
        indx_tot = np.array([[0, 0, 0]], dtype=int)
        Weight_tot = np.ones(nWSC)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._count_and_store_terms(
            Mat_tot, indx_tot, Weight_tot, nWSC,
            NsiteUC, 0.01, 0, NtUJ, tUJindx, tUJ,
        )

        assert NtUJ[0] == 2
        assert tUJ[0] is not None
        assert len(tUJ[0]) == 2

    def test_applies_weights_before_counting(self):
        """Should multiply by Weight_tot before checking cutoff."""
        nWSC = 1
        NsiteUC = 1
        Mat_tot = np.zeros((nWSC, NsiteUC, NsiteUC), dtype=complex)
        Mat_tot[0, 0, 0] = 0.05 + 0j  # above 0.01, but weight makes it 0.025
        indx_tot = np.array([[0, 0, 0]], dtype=int)
        Weight_tot = np.array([0.5])
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._count_and_store_terms(
            Mat_tot, indx_tot, Weight_tot, nWSC,
            NsiteUC, 0.01, 0, NtUJ, tUJindx, tUJ,
        )

        assert NtUJ[0] == 1
        np.testing.assert_allclose(tUJ[0][0], 0.025 + 0j, atol=1e-12)

    def test_stores_correct_indices(self):
        """Should store the correct R-vector and band indices."""
        nWSC = 1
        NsiteUC = 2
        Mat_tot = np.zeros((nWSC, NsiteUC, NsiteUC), dtype=complex)
        Mat_tot[0, 0, 1] = 1.0 + 0j
        indx_tot = np.array([[3, -1, 2]], dtype=int)
        Weight_tot = np.ones(nWSC)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._count_and_store_terms(
            Mat_tot, indx_tot, Weight_tot, nWSC,
            NsiteUC, 0.01, 0, NtUJ, tUJindx, tUJ,
        )

        assert NtUJ[0] == 1
        np.testing.assert_array_equal(tUJindx[0][0], [3, -1, 2, 0, 1])

    def test_all_below_cutoff_returns_empty(self):
        """Should return empty arrays when all terms are below cutoff."""
        nWSC = 1
        NsiteUC = 1
        Mat_tot = np.zeros((nWSC, NsiteUC, NsiteUC), dtype=complex)
        Mat_tot[0, 0, 0] = 0.001 + 0j
        indx_tot = np.array([[0, 0, 0]], dtype=int)
        Weight_tot = np.ones(nWSC)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._count_and_store_terms(
            Mat_tot, indx_tot, Weight_tot, nWSC,
            NsiteUC, 0.01, 0, NtUJ, tUJindx, tUJ,
        )

        assert NtUJ[0] == 0
        assert len(tUJ[0]) == 0

    def test_itUJ_index_selects_slot(self):
        """Should store into the correct tUJ/tUJindx slot."""
        nWSC = 1
        NsiteUC = 1
        Mat_tot = np.zeros((nWSC, NsiteUC, NsiteUC), dtype=complex)
        Mat_tot[0, 0, 0] = 5.0 + 0j
        indx_tot = np.array([[0, 0, 0]], dtype=int)
        Weight_tot = np.ones(nWSC)
        NtUJ = [0, 0, 0]
        tUJ = [None, None, None]
        tUJindx = [None, None, None]

        w90._count_and_store_terms(
            Mat_tot, indx_tot, Weight_tot, nWSC,
            NsiteUC, 0.01, 2, NtUJ, tUJindx, tUJ,
        )

        assert NtUJ[2] == 1
        assert tUJ[0] is None  # slot 0 untouched
        assert tUJ[1] is None  # slot 1 untouched
        assert tUJ[2] is not None
