"""Unit tests for wannier90 module.

Tests for the Python translation of Wannier90.c.
"""
from __future__ import annotations

import math
import os
import textwrap

import numpy as np
import pytest

from stdface_vals import StdIntList
import wannier90 as w90


# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------

NaN_d = float("nan")
NaN_i = 2147483647
NaN_c = complex(float("nan"), 0.0)


# ---------------------------------------------------------------------------
#  Helpers: inverse matrix / in-box
# ---------------------------------------------------------------------------


class TestCalcInverseMatrix:
    """Tests for the internal _calc_inverse_matrix helper."""

    def test_identity(self):
        """Inverse of the identity matrix should be the identity matrix."""
        mat = np.eye(3)
        inv = w90._calc_inverse_matrix(mat)
        np.testing.assert_allclose(inv, np.eye(3), atol=1e-12)

    def test_diagonal(self):
        """Inverse of a diagonal matrix should have reciprocal entries."""
        mat = np.diag([2.0, 4.0, 5.0])
        inv = w90._calc_inverse_matrix(mat)
        expected = np.diag([0.5, 0.25, 0.2])
        np.testing.assert_allclose(inv, expected, atol=1e-12)

    def test_general_3x3(self):
        """Inverse of a general non-singular 3x3 matrix."""
        mat = np.array([[1.0, 2.0, 3.0],
                        [0.0, 1.0, 4.0],
                        [5.0, 6.0, 0.0]])
        inv = w90._calc_inverse_matrix(mat)
        product = mat @ inv
        np.testing.assert_allclose(product, np.eye(3), atol=1e-12)

    def test_roundtrip(self):
        """M @ inv(M) should be the identity for a random matrix."""
        rng = np.random.default_rng(42)
        mat = rng.random((3, 3)) + np.eye(3)
        inv = w90._calc_inverse_matrix(mat)
        np.testing.assert_allclose(mat @ inv, np.eye(3), atol=1e-10)


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
    s.S2 = NaN_i

    return s


# ---------------------------------------------------------------------------
#  Tests: module structure
# ---------------------------------------------------------------------------


class TestModuleStructure:
    """Verify module can be imported and has expected attributes."""

    def test_import(self):
        """wannier90 module should import without error."""
        import wannier90
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

        with pytest.raises(SystemExit):
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
        s.NaN_i = NaN_i
        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = NaN_i
        s.L = NaN_i
        s.Height = NaN_i

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
        s.NaN_i = NaN_i

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
        s.NaN_i = NaN_i
        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = NaN_i
        s.L = NaN_i
        s.Height = NaN_i

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
        s.NaN_i = NaN_i
        s.NsiteUC = NsiteUC
        s.direct = np.eye(3)
        s.tau = np.zeros((NsiteUC, 3))
        s.W = NaN_i
        s.L = NaN_i
        s.Height = NaN_i

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

        with pytest.raises(SystemExit):
            w90._read_density_matrix(s, str(tmp_path / "nonexistent_dr.dat"))


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

        assert (tmp_path / "lattice.xsf").exists()
        assert (tmp_path / "wan2site.dat").exists()
        assert (tmp_path / "geometry.dat").exists()

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

        assert s.transindx is not None
        assert s.trans is not None

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

        with pytest.raises(SystemExit):
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

        with pytest.raises(SystemExit):
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

        with pytest.raises(SystemExit):
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

        with pytest.raises(SystemExit):
            w90.wannier90(s)
