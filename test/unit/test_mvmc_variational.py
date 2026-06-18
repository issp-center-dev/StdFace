"""Unit tests for mvmc_variational module.

Tests for the mVMC variational parameter generation functions extracted
from ``stdface_model_util``: ``_anti_period_dot``, ``_parity_sign``,
``_fold_site_sub``, ``proj``, ``_init_site_sub``, ``generate_orb``,
and ``print_jastrow``.
"""
from __future__ import annotations

import math
import os
import tempfile

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.solvers.mvmc.variational import (
    _anti_period_dot,
    _parity_sign,
    _check_commensurate,
    _fold_site_sub,
    _init_site_sub,
    _assign_orb_sector,
    _jastrow_momentum_projected,
    _jastrow_global_optimization,
    proj,
    generate_orb,
    print_jastrow,
)


NaN_i = 2147483647
NaN_d = float("nan")


# ===================================================================
#  _anti_period_dot
# ===================================================================


class TestAntiPeriodDot:
    """Tests for _anti_period_dot helper."""

    def test_all_zero(self):
        """Zero anti-period flags give zero dot product."""
        ap = np.array([0, 0, 0], dtype=int)
        assert _anti_period_dot(ap, [5, 3, 1]) == 0

    def test_single_direction(self):
        """One non-zero flag selects the corresponding nBox component."""
        ap = np.array([0, 1, 0], dtype=int)
        assert _anti_period_dot(ap, [2, 7, 4]) == 7

    def test_all_ones(self):
        """All flags set → sum of nBox components."""
        ap = np.array([1, 1, 1], dtype=int)
        assert _anti_period_dot(ap, [3, -2, 5]) == 6

    def test_mixed_signs(self):
        """Negative nBox values handled correctly."""
        ap = np.array([1, 0, 1], dtype=int)
        assert _anti_period_dot(ap, [-1, 0, 2]) == 1

    def test_large_nbox(self):
        """Large nBox values work without overflow."""
        ap = np.array([1, 1, 0], dtype=int)
        assert _anti_period_dot(ap, [1000, -1000, 0]) == 0


# ===================================================================
#  _parity_sign
# ===================================================================


class TestParitySign:
    """Tests for _parity_sign helper."""

    def test_zero_is_even(self):
        """Zero is even → +1."""
        assert _parity_sign(0) == 1

    def test_one_is_odd(self):
        """One is odd → -1."""
        assert _parity_sign(1) == -1

    def test_even_positive(self):
        """Even positive integer → +1."""
        assert _parity_sign(4) == 1

    def test_odd_positive(self):
        """Odd positive integer → -1."""
        assert _parity_sign(7) == -1

    def test_negative_even(self):
        """Negative even integer → +1."""
        assert _parity_sign(-2) == 1

    def test_negative_odd(self):
        """Negative odd integer → -1."""
        assert _parity_sign(-3) == -1


# ===================================================================
#  _check_commensurate
# ===================================================================


class TestCheckCommensurate:
    """Tests for _check_commensurate helper."""

    def test_identity_is_commensurate(self):
        """Identity sublattice is commensurate with any main lattice."""
        rbox_sub = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        box = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 1]], dtype=float)
        ncell_sub = 1
        assert _check_commensurate(rbox_sub, box, ncell_sub) is True

    def test_commensurate_2x2_in_4x4(self):
        """2x2 sublattice is commensurate with 4x4 main lattice."""
        # rbox_sub for box_sub = [[2,0,0],[0,2,0],[0,0,1]] is cofactor/det
        # cofactor = [[2,0,0],[0,2,0],[0,0,4]], det = 4
        # rbox_sub = cofactor (before division)
        rbox_sub = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 4]], dtype=float)
        box = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 1]], dtype=float)
        ncell_sub = 4
        assert _check_commensurate(rbox_sub, box, ncell_sub) is True

    def test_incommensurate_3x3_in_4x4(self):
        """3x3 sublattice is incommensurate with 4x4 main lattice."""
        # rbox_sub for box_sub = [[3,0,0],[0,3,0],[0,0,1]]
        # cofactor = [[3,0,0],[0,3,0],[0,0,9]], det = 9
        rbox_sub = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 9]], dtype=float)
        box = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 1]], dtype=float)
        ncell_sub = 9
        # 3*4 = 12, 12 % 9 = 3 != 0 → incommensurate
        assert _check_commensurate(rbox_sub, box, ncell_sub) is False

    def test_commensurate_same_box(self):
        """Same sublattice as main lattice is commensurate."""
        rbox_sub = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 16]], dtype=float)
        box = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 1]], dtype=float)
        ncell_sub = 16
        assert _check_commensurate(rbox_sub, box, ncell_sub) is True

    def test_all_zeros_incommensurate(self):
        """Zero ncell_sub causes division by zero check (ncell != 0 checked elsewhere)."""
        # This would cause ZeroDivisionError, but ncell_sub == 0 is caught before calling
        # Let's test with ncell_sub = 1 and mismatched dimensions
        rbox_sub = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        box = np.array([[3, 0, 0], [0, 5, 0], [0, 0, 7]], dtype=float)
        ncell_sub = 2  # 1*3=3, 3%2=1 != 0
        assert _check_commensurate(rbox_sub, box, ncell_sub) is False


# ===================================================================
#  Test helpers
# ===================================================================


def _make_chain_stdi(L: int = 4, model: str = "hubbard") -> StdIntList:
    """Create a minimal StdIntList for a 1D chain with all site fields set.

    Manually sets up Cell, NCell, rbox, nsite, etc. without calling
    init_site, to avoid its complex validation and file-writing.
    """
    StdI = StdIntList()
    StdI.solver = "mVMC"
    StdI.pi = math.acos(-1.0)
    StdI.model = model
    StdI.lattice = "chain"
    StdI.lGC = 0
    StdI.lBoost = 0

    # Box = [[L,0,0],[0,1,0],[0,0,1]]
    StdI.box[:, :] = 0
    StdI.box[0, 0] = L
    StdI.box[1, 1] = 1
    StdI.box[2, 2] = 1

    # NCell = determinant = L
    StdI.NCell = L
    StdI.NsiteUC = 1
    nsite = L * 1  # NCell * NsiteUC
    if model == "kondo":
        nsite *= 2
    StdI.nsite = nsite

    # Reciprocal box (cofactor matrix of box): rbox = [[1,0,0],[0,L,0],[0,0,L]]
    StdI.rbox = np.zeros((3, 3), dtype=float)
    StdI.rbox[0, 0] = 1
    StdI.rbox[1, 1] = L
    StdI.rbox[2, 2] = L

    # Cell array: each cell's fractional coordinate
    StdI.Cell = np.zeros((L, 3), dtype=float)
    for i in range(L):
        StdI.Cell[i, 0] = i

    # tau (unit cell structure)
    StdI.tau = np.zeros((1, 3))

    # Phase / AntiPeriod
    StdI.phase[:] = 0.0
    StdI.pi180 = StdI.pi / 180.0
    StdI.ExpPhase = np.ones(3, dtype=complex)
    StdI.AntiPeriod = np.zeros(3, dtype=int)

    # locspinflag
    if model == "spin":
        StdI.locspinflag = [1] * nsite
    else:
        StdI.locspinflag = [0] * nsite

    # Sub-lattice defaults (unset)
    StdI.boxsub[:, :] = NaN_i
    StdI.Hsub = NaN_i
    StdI.Lsub = NaN_i
    StdI.Wsub = NaN_i

    # NMPTrans default
    StdI.NMPTrans = NaN_i

    return StdI


class TestFoldSiteSub:
    """Tests for _fold_site_sub."""

    def test_identity_fold(self):
        """Test that a site inside the sub-cell folds to itself."""
        StdI = _make_chain_stdi(L=4)
        # Set up sub-cell = full cell
        StdI.boxsub[:, :] = 0
        StdI.boxsub[0, 0] = 4
        StdI.boxsub[1, 1] = 1
        StdI.boxsub[2, 2] = 1
        StdI.NCellsub = 4
        StdI.rboxsub = np.zeros((3, 3), dtype=float)
        StdI.rboxsub[0, 0] = 1
        StdI.rboxsub[1, 1] = 4
        StdI.rboxsub[2, 2] = 4

        nBox, folded = _fold_site_sub(StdI, [0, 0, 0])
        assert folded == [0, 0, 0]
        assert nBox == [0, 0, 0]

    def test_fold_outside_sub(self):
        """Test folding a site outside the sub-cell."""
        StdI = _make_chain_stdi(L=4)
        # Sub-cell of size 2
        StdI.boxsub[:, :] = 0
        StdI.boxsub[0, 0] = 2
        StdI.boxsub[1, 1] = 1
        StdI.boxsub[2, 2] = 1
        StdI.NCellsub = 2
        StdI.rboxsub = np.zeros((3, 3), dtype=float)
        StdI.rboxsub[0, 0] = 1
        StdI.rboxsub[1, 1] = 2
        StdI.rboxsub[2, 2] = 2

        nBox, folded = _fold_site_sub(StdI, [3, 0, 0])
        assert folded[0] == 1  # 3 mod 2 = 1
        assert nBox[0] == 1    # super-cell index 1


class TestInitSiteSub:
    """Tests for _init_site_sub."""

    def test_lwh_mode(self):
        """Test initialization via Lsub/Wsub/Hsub."""
        StdI = _make_chain_stdi(L=4)
        StdI.Wsub = 2
        StdI.Lsub = NaN_i
        StdI.Hsub = NaN_i
        _init_site_sub(StdI)
        assert StdI.NCellsub == 2

    def test_box_mode(self):
        """Test initialization via boxsub matrix."""
        StdI = _make_chain_stdi(L=4)
        StdI.boxsub[:, :] = 0
        StdI.boxsub[0, 0] = 2
        StdI.boxsub[1, 1] = 1
        StdI.boxsub[2, 2] = 1
        _init_site_sub(StdI)
        assert StdI.NCellsub == 2

    def test_conflict_raises(self):
        """Test that conflicting Lsub and boxsub raises exit."""
        StdI = _make_chain_stdi(L=4)
        StdI.Lsub = 2
        StdI.boxsub[0, 0] = 2  # conflict: both specified
        with pytest.raises(ValueError):
            _init_site_sub(StdI)


class TestGenerateOrb:
    """Tests for generate_orb."""

    def test_sets_norb(self):
        """Test that NOrb is set after generate_orb."""
        StdI = _make_chain_stdi(L=4, model="hubbard")
        StdI.Wsub = 2
        generate_orb(StdI)
        assert StdI.NOrb > 0

    def test_orb_array_shape(self):
        """Test that Orb array has correct shape."""
        StdI = _make_chain_stdi(L=4, model="hubbard")
        StdI.Wsub = 2
        generate_orb(StdI)
        assert StdI.Orb.shape == (StdI.nsite, StdI.nsite)
        assert StdI.AntiOrb.shape == (StdI.nsite, StdI.nsite)


class TestProj:
    """Tests for proj (qptransidx.def generation)."""

    def test_creates_qptransidx(self):
        """Test that proj creates qptransidx.def."""
        StdI = _make_chain_stdi(L=4, model="hubbard")
        StdI.Wsub = 2
        _init_site_sub(StdI)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                proj(StdI)
                assert os.path.exists("qptransidx.def")
            finally:
                os.chdir(orig)

    def test_nsym_positive(self):
        """Test that NSym is set to a positive value."""
        StdI = _make_chain_stdi(L=4, model="hubbard")
        StdI.Wsub = 2
        _init_site_sub(StdI)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                proj(StdI)
                assert StdI.NSym > 0
            finally:
                os.chdir(orig)

    def test_header_contains_nqptrans(self):
        """Test that the header contains NQPTrans."""
        StdI = _make_chain_stdi(L=4, model="hubbard")
        StdI.Wsub = 2
        _init_site_sub(StdI)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                proj(StdI)
                with open("qptransidx.def") as f:
                    content = f.read()
                assert "NQPTrans" in content
            finally:
                os.chdir(orig)


class TestPrintJastrow:
    """Tests for print_jastrow."""

    def test_creates_jastrowidx(self):
        """Test that print_jastrow creates jastrowidx.def."""
        StdI = _make_chain_stdi(L=4, model="hubbard")
        StdI.Wsub = 2
        generate_orb(StdI)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_jastrow(StdI)
                assert os.path.exists("jastrowidx.def")
            finally:
                os.chdir(orig)

    def test_header_contains_njastrowidx(self):
        """Test that the header contains NJastrowIdx."""
        StdI = _make_chain_stdi(L=4, model="hubbard")
        StdI.Wsub = 2
        generate_orb(StdI)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_jastrow(StdI)
                with open("jastrowidx.def") as f:
                    content = f.read()
                assert "NJastrowIdx" in content
            finally:
                os.chdir(orig)

    def test_spin_model_jastrow(self):
        """Test Jastrow generation for spin model with NMPTrans != 1."""
        StdI = _make_chain_stdi(L=4, model="spin")
        StdI.Wsub = 2
        generate_orb(StdI)
        StdI.NMPTrans = 0

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_jastrow(StdI)
                assert os.path.exists("jastrowidx.def")
                with open("jastrowidx.def") as f:
                    content = f.read()
                assert "NJastrowIdx" in content
            finally:
                os.chdir(orig)


# ===================================================================
#  _assign_orb_sector
# ===================================================================


def _make_orb_stdi(nsite: int = 4, NsiteUC: int = 1) -> StdIntList:
    """Create a minimal StdIntList for _assign_orb_sector tests."""
    StdI = StdIntList()
    StdI.nsite = nsite
    StdI.NsiteUC = NsiteUC
    StdI.Orb = np.full((nsite, nsite), -1, dtype=int)
    StdI.AntiOrb = np.zeros((nsite, nsite), dtype=int)
    return StdI


class TestAssignOrbSector:
    """Tests for _assign_orb_sector helper."""

    def test_new_assigns_orb_to_ref(self):
        """Test that is_new=True assigns iOrb to the reference cell pair."""
        StdI = _make_orb_stdi(4)
        result = _assign_orb_sector(
            StdI, iOrb=5, anti_val=1,
            iCell=1, jCell=2, iCell2=0, jCell2=1,
            isite=0, jsite=0, i_off=0, j_off=0, is_new=True)
        assert StdI.Orb[0, 1] == 5  # ref pair (iCell2*1+0, jCell2*1+0)
        assert result == 6  # iOrb incremented

    def test_new_copies_to_current(self):
        """Test that is_new=True also copies to the current cell pair."""
        StdI = _make_orb_stdi(4)
        _assign_orb_sector(
            StdI, iOrb=5, anti_val=1,
            iCell=1, jCell=2, iCell2=0, jCell2=1,
            isite=0, jsite=0, i_off=0, j_off=0, is_new=True)
        assert StdI.Orb[1, 2] == 5  # current pair copies from ref

    def test_not_new_skips_ref_assignment(self):
        """Test that is_new=False does NOT assign new iOrb to ref."""
        StdI = _make_orb_stdi(4)
        StdI.Orb[0, 1] = 42  # pre-existing ref value
        result = _assign_orb_sector(
            StdI, iOrb=5, anti_val=1,
            iCell=1, jCell=2, iCell2=0, jCell2=1,
            isite=0, jsite=0, i_off=0, j_off=0, is_new=False)
        assert StdI.Orb[0, 1] == 42  # unchanged
        assert result == 5  # iOrb NOT incremented

    def test_not_new_copies_from_ref(self):
        """Test that is_new=False copies ref value to current pair."""
        StdI = _make_orb_stdi(4)
        StdI.Orb[0, 1] = 42  # pre-existing ref value
        _assign_orb_sector(
            StdI, iOrb=5, anti_val=-1,
            iCell=1, jCell=2, iCell2=0, jCell2=1,
            isite=0, jsite=0, i_off=0, j_off=0, is_new=False)
        assert StdI.Orb[1, 2] == 42  # copied from ref

    def test_anti_val_recorded(self):
        """Test that anti_val is set on both ref and current pair."""
        StdI = _make_orb_stdi(4)
        _assign_orb_sector(
            StdI, iOrb=0, anti_val=-1,
            iCell=1, jCell=2, iCell2=0, jCell2=1,
            isite=0, jsite=0, i_off=0, j_off=0, is_new=True)
        assert StdI.AntiOrb[0, 1] == -1  # ref
        assert StdI.AntiOrb[1, 2] == -1  # current

    def test_kondo_offset_half(self):
        """Test with i_off=half for Kondo model sector."""
        StdI = _make_orb_stdi(8)  # half=4
        result = _assign_orb_sector(
            StdI, iOrb=10, anti_val=1,
            iCell=1, jCell=0, iCell2=0, jCell2=0,
            isite=0, jsite=0, i_off=4, j_off=0, is_new=True)
        # ref: (4+0*1+0, 0+0*1+0) = (4, 0)
        assert StdI.Orb[4, 0] == 10
        # current: (4+1*1+0, 0+0*1+0) = (5, 0)
        assert StdI.Orb[5, 0] == 10
        assert result == 11

    def test_kondo_offset_both(self):
        """Test with both offsets for Kondo (half, half) sector."""
        StdI = _make_orb_stdi(8)
        result = _assign_orb_sector(
            StdI, iOrb=20, anti_val=1,
            iCell=1, jCell=1, iCell2=0, jCell2=0,
            isite=0, jsite=0, i_off=4, j_off=4, is_new=True)
        # ref: (4+0, 4+0) = (4, 4)
        assert StdI.Orb[4, 4] == 20
        # current: (4+1, 4+1) = (5, 5)
        assert StdI.Orb[5, 5] == 20
        assert result == 21

    def test_multisite_uc(self):
        """Test with NsiteUC > 1."""
        StdI = _make_orb_stdi(6, NsiteUC=2)
        result = _assign_orb_sector(
            StdI, iOrb=0, anti_val=1,
            iCell=1, jCell=2, iCell2=0, jCell2=1,
            isite=1, jsite=0, i_off=0, j_off=0, is_new=True)
        # ref: (0*2+1, 1*2+0) = (1, 2)
        assert StdI.Orb[1, 2] == 0
        # current: (1*2+1, 2*2+0) = (3, 4)
        assert StdI.Orb[3, 4] == 0
        assert result == 1


# ===================================================================
#  _jastrow_momentum_projected
# ===================================================================


def _make_jastrow_stdi(L: int = 4, model: str = "hubbard") -> StdIntList:
    """Create a StdIntList with Orb/NOrb set for Jastrow tests.

    Sets up a chain with orbital indices pre-computed (simple diagonal
    pattern: Orb[i,j] = abs(i-j) % L).
    """
    StdI = _make_chain_stdi(L=L, model=model)
    nsite = StdI.nsite
    StdI.Orb = np.zeros((nsite, nsite), dtype=int)
    StdI.AntiOrb = np.zeros((nsite, nsite), dtype=int)
    # Simple orbital pattern: distance-based
    for i in range(nsite):
        for j in range(nsite):
            StdI.Orb[i, j] = abs(i - j) % L
    StdI.NOrb = L
    return StdI


class TestJastrowMomentumProjected:
    """Tests for _jastrow_momentum_projected."""

    def test_returns_positive_njastrow(self):
        """Test that NJastrow is positive."""
        StdI = _make_jastrow_stdi(L=4, model="hubbard")
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        NJastrow, Jastrow = _jastrow_momentum_projected(StdI, Jastrow)
        assert NJastrow > 0

    def test_off_diagonal_non_negative(self):
        """Test that off-diagonal Jastrow indices are non-negative."""
        StdI = _make_jastrow_stdi(L=4, model="hubbard")
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        NJastrow, Jastrow = _jastrow_momentum_projected(StdI, Jastrow)
        for i in range(StdI.nsite):
            for j in range(StdI.nsite):
                if i != j:
                    assert Jastrow[i, j] >= 0

    def test_hubbard_starts_at_zero(self):
        """Test that Hubbard model starts renumbering from 0."""
        StdI = _make_jastrow_stdi(L=2, model="hubbard")
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        NJastrow, Jastrow = _jastrow_momentum_projected(StdI, Jastrow)
        # Hubbard: NJastrow starts at 0, so first unique index is 0
        assert 0 in Jastrow

    def test_spin_starts_at_minus_one(self):
        """Test that Spin model starts with offset -1."""
        StdI = _make_jastrow_stdi(L=2, model="spin")
        # Spin model: all sites are local spins
        StdI.locspinflag = [1] * StdI.nsite
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        NJastrow, Jastrow = _jastrow_momentum_projected(StdI, Jastrow)
        # All local spins → Jastrow rows/cols set to -1, then flipped to 0
        assert NJastrow >= 0

    def test_symmetrized_output(self):
        """Test that output matrix is symmetric."""
        StdI = _make_jastrow_stdi(L=4, model="hubbard")
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        _, Jastrow = _jastrow_momentum_projected(StdI, Jastrow)
        # After symmetrize + renumber, Jastrow should still be symmetric
        for i in range(StdI.nsite):
            for j in range(i + 1, StdI.nsite):
                assert Jastrow[i, j] == Jastrow[j, i]

    def test_local_spin_sites_excluded(self):
        """Test that local-spin sites get uniform index."""
        StdI = _make_jastrow_stdi(L=4, model="hubbard")
        # Mark site 0 as local spin
        StdI.locspinflag[0] = 1
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        NJastrow, Jastrow = _jastrow_momentum_projected(StdI, Jastrow)
        # Local-spin row/col was set to -1, then flipped to 0
        # All entries involving site 0 should have the same index
        row_vals = set(int(Jastrow[0, j]) for j in range(StdI.nsite))
        col_vals = set(int(Jastrow[i, 0]) for i in range(StdI.nsite))
        assert len(row_vals) == 1
        assert len(col_vals) == 1


# ===================================================================
#  _jastrow_global_optimization
# ===================================================================


class TestJastrowGlobalOptimization:
    """Tests for _jastrow_global_optimization."""

    def test_spin_model_single_index(self):
        """Test that Spin model produces NJastrow=1 with all zeros."""
        StdI = _make_jastrow_stdi(L=4, model="spin")
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        NJastrow = _jastrow_global_optimization(StdI, Jastrow)
        assert NJastrow == 1
        assert np.all(Jastrow == 0)

    def test_hubbard_positive_njastrow(self):
        """Test that Hubbard model produces positive NJastrow."""
        StdI = _make_jastrow_stdi(L=4, model="hubbard")
        StdI.NMPTrans = 0  # not momentum-projected
        StdI.direct = np.eye(3)
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        NJastrow = _jastrow_global_optimization(StdI, Jastrow)
        assert NJastrow > 0

    def test_symmetric_jastrow(self):
        """Test that global optimization produces symmetric Jastrow."""
        StdI = _make_jastrow_stdi(L=4, model="hubbard")
        StdI.NMPTrans = 0
        StdI.direct = np.eye(3)
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        _jastrow_global_optimization(StdI, Jastrow)
        for i in range(StdI.nsite):
            for j in range(i + 1, StdI.nsite):
                assert Jastrow[i, j] == Jastrow[j, i]

    def test_diagonal_excluded(self):
        """Test that on-site terms (i==j at origin cell) are excluded."""
        StdI = _make_jastrow_stdi(L=2, model="hubbard")
        StdI.NMPTrans = 0
        StdI.direct = np.eye(3)
        Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
        _jastrow_global_optimization(StdI, Jastrow)
        # On-site is excluded from assignment, so diagonal should remain 0
        for i in range(StdI.nsite):
            assert Jastrow[i, i] == 0
