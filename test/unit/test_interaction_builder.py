"""Unit tests for interaction_builder module.

Tests for ``trans``, ``hopping``, ``hubbard_local``, ``mag_field``,
``intr``, ``general_j``, ``coulomb``, and ``malloc_interactions``.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList, MethodType, SolverType
from stdface.lattice.interaction_builder import (
    trans,
    hopping,
    hubbard_local,
    mag_field,
    intr,
    general_j,
    coulomb,
    compute_max_interactions,
    malloc_interactions,
    add_neighbor_interaction,
    add_neighbor_interaction_3d,
    add_local_terms,
    _spin_ladder_factor,
    _add_spin_half_terms,
)


def _make_stdi(
    solver: str = "HPhi",
    model: str = "hubbard",
    ntrans: int = 100,
    nintr: int = 100,
) -> StdIntList:
    """Create a minimal StdIntList with pre-allocated interaction arrays."""
    StdI = StdIntList()
    StdI.solver = solver
    StdI.model = model
    StdI.method = "lanczos"
    StdI.PumpBody = 0
    malloc_interactions(StdI, ntrans, nintr)
    return StdI


# ===================================================================
#  malloc_interactions
# ===================================================================


class TestMallocInteractions:
    """Tests for malloc_interactions."""

    def test_trans_arrays_allocated(self):
        """Test that transfer arrays are allocated with correct shapes."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        StdI.model = "hubbard"
        StdI.method = "lanczos"
        StdI.PumpBody = 0
        malloc_interactions(StdI, 50, 80)
        assert len(StdI.trans_list) == 0

    def test_intr_arrays_allocated(self):
        """Test that interaction arrays are allocated with correct shapes."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        StdI.model = "hubbard"
        StdI.method = "lanczos"
        StdI.PumpBody = 0
        malloc_interactions(StdI, 50, 80)
        assert len(StdI.intr_list) == 0

    def test_coulomb_arrays_allocated(self):
        """Test that Coulomb arrays are allocated."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        StdI.model = "hubbard"
        StdI.method = "lanczos"
        StdI.PumpBody = 0
        malloc_interactions(StdI, 50, 80)
        assert len(StdI.Cintra_list) == 0
        assert len(StdI.Cinter_list) == 0

    def test_exchange_arrays_allocated(self):
        """Test that exchange and pair arrays are allocated."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        StdI.model = "hubbard"
        StdI.method = "lanczos"
        StdI.PumpBody = 0
        malloc_interactions(StdI, 50, 80)
        assert len(StdI.Hund_list) == 0
        assert len(StdI.Ex_list) == 0
        assert len(StdI.PairLift_list) == 0
        assert len(StdI.PairHopp_list) == 0

    def test_pump_arrays_for_timeevolution(self):
        """Test that pump arrays are allocated for HPhi time-evolution."""
        StdI = StdIntList()
        StdI.solver = "HPhi"
        StdI.model = "hubbard"
        StdI.method = "timeevolution"
        StdI.PumpBody = 1
        StdI.Lanczos_max = 10
        malloc_interactions(StdI, 50, 80)
        assert StdI.npump.shape == (10,)
        assert StdI.pumpindx.shape == (10, 50, 4)
        assert StdI.pump.shape == (10, 50)


# ===================================================================
#  trans
# ===================================================================


class TestTrans:
    """Tests for trans."""

    def test_adds_transfer_term(self):
        """Test that a transfer term is added correctly."""
        StdI = _make_stdi()
        trans(StdI, 1.0 + 0j, 0, 0, 1, 0)
        assert len(StdI.trans_list) == 1
        assert StdI.trans_list[0][0] == 1.0 + 0j
        assert list(StdI.trans_list[0][1:]) == [0, 0, 1, 0]

    def test_skips_small_value(self):
        """Test that values below threshold are skipped."""
        StdI = _make_stdi()
        trans(StdI, 1e-13, 0, 0, 1, 0)
        assert len(StdI.trans_list) == 0

    def test_increments_counter(self):
        """Test that multiple calls increment the counter."""
        StdI = _make_stdi()
        trans(StdI, 1.0, 0, 0, 1, 0)
        trans(StdI, 2.0, 1, 1, 0, 1)
        assert len(StdI.trans_list) == 2
        assert StdI.trans_list[1][0] == 2.0

    def test_complex_value(self):
        """Test that complex transfer values are stored."""
        StdI = _make_stdi()
        trans(StdI, 1.0 + 0.5j, 0, 0, 1, 1)
        assert StdI.trans_list[0][0] == pytest.approx(1.0 + 0.5j)


class TestHoppingPump:
    """Covers ``hopping`` when pump arrays are used instead of ``trans``."""

    def test_pump_mode_populates_pump_arrays(self):
        """``interaction_builder.hopping`` pump branch (time evolution + PumpBody)."""
        StdI = StdIntList()
        StdI.solver = SolverType.HPhi
        StdI.model = "hubbard"
        StdI.method = MethodType.TIME_EVOLUTION
        StdI.PumpBody = 1
        StdI.Lanczos_max = 1
        StdI.At = np.array([[0.1, 0.0, 0.0]], dtype=float)
        malloc_interactions(StdI, 50, 80)
        StdI.npump[:] = 0
        dR = np.array([1.0, 0.0, 0.0])
        hopping(StdI, 0.5 + 0.25j, 0, 1, dR)
        # Two spins × two entries (forward + conjugate) per time slice
        assert StdI.npump[0] == 4
        assert StdI.pump[0][0] != 0.0
        assert list(StdI.pumpindx[0][0]) == [0, 0, 1, 0]
        assert list(StdI.pumpindx[0][1]) == [1, 0, 0, 0]


# ===================================================================
#  hopping (default branch)
# ===================================================================


class TestHopping:
    """Tests for hopping."""

    def test_adds_both_spins(self):
        """Test that hopping adds terms for both spin channels."""
        StdI = _make_stdi()
        dR = np.array([1.0, 0.0, 0.0])
        hopping(StdI, 1.0 + 0j, 0, 1, dR)
        # 2 spins × 2 directions = 4 terms
        assert len(StdI.trans_list) == 4

    def test_hermitian_conjugate(self):
        """Test that hopping adds conjugate pairs."""
        StdI = _make_stdi()
        dR = np.array([1.0, 0.0, 0.0])
        hopping(StdI, 1.0 + 0.5j, 0, 1, dR)
        # For spin 0: trans(t, j, 0, i, 0) and trans(conj(t), i, 0, j, 0)
        assert StdI.trans_list[0][0] == pytest.approx(1.0 + 0.5j)
        assert StdI.trans_list[1][0] == pytest.approx(1.0 - 0.5j)


# ===================================================================
#  hubbard_local
# ===================================================================


class TestHubbardLocal:
    """Tests for hubbard_local."""

    def test_chemical_potential(self):
        """Test that chemical potential adds transfer terms."""
        StdI = _make_stdi()
        hubbard_local(StdI, mu0=1.0, h0=0.0, Gamma0=0.0, Gamma0_y=0.0,
                      U0=0.0, isite=0)
        # mu0 contributes 2 terms (spin up and down)
        assert len(StdI.trans_list) == 2
        assert StdI.trans_list[0][0] == pytest.approx(1.0)  # mu - 0.5*h
        assert StdI.trans_list[1][0] == pytest.approx(1.0)  # mu + 0.5*h

    def test_magnetic_field(self):
        """Test that longitudinal field splits spin channels."""
        StdI = _make_stdi()
        hubbard_local(StdI, mu0=0.0, h0=2.0, Gamma0=0.0, Gamma0_y=0.0,
                      U0=0.0, isite=0)
        assert len(StdI.trans_list) == 2
        assert StdI.trans_list[0][0] == pytest.approx(-1.0)  # 0 - 0.5*2
        assert StdI.trans_list[1][0] == pytest.approx(1.0)   # 0 + 0.5*2

    def test_intra_coulomb(self):
        """Test that intra-Coulomb U is added."""
        StdI = _make_stdi()
        hubbard_local(StdI, mu0=0.0, h0=0.0, Gamma0=0.0, Gamma0_y=0.0,
                      U0=4.0, isite=2)
        assert len(StdI.Cintra_list) == 1
        assert StdI.Cintra_list[0][0] == pytest.approx(4.0)
        assert StdI.Cintra_list[0][1] == 2

    def test_transverse_field(self):
        """Test that transverse field Gamma adds off-diagonal terms."""
        StdI = _make_stdi()
        hubbard_local(StdI, mu0=0.0, h0=0.0, Gamma0=1.0, Gamma0_y=0.0,
                      U0=0.0, isite=0)
        # Gamma: 2 terms (spin-flip)
        assert len(StdI.trans_list) == 2
        assert StdI.trans_list[0][0] == pytest.approx(-0.5)
        assert StdI.trans_list[1][0] == pytest.approx(-0.5)


# ===================================================================
#  mag_field
# ===================================================================


class TestMagField:
    """Tests for mag_field."""

    def test_spin_half_longitudinal(self):
        """Test longitudinal field for S=1/2."""
        StdI = _make_stdi()
        mag_field(StdI, S2=1, h=1.0, Gamma=0.0, Gamma_y=0.0, isite=0)
        # S=1/2: 2 states → 2 diagonal terms
        assert len(StdI.trans_list) == 2
        # ispin=0: Sz=+0.5, trans = -h*0.5 = -0.5
        assert StdI.trans_list[0][0] == pytest.approx(-0.5)
        # ispin=1: Sz=-0.5, trans = -h*(-0.5) = 0.5
        assert StdI.trans_list[1][0] == pytest.approx(0.5)

    def test_spin_half_transverse(self):
        """Test transverse field for S=1/2."""
        StdI = _make_stdi()
        mag_field(StdI, S2=1, h=0.0, Gamma=2.0, Gamma_y=0.0, isite=0)
        # 2 diagonal (h=0 → both zero, skipped) + 2 off-diagonal from Gamma
        assert len(StdI.trans_list) == 2

    def test_spin_one(self):
        """Test field for S=1 (3 states)."""
        StdI = _make_stdi()
        mag_field(StdI, S2=2, h=1.0, Gamma=0.0, Gamma_y=0.0, isite=0)
        # S=1: Sz=1,0,-1 → 3 diagonal terms, but Sz=0 gives 0 which is skipped
        assert len(StdI.trans_list) == 2


# ===================================================================
#  intr
# ===================================================================


class TestIntr:
    """Tests for intr."""

    def test_adds_interaction_term(self):
        """Test that a two-body term is added correctly."""
        StdI = _make_stdi()
        intr(StdI, 1.5 + 0j, 0, 0, 0, 0, 1, 0, 1, 0)
        assert len(StdI.intr_list) == 1
        assert StdI.intr_list[0][0] == pytest.approx(1.5 + 0j)
        assert list(StdI.intr_list[0][1:]) == [0, 0, 0, 0, 1, 0, 1, 0]

    def test_skips_small_value(self):
        """Test that values below threshold are skipped."""
        StdI = _make_stdi()
        intr(StdI, 1e-13, 0, 0, 0, 0, 1, 0, 1, 0)
        assert len(StdI.intr_list) == 0


# ===================================================================
#  general_j
# ===================================================================


class TestGeneralJ:
    """Tests for general_j."""

    def test_isotropic_spin_half(self):
        """Test isotropic J for S=1/2 sites."""
        StdI = _make_stdi()
        J = np.eye(3) * 1.0
        general_j(StdI, J, Si2=1, Sj2=1, isite=0, jsite=1)
        # For S=1/2, should generate Hund, Cinter, Ex, PairLift terms
        assert len(StdI.Hund_list) == 1
        assert len(StdI.Cinter_list) == 1
        assert len(StdI.Ex_list) == 1
        assert len(StdI.PairLift_list) == 1
        # Hund = -0.5 * Jzz = -0.5
        assert StdI.Hund_list[0][0] == pytest.approx(-0.5)

    def test_anisotropic_jzz_only(self):
        """Test Jzz-only interaction for S=1/2."""
        StdI = _make_stdi()
        J = np.zeros((3, 3))
        J[2, 2] = 2.0
        general_j(StdI, J, Si2=1, Sj2=1, isite=0, jsite=1)
        assert len(StdI.Hund_list) == 1
        assert StdI.Hund_list[0][0] == pytest.approx(-1.0)

    def test_spin_one_uses_general_path(self):
        """Test that S=1 (Si2=2) uses the general ZGeneral=1 path."""
        StdI = _make_stdi()
        J = np.eye(3)
        general_j(StdI, J, Si2=2, Sj2=2, isite=0, jsite=1)
        # S=1: ZGeneral stays 1, so all terms go through intr()
        assert len(StdI.Hund_list) == 0  # No Hund shortcut for S>1/2
        assert len(StdI.intr_list) > 0

    def test_kondo_exchange_sign(self):
        """Test that kondo model uses negative exchange sign."""
        StdI = _make_stdi(model="kondo")
        J = np.eye(3)
        general_j(StdI, J, Si2=1, Sj2=1, isite=0, jsite=1)
        # For kondo: Ex = -0.25 * (Jxx + Jyy) = -0.5
        assert StdI.Ex_list[0][0] == pytest.approx(-0.5)

    def test_non_kondo_exchange_sign(self):
        """Test that non-kondo model uses positive exchange sign."""
        StdI = _make_stdi(model="hubbard")
        J = np.eye(3)
        general_j(StdI, J, Si2=1, Sj2=1, isite=0, jsite=1)
        # For hubbard: Ex = +0.25 * (Jxx + Jyy) = +0.5
        assert StdI.Ex_list[0][0] == pytest.approx(0.5)


# ===================================================================
#  coulomb
# ===================================================================


class TestCoulomb:
    """Tests for coulomb."""

    def test_adds_coulomb_term(self):
        """Test that a Coulomb interaction is added."""
        StdI = _make_stdi()
        coulomb(StdI, 2.0, 0, 1)
        assert len(StdI.Cinter_list) == 1
        assert StdI.Cinter_list[0][0] == pytest.approx(2.0)
        assert StdI.Cinter_list[0][1] == 0
        assert StdI.Cinter_list[0][2] == 1

    def test_increments_counter(self):
        """Test that multiple calls increment the counter."""
        StdI = _make_stdi()
        coulomb(StdI, 1.0, 0, 1)
        coulomb(StdI, 2.0, 1, 2)
        assert len(StdI.Cinter_list) == 2


# ===================================================================
#  Backward compatibility
# ===================================================================


class TestBackwardCompatibility:
    """Test that functions are still importable from stdface_model_util."""

    def test_import_from_stdface_model_util(self):
        """Test that all 8 functions are re-exported."""
        from stdface.core.stdface_model_util import (
            trans as t,
            hopping as h,
            hubbard_local as hl,
            mag_field as mf,
            intr as i,
            general_j as gj,
            coulomb as c,
            malloc_interactions as mi,
        )
        from stdface.lattice.interaction_builder import (
            trans,
            hopping,
            hubbard_local,
            mag_field,
            intr,
            general_j,
            coulomb,
            malloc_interactions,
        )
        assert t is trans
        assert h is hopping
        assert hl is hubbard_local
        assert mf is mag_field
        assert i is intr
        assert gj is general_j
        assert c is coulomb
        assert mi is malloc_interactions


# ===================================================================
#  _spin_ladder_factor
# ===================================================================


class TestSpinLadderFactor:
    """Tests for _spin_ladder_factor helper."""

    def test_spin_half_sz_minus_half(self):
        """Test S=1/2, Sz=-1/2: sqrt(0.5*1.5 - (-0.5)*0.5) = 1."""
        assert _spin_ladder_factor(0.5, -0.5) == pytest.approx(1.0)

    def test_spin_half_sz_plus_half(self):
        """Test S=1/2, Sz=+1/2: sqrt(0.5*1.5 - 0.5*1.5) = 0."""
        assert _spin_ladder_factor(0.5, 0.5) == pytest.approx(0.0)

    def test_spin_one_sz_minus_one(self):
        """Test S=1, Sz=-1: sqrt(1*2 - (-1)*0) = sqrt(2)."""
        assert _spin_ladder_factor(1.0, -1.0) == pytest.approx(math.sqrt(2.0))

    def test_spin_one_sz_zero(self):
        """Test S=1, Sz=0: sqrt(1*2 - 0*1) = sqrt(2)."""
        assert _spin_ladder_factor(1.0, 0.0) == pytest.approx(math.sqrt(2.0))

    def test_spin_one_sz_plus_one(self):
        """Test S=1, Sz=+1: sqrt(1*2 - 1*2) = 0."""
        assert _spin_ladder_factor(1.0, 1.0) == pytest.approx(0.0)

    def test_spin_three_halves(self):
        """Test S=3/2, Sz=-1/2: sqrt(3/2*5/2 - (-1/2)*(1/2)) = sqrt(4) = 2."""
        result = _spin_ladder_factor(1.5, -0.5)
        # S(S+1) = 3.75, Sz(Sz+1) = -0.25 → sqrt(4.0) = 2.0
        assert result == pytest.approx(2.0)

    def test_matches_inline_computation(self):
        """Test that helper matches the original inline formula."""
        S, Sz = 2.0, -1.0
        expected = math.sqrt(S * (S + 1.0) - Sz * (Sz + 1.0))
        assert _spin_ladder_factor(S, Sz) == pytest.approx(expected)


# ===================================================================
#  _add_spin_half_terms
# ===================================================================


class TestAddSpinHalfTerms:
    """Tests for _add_spin_half_terms helper."""

    def test_adds_hund_term(self):
        """Test that Hund term is always added."""
        StdI = _make_stdi()
        J = np.eye(3) * 2.0
        _add_spin_half_terms(StdI, J, isite=0, jsite=1)
        assert len(StdI.Hund_list) == 1
        assert StdI.Hund_list[0][0] == pytest.approx(-1.0)  # -0.5 * 2.0

    def test_adds_cinter_term(self):
        """Test that Cinter term is always added."""
        StdI = _make_stdi()
        J = np.eye(3) * 4.0
        _add_spin_half_terms(StdI, J, isite=0, jsite=1)
        assert len(StdI.Cinter_list) == 1
        assert StdI.Cinter_list[0][0] == pytest.approx(-1.0)  # -0.25 * 4.0

    def test_hund_site_indices(self):
        """Test that Hund site indices are recorded."""
        StdI = _make_stdi()
        J = np.eye(3)
        _add_spin_half_terms(StdI, J, isite=3, jsite=7)
        assert StdI.Hund_list[0][1] == 3
        assert StdI.Hund_list[0][2] == 7

    def test_returns_false_false_for_diagonal_j(self):
        """Test return flags for diagonal J (off-diagonal negligible)."""
        StdI = _make_stdi()
        J = np.diag([1.0, 1.0, 1.0])
        use_z, use_ex = _add_spin_half_terms(StdI, J, isite=0, jsite=1)
        assert use_z is False
        assert use_ex is False

    def test_returns_false_true_for_offdiag_j(self):
        """Test return flags when J has off-diagonal xy components."""
        StdI = _make_stdi()
        J = np.eye(3)
        J[0, 1] = 1.0  # non-negligible off-diagonal
        use_z, use_ex = _add_spin_half_terms(StdI, J, isite=0, jsite=1)
        assert use_z is False
        assert use_ex is True

    def test_adds_ex_and_pairlift_for_diagonal(self):
        """Test that Ex and PairLift are added for diagonal J."""
        StdI = _make_stdi()
        J = np.diag([2.0, 2.0, 1.0])
        _add_spin_half_terms(StdI, J, isite=0, jsite=1)
        assert len(StdI.Ex_list) == 1
        assert len(StdI.PairLift_list) == 1
        # Ex = +0.25 * (2 + 2) = 1.0 (hubbard, non-kondo)
        assert StdI.Ex_list[0][0] == pytest.approx(1.0)
        # PairLift = 0.25 * (2 - 2) = 0.0
        assert StdI.PairLift_list[0][0] == pytest.approx(0.0)

    def test_kondo_exchange_sign(self):
        """Test that kondo model uses negative exchange sign."""
        StdI = _make_stdi(model="kondo")
        J = np.diag([2.0, 2.0, 1.0])
        _add_spin_half_terms(StdI, J, isite=0, jsite=1)
        # Kondo: Ex = -0.25 * (2 + 2) = -1.0
        assert StdI.Ex_list[0][0] == pytest.approx(-1.0)

    def test_mvmc_exchange_sign(self):
        """Test that mVMC solver uses negative exchange sign."""
        StdI = _make_stdi(solver="mVMC")
        J = np.diag([2.0, 2.0, 1.0])
        _add_spin_half_terms(StdI, J, isite=0, jsite=1)
        # mVMC: Ex = -0.25 * (2 + 2) = -1.0
        assert StdI.Ex_list[0][0] == pytest.approx(-1.0)

    def test_mvmc_requires_jxx_eq_jyy(self):
        """Test that mVMC needs Jxx == Jyy for shortcut path."""
        StdI = _make_stdi(solver="mVMC")
        J = np.diag([3.0, 1.0, 1.0])  # Jxx != Jyy
        use_z, use_ex = _add_spin_half_terms(StdI, J, isite=0, jsite=1)
        assert use_z is False
        assert use_ex is True  # falls back to general path
        assert len(StdI.Ex_list) == 0  # no shortcut Ex added

    def test_no_ex_for_offdiag(self):
        """Test that Ex/PairLift are NOT added for off-diagonal J."""
        StdI = _make_stdi()
        J = np.eye(3)
        J[0, 1] = 1.0
        _add_spin_half_terms(StdI, J, isite=0, jsite=1)
        assert len(StdI.Ex_list) == 0
        assert len(StdI.PairLift_list) == 0


# ===================================================================
#  add_neighbor_interaction
# ===================================================================


def _make_stdi_square_for_neighbor(
    W: int = 4, L: int = 4, model: str = "hubbard",
) -> StdIntList:
    """Create a StdIntList with a WxL square lattice and pre-allocated arrays.

    Manually sets up the super-cell (bypassing ``init_site``) to avoid
    the L/W/Height vs box conflict check.
    """
    import math

    StdI = StdIntList()
    StdI.solver = "HPhi"
    StdI.model = model
    StdI.method = "lanczos"
    StdI.PumpBody = 0
    StdI.pi = math.acos(-1.0)
    StdI.pi180 = StdI.pi / 180.0
    StdI.S2 = 1
    StdI.NsiteUC = 1

    # box = diag(W, L, 1)
    StdI.box[:, :] = 0
    StdI.box[0, 0] = W
    StdI.box[1, 1] = L
    StdI.box[2, 2] = 1

    StdI.NCell = W * L

    # direct lattice = identity
    StdI.direct[:, :] = 0.0
    StdI.direct[0, 0] = 1.0
    StdI.direct[1, 1] = 1.0
    StdI.direct[2, 2] = 1.0

    # Reciprocal box (cofactor of diag(W,L,1))
    StdI.rbox[:, :] = 0
    StdI.rbox[0, 0] = L       # cofactor of W => L*1
    StdI.rbox[1, 1] = W       # cofactor of L => W*1
    StdI.rbox[2, 2] = W * L   # cofactor of 1 => W*L

    # Cell array: all (iW, iL, 0) for 0<=iW<W, 0<=iL<L
    StdI.Cell = np.zeros((StdI.NCell, 3), dtype=int)
    idx = 0
    for iL_idx in range(L):
        for iW_idx in range(W):
            StdI.Cell[idx, 0] = iW_idx
            StdI.Cell[idx, 1] = iL_idx
            idx += 1

    # tau
    StdI.tau = np.zeros((1, 3))

    # Phase factors
    StdI.phase[:] = 0.0
    StdI.ExpPhase[:] = 1.0 + 0j
    StdI.AntiPeriod[:] = 0

    StdI.nsite = StdI.NCell
    if model == "kondo":
        StdI.nsite *= 2

    malloc_interactions(StdI, 500, 500)
    return StdI


class TestAddNeighborInteraction:
    """Tests for add_neighbor_interaction."""

    def test_returns_site_indices(self):
        """Returns (isite, jsite, Cphase, dR) from set_label."""
        StdI = _make_stdi_square_for_neighbor(4, 4, "hubbard")
        J = np.zeros((3, 3))
        # Neighbor along W: iW=0, iL=0, diW=1, diL=0 → sites 0 and 1
        isite, jsite, Cphase, dR = add_neighbor_interaction(
            StdI, None, 0, 0, 1, 0, 0, 0, 1, J, 1.0, 0.0)
        assert isite == 0
        assert jsite == 1

    def test_spin_model_calls_general_j(self):
        """For spin model, general_j adds exchange/Hund terms (S2=1)."""
        StdI = _make_stdi_square_for_neighbor(4, 4, "spin")
        J = np.eye(3)
        nex_before = len(StdI.Ex_list)
        nhund_before = len(StdI.Hund_list)
        add_neighbor_interaction(
            StdI, None, 0, 0, 1, 0, 0, 0, 1, J, 0.0, 0.0)
        # For S2=1 (spin-1/2), general_j uses shortcut path: Ex + Hund
        assert len(StdI.Ex_list) > nex_before or len(StdI.Hund_list) > nhund_before

    def test_hubbard_model_calls_hopping(self):
        """For hubbard model, hopping is called (ntrans increases)."""
        StdI = _make_stdi_square_for_neighbor(4, 4, "hubbard")
        J = np.zeros((3, 3))
        ntrans_before = len(StdI.trans_list)
        add_neighbor_interaction(
            StdI, None, 0, 0, 1, 0, 0, 0, 1, J, 1.0, 0.0)
        assert len(StdI.trans_list) > ntrans_before

    def test_hubbard_model_calls_coulomb(self):
        """For hubbard model with V != 0, coulomb adds interaction."""
        StdI = _make_stdi_square_for_neighbor(4, 4, "hubbard")
        J = np.zeros((3, 3))
        ncinter_before = len(StdI.Cinter_list)
        add_neighbor_interaction(
            StdI, None, 0, 0, 1, 0, 0, 0, 1, J, 0.0, 1.5)
        assert len(StdI.Cinter_list) > ncinter_before

    def test_spin_model_no_hopping(self):
        """For spin model, ntrans should not change."""
        StdI = _make_stdi_square_for_neighbor(4, 4, "spin")
        J = np.eye(3)
        ntrans_before = len(StdI.trans_list)
        add_neighbor_interaction(
            StdI, None, 0, 0, 1, 0, 0, 0, 1, J, 1.0, 1.0)
        assert len(StdI.trans_list) == ntrans_before

    def test_hubbard_model_no_general_j(self):
        """For hubbard model, general_j (spin-spin) should not be called."""
        StdI = _make_stdi_square_for_neighbor(4, 4, "hubbard")
        J = np.eye(3)  # non-zero J, but should be ignored
        nintr_before = len(StdI.intr_list)
        add_neighbor_interaction(
            StdI, None, 0, 0, 1, 0, 0, 0, 1, J, 0.0, 0.0)
        # For hubbard with t=0 and V=0, no interactions added
        assert len(StdI.intr_list) == nintr_before

    def test_fp_none_no_error(self):
        """fp=None does not cause an error."""
        StdI = _make_stdi_square_for_neighbor(4, 4, "hubbard")
        J = np.zeros((3, 3))
        # Should not raise
        add_neighbor_interaction(
            StdI, None, 0, 0, 1, 0, 0, 0, 1, J, 1.0, 0.0)

    def test_gnuplot_output_written(self):
        """When a buffer is provided, the bond is recorded."""
        from stdface.lattice.site_util import GnuplotBuffer
        StdI = _make_stdi_square_for_neighbor(4, 4, "hubbard")
        J = np.zeros((3, 3))
        buf = GnuplotBuffer()
        add_neighbor_interaction(
            StdI, buf, 0, 0, 1, 0, 0, 0, 1, J, 1.0, 0.0)
        content = buf.build(StdI).content
        assert "set label" in content


# ===================================================================
#  add_neighbor_interaction_3d
# ===================================================================


def _make_stdi_ortho_for_neighbor(
    W: int = 2, L: int = 2, H: int = 2, model: str = "hubbard",
) -> StdIntList:
    """Create a StdIntList with a WxLxH orthorhombic lattice.

    Manually sets up the 3D super-cell for testing
    ``add_neighbor_interaction_3d``.
    """
    import math

    StdI = StdIntList()
    StdI.solver = "HPhi"
    StdI.model = model
    StdI.method = "lanczos"
    StdI.PumpBody = 0
    StdI.pi = math.acos(-1.0)
    StdI.pi180 = StdI.pi / 180.0
    StdI.S2 = 1
    StdI.NsiteUC = 1

    # box = diag(W, L, H)
    StdI.box[:, :] = 0
    StdI.box[0, 0] = W
    StdI.box[1, 1] = L
    StdI.box[2, 2] = H

    StdI.NCell = W * L * H

    # direct lattice = identity
    StdI.direct[:, :] = 0.0
    StdI.direct[0, 0] = 1.0
    StdI.direct[1, 1] = 1.0
    StdI.direct[2, 2] = 1.0

    # Reciprocal box (cofactor of diag(W,L,H))
    StdI.rbox[:, :] = 0
    StdI.rbox[0, 0] = L * H
    StdI.rbox[1, 1] = W * H
    StdI.rbox[2, 2] = W * L

    # Cell array: all (iW, iL, iH)
    StdI.Cell = np.zeros((StdI.NCell, 3), dtype=int)
    idx = 0
    for iH_idx in range(H):
        for iL_idx in range(L):
            for iW_idx in range(W):
                StdI.Cell[idx, 0] = iW_idx
                StdI.Cell[idx, 1] = iL_idx
                StdI.Cell[idx, 2] = iH_idx
                idx += 1

    # tau
    StdI.tau = np.zeros((1, 3))

    # Phase factors
    StdI.phase[:] = 0.0
    StdI.ExpPhase[:] = 1.0 + 0j
    StdI.AntiPeriod[:] = 0

    StdI.nsite = StdI.NCell
    if model == "kondo":
        StdI.nsite *= 2

    malloc_interactions(StdI, 500, 500)
    return StdI


class TestAddNeighborInteraction3D:
    """Tests for add_neighbor_interaction_3d."""

    def test_returns_site_indices(self):
        """Returns (isite, jsite, Cphase, dR) from find_site."""
        StdI = _make_stdi_ortho_for_neighbor(2, 2, 2, "hubbard")
        J = np.zeros((3, 3))
        isite, jsite, Cphase, dR = add_neighbor_interaction_3d(
            StdI, 0, 0, 0, 1, 0, 0, 0, 0, J, 1.0, 0.0)
        assert isite == 0
        assert jsite == 1

    def test_spin_model_calls_general_j(self):
        """For spin model, general_j adds exchange terms."""
        StdI = _make_stdi_ortho_for_neighbor(2, 2, 2, "spin")
        J = np.eye(3)
        nex_before = len(StdI.Ex_list)
        nhund_before = len(StdI.Hund_list)
        add_neighbor_interaction_3d(
            StdI, 0, 0, 0, 1, 0, 0, 0, 0, J, 0.0, 0.0)
        assert len(StdI.Ex_list) > nex_before or len(StdI.Hund_list) > nhund_before

    def test_hubbard_model_calls_hopping(self):
        """For hubbard model, hopping terms are added."""
        StdI = _make_stdi_ortho_for_neighbor(2, 2, 2, "hubbard")
        J = np.zeros((3, 3))
        ntrans_before = len(StdI.trans_list)
        add_neighbor_interaction_3d(
            StdI, 0, 0, 0, 1, 0, 0, 0, 0, J, 1.0, 0.0)
        assert len(StdI.trans_list) > ntrans_before

    def test_hubbard_model_calls_coulomb(self):
        """For hubbard model with V != 0, coulomb adds interaction."""
        StdI = _make_stdi_ortho_for_neighbor(2, 2, 2, "hubbard")
        J = np.zeros((3, 3))
        ncinter_before = len(StdI.Cinter_list)
        add_neighbor_interaction_3d(
            StdI, 0, 0, 0, 1, 0, 0, 0, 0, J, 0.0, 1.5)
        assert len(StdI.Cinter_list) > ncinter_before

    def test_spin_model_no_hopping(self):
        """For spin model, ntrans should not change from neighbor interaction."""
        StdI = _make_stdi_ortho_for_neighbor(2, 2, 2, "spin")
        J = np.eye(3)
        ntrans_before = len(StdI.trans_list)
        add_neighbor_interaction_3d(
            StdI, 0, 0, 0, 1, 0, 0, 0, 0, J, 1.0, 1.0)
        assert len(StdI.trans_list) == ntrans_before

    def test_hubbard_model_no_general_j(self):
        """For hubbard model, general_j (spin-spin) should not be called."""
        StdI = _make_stdi_ortho_for_neighbor(2, 2, 2, "hubbard")
        J = np.eye(3)
        nintr_before = len(StdI.intr_list)
        add_neighbor_interaction_3d(
            StdI, 0, 0, 0, 1, 0, 0, 0, 0, J, 0.0, 0.0)
        assert len(StdI.intr_list) == nintr_before

    def test_neighbor_along_height(self):
        """3D neighbor along the H direction returns correct sites."""
        StdI = _make_stdi_ortho_for_neighbor(2, 2, 2, "hubbard")
        J = np.zeros((3, 3))
        # iW=0, iL=0, iH=0, diH=1 → should go to cell (0,0,1)
        isite, jsite, Cphase, dR = add_neighbor_interaction_3d(
            StdI, 0, 0, 0, 0, 0, 1, 0, 0, J, 1.0, 0.0)
        assert isite != jsite
        assert isite == 0


# ===================================================================
#  compute_max_interactions
# ===================================================================


def _make_stdi_for_max(
    model: str = "hubbard",
    NCell: int = 8,
    NsiteUC: int = 1,
    S2: int = 1,
) -> StdIntList:
    """Create a minimal StdIntList for testing compute_max_interactions."""
    StdI = StdIntList()
    StdI.model = model
    StdI.NCell = NCell
    StdI.NsiteUC = NsiteUC
    StdI.S2 = S2
    StdI.nsite = NCell * NsiteUC
    if model == "kondo":
        StdI.nsite *= 2
    return StdI


class TestComputeMaxInteractions:
    """Tests for compute_max_interactions."""

    def test_spin_ntrans(self):
        """Test ntransMax formula for spin model."""
        StdI = _make_stdi_for_max(model="spin", NCell=4, NsiteUC=1, S2=1)
        # nsite = 4, S2 = 1
        # ntransMax = nsite * (S2 + 1 + 2*S2) = 4 * (1 + 1 + 2) = 16
        ntrans, nintr = compute_max_interactions(StdI, n_bonds=6)
        assert ntrans == 16

    def test_spin_nintr(self):
        """Test nintrMax formula for spin model."""
        StdI = _make_stdi_for_max(model="spin", NCell=4, NsiteUC=1, S2=1)
        # nintrMax = NCell * (NsiteUC + n_bonds) * (3*S2+1)^2
        #          = 4 * (1 + 6) * 4^2 = 4 * 7 * 16 = 448
        ntrans, nintr = compute_max_interactions(StdI, n_bonds=6)
        assert nintr == 448

    def test_spin_s1(self):
        """Test spin model with S=1 (S2=2)."""
        StdI = _make_stdi_for_max(model="spin", NCell=4, NsiteUC=1, S2=2)
        # ntransMax = 4 * (2 + 1 + 4) = 28
        # nintrMax = 4 * (1 + 6) * (7)^2 = 4 * 7 * 49 = 1372
        ntrans, nintr = compute_max_interactions(StdI, n_bonds=6)
        assert ntrans == 28
        assert nintr == 1372

    def test_hubbard_ntrans(self):
        """Test ntransMax formula for hubbard model."""
        StdI = _make_stdi_for_max(model="hubbard", NCell=4, NsiteUC=1, S2=1)
        # ntransMax = NCell * 2 * (2*NsiteUC + 2*n_bonds)
        #           = 4 * 2 * (2 + 12) = 4 * 2 * 14 = 112
        ntrans, nintr = compute_max_interactions(StdI, n_bonds=6)
        assert ntrans == 112

    def test_hubbard_nintr(self):
        """Test nintrMax formula for hubbard model."""
        StdI = _make_stdi_for_max(model="hubbard", NCell=4, NsiteUC=1, S2=1)
        # nintrMax = NCell * (NsiteUC + 4*n_bonds) = 4 * (1 + 24) = 100
        ntrans, nintr = compute_max_interactions(StdI, n_bonds=6)
        assert nintr == 100

    def test_kondo_adds_spin_terms(self):
        """Test kondo model adds extra spin contribution."""
        StdI = _make_stdi_for_max(model="kondo", NCell=4, NsiteUC=1, S2=1)
        # nsite = 8 (doubled for kondo)
        # Base: ntransMax = 4 * 2 * (2 + 12) = 112
        # Kondo extra: nsite//2 * (1+1+2) = 4 * 4 = 16
        # Total ntransMax = 128
        ntrans, nintr = compute_max_interactions(StdI, n_bonds=6)
        assert ntrans == 128
        # Base nintrMax = 4 * (1 + 24) = 100
        # Kondo extra: 4 * (3*1+1)^2 = 4 * 16 = 64
        # Total nintrMax = 164
        assert nintr == 164

    def test_multi_site_uc(self):
        """Test with NsiteUC > 1 (e.g. pyrochlore with 4 sites per UC)."""
        StdI = _make_stdi_for_max(model="spin", NCell=2, NsiteUC=4, S2=1)
        # nsite = 8, S2 = 1
        # ntransMax = 8 * (1+1+2) = 32
        # nintrMax = 2 * (4 + 12) * 16 = 2 * 16 * 16 = 512
        ntrans, nintr = compute_max_interactions(StdI, n_bonds=12)
        assert ntrans == 32
        assert nintr == 512

    def test_returns_tuple(self):
        """Test that function returns a 2-tuple."""
        StdI = _make_stdi_for_max()
        result = compute_max_interactions(StdI, n_bonds=6)
        assert isinstance(result, tuple)
        assert len(result) == 2


# ===================================================================
#  add_local_terms
# ===================================================================


class TestAddLocalTerms:
    """Tests for add_local_terms."""

    def test_spin_model_adds_mag_field(self):
        """For spin model with h != 0, magnetic field terms are added."""
        StdI = _make_stdi(model="spin")
        StdI.S2 = 1
        StdI.h = 1.0
        StdI.Gamma = 0.0
        StdI.Gamma_y = 0.0
        StdI.D = np.zeros((3, 3))
        ntrans_before = len(StdI.trans_list)
        add_local_terms(StdI, 0, 0)
        # h != 0 means mag_field adds transfer terms
        assert len(StdI.trans_list) > ntrans_before

    def test_spin_model_adds_anisotropy(self):
        """For spin model with D != 0, anisotropy terms are added via general_j."""
        StdI = _make_stdi(model="spin")
        StdI.S2 = 3  # S=3/2 so general_j goes through full intr path
        StdI.h = 0.0
        StdI.Gamma = 0.0
        StdI.Gamma_y = 0.0
        StdI.D = np.zeros((3, 3))
        StdI.D[2, 2] = 1.0
        nintr_before = len(StdI.intr_list)
        add_local_terms(StdI, 0, 0)
        assert len(StdI.intr_list) > nintr_before

    def test_hubbard_model_adds_local_terms(self):
        """For hubbard model, hubbard_local adds transfer + Cintra terms."""
        StdI = _make_stdi(model="hubbard")
        StdI.mu = 1.0
        StdI.h = 0.0
        StdI.Gamma = 0.0
        StdI.Gamma_y = 0.0
        StdI.U = 2.0
        ntrans_before = len(StdI.trans_list)
        ncintra_before = len(StdI.Cintra_list)
        add_local_terms(StdI, 0, 0)
        assert len(StdI.trans_list) > ntrans_before
        assert len(StdI.Cintra_list) > ncintra_before

    def test_hubbard_model_no_kondo_coupling(self):
        """For hubbard model, no Kondo J-coupling is added."""
        StdI = _make_stdi(model="hubbard")
        StdI.mu = 0.0
        StdI.h = 0.0
        StdI.Gamma = 0.0
        StdI.Gamma_y = 0.0
        StdI.U = 0.0
        StdI.J = np.eye(3)  # non-zero J, but should be ignored
        nintr_before = len(StdI.intr_list)
        nhund_before = len(StdI.Hund_list)
        nex_before = len(StdI.Ex_list)
        add_local_terms(StdI, 0, 0)
        # No interactions from J coupling
        assert len(StdI.intr_list) == nintr_before
        assert len(StdI.Hund_list) == nhund_before
        assert len(StdI.Ex_list) == nex_before

    def test_kondo_model_adds_j_coupling(self):
        """For kondo model, J-coupling terms are added between isite and jsite_kondo."""
        StdI = _make_stdi(model="kondo")
        StdI.S2 = 1
        StdI.mu = 0.0
        StdI.h = 0.0
        StdI.Gamma = 0.0
        StdI.Gamma_y = 0.0
        StdI.U = 0.0
        StdI.J = np.eye(3)
        nhund_before = len(StdI.Hund_list)
        nex_before = len(StdI.Ex_list)
        add_local_terms(StdI, 10, 5)
        # S2=1 uses shortcut: Hund + Ex
        assert len(StdI.Hund_list) > nhund_before or len(StdI.Ex_list) > nex_before

    def test_kondo_model_adds_localized_mag_field(self):
        """For kondo model with h != 0, mag_field is added on jsite_kondo."""
        StdI = _make_stdi(model="kondo")
        StdI.S2 = 1
        StdI.mu = 0.0
        StdI.h = 1.0
        StdI.Gamma = 0.0
        StdI.Gamma_y = 0.0
        StdI.U = 0.0
        StdI.J = np.zeros((3, 3))
        ntrans_before = len(StdI.trans_list)
        add_local_terms(StdI, 10, 5)
        # hubbard_local adds h on isite=10, mag_field adds h on jsite_kondo=5
        assert len(StdI.trans_list) > ntrans_before

    def test_spin_model_no_hubbard_terms(self):
        """For spin model, no Cintra or Cinter terms are added."""
        StdI = _make_stdi(model="spin")
        StdI.S2 = 1
        StdI.h = 1.0
        StdI.Gamma = 0.0
        StdI.Gamma_y = 0.0
        StdI.D = np.zeros((3, 3))
        ncintra_before = len(StdI.Cintra_list)
        add_local_terms(StdI, 0, 0)
        assert len(StdI.Cintra_list) == ncintra_before

    def test_spin_model_transverse_field(self):
        """For spin model with Gamma != 0, transverse field terms are added."""
        StdI = _make_stdi(model="spin")
        StdI.S2 = 1
        StdI.h = 0.0
        StdI.Gamma = 1.0
        StdI.Gamma_y = 0.0
        StdI.D = np.zeros((3, 3))
        ntrans_before = len(StdI.trans_list)
        add_local_terms(StdI, 0, 0)
        assert len(StdI.trans_list) > ntrans_before
