"""Unit tests for chain_lattice module.

Tests for the Python translation of ChainLattice.c.
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.lattice import chain_lattice as cl


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

NaN_d = float("nan")
NaN_i = 2147483647
NaN_c = complex(float("nan"), 0.0)


def _make_spin_chain(L: int = 4) -> StdIntList:
    """Return an StdIntList pre-configured for a spin chain."""
    s = StdIntList()
    # Set sentinel / NaN values for unused parameters
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "spin"
    s.solver = "HPhi"
    s.lattice = "chain"

    # Lattice parameters -- set NaN for those that should get defaults
    s.a = NaN_d
    s.length[0] = NaN_d
    s.length[1] = NaN_d
    s.direct[0, 0] = NaN_d
    s.direct[0, 1] = NaN_d
    s.direct[1, 0] = NaN_d
    s.direct[1, 1] = NaN_d
    s.phase[0] = NaN_d
    s.phase[1] = NaN_d

    s.L = L
    s.W = NaN_i
    s.Height = NaN_i
    # Set box to NaN_i so init_site uses L/W/Height path
    s.box[:, :] = NaN_i

    # Spin parameters
    s.S2 = NaN_i
    s.h = NaN_d
    s.Gamma = NaN_d
    s.Gamma_y = NaN_d
    s.D[2, 2] = NaN_d

    # J couplings -- set to NaN
    s.JAll = NaN_d
    s.JpAll = NaN_d
    s.JppAll = NaN_d
    s.J0All = NaN_d
    s.J0pAll = NaN_d
    s.J0ppAll = NaN_d
    s.J1All = NaN_d
    s.J1pAll = NaN_d
    s.J2All = NaN_d
    s.J2pAll = NaN_d
    s.J[:, :] = NaN_d
    s.Jp[:, :] = NaN_d
    s.Jpp[:, :] = NaN_d
    s.J0[:, :] = NaN_d
    s.J0p[:, :] = NaN_d
    s.J0pp[:, :] = NaN_d
    s.J1[:, :] = NaN_d
    s.J1p[:, :] = NaN_d
    s.J2[:, :] = NaN_d
    s.J2p[:, :] = NaN_d

    # Unused parameters -- set to NaN so not_used checks pass
    s.mu = NaN_d
    s.U = NaN_d
    s.t = NaN_c
    s.t0 = NaN_c
    s.tp = NaN_c
    s.t1 = NaN_c
    s.t2 = NaN_c
    s.t1p = NaN_d
    s.t2p = NaN_d
    s.V = NaN_d
    s.V0 = NaN_d
    s.Vp = NaN_d
    s.V1 = NaN_d
    s.V2 = NaN_d
    s.V1p = NaN_d
    s.V2p = NaN_d
    s.K = NaN_d

    # Set JAll to 1.0 so the spin model gets a non-trivial coupling
    s.JAll = 1.0

    return s


# ---------------------------------------------------------------------------
#  Tests: chain function exists and is callable
# ---------------------------------------------------------------------------


class TestChainCallable:
    """Verify that the chain function is importable and callable."""

    def test_chain_exists(self):
        """chain should be a callable attribute of chain_lattice."""
        assert hasattr(cl, "chain")
        assert callable(cl.chain)

    def test_chain_boost_exists(self):
        """chain_boost should be a callable attribute of chain_lattice."""
        assert hasattr(cl, "chain_boost")
        assert callable(cl.chain_boost)


# ---------------------------------------------------------------------------
#  Tests: spin model
# ---------------------------------------------------------------------------


class TestSpinModel:
    """Test that chain() correctly sets up a spin model."""

    def test_spin_chain_basic(self, tmp_path):
        """chain() with a spin model should populate locspinflag and arrays."""
        os.chdir(tmp_path)
        s = _make_spin_chain(L=4)
        cl.chain(s)

        # NsiteUC should be 1
        assert s.NsiteUC == 1
        # nsite should equal L for spin model
        assert s.nsite == 4
        # locspinflag should be set to S2=1 for all sites
        assert s.locspinflag is not None
        assert len(s.locspinflag) == 4
        assert all(s.locspinflag[i] == 1 for i in range(4))

    def test_spin_chain_interactions_allocated(self, tmp_path):
        """After chain(), interaction arrays should be allocated."""
        os.chdir(tmp_path)
        s = _make_spin_chain(L=4)
        cl.chain(s)

        assert s.transindx is not None
        assert s.intrindx is not None
        assert s.trans is not None
        assert s.intr is not None

    def test_spin_chain_tau(self, tmp_path):
        """tau should be set to zeros for a chain lattice."""
        os.chdir(tmp_path)
        s = _make_spin_chain(L=4)
        cl.chain(s)

        assert s.tau is not None
        np.testing.assert_allclose(s.tau[0], [0.0, 0.0, 0.0])

    def test_spin_chain_defaults(self, tmp_path):
        """Default values should be applied for h, Gamma, etc."""
        os.chdir(tmp_path)
        s = _make_spin_chain(L=4)
        cl.chain(s)

        assert s.h == 0.0
        assert s.Gamma == 0.0
        assert s.Gamma_y == 0.0
        assert s.S2 == 1


# ---------------------------------------------------------------------------
#  Tests: lattice.gp output
# ---------------------------------------------------------------------------


class TestLatticeGP:
    """Test that lattice.gp file is created."""

    def test_lattice_gp_created(self, tmp_path):
        """chain() should produce a lattice.gp file."""
        os.chdir(tmp_path)
        s = _make_spin_chain(L=4)
        cl.chain(s)

        gp_file = tmp_path / "lattice.gp"
        assert gp_file.exists()

    def test_lattice_gp_content(self, tmp_path):
        """lattice.gp should end with the expected gnuplot command."""
        os.chdir(tmp_path)
        s = _make_spin_chain(L=4)
        cl.chain(s)

        gp_file = tmp_path / "lattice.gp"
        content = gp_file.read_text()
        assert "plot '-' w d lc 7" in content
        assert "pause -1" in content

    def test_lattice_gp_not_created_hwave(self, tmp_path):
        """When solver is HWAVE and lattice_gp is 0, no lattice.gp file."""
        os.chdir(tmp_path)
        s = _make_spin_chain(L=4)
        s.solver = "HWAVE"
        s.lattice_gp = 0
        cl.chain(s)

        gp_file = tmp_path / "lattice.gp"
        assert not gp_file.exists()

    def test_lattice_gp_created_hwave_flag(self, tmp_path):
        """When solver is HWAVE but lattice_gp is 1, lattice.gp is created."""
        os.chdir(tmp_path)
        s = _make_spin_chain(L=4)
        s.solver = "HWAVE"
        s.lattice_gp = 1
        cl.chain(s)

        gp_file = tmp_path / "lattice.gp"
        assert gp_file.exists()


# ---------------------------------------------------------------------------
#  Tests: chain_boost
# ---------------------------------------------------------------------------


class TestChainBoost:
    """Test chain_boost function."""

    def test_boost_skipped_non_hphi(self, tmp_path):
        """chain_boost should return early if solver is not HPhi."""
        os.chdir(tmp_path)
        s = StdIntList()
        s.solver = "mVMC"
        cl.chain_boost(s)
        # Should not create boost.def
        assert not (tmp_path / "boost.def").exists()

    def test_boost_creates_file(self, tmp_path):
        """chain_boost should create boost.def for HPhi solver."""
        os.chdir(tmp_path)
        s = _make_spin_chain(L=16)
        # Run chain first to set up J0, Jp, etc.
        cl.chain(s)
        # Now run boost
        s.L = 16
        s.S2 = 1
        cl.chain_boost(s)

        boost_file = tmp_path / "boost.def"
        assert boost_file.exists()

    def test_boost_s2_error(self, tmp_path):
        """chain_boost should exit if S2 != 1."""
        os.chdir(tmp_path)
        s = StdIntList()
        s.solver = "HPhi"
        s.S2 = 2
        s.L = 16
        s.Gamma = 0.0
        s.Gamma_y = 0.0
        s.h = 0.0
        s.J0 = np.zeros((3, 3))
        s.Jp = np.zeros((3, 3))
        with pytest.raises(SystemExit):
            cl.chain_boost(s)

    def test_boost_L_mod8_error(self, tmp_path):
        """chain_boost should exit if L %% 8 != 0."""
        os.chdir(tmp_path)
        s = StdIntList()
        s.solver = "HPhi"
        s.S2 = 1
        s.L = 10  # not divisible by 8
        s.Gamma = 0.0
        s.Gamma_y = 0.0
        s.h = 0.0
        s.J0 = np.zeros((3, 3))
        s.Jp = np.zeros((3, 3))
        with pytest.raises(SystemExit):
            cl.chain_boost(s)
