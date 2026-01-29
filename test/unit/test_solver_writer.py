"""Unit tests for solver plugin system.

Tests for the SolverPlugin classes and the plugin registry.
"""
from __future__ import annotations

import os
import math
import tempfile

import pytest

from stdface.plugin import SolverPlugin, get_plugin
from stdface.solvers.hphi import HPhiPlugin
from stdface.solvers.mvmc import MVMCPlugin
from stdface.solvers.uhf import UHFPlugin
from stdface.solvers.hwave import HWavePlugin
from stdface.core.stdface_vals import StdIntList

# Sentinel values matching the C code
NaN_i = 2147483647
NaN_d = float("nan")
NaN_c = complex(float("nan"), 0.0)


def _make_stdi_for_hphi(nsite: int = 4) -> StdIntList:
    """Create a minimal StdIntList ready for HPhi writer."""
    StdI = StdIntList()
    StdI.solver = "HPhi"
    StdI.pi = math.acos(-1.0)
    StdI.nsite = nsite
    StdI.model = "hubbard"
    StdI.lattice = "chain"
    StdI.method = "lanczos"
    StdI.lGC = 0
    StdI.lBoost = 0
    StdI.Sz2 = 0
    StdI.ncond = 2

    # Green function / output settings
    StdI.ioutputmode = 1  # correlation
    StdI.lGC = 0

    # Namelist file head
    StdI.CDataFileHead = "zvo"

    # Fields needed by print_mod_para
    StdI.Lanczos_max = 1000
    StdI.initial_iv = 0
    StdI.nvec = 1
    StdI.exct = 1
    StdI.LanczosEps = 14
    StdI.LanczosTarget = 1
    StdI.NumAve = 1
    StdI.ExpecInterval = 1
    StdI.LargeValue = 10.0
    StdI.Nomega = 200
    StdI.OmegaMax = 10.0
    StdI.OmegaMin = -10.0
    StdI.OmegaOrg = 0.0
    StdI.OmegaIm = 0.01
    StdI.FlgTemp = 1

    # CalcSpec / SpectrumType
    StdI.CalcSpec = "none"
    StdI.SpectrumType = "****"
    StdI.SpectrumQ[:] = NaN_d
    StdI.Restart = "none"
    StdI.EigenVecIO = "none"
    StdI.InitialVecType = "c"
    StdI.HamIO = "none"
    StdI.OutputExVec = "none"
    StdI.PumpType = "****"
    StdI.NGPU = NaN_i
    StdI.Scalapack = NaN_i

    # Interactions (empty)
    StdI.NCintra = 0
    StdI.NCinter = 0
    StdI.NHund = 0
    StdI.NEx = 0
    StdI.NPairLift = 0
    StdI.NISpin = 0
    StdI.NCisAjtCkuAlv = 0
    StdI.NCisAjt = 0
    StdI.NTransfer = 0

    # Transfer array
    StdI.transindx = []
    StdI.trans = []

    # locspinflag
    StdI.locspinflag = [0] * nsite

    # outputmode
    StdI.outputmode = "****"

    return StdI


# =====================================================================
#  Tests for get_plugin registry
# =====================================================================


class TestGetPlugin:
    """Tests for the get_plugin registry function."""

    def test_returns_hphi_plugin(self):
        """Test that 'HPhi' returns an HPhiPlugin."""
        plugin = get_plugin("HPhi")
        assert isinstance(plugin, HPhiPlugin)
        assert plugin.name == "HPhi"

    def test_returns_mvmc_plugin(self):
        """Test that 'mVMC' returns an MVMCPlugin."""
        plugin = get_plugin("mVMC")
        assert isinstance(plugin, MVMCPlugin)
        assert plugin.name == "mVMC"

    def test_returns_uhf_plugin(self):
        """Test that 'UHF' returns a UHFPlugin."""
        plugin = get_plugin("UHF")
        assert isinstance(plugin, UHFPlugin)
        assert plugin.name == "UHF"

    def test_returns_hwave_plugin(self):
        """Test that 'HWAVE' returns an HWavePlugin."""
        plugin = get_plugin("HWAVE")
        assert isinstance(plugin, HWavePlugin)
        assert plugin.name == "HWAVE"

    def test_raises_on_unknown_solver(self):
        """Test that an unknown solver raises KeyError."""
        with pytest.raises(KeyError, match="No solver plugin"):
            get_plugin("unknown")


# =====================================================================
#  Tests for SolverPlugin subclasses
# =====================================================================


class TestSolverPluginIsAbstract:
    """Tests for the abstract base class."""

    def test_cannot_instantiate_directly(self):
        """Test that SolverPlugin cannot be instantiated directly."""
        with pytest.raises(TypeError):
            SolverPlugin()


class TestHPhiPlugin:
    """Tests for the HPhiPlugin class."""

    def test_writes_namelist_def(self):
        """Test that HPhiPlugin creates namelist.def."""
        StdI = _make_stdi_for_hphi(nsite=4)
        plugin = get_plugin("HPhi")

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                plugin.write(StdI)
                assert os.path.exists("namelist.def")
                assert os.path.exists("calcmod.def")
                assert os.path.exists("modpara.def")
                assert os.path.exists("locspn.def")
            finally:
                os.chdir(orig)

    def test_writes_trans_def(self):
        """Test that HPhiPlugin creates trans.def."""
        StdI = _make_stdi_for_hphi(nsite=4)
        plugin = get_plugin("HPhi")

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                plugin.write(StdI)
                assert os.path.exists("trans.def")
            finally:
                os.chdir(orig)

    def test_writes_green_files(self):
        """Test that HPhiPlugin creates greenone.def and greentwo.def."""
        StdI = _make_stdi_for_hphi(nsite=4)
        plugin = get_plugin("HPhi")

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                plugin.write(StdI)
                assert os.path.exists("greenone.def")
                assert os.path.exists("greentwo.def")
            finally:
                os.chdir(orig)


class TestUHFPlugin:
    """Tests for the UHFPlugin class."""

    def test_writes_expected_files(self):
        """Test that UHFPlugin creates the expected set of files."""
        StdI = _make_stdi_for_hphi(nsite=4)
        StdI.solver = "UHF"
        StdI.outputmode = "****"
        plugin = get_plugin("UHF")

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                plugin.write(StdI)
                assert os.path.exists("locspn.def")
                assert os.path.exists("trans.def")
                assert os.path.exists("modpara.def")
                assert os.path.exists("namelist.def")
                assert os.path.exists("greenone.def")
            finally:
                os.chdir(orig)


class TestHWavePlugin:
    """Tests for the HWavePlugin class."""

    def test_uhfr_mode_writes_trans(self):
        """Test that HWAVE in uhfr mode writes trans.def."""
        StdI = _make_stdi_for_hphi(nsite=4)
        StdI.solver = "HWAVE"
        StdI.calcmode = "uhfr"
        StdI.outputmode = "****"
        plugin = get_plugin("HWAVE")

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                plugin.write(StdI)
                assert os.path.exists("trans.def")
                assert os.path.exists("greenone.def")
            finally:
                os.chdir(orig)
