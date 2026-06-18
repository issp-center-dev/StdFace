"""Unit tests for solver plugin system.

Tests for the SolverPlugin classes and the plugin registry.
"""
from __future__ import annotations

import os
import math
import tempfile

import numpy as np
import pytest

from stdface.plugin import SolverPlugin, get_plugin
from stdface.solvers.hphi import HPhiPlugin
from stdface.solvers.mvmc import MVMCPlugin
from stdface.solvers.uhf import UHFPlugin
from stdface.solvers.hwave import HWavePlugin
from stdface.core.stdface_vals import StdIntList
from stdface.lattice import chain_lattice as cl

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


def _make_spin_chain(L: int = 16) -> StdIntList:
    """Return a spin-chain StdIntList (aligned with ``test_chain_lattice``)."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "spin"
    s.solver = "HPhi"
    s.lattice = "chain"
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
    s.box[:, :] = NaN_i
    s.S2 = NaN_i
    s.h = NaN_d
    s.Gamma = NaN_d
    s.Gamma_y = NaN_d
    s.D[2, 2] = NaN_d
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
    s.JAll = 1.0
    return s


def _make_stdi_for_mvmc_write(
    nsite: int = 4,
    *,
    lgc_orb_para: bool = False,
) -> StdIntList:
    """Build StdIntList for :meth:`MVMCPlugin.write` (1D Hubbard chain + subcell)."""
    StdI = _make_stdi_for_hphi(nsite=nsite)
    StdI.solver = "mVMC"
    StdI.lattice = "chain"
    StdI.model = "hubbard"
    L = nsite
    StdI.box[:, :] = 0
    StdI.box[0, 0] = L
    StdI.box[1, 1] = 1
    StdI.box[2, 2] = 1
    StdI.NCell = L
    StdI.NsiteUC = 1
    StdI.nsite = nsite
    StdI.rbox = np.zeros((3, 3), dtype=float)
    StdI.rbox[0, 0] = 1
    StdI.rbox[1, 1] = L
    StdI.rbox[2, 2] = L
    StdI.Cell = np.zeros((L, 3), dtype=float)
    for i in range(L):
        StdI.Cell[i, 0] = float(i)
    StdI.tau = np.zeros((1, 3))
    StdI.phase[:] = 0.0
    StdI.pi180 = StdI.pi / 180.0
    StdI.ExpPhase = np.ones(3, dtype=complex)
    StdI.AntiPeriod = np.zeros(3, dtype=int)
    StdI.locspinflag = [0] * nsite
    StdI.boxsub[:, :] = NaN_i
    StdI.Hsub = NaN_i
    StdI.Lsub = NaN_i
    StdI.Wsub = 2
    StdI.NMPTrans = NaN_i
    if lgc_orb_para:
        StdI.lGC = 1
        StdI.Sz2 = NaN_i
    else:
        StdI.lGC = 0
        StdI.Sz2 = 0
    # ``_check_mod_para_mvmc``: when ``NVMCCalMode == 0``, ``NDataQtySmp`` must stay unset
    StdI.NDataQtySmp = NaN_i
    StdI.NSPGaussLeg = NaN_i
    StdI.NSPStot = NaN_i
    return StdI


def _make_hwave_wannier_export_stdi(nsiteUC: int = 2, ncell: int = 2) -> StdIntList:
    """StdIntList for HWave Wannier export (geometry + empty interactions)."""
    s = StdIntList()
    s.solver = "HWAVE"
    s.NsiteUC = nsiteUC
    s.NCell = ncell
    s.fileprefix = ""
    s.export_all = NaN_i
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
    s.direct = np.eye(3)
    s.tau = np.zeros((nsiteUC, 3))
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

    def test_post_lattice_chain_boost(self, tmp_path, monkeypatch):
        """``post_lattice`` runs lattice Boost when ``lBoost == 1``."""
        monkeypatch.chdir(tmp_path)
        s = _make_spin_chain(L=16)
        cl.chain(s)
        s.lBoost = 1
        plugin = get_plugin("HPhi")
        plugin.post_lattice(s)
        assert (tmp_path / "boost.def").exists()

    def test_time_evolution_pump_writes_teone(self, tmp_path, monkeypatch):
        """``print_pump`` writes ``teone.def`` for time evolution + pump."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_for_hphi(nsite=4)
        StdI.method = "timeevolution"
        StdI.PumpBody = 1
        StdI.Lanczos_max = 2
        StdI.dt = 0.1
        StdI.npump = [1, 1]
        StdI.pumpindx = [
            [[0, 0, 1, 0]],
            [[0, 0, 1, 0]],
        ]
        StdI.pump = [
            [complex(1.0, 0.0)],
            [complex(0.5, 0.0)],
        ]
        from stdface.solvers.hphi.writer import print_pump
        print_pump(StdI)
        assert os.path.exists("teone.def")


class TestMVMCPlugin:
    """Tests for :class:`MVMCPlugin` full :meth:`~MVMCPlugin.write` path."""

    def test_write_creates_variational_files(self, tmp_path, monkeypatch):
        """mVMC ``write`` produces variational files (orbital / Jastrow / Gutzwiller)."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_for_mvmc_write(nsite=4, lgc_orb_para=False)
        plugin = get_plugin("mVMC")
        plugin.write(StdI)
        assert os.path.exists("orbitalidx.def")
        assert os.path.exists("jastrowidx.def")
        assert os.path.exists("gutzwilleridx.def")
        assert os.path.exists("qptransidx.def")

    def test_write_lgc_prints_orb_para(self, tmp_path, monkeypatch):
        """``lGC == 1`` triggers ``print_orb_para`` (orbitalidxpara / gen)."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_for_mvmc_write(nsite=4, lgc_orb_para=True)
        plugin = get_plugin("mVMC")
        plugin.write(StdI)
        assert os.path.exists("orbitalidxpara.def")
        assert os.path.exists("orbitalidxgen.def")


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

    def test_wannier90_export_writes_geom(self, tmp_path, monkeypatch):
        """Non-``uhfr`` ``calcmode`` runs ``export_geometry`` / ``export_interaction``."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_hwave_wannier_export_stdi(nsiteUC=2, ncell=2)
        StdI.calcmode = "wannier90"
        plugin = get_plugin("HWAVE")
        plugin.write(StdI)
        assert os.path.exists("geom.dat")
        # ``ntrans == 0`` → transfer export is skipped (no ``transfer.dat``)
