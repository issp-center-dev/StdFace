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
from stdface.core.stdface_main import _attach_solver_config
from stdface.lattice import chain_lattice as cl

# Sentinel values matching the C code
NaN_i = 2147483647
NaN_d = float("nan")
NaN_c = complex(float("nan"), 0.0)


def _make_stdi_for_hphi(nsite: int = 4) -> StdIntList:
    """Create a minimal StdIntList ready for HPhi writer."""
    from stdface.core.stdface_main import _attach_solver_config
    StdI = StdIntList()
    StdI.solver = "HPhi"
    _attach_solver_config(StdI)  # C3: solver-specific fields live on the config
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
    StdI.SpectrumType = None
    StdI.SpectrumQ[:] = NaN_d
    StdI.Restart = "none"
    StdI.EigenVecIO = "none"
    StdI.InitialVecType = "c"
    StdI.HamIO = "none"
    StdI.OutputExVec = "none"
    StdI.PumpType = None
    StdI.NGPU = None
    StdI.Scalapack = None

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

    # Transfer terms
    StdI.trans_list = []

    # locspinflag
    StdI.locspinflag = [0] * nsite

    # outputmode
    StdI.outputmode = None

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
    s.W = None
    s.Height = None
    s.box[:, :] = NaN_i
    s.S2 = None
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
    # Attach MVMCConfig (C3) so mVMC-only fields resolve off the config.
    from stdface.core.stdface_main import _attach_solver_config
    _attach_solver_config(StdI)
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
    StdI.Hsub = None
    StdI.Lsub = None
    StdI.Wsub = 2
    StdI.NMPTrans = None
    if lgc_orb_para:
        StdI.lGC = 1
        StdI.Sz2 = None
    else:
        StdI.lGC = 0
        StdI.Sz2 = 0
    # ``_check_mod_para_mvmc``: when ``NVMCCalMode == 0``, ``NDataQtySmp`` must stay unset
    StdI.NDataQtySmp = None
    StdI.NSPGaussLeg = None
    StdI.NSPStot = None
    return StdI


def _make_hwave_wannier_export_stdi(nsiteUC: int = 2, ncell: int = 2) -> StdIntList:
    """StdIntList for HWave Wannier export (geometry + empty interactions)."""
    from stdface.core.stdface_main import _attach_solver_config
    s = StdIntList()
    s.solver = "HWAVE"
    _attach_solver_config(s)  # C3: H-wave fields (fileprefix/...) live on the config
    s.NsiteUC = nsiteUC
    s.NCell = ncell
    s.fileprefix = ""
    s.export_all = None
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
    s.trans_list = []
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

    def test_uhfr_uhfk_resolve_to_hwave_plugin(self):
        """C4-2b: UHFR/UHFK/HWAVE all resolve to the single HWavePlugin."""
        from stdface.solvers.hwave._plugin import HWavePlugin
        assert isinstance(get_plugin("UHFR"), HWavePlugin)
        assert isinstance(get_plugin("UHFK"), HWavePlugin)
        assert isinstance(get_plugin("HWAVE"), HWavePlugin)

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
        _attach_solver_config(StdI)
        StdI.outputmode = None
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

    def test_build_output_returns_container(self):
        """UHF build_output returns a fully data-backed ExpertModeOutput."""
        from stdface.core.output import ExpertModeOutput
        StdI = _make_stdi_for_hphi(nsite=4)
        StdI.solver = "UHF"
        _attach_solver_config(StdI)
        StdI.outputmode = None
        out = get_plugin("UHF").build_output(StdI)
        assert isinstance(out, ExpertModeOutput)
        d = out.to_dict()
        # UHF: no two-body green, but everything else present and serialisable
        assert d["green_two"] is None
        assert d["locspn"] and d["trans"] and d["modpara"] and d["namelist"]
        assert "params" in d["modpara"]
        assert d["namelist"]["entries"][0] == ["ModPara", "modpara.def"]

    def test_hphi_build_output_emits_solver_files(self, tmp_path, monkeypatch):
        """HPhi build_output writes solver-specific files via the hook."""
        from stdface.core.output import ExpertModeOutput
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_for_hphi(nsite=4)
        StdI.solver = "HPhi"
        StdI.outputmode = None
        out = get_plugin("HPhi").build_output(StdI)
        assert isinstance(out, ExpertModeOutput)
        # excitation / calcmod written eagerly by write_solver_files
        assert os.path.exists("calcmod.def")
        assert os.path.exists("pair.def") or os.path.exists("single.def")


class TestHWaveFamilyPlugin:
    """C4-2a: HWavePlugin owns H-wave input; output mode resolves by calcmode."""

    def test_hwave_registered_and_owns_input(self):
        from stdface.solvers.hwave._plugin import HWavePlugin
        plugin = get_plugin("HWAVE")
        assert isinstance(plugin, HWavePlugin)
        kw = plugin.keyword_table
        # union: UHF base + calcmode + UHFK Wannier-export keys
        assert "calcmode" in kw and "fileprefix" in kw and "mix" in kw
        reset_names = {n for n, _ in plugin.reset_scalars}
        assert {"export_all", "lattice_gp", "mix", "NMPTrans"} <= reset_names

    def test_output_mode_resolution(self):
        from stdface.solvers.hwave._plugin import _hwave_output_mode
        from stdface.core.stdface_vals import SolverType

        def mode(solver, calcmode):
            s = StdIntList()
            s.solver = solver
            s.calcmode = calcmode
            return _hwave_output_mode(s)

        assert mode("HWAVE", "uhfr") == SolverType.UHFR
        assert mode("HWAVE", "uhfk") == SolverType.UHFK
        assert mode("HWAVE", "rpa") == SolverType.UHFK
        assert mode("HWAVE", None) == SolverType.UHFK
        assert mode("UHFR", None) == SolverType.UHFR
        assert mode("UHFK", None) == SolverType.UHFK


class TestHWaveSplit:
    """Tests for the UHFR / UHFK plugins (the H-wave split)."""

    def test_uhfr_writes_trans(self):
        """UHFR writes trans.def / greenone.def."""
        StdI = _make_stdi_for_hphi(nsite=4)
        StdI.solver = "UHFR"
        _attach_solver_config(StdI)
        StdI.calcmode = "uhfr"
        StdI.outputmode = None
        plugin = get_plugin("UHFR")

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                plugin.write(StdI)
                assert os.path.exists("trans.def")
                assert os.path.exists("greenone.def")
            finally:
                os.chdir(orig)

    def test_uhfk_export_writes_geom(self, tmp_path, monkeypatch):
        """UHFK runs ``export_geometry`` / ``export_interaction``."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_hwave_wannier_export_stdi(nsiteUC=2, ncell=2)
        StdI.solver = "UHFK"
        plugin = get_plugin("UHFK")
        plugin.write(StdI)
        assert os.path.exists("geom.dat")
        # ``ntrans == 0`` → transfer export is skipped (no ``transfer.dat``)
