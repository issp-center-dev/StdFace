"""Unit tests for common_writer module.

Tests for the common output writer functions shared across all solvers,
extracted from ``stdface_main.py``.
"""
from __future__ import annotations

import os
import tempfile

import pytest

from stdface.core.stdface_vals import StdIntList, ModelType, SolverType
from stdface.writer.common_writer import (
    print_loc_spin,
    print_trans,
    print_namelist,
    print_mod_para,
    print_1_green,
    print_2_green,
    unsupported_system,
    check_output_mode,
    check_mod_para,
    OUTPUT_MODE_TO_INT,
    MODEL_GC_TO_EX_UPDATE_PATH,
    _check_mod_para_hphi,
    _check_mod_para_mvmc,
    _check_mod_para_uhf,
    _check_conserved_quantities,
    _CONSERVED_QTY_RULES,
    _write_modpara_hphi,
    _write_modpara_mvmc,
    _write_modpara_uhf_hwave,
    _MODPARA_BANNER,
    _write_namelist_hphi,
    _write_namelist_mvmc,
    _INTERACTION_FLAGS,
    GreenFunctionIndices,
    _merge_duplicate_terms,
)
from stdface.writer.interaction_writer import print_interactions

# Sentinel values matching what _reset_vals sets at runtime
NaN_i = 2147483647


def _make_stdi_base(**overrides) -> StdIntList:
    """Create a minimal StdIntList with common fields."""
    StdI = StdIntList()
    StdI.nsite = overrides.get("nsite", 4)
    StdI.solver = overrides.get("solver", "HPhi")
    StdI.model = overrides.get("model", "hubbard")
    StdI.lGC = overrides.get("lGC", 0)
    StdI.lBoost = overrides.get("lBoost", 0)
    StdI.locspinflag = overrides.get("locspinflag", [0] * StdI.nsite)
    StdI.ioutputmode = overrides.get("ioutputmode", 1)
    StdI.outputmode = overrides.get("outputmode", "****")
    StdI.NsiteUC = overrides.get("NsiteUC", 1)
    StdI.Sz2 = overrides.get("Sz2", NaN_i)
    StdI.ncond = overrides.get("ncond", NaN_i)
    return StdI


class TestPrintLocSpin:
    """Tests for the print_loc_spin function."""

    def test_writes_locspn_def(self):
        """Test that locspn.def is created with correct content."""
        StdI = _make_stdi_base(nsite=4)
        StdI.locspinflag = [0, 0, 1, 1]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_loc_spin(StdI)
                assert os.path.exists("locspn.def")
                content = open("locspn.def").read()
                assert "NlocalSpin     2" in content
            finally:
                os.chdir(orig)

    def test_all_itinerant(self):
        """Test with all itinerant electrons (nlocspin=0)."""
        StdI = _make_stdi_base(nsite=2)
        StdI.locspinflag = [0, 0]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_loc_spin(StdI)
                content = open("locspn.def").read()
                assert "NlocalSpin     0" in content
            finally:
                os.chdir(orig)


class TestPrintTrans:
    """Tests for the print_trans function."""

    def test_writes_trans_def(self):
        """Test that trans.def is created."""
        StdI = _make_stdi_base(nsite=2)
        StdI.trans_list = [(1.0 + 0j, 0, 0, 1, 0), (1.0 + 0j, 1, 0, 0, 0)]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_trans(StdI)
                assert os.path.exists("trans.def")
                content = open("trans.def").read()
                assert "NTransfer" in content
            finally:
                os.chdir(orig)

    def test_merges_duplicates(self):
        """Test that duplicate transfer terms are merged."""
        StdI = _make_stdi_base(nsite=2)
        StdI.trans_list = [
            (1.0 + 0j, 0, 0, 1, 0),
            (2.0 + 0j, 0, 0, 1, 0),
            (0.5 + 0j, 1, 0, 0, 0),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_trans(StdI)
                content = open("trans.def").read()
                # First two are merged: 1+2=3, third is 0.5
                # So NTransfer should be 2
                assert "NTransfer       2" in content
            finally:
                os.chdir(orig)

    def test_suppresses_small_values(self):
        """Test that entries below threshold are suppressed."""
        StdI = _make_stdi_base(nsite=2)
        StdI.trans_list = [(1e-8 + 0j, 0, 0, 1, 0), (1.0 + 0j, 1, 0, 0, 0)]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_trans(StdI)
                content = open("trans.def").read()
                assert "NTransfer       1" in content
            finally:
                os.chdir(orig)


class TestCheckOutputMode:
    """Tests for the check_output_mode function."""

    def test_none_mode(self):
        """Test 'none' sets ioutputmode=0."""
        StdI = _make_stdi_base()
        StdI.outputmode = "none"
        check_output_mode(StdI)
        assert StdI.ioutputmode == 0

    def test_correlation_mode(self):
        """Test 'correlation' sets ioutputmode=1."""
        StdI = _make_stdi_base()
        StdI.outputmode = "correlation"
        check_output_mode(StdI)
        assert StdI.ioutputmode == 1

    def test_default_sentinel(self):
        """Test '****' (default) sets ioutputmode=1."""
        StdI = _make_stdi_base()
        StdI.outputmode = "****"
        check_output_mode(StdI)
        assert StdI.ioutputmode == 1

    def test_raw_mode(self):
        """Test 'raw' sets ioutputmode=2."""
        StdI = _make_stdi_base()
        StdI.outputmode = "raw"
        check_output_mode(StdI)
        assert StdI.ioutputmode == 2

    def test_full_mode(self):
        """Test 'full' sets ioutputmode=2."""
        StdI = _make_stdi_base()
        StdI.outputmode = "full"
        check_output_mode(StdI)
        assert StdI.ioutputmode == 2

    def test_invalid_mode_exits(self):
        """Test that invalid mode raises ValueError."""
        StdI = _make_stdi_base()
        StdI.outputmode = "bogus"
        with pytest.raises(ValueError):
            check_output_mode(StdI)


class TestUnsupportedSystem:
    """Tests for the unsupported_system function."""

    def test_exits(self):
        """Test that unsupported_system always raises ValueError."""
        with pytest.raises(ValueError):
            unsupported_system("hubbard", "dodecahedron")


class TestPrintNamelist:
    """Tests for the print_namelist function."""

    def test_writes_namelist_def(self):
        """Test that namelist.def is created with basic entries."""
        StdI = _make_stdi_base(solver="HPhi")
        StdI.LCintra = 1
        StdI.LCinter = 0
        StdI.LHund = 0
        StdI.LEx = 0
        StdI.LPairLift = 0
        StdI.LPairHopp = 0
        StdI.Lintr = 0
        StdI.ioutputmode = 1
        StdI.SpectrumBody = 1
        StdI.method = "lanczos"
        StdI.PumpBody = 0
        StdI.CDataFileHead = "zvo"
        StdI.lBoost = 0

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_namelist(StdI)
                assert os.path.exists("namelist.def")
                content = open("namelist.def").read()
                assert "ModPara  modpara.def" in content
                assert "LocSpin  locspn.def" in content
                assert "Trans  trans.def" in content
                assert "CoulombIntra  coulombintra.def" in content
                assert "CalcMod  calcmod.def" in content
            finally:
                os.chdir(orig)

    def test_mvmc_entries(self):
        """Test mVMC-specific namelist entries."""
        StdI = _make_stdi_base(solver="mVMC")
        StdI.LCintra = 0
        StdI.LCinter = 0
        StdI.LHund = 0
        StdI.LEx = 0
        StdI.LPairLift = 0
        StdI.LPairHopp = 0
        StdI.Lintr = 0
        StdI.ioutputmode = 1
        StdI.lGC = 0
        StdI.Sz2 = 0

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_namelist(StdI)
                content = open("namelist.def").read()
                assert "Gutzwiller  gutzwilleridx.def" in content
                assert "Jastrow  jastrowidx.def" in content
                assert "Orbital  orbitalidx.def" in content
                assert "TransSym  qptransidx.def" in content
            finally:
                os.chdir(orig)


class TestNamelistBodyDispatch:
    """Tests for the namelist body-writer functions and helpers."""

    def test_interaction_flags_table(self):
        """Test that _INTERACTION_FLAGS has 7 entries with valid attributes."""
        assert len(_INTERACTION_FLAGS) == 7
        stdi = StdIntList()
        for attr, line in _INTERACTION_FLAGS:
            assert hasattr(stdi, attr), f"StdIntList missing attribute {attr}"
            assert line.endswith("\n")

    def test_hphi_body_calcmod(self):
        """Test that HPhi body writes CalcMod entry."""
        import io
        StdI = _make_stdi_base(solver="HPhi")
        StdI.SpectrumBody = 1
        StdI.method = "lanczos"
        StdI.PumpBody = 0
        StdI.CDataFileHead = "zvo"
        StdI.lBoost = 0
        buf = io.StringIO()
        _write_namelist_hphi(buf, StdI)
        content = buf.getvalue()
        assert "CalcMod  calcmod.def" in content
        assert "SingleExcitation  single.def" in content
        assert "SpectrumVec  zvo_eigenvec_0" in content

    def test_hphi_body_pair_excitation(self):
        """Test that HPhi body writes PairExcitation when SpectrumBody != 1."""
        import io
        StdI = _make_stdi_base(solver="HPhi")
        StdI.SpectrumBody = 0
        StdI.method = "lanczos"
        StdI.PumpBody = 0
        StdI.CDataFileHead = "zvo"
        StdI.lBoost = 0
        buf = io.StringIO()
        _write_namelist_hphi(buf, StdI)
        content = buf.getvalue()
        assert "PairExcitation  pair.def" in content
        assert "SingleExcitation" not in content

    def test_hphi_body_time_evolution_onebody(self):
        """Test HPhi body writes TEOneBody for PumpBody=1 in TE mode."""
        import io
        StdI = _make_stdi_base(solver="HPhi")
        StdI.SpectrumBody = 1
        StdI.method = "timeevolution"
        StdI.PumpBody = 1
        StdI.CDataFileHead = "zvo"
        StdI.lBoost = 0
        buf = io.StringIO()
        _write_namelist_hphi(buf, StdI)
        content = buf.getvalue()
        assert "TEOneBody  teone.def" in content
        assert "TETwoBody" not in content

    def test_hphi_body_time_evolution_twobody(self):
        """Test HPhi body writes TETwoBody for PumpBody=2 in TE mode."""
        import io
        StdI = _make_stdi_base(solver="HPhi")
        StdI.SpectrumBody = 1
        StdI.method = "timeevolution"
        StdI.PumpBody = 2
        StdI.CDataFileHead = "zvo"
        StdI.lBoost = 0
        buf = io.StringIO()
        _write_namelist_hphi(buf, StdI)
        content = buf.getvalue()
        assert "TETwoBody  tetwo.def" in content
        assert "TEOneBody" not in content

    def test_hphi_body_boost(self):
        """Test HPhi body writes Boost entry when lBoost=1."""
        import io
        StdI = _make_stdi_base(solver="HPhi")
        StdI.SpectrumBody = 1
        StdI.method = "lanczos"
        StdI.PumpBody = 0
        StdI.CDataFileHead = "zvo"
        StdI.lBoost = 1
        buf = io.StringIO()
        _write_namelist_hphi(buf, StdI)
        content = buf.getvalue()
        assert "Boost  boost.def" in content

    def test_mvmc_body_basic(self):
        """Test that mVMC body writes required entries."""
        import io
        StdI = _make_stdi_base(solver="mVMC")
        StdI.lGC = 0
        StdI.Sz2 = 0
        buf = io.StringIO()
        _write_namelist_mvmc(buf, StdI)
        content = buf.getvalue()
        assert "Gutzwiller  gutzwilleridx.def" in content
        assert "Jastrow  jastrowidx.def" in content
        assert "Orbital  orbitalidx.def" in content
        assert "TransSym  qptransidx.def" in content
        assert "OrbitalParallel" not in content

    def test_mvmc_body_gc_orbital_parallel(self):
        """Test that mVMC body writes OrbitalParallel when lGC=1."""
        import io
        StdI = _make_stdi_base(solver="mVMC")
        StdI.lGC = 1
        StdI.Sz2 = 0
        buf = io.StringIO()
        _write_namelist_mvmc(buf, StdI)
        content = buf.getvalue()
        assert "OrbitalParallel  orbitalidxpara.def" in content
        assert "OrbitalGeneral  orbitalidxgen.def" in content

    def test_mvmc_body_nonzero_sz2_orbital_parallel(self):
        """Test that mVMC body writes OrbitalParallel when Sz2 nonzero."""
        import io
        StdI = _make_stdi_base(solver="mVMC")
        StdI.lGC = 0
        StdI.Sz2 = 2
        buf = io.StringIO()
        _write_namelist_mvmc(buf, StdI)
        content = buf.getvalue()
        assert "OrbitalParallel  orbitalidxpara.def" in content

    def test_uhf_no_solver_specific_entries(self):
        """Test that UHF namelist has no solver-specific entries."""
        StdI = _make_stdi_base(solver="UHF")
        StdI.LCintra = 0
        StdI.LCinter = 0
        StdI.LHund = 0
        StdI.LEx = 0
        StdI.LPairLift = 0
        StdI.LPairHopp = 0
        StdI.Lintr = 0
        StdI.ioutputmode = 0

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_namelist(StdI)
                content = open("namelist.def").read()
                assert "ModPara" in content
                assert "LocSpin" in content
                assert "Trans" in content
                assert "CalcMod" not in content
                assert "Gutzwiller" not in content
            finally:
                os.chdir(orig)

    def test_interaction_flags_conditional(self):
        """Test that interaction flags correctly gate namelist entries."""
        StdI = _make_stdi_base(solver="HPhi")
        StdI.LCintra = 0
        StdI.LCinter = 1
        StdI.LHund = 0
        StdI.LEx = 1
        StdI.LPairLift = 0
        StdI.LPairHopp = 0
        StdI.Lintr = 1
        StdI.ioutputmode = 0
        StdI.SpectrumBody = 1
        StdI.method = "lanczos"
        StdI.PumpBody = 0
        StdI.CDataFileHead = "zvo"
        StdI.lBoost = 0

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_namelist(StdI)
                content = open("namelist.def").read()
                assert "CoulombIntra" not in content
                assert "CoulombInter  coulombinter.def" in content
                assert "Hund" not in content
                assert "Exchange  exchange.def" in content
                assert "PairLift" not in content
                assert "PairHop" not in content
                assert "InterAll  interall.def" in content
            finally:
                os.chdir(orig)


class TestPrint1Green:
    """Tests for the print_1_green function."""

    def test_no_output_when_mode_zero(self):
        """Test that no file is written when ioutputmode=0."""
        StdI = _make_stdi_base(nsite=2)
        StdI.ioutputmode = 0

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_1_green(StdI)
                assert not os.path.exists("greenone.def")
            finally:
                os.chdir(orig)

    def test_writes_greenone_correlation(self):
        """Test that greenone.def is written in correlation mode."""
        StdI = _make_stdi_base(nsite=2, model="hubbard")
        StdI.ioutputmode = 1
        StdI.NsiteUC = 1
        StdI.locspinflag = [0, 0]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_1_green(StdI)
                assert os.path.exists("greenone.def")
                content = open("greenone.def").read()
                assert "NCisAjs" in content
            finally:
                os.chdir(orig)

    def test_writes_greenone_raw(self):
        """Test that greenone.def is written in raw mode."""
        StdI = _make_stdi_base(nsite=2, model="hubbard")
        StdI.ioutputmode = 2
        StdI.locspinflag = [0, 0]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_1_green(StdI)
                assert os.path.exists("greenone.def")
                content = open("greenone.def").read()
                assert "NCisAjs" in content
            finally:
                os.chdir(orig)


class TestPrint2Green:
    """Tests for the print_2_green function."""

    def test_writes_greentwo_correlation(self):
        """Test that greentwo.def is written in correlation mode."""
        StdI = _make_stdi_base(nsite=2, model="hubbard", solver="HPhi")
        StdI.ioutputmode = 1
        StdI.NsiteUC = 1
        StdI.locspinflag = [0, 0]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_2_green(StdI)
                assert os.path.exists("greentwo.def")
                content = open("greentwo.def").read()
                assert "NCisAjsCktAltDC" in content
            finally:
                os.chdir(orig)

    def test_no_file_when_ioutputmode_off(self):
        """``ioutputmode`` 0 skips writing ``greentwo.def``."""
        StdI = _make_stdi_base(nsite=2, model="hubbard", solver="HPhi")
        StdI.ioutputmode = 0
        StdI.NsiteUC = 1
        StdI.locspinflag = [0, 0]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_2_green(StdI)
                assert not os.path.exists("greentwo.def")
            finally:
                os.chdir(orig)

    def test_writes_greentwo_raw_mode(self):
        """Raw / full output mode uses ``green2_raw`` index list."""
        StdI = _make_stdi_base(nsite=2, model="hubbard", solver="HPhi")
        StdI.ioutputmode = 2
        StdI.NsiteUC = 1
        StdI.locspinflag = [0, 0]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_2_green(StdI)
                assert os.path.exists("greentwo.def")
                content = open("greentwo.def").read()
                assert "NCisAjsCktAltDC" in content
            finally:
                os.chdir(orig)


class TestPrintInteractions:
    """Tests for the print_interactions function."""

    def test_no_interactions_no_files(self):
        """Test that no files are written when all interaction counts are zero."""
        StdI = _make_stdi_base(nsite=2)
        StdI.Cintra_list = []
        StdI.Cinter_list = []
        StdI.Hund_list = []
        StdI.Ex_list = []
        StdI.PairLift_list = []
        StdI.PairHopp_list = []
        StdI.intr_list = []

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_interactions(StdI)
                assert not os.path.exists("coulombintra.def")
                assert not os.path.exists("coulombinter.def")
                assert not os.path.exists("interall.def")
            finally:
                os.chdir(orig)

    def test_coulomb_intra_written(self):
        """Test that coulombintra.def is written with non-zero terms."""
        StdI = _make_stdi_base(nsite=2)
        StdI.Cintra_list = [(4.0, 0), (4.0, 1)]
        StdI.Cinter_list = []
        StdI.Hund_list = []
        StdI.Ex_list = []
        StdI.PairLift_list = []
        StdI.PairHopp_list = []
        StdI.intr_list = []

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_interactions(StdI)
                assert os.path.exists("coulombintra.def")
                content = open("coulombintra.def").read()
                assert "NCoulombIntra          2" in content
                assert StdI.LCintra == 1
            finally:
                os.chdir(orig)


class TestPrintModPara:
    """Tests for the print_mod_para function."""

    def test_writes_modpara_hphi(self):
        """Test that modpara.def is created for HPhi solver."""
        StdI = _make_stdi_base(solver="HPhi", nsite=4)
        StdI.CDataFileHead = "zvo"
        StdI.Lanczos_max = 2000
        StdI.initial_iv = -1
        StdI.nvec = NaN_i
        StdI.exct = 1
        StdI.LanczosEps = 14
        StdI.LanczosTarget = 2
        StdI.LargeValue = 10.0
        StdI.NumAve = 5
        StdI.ExpecInterval = 20
        StdI.Nomega = 200
        StdI.OmegaMax = 40.0
        StdI.OmegaMin = -40.0
        StdI.OmegaOrg = 0.0
        StdI.OmegaIm = 0.1
        StdI.method = "lanczos"
        StdI.ExpandCoef = 10

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_mod_para(StdI)
                assert os.path.exists("modpara.def")
                content = open("modpara.def").read()
                assert "HPhi_Cal_Parameters" in content
                assert "CDataFileHead  zvo" in content
                assert "Nsite" in content
            finally:
                os.chdir(orig)

    def test_writes_modpara_uhf(self):
        """Test that modpara.def is created for UHF solver."""
        StdI = _make_stdi_base(solver="UHF", nsite=4)
        StdI.CDataFileHead = "zvo"
        StdI.Iteration_max = 1000
        StdI.eps = 8
        StdI.mix = 0.5
        StdI.RndSeed = 123456789
        StdI.eps_slater = 6
        StdI.NMPTrans = 0

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_mod_para(StdI)
                content = open("modpara.def").read()
                assert "UHF_Cal_Parameters" in content
            finally:
                os.chdir(orig)

    # Note: HWAVE no longer writes modpara.def via print_mod_para (it is not an
    # ExpertModeSolverPlugin and overrides write()); the former
    # test_writes_modpara_hwave / unknown-solver tests were removed in B1.


class TestModparaBodyDispatch:
    """Tests for the modpara body-writer functions and banner."""

    def test_banner_has_uhf_and_hwave(self):
        assert SolverType.UHF in _MODPARA_BANNER
        assert SolverType.HWAVE in _MODPARA_BANNER

    def test_banner_values(self):
        assert _MODPARA_BANNER[SolverType.UHF] == "UHF_Cal_Parameters"
        assert _MODPARA_BANNER[SolverType.HWAVE] == "HWAVE_Cal_Parameters"

    def test_write_modpara_hphi_content(self):
        """Verify HPhi body writer produces expected fields."""
        StdI = _make_stdi_base(solver="HPhi", nsite=4)
        StdI.CDataFileHead = "zvo"
        StdI.Lanczos_max = 2000
        StdI.initial_iv = -1
        StdI.nvec = NaN_i
        StdI.exct = 1
        StdI.LanczosEps = 14
        StdI.LanczosTarget = 2
        StdI.LargeValue = 10.0
        StdI.NumAve = 5
        StdI.ExpecInterval = 20
        StdI.Nomega = 200
        StdI.OmegaMax = 40.0
        StdI.OmegaMin = -40.0
        StdI.OmegaOrg = 0.0
        StdI.OmegaIm = 0.1
        StdI.method = "lanczos"
        StdI.ExpandCoef = 10

        import io
        fp = io.StringIO()
        _write_modpara_hphi(fp, StdI)
        content = fp.getvalue()
        assert "HPhi_Cal_Parameters" in content
        assert "CDataFileHead  zvo" in content
        assert "Lanczos_max" in content
        assert "LargeValue" in content

    def test_write_modpara_mvmc_content(self):
        """Verify mVMC body writer produces expected fields."""
        StdI = _make_stdi_base(solver="mVMC", nsite=4)
        StdI.CDataFileHead = "zvo"
        StdI.CParaFileHead = "zqp"
        StdI.NVMCCalMode = 0
        StdI.NLanczosMode = 0
        StdI.NDataIdxStart = 1
        StdI.NDataQtySmp = 1
        StdI.ncond = 2
        StdI.NMPTrans = -1
        StdI.NSROptItrStep = 1000
        StdI.NSROptItrSmp = 100
        StdI.DSROptRedCut = 0.001
        StdI.DSROptStaDel = 0.02
        StdI.DSROptStepDt = 0.02
        StdI.NVMCWarmUp = 10
        StdI.NVMCInterval = 1
        StdI.NVMCSample = 1000
        StdI.NExUpdatePath = 0
        StdI.RndSeed = 123456789
        StdI.NSplitSize = 1
        StdI.NStore = 1
        StdI.NSRCG = 0

        import io
        fp = io.StringIO()
        _write_modpara_mvmc(fp, StdI)
        content = fp.getvalue()
        assert "VMC_Cal_Parameters" in content
        assert "NVMCSample" in content
        assert "DSROptRedCut" in content

    def test_write_modpara_uhf_hwave_uhf(self):
        """Verify UHF banner in shared writer."""
        StdI = _make_stdi_base(solver="UHF", nsite=4)
        StdI.CDataFileHead = "zvo"
        StdI.Iteration_max = 1000
        StdI.eps = 8
        StdI.mix = 0.5
        StdI.RndSeed = 123456789
        StdI.eps_slater = 6
        StdI.NMPTrans = 0

        import io
        fp = io.StringIO()
        _write_modpara_uhf_hwave(fp, StdI)
        content = fp.getvalue()
        assert "UHF_Cal_Parameters" in content
        assert "IterationMax" in content
        assert "EpsSlater" in content

    def test_write_modpara_uhf_hwave_hwave(self):
        """Verify HWAVE banner in shared writer."""
        StdI = _make_stdi_base(solver="HWAVE", nsite=4)
        StdI.CDataFileHead = "zvo"
        StdI.Iteration_max = 1000
        StdI.eps = 8
        StdI.mix = 0.5
        StdI.RndSeed = 123456789
        StdI.eps_slater = 6
        StdI.NMPTrans = 0

        import io
        fp = io.StringIO()
        _write_modpara_uhf_hwave(fp, StdI)
        content = fp.getvalue()
        assert "HWAVE_Cal_Parameters" in content
        assert "UHF_Cal_Parameters" not in content

    def test_hphi_omits_expand_coef_for_non_te(self):
        """ExpandCoef only written for time-evolution method."""
        StdI = _make_stdi_base(solver="HPhi", nsite=4)
        StdI.CDataFileHead = "zvo"
        StdI.Lanczos_max = 2000
        StdI.initial_iv = -1
        StdI.nvec = NaN_i
        StdI.exct = 1
        StdI.LanczosEps = 14
        StdI.LanczosTarget = 2
        StdI.LargeValue = 10.0
        StdI.NumAve = 5
        StdI.ExpecInterval = 20
        StdI.Nomega = 200
        StdI.OmegaMax = 40.0
        StdI.OmegaMin = -40.0
        StdI.OmegaOrg = 0.0
        StdI.OmegaIm = 0.1
        StdI.method = "lanczos"
        StdI.ExpandCoef = 10

        import io
        fp = io.StringIO()
        _write_modpara_hphi(fp, StdI)
        content = fp.getvalue()
        assert "ExpandCoef" not in content

    def test_hphi_includes_expand_coef_for_te(self):
        """ExpandCoef written for time-evolution method."""
        StdI = _make_stdi_base(solver="HPhi", nsite=4)
        StdI.CDataFileHead = "zvo"
        StdI.Lanczos_max = 2000
        StdI.initial_iv = -1
        StdI.nvec = NaN_i
        StdI.exct = 1
        StdI.LanczosEps = 14
        StdI.LanczosTarget = 2
        StdI.LargeValue = 10.0
        StdI.NumAve = 5
        StdI.ExpecInterval = 20
        StdI.Nomega = 200
        StdI.OmegaMax = 40.0
        StdI.OmegaMin = -40.0
        StdI.OmegaOrg = 0.0
        StdI.OmegaIm = 0.1
        StdI.method = "timeevolution"
        StdI.ExpandCoef = 10

        import io
        fp = io.StringIO()
        _write_modpara_hphi(fp, StdI)
        content = fp.getvalue()
        assert "ExpandCoef" in content


class TestOutputModeToInt:
    """Tests for the OUTPUT_MODE_TO_INT dispatch dict."""

    # -- "off" group → 0 --

    def test_non(self):
        assert OUTPUT_MODE_TO_INT["non"] == 0

    def test_none(self):
        assert OUTPUT_MODE_TO_INT["none"] == 0

    def test_off(self):
        assert OUTPUT_MODE_TO_INT["off"] == 0

    # -- "correlation" group → 1 --

    def test_cor(self):
        assert OUTPUT_MODE_TO_INT["cor"] == 1

    def test_corr(self):
        assert OUTPUT_MODE_TO_INT["corr"] == 1

    def test_correlation(self):
        assert OUTPUT_MODE_TO_INT["correlation"] == 1

    # -- "raw" group → 2 --

    def test_raw(self):
        assert OUTPUT_MODE_TO_INT["raw"] == 2

    def test_all(self):
        assert OUTPUT_MODE_TO_INT["all"] == 2

    def test_full(self):
        assert OUTPUT_MODE_TO_INT["full"] == 2

    # -- unknown --

    def test_unknown_returns_none(self):
        assert OUTPUT_MODE_TO_INT.get("bogus") is None


class TestModelGcToExUpdatePath:
    """Tests for the MODEL_GC_TO_EX_UPDATE_PATH dispatch dict."""

    # -- Hubbard: always 0 regardless of lGC --

    def test_hubbard_canonical(self):
        assert MODEL_GC_TO_EX_UPDATE_PATH[(ModelType.HUBBARD, 0)] == 0

    def test_hubbard_gc(self):
        assert MODEL_GC_TO_EX_UPDATE_PATH[(ModelType.HUBBARD, 1)] == 0

    # -- Spin: always 2 regardless of lGC --

    def test_spin_canonical(self):
        assert MODEL_GC_TO_EX_UPDATE_PATH[(ModelType.SPIN, 0)] == 2

    def test_spin_gc(self):
        assert MODEL_GC_TO_EX_UPDATE_PATH[(ModelType.SPIN, 1)] == 2

    # -- Kondo: 1 for canonical, 3 for grand-canonical --

    def test_kondo_canonical(self):
        assert MODEL_GC_TO_EX_UPDATE_PATH[(ModelType.KONDO, 0)] == 1

    def test_kondo_gc(self):
        assert MODEL_GC_TO_EX_UPDATE_PATH[(ModelType.KONDO, 1)] == 3

    # -- completeness --

    def test_dict_has_six_entries(self):
        assert len(MODEL_GC_TO_EX_UPDATE_PATH) == 6

    def test_unknown_returns_none(self):
        assert MODEL_GC_TO_EX_UPDATE_PATH.get(("bogus", 0)) is None


class TestSolverDefaultsDispatch:
    """Tests for the solver-defaults helper functions (set via plugins)."""

    # -- HPhi handler sets defaults --

    def test_hphi_sets_lanczos_max(self):
        StdI = _make_stdi_base(solver="HPhi")
        StdI.LargeValue = 10.0
        StdI.Lanczos_max = NaN_i
        StdI.initial_iv = NaN_i
        StdI.exct = NaN_i
        StdI.LanczosEps = NaN_i
        StdI.LanczosTarget = NaN_i
        StdI.NumAve = NaN_i
        StdI.ExpecInterval = NaN_i
        StdI.Nomega = NaN_i
        _check_mod_para_hphi(StdI)
        assert StdI.Lanczos_max == 2000

    def test_hphi_sets_exct(self):
        StdI = _make_stdi_base(solver="HPhi")
        StdI.LargeValue = 10.0
        StdI.Lanczos_max = NaN_i
        StdI.initial_iv = NaN_i
        StdI.exct = NaN_i
        StdI.LanczosEps = NaN_i
        StdI.LanczosTarget = NaN_i
        StdI.NumAve = NaN_i
        StdI.ExpecInterval = NaN_i
        StdI.Nomega = NaN_i
        _check_mod_para_hphi(StdI)
        assert StdI.exct == 1

    def test_hphi_sets_omega_defaults(self):
        StdI = _make_stdi_base(solver="HPhi", nsite=4)
        StdI.LargeValue = 10.0
        StdI.Lanczos_max = NaN_i
        StdI.initial_iv = NaN_i
        StdI.exct = NaN_i
        StdI.LanczosEps = NaN_i
        StdI.LanczosTarget = NaN_i
        StdI.NumAve = NaN_i
        StdI.ExpecInterval = NaN_i
        StdI.Nomega = NaN_i
        StdI.OmegaMax = float("nan")
        StdI.OmegaMin = float("nan")
        StdI.OmegaOrg = float("nan")
        StdI.OmegaIm = float("nan")
        _check_mod_para_hphi(StdI)
        assert StdI.Nomega == 200
        assert StdI.OmegaMax == 40.0
        assert StdI.OmegaMin == -40.0

    def test_hphi_clamps_lanczos_target_when_below_exct(self):
        """``LanczosTarget`` is raised to at least ``exct`` when smaller."""
        StdI = _make_stdi_base(solver="HPhi", nsite=4)
        StdI.LargeValue = 10.0
        StdI.Lanczos_max = NaN_i
        StdI.initial_iv = NaN_i
        StdI.exct = 5
        StdI.LanczosEps = NaN_i
        StdI.LanczosTarget = 2
        StdI.NumAve = NaN_i
        StdI.ExpecInterval = NaN_i
        StdI.Nomega = NaN_i
        _check_mod_para_hphi(StdI)
        assert StdI.LanczosTarget == 5

    # -- mVMC handler sets defaults --

    def test_mvmc_sets_cpara_default(self):
        StdI = _make_stdi_base(solver="mVMC", model="hubbard")
        StdI.NVMCCalMode = NaN_i
        StdI.NLanczosMode = NaN_i
        StdI.NDataIdxStart = NaN_i
        StdI.NDataQtySmp = NaN_i
        StdI.NSPGaussLeg = NaN_i
        StdI.NSPStot = NaN_i
        StdI.NMPTrans = NaN_i
        StdI.NSROptItrStep = NaN_i
        StdI.NSROptItrSmp = NaN_i
        StdI.NVMCWarmUp = NaN_i
        StdI.NVMCInterval = NaN_i
        StdI.NVMCSample = NaN_i
        StdI.RndSeed = NaN_i
        StdI.NSplitSize = NaN_i
        StdI.NStore = NaN_i
        StdI.NSRCG = NaN_i
        _check_mod_para_mvmc(StdI)
        assert StdI.CParaFileHead == "zqp"

    def test_mvmc_respects_explicit_cparafilehead(self):
        """Non-sentinel ``CParaFileHead`` takes the else branch (no default)."""
        StdI = _make_stdi_base(solver="mVMC", model="hubbard")
        StdI.CParaFileHead = "my_cpara"
        StdI.NVMCCalMode = NaN_i
        StdI.NLanczosMode = NaN_i
        StdI.NDataIdxStart = NaN_i
        StdI.NDataQtySmp = NaN_i
        StdI.NSPGaussLeg = NaN_i
        StdI.NSPStot = NaN_i
        StdI.NMPTrans = NaN_i
        StdI.NSROptItrStep = NaN_i
        StdI.NSROptItrSmp = NaN_i
        StdI.NVMCWarmUp = NaN_i
        StdI.NVMCInterval = NaN_i
        StdI.NVMCSample = NaN_i
        StdI.RndSeed = NaN_i
        StdI.NSplitSize = NaN_i
        StdI.NStore = NaN_i
        StdI.NSRCG = NaN_i
        _check_mod_para_mvmc(StdI)
        assert StdI.CParaFileHead == "my_cpara"

    def test_mvmc_nvmccalmode_one_marks_nsroptitrsmp_not_used(self):
        """``NVMCCalMode == 1`` runs the ``not_used_i`` branch for ``NSROptItrSmp``."""
        StdI = _make_stdi_base(solver="mVMC", model="hubbard")
        StdI.NVMCCalMode = 1
        StdI.NLanczosMode = NaN_i
        StdI.NDataIdxStart = NaN_i
        StdI.NDataQtySmp = NaN_i
        StdI.NSPGaussLeg = NaN_i
        StdI.NSPStot = NaN_i
        StdI.NMPTrans = NaN_i
        StdI.NSROptItrStep = 1000
        StdI.NSROptItrSmp = NaN_i
        StdI.NVMCWarmUp = NaN_i
        StdI.NVMCInterval = NaN_i
        StdI.NVMCSample = NaN_i
        StdI.RndSeed = NaN_i
        StdI.NSplitSize = NaN_i
        StdI.NStore = NaN_i
        StdI.NSRCG = NaN_i
        _check_mod_para_mvmc(StdI)
        assert StdI.NSROptItrSmp == 100

    def test_mvmc_sets_nvmccalmode(self):
        StdI = _make_stdi_base(solver="mVMC", model="hubbard")
        StdI.NVMCCalMode = NaN_i
        StdI.NLanczosMode = NaN_i
        StdI.NDataIdxStart = NaN_i
        StdI.NDataQtySmp = NaN_i
        StdI.NSPGaussLeg = NaN_i
        StdI.NSPStot = NaN_i
        StdI.NMPTrans = NaN_i
        StdI.NSROptItrStep = NaN_i
        StdI.NSROptItrSmp = NaN_i
        StdI.NVMCWarmUp = NaN_i
        StdI.NVMCInterval = NaN_i
        StdI.NVMCSample = NaN_i
        StdI.RndSeed = NaN_i
        StdI.NSplitSize = NaN_i
        StdI.NStore = NaN_i
        StdI.NSRCG = NaN_i
        _check_mod_para_mvmc(StdI)
        assert StdI.NVMCCalMode == 0

    def test_mvmc_sets_rndseed(self):
        StdI = _make_stdi_base(solver="mVMC", model="hubbard")
        StdI.NVMCCalMode = NaN_i
        StdI.NLanczosMode = NaN_i
        StdI.NDataIdxStart = NaN_i
        StdI.NDataQtySmp = NaN_i
        StdI.NSPGaussLeg = NaN_i
        StdI.NSPStot = NaN_i
        StdI.NMPTrans = NaN_i
        StdI.NSROptItrStep = NaN_i
        StdI.NSROptItrSmp = NaN_i
        StdI.NVMCWarmUp = NaN_i
        StdI.NVMCInterval = NaN_i
        StdI.NVMCSample = NaN_i
        StdI.RndSeed = NaN_i
        StdI.NSplitSize = NaN_i
        StdI.NStore = NaN_i
        StdI.NSRCG = NaN_i
        _check_mod_para_mvmc(StdI)
        assert StdI.RndSeed == 123456789

    # -- UHF handler sets defaults --

    def test_uhf_sets_rndseed(self):
        StdI = _make_stdi_base(solver="UHF")
        StdI.RndSeed = NaN_i
        StdI.Iteration_max = NaN_i
        StdI.eps = NaN_i
        StdI.eps_slater = NaN_i
        StdI.NMPTrans = NaN_i
        _check_mod_para_uhf(StdI)
        assert StdI.RndSeed == 123456789

    def test_uhf_sets_iteration_max(self):
        StdI = _make_stdi_base(solver="UHF")
        StdI.RndSeed = NaN_i
        StdI.Iteration_max = NaN_i
        StdI.eps = NaN_i
        StdI.eps_slater = NaN_i
        StdI.NMPTrans = NaN_i
        _check_mod_para_uhf(StdI)
        assert StdI.Iteration_max == 1000

    def test_uhf_sets_mix(self):
        StdI = _make_stdi_base(solver="UHF")
        StdI.RndSeed = NaN_i
        StdI.Iteration_max = NaN_i
        StdI.mix = float("nan")
        StdI.eps = NaN_i
        StdI.eps_slater = NaN_i
        StdI.NMPTrans = NaN_i
        _check_mod_para_uhf(StdI)
        assert StdI.mix == 0.5


class TestCheckConservedQuantities:
    """Tests for _check_conserved_quantities helper."""

    # -- Hubbard + non-HPhi + canonical: 2Sz defaults to 0 --

    def test_hubbard_mvmc_canonical_sz2_defaults(self):
        StdI = _make_stdi_base(model="hubbard", solver="mVMC", lGC=0)
        StdI.ncond = 4
        StdI.Sz2 = NaN_i
        _check_conserved_quantities(StdI)
        assert StdI.Sz2 == 0

    def test_hubbard_mvmc_canonical_ncond_unchanged(self):
        StdI = _make_stdi_base(model="hubbard", solver="mVMC", lGC=0)
        StdI.ncond = 4
        StdI.Sz2 = NaN_i
        _check_conserved_quantities(StdI)
        assert StdI.ncond == 4

    # -- Hubbard + non-HPhi + GC: 2Sz not used, no changes --

    def test_hubbard_mvmc_gc_no_sz2_change(self):
        StdI = _make_stdi_base(model="hubbard", solver="mVMC", lGC=1)
        StdI.ncond = 4
        StdI.Sz2 = NaN_i
        _check_conserved_quantities(StdI)
        assert StdI.Sz2 == NaN_i

    # -- Hubbard + HPhi + canonical: nelec required (passes if set) --

    def test_hubbard_hphi_canonical_passes_with_ncond(self):
        StdI = _make_stdi_base(model="hubbard", solver="HPhi", lGC=0)
        StdI.ncond = 4
        StdI.Sz2 = NaN_i
        _check_conserved_quantities(StdI)  # should not exit
        assert StdI.ncond == 4

    # -- Hubbard + HPhi + GC: nelec and 2Sz not used --

    def test_hubbard_hphi_gc_passes(self):
        StdI = _make_stdi_base(model="hubbard", solver="HPhi", lGC=1)
        StdI.ncond = NaN_i
        StdI.Sz2 = NaN_i
        _check_conserved_quantities(StdI)  # should not exit

    # -- Spin + mVMC: ncond set to 0 --

    def test_spin_mvmc_sets_ncond_zero(self):
        StdI = _make_stdi_base(model="spin", solver="mVMC", lGC=0)
        StdI.ncond = NaN_i
        StdI.Sz2 = 0
        _check_conserved_quantities(StdI)
        assert StdI.ncond == 0

    # -- Spin + HPhi + canonical: 2Sz required --

    def test_spin_hphi_canonical_passes_with_sz2(self):
        StdI = _make_stdi_base(model="spin", solver="HPhi", lGC=0)
        StdI.ncond = NaN_i
        StdI.Sz2 = 0
        _check_conserved_quantities(StdI)  # should not exit
        assert StdI.Sz2 == 0

    # -- Spin + HPhi + GC: 2Sz not used --

    def test_spin_hphi_gc_passes(self):
        StdI = _make_stdi_base(model="spin", solver="HPhi", lGC=1)
        StdI.ncond = NaN_i
        StdI.Sz2 = NaN_i
        _check_conserved_quantities(StdI)  # should not exit

    # -- Kondo + non-HPhi + canonical: 2Sz defaults to 0 --

    def test_kondo_mvmc_canonical_sz2_defaults(self):
        StdI = _make_stdi_base(model="kondo", solver="mVMC", lGC=0)
        StdI.ncond = 4
        StdI.Sz2 = NaN_i
        _check_conserved_quantities(StdI)
        assert StdI.Sz2 == 0

    # -- Kondo + HPhi + canonical: ncond required --

    def test_kondo_hphi_canonical_passes_with_ncond(self):
        StdI = _make_stdi_base(model="kondo", solver="HPhi", lGC=0)
        StdI.ncond = 4
        StdI.Sz2 = NaN_i
        _check_conserved_quantities(StdI)  # should not exit
        assert StdI.ncond == 4

    # -- Kondo + HPhi + GC: nelec and 2Sz not used --

    def test_kondo_hphi_gc_passes(self):
        StdI = _make_stdi_base(model="kondo", solver="HPhi", lGC=1)
        StdI.ncond = NaN_i
        StdI.Sz2 = NaN_i
        _check_conserved_quantities(StdI)  # should not exit

    def test_invalid_lgc_has_no_rule_returns_early(self):
        """``(model, is_hphi, lGC)`` not in the table → no validation."""
        StdI = _make_stdi_base(model="hubbard", solver="HPhi")
        StdI.lGC = 2
        _check_conserved_quantities(StdI)


class TestConservedQtyRulesTable:
    """Tests for the _CONSERVED_QTY_RULES dispatch table structure."""

    def test_has_12_entries(self):
        """Test that table has exactly 12 entries (3 models × 2 hphi × 2 lGC)."""
        assert len(_CONSERVED_QTY_RULES) == 12

    def test_all_models_present(self):
        """Test that all three model types appear in the rules."""
        models = {key[0] for key in _CONSERVED_QTY_RULES}
        assert ModelType.HUBBARD in models
        assert ModelType.SPIN in models
        assert ModelType.KONDO in models

    def test_both_hphi_flags_per_model(self):
        """Test that each model has both is_hphi=True and is_hphi=False."""
        for model in (ModelType.HUBBARD, ModelType.SPIN, ModelType.KONDO):
            hphi_flags = {key[1] for key in _CONSERVED_QTY_RULES if key[0] == model}
            assert True in hphi_flags, f"{model}: missing is_hphi=True"
            assert False in hphi_flags, f"{model}: missing is_hphi=False"

    def test_both_lgc_values_per_model(self):
        """Test that each model has both lGC=0 and lGC=1."""
        for model in (ModelType.HUBBARD, ModelType.SPIN, ModelType.KONDO):
            lgc_vals = {key[2] for key in _CONSERVED_QTY_RULES if key[0] == model}
            assert 0 in lgc_vals, f"{model}: missing lGC=0"
            assert 1 in lgc_vals, f"{model}: missing lGC=1"

    def test_entries_are_3_tuples(self):
        """Test that every rule value is a (label, ncond_action, sz2_action) 3-tuple."""
        for key, rule in _CONSERVED_QTY_RULES.items():
            assert isinstance(rule, tuple), f"{key}: not tuple"
            assert len(rule) == 3, f"{key}: len={len(rule)}"

    def test_ncond_label_is_string(self):
        """Test that the ncond label is always a string."""
        for key, (label, _, _) in _CONSERVED_QTY_RULES.items():
            assert isinstance(label, str), f"{key}: label {label!r} not str"
            assert label in ("nelec", "ncond"), f"{key}: unexpected label {label!r}"

    def test_actions_are_valid(self):
        """Test that actions are valid action strings or None."""
        valid_actions = {"required", "not_used", "default_0", None}
        for key, (_, ncond_act, sz2_act) in _CONSERVED_QTY_RULES.items():
            assert ncond_act in valid_actions, f"{key}: bad ncond_action {ncond_act!r}"
            assert sz2_act in valid_actions, f"{key}: bad sz2_action {sz2_act!r}"

    def test_hubbard_hphi_canonical_uses_nelec(self):
        """Test that Hubbard+HPhi+canonical uses 'nelec' label."""
        rule = _CONSERVED_QTY_RULES[(ModelType.HUBBARD, True, 0)]
        assert rule[0] == "nelec"
        assert rule[1] == "required"

    def test_hubbard_non_hphi_canonical_defaults_sz2(self):
        """Test that Hubbard+non-HPhi+canonical defaults Sz2."""
        rule = _CONSERVED_QTY_RULES[(ModelType.HUBBARD, False, 0)]
        assert rule[0] == "ncond"
        assert rule[1] == "required"
        assert rule[2] == "default_0"

    def test_spin_rules_identical_for_hphi_flag(self):
        """Test that Spin model rules are the same regardless of is_hphi."""
        for lgc in (0, 1):
            assert (_CONSERVED_QTY_RULES[(ModelType.SPIN, True, lgc)]
                    == _CONSERVED_QTY_RULES[(ModelType.SPIN, False, lgc)])


# ===================================================================
#  _merge_duplicate_terms
# ===================================================================


class TestMergeDuplicateTerms:
    """Tests for _merge_duplicate_terms helper."""

    def test_no_duplicates(self):
        """Test that non-duplicate terms are unchanged."""
        indx = [[0, 0, 1, 0], [0, 0, 2, 0], [0, 0, 3, 0]]
        vals = [1.0 + 0j, 2.0 + 0j, 3.0 + 0j]
        count = _merge_duplicate_terms(indx, vals, 3)
        assert count == 3
        assert vals[0] == pytest.approx(1.0)
        assert vals[1] == pytest.approx(2.0)
        assert vals[2] == pytest.approx(3.0)

    def test_two_duplicates_merged(self):
        """Test that duplicate index quadruples are merged."""
        indx = [[0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 2, 0]]
        vals = [1.0 + 0j, 2.0 + 0j, 3.0 + 0j]
        count = _merge_duplicate_terms(indx, vals, 3)
        assert count == 2  # merged pair + one unique
        assert vals[0] == pytest.approx(3.0)  # 1 + 2
        assert vals[1] == pytest.approx(0.0)  # zeroed
        assert vals[2] == pytest.approx(3.0)  # unchanged

    def test_three_duplicates_merged(self):
        """Test that three duplicates are merged into one."""
        indx = [[1, 0, 2, 1]] * 3
        vals = [1.0 + 0j, 1.0 + 0j, 1.0 + 0j]
        count = _merge_duplicate_terms(indx, vals, 3)
        assert count == 1
        assert vals[0] == pytest.approx(3.0)
        assert vals[1] == pytest.approx(0.0)
        assert vals[2] == pytest.approx(0.0)

    def test_cancellation_yields_zero_count(self):
        """Test that cancelling terms result in zero count."""
        indx = [[0, 0, 1, 0], [0, 0, 1, 0]]
        vals = [1.0 + 0j, -1.0 + 0j]
        count = _merge_duplicate_terms(indx, vals, 2)
        assert count == 0

    def test_negligible_entries_excluded(self):
        """Test that entries below 1e-6 threshold are excluded from count."""
        indx = [[0, 0, 1, 0], [0, 0, 2, 0]]
        vals = [1e-7 + 0j, 1.0 + 0j]
        count = _merge_duplicate_terms(indx, vals, 2)
        assert count == 1

    def test_empty_input(self):
        """Test that zero-length input returns zero."""
        count = _merge_duplicate_terms([], [], 0)
        assert count == 0

    def test_single_entry(self):
        """Test with a single entry."""
        indx = [[0, 0, 1, 0]]
        vals = [5.0 + 0j]
        count = _merge_duplicate_terms(indx, vals, 1)
        assert count == 1
        assert vals[0] == pytest.approx(5.0)

    def test_complex_amplitudes(self):
        """Test that complex amplitudes are summed correctly."""
        indx = [[0, 0, 1, 0], [0, 0, 1, 0]]
        vals = [1.0 + 2.0j, 3.0 + 4.0j]
        count = _merge_duplicate_terms(indx, vals, 2)
        assert count == 1
        assert vals[0] == pytest.approx(4.0 + 6.0j)

    def test_partial_index_match_not_merged(self):
        """Test that partial index matches are NOT merged."""
        indx = [[0, 0, 1, 0], [0, 0, 1, 1]]  # differ in 4th element
        vals = [1.0 + 0j, 2.0 + 0j]
        count = _merge_duplicate_terms(indx, vals, 2)
        assert count == 2
        assert vals[0] == pytest.approx(1.0)
        assert vals[1] == pytest.approx(2.0)

    def test_n_less_than_array_length(self):
        """Test that only first n entries are processed."""
        indx = [[0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0]]
        vals = [1.0 + 0j, 2.0 + 0j, 3.0 + 0j]
        count = _merge_duplicate_terms(indx, vals, 2)
        assert count == 1  # only first 2 processed
        assert vals[0] == pytest.approx(3.0)
        assert vals[1] == pytest.approx(0.0)
        assert vals[2] == pytest.approx(3.0)  # untouched


# ===================================================================
#  GreenFunctionIndices class
# ===================================================================


class TestGreenFunctionIndicesConstruction:
    """Tests for GreenFunctionIndices construction and attributes."""

    def test_basic_construction(self):
        """Test that the class stores its attributes."""
        gf = GreenFunctionIndices(
            nsite=4, NsiteUC=2, locspinflag=[0, 0, 0, 0],
            is_kondo=False, is_mvmc=False,
        )
        assert gf.nsite == 4
        assert gf.NsiteUC == 2
        assert gf.locspinflag == [0, 0, 0, 0]
        assert gf.is_kondo is False
        assert gf.is_mvmc is False

    def test_is_mvmc_defaults_false(self):
        """Test that is_mvmc defaults to False."""
        gf = GreenFunctionIndices(
            nsite=2, NsiteUC=1, locspinflag=[0, 0], is_kondo=False,
        )
        assert gf.is_mvmc is False

    def test_kondo_mode(self):
        """Test that Kondo flag is stored."""
        gf = GreenFunctionIndices(
            nsite=4, NsiteUC=1, locspinflag=[0, 0, 1, 1], is_kondo=True,
        )
        assert gf.is_kondo is True


class TestGreenFunctionIndicesHelpers:
    """Tests for GreenFunctionIndices low-level helper methods."""

    def test_spin_max_itinerant(self):
        """Test spin_max returns 1 for itinerant site."""
        gf = GreenFunctionIndices(3, 1, [0, 0, 0], is_kondo=False)
        assert gf.spin_max(0) == 1
        assert gf.spin_max(2) == 1

    def test_spin_max_local_spin_half(self):
        """Test spin_max returns 1 for S=1/2 local-spin site (flag=1)."""
        gf = GreenFunctionIndices(3, 1, [1, 0, 0], is_kondo=False)
        assert gf.spin_max(0) == 1

    def test_spin_max_local_spin_one(self):
        """Test spin_max returns 2 for S=1 local-spin site (flag=2)."""
        gf = GreenFunctionIndices(3, 1, [0, 2, 0], is_kondo=False)
        assert gf.spin_max(1) == 2

    def test_spin_max_local_spin_three_half(self):
        """Test spin_max returns locspinflag value for S=3/2 site (flag=3)."""
        gf = GreenFunctionIndices(3, 1, [0, 0, 3], is_kondo=False)
        assert gf.spin_max(2) == 3

    def test_skip_local_spin_pair_both_itinerant(self):
        """Test that two itinerant sites are not skipped."""
        gf = GreenFunctionIndices(2, 1, [0, 0], is_kondo=False)
        assert gf.skip_local_spin_pair(0, 1) is False

    def test_skip_local_spin_pair_both_local_distinct(self):
        """Test that two distinct local-spin sites are skipped."""
        gf = GreenFunctionIndices(2, 1, [1, 1], is_kondo=False)
        assert gf.skip_local_spin_pair(0, 1) is True
        gf2 = GreenFunctionIndices(2, 1, [2, 3], is_kondo=False)
        assert gf2.skip_local_spin_pair(0, 1) is True

    def test_skip_local_spin_pair_same_site(self):
        """Test that same-site pair is never skipped (even if local-spin)."""
        gf = GreenFunctionIndices(2, 1, [1, 1], is_kondo=False)
        assert gf.skip_local_spin_pair(0, 0) is False
        gf2 = GreenFunctionIndices(2, 1, [2, 2], is_kondo=False)
        assert gf2.skip_local_spin_pair(1, 1) is False

    def test_skip_local_spin_pair_one_itinerant(self):
        """Test that mixed itinerant/local-spin pair is not skipped."""
        gf = GreenFunctionIndices(2, 1, [0, 1], is_kondo=False)
        assert gf.skip_local_spin_pair(0, 1) is False
        gf2 = GreenFunctionIndices(2, 1, [1, 0], is_kondo=False)
        assert gf2.skip_local_spin_pair(0, 1) is False

    def test_kondo_site_within_uc(self):
        """Test kondo_site maps within-UC site to itself."""
        gf = GreenFunctionIndices(8, 2, [0]*8, is_kondo=True)
        assert gf.kondo_site(0) == 0
        assert gf.kondo_site(1) == 1

    def test_kondo_site_beyond_uc(self):
        """Test kondo_site maps beyond-UC site to second half."""
        gf = GreenFunctionIndices(8, 2, [0]*8, is_kondo=True)
        assert gf.kondo_site(2) == 4
        assert gf.kondo_site(3) == 5

    def test_kondo_site_single_site_uc(self):
        """Test kondo_site with NsiteUC=1."""
        gf = GreenFunctionIndices(4, 1, [0]*4, is_kondo=True)
        assert gf.kondo_site(0) == 0
        assert gf.kondo_site(1) == 2  # 1 - 1 + 2 = 2


class TestGreenFunctionIndicesGreen1:
    """Tests for GreenFunctionIndices one-body index generators."""

    def test_green1_corr_2site_hubbard_basic(self):
        """Test 2-site Hubbard: NsiteUC=1, nsite=2, itinerant."""
        gf = GreenFunctionIndices(
            nsite=2, NsiteUC=1, locspinflag=[0, 0], is_kondo=False,
        )
        indices = gf.green1_corr()
        assert len(indices) == 4
        assert (0, 0, 0, 0) in indices
        assert (0, 0, 1, 0) in indices
        assert (0, 1, 0, 1) in indices
        assert (0, 1, 1, 1) in indices

    def test_green1_corr_only_same_spin(self):
        """Test that green1_corr only keeps same-spin pairs."""
        gf = GreenFunctionIndices(
            nsite=1, NsiteUC=1, locspinflag=[0], is_kondo=False,
        )
        for isite, ispin, jsite, jspin in gf.green1_corr():
            assert ispin == jspin

    def test_green1_corr_kondo(self):
        """Test that green1_corr doubles UC range in Kondo mode."""
        gf_normal = GreenFunctionIndices(
            nsite=4, NsiteUC=1, locspinflag=[0, 0, 0, 0], is_kondo=False,
        )
        gf_kondo = GreenFunctionIndices(
            nsite=4, NsiteUC=1, locspinflag=[0, 0, 0, 0], is_kondo=True,
        )
        assert len(gf_kondo.green1_corr()) > len(gf_normal.green1_corr())

    def test_green1_corr_local_spin_pair_skipped(self):
        """Test that distinct local-spin site pairs are skipped in corr mode."""
        gf = GreenFunctionIndices(
            nsite=2, NsiteUC=1, locspinflag=[1, 1], is_kondo=False,
        )
        for isite, ispin, jsite, jspin in gf.green1_corr():
            if isite != jsite:
                assert False, "distinct local-spin pair should be skipped"

    def test_green1_raw_2site_itinerant(self):
        """Test 2-site itinerant raw: all combinations."""
        gf = GreenFunctionIndices(
            nsite=2, NsiteUC=1, locspinflag=[0, 0], is_kondo=False,
        )
        indices = gf.green1_raw()
        # 2 sites × 2 spins × 2 sites × 2 spins = 16
        assert len(indices) == 16

    def test_green1_raw_1site_itinerant(self):
        """Test 1-site itinerant raw: 4 combinations."""
        gf = GreenFunctionIndices(
            nsite=1, NsiteUC=1, locspinflag=[0], is_kondo=False,
        )
        indices = gf.green1_raw()
        assert len(indices) == 4
        assert (0, 0, 0, 0) in indices
        assert (0, 0, 0, 1) in indices
        assert (0, 1, 0, 0) in indices
        assert (0, 1, 0, 1) in indices

    def test_green1_raw_local_spin_pair_skipped(self):
        """Test that distinct local-spin pairs are skipped in raw mode."""
        gf = GreenFunctionIndices(
            nsite=2, NsiteUC=1, locspinflag=[1, 1], is_kondo=False,
        )
        for isite, ispin, jsite, jspin in gf.green1_raw():
            if isite != jsite:
                assert False, "distinct local-spin pair should be skipped"


class TestGreenFunctionIndicesGreen2:
    """Tests for GreenFunctionIndices two-body index generators."""

    def test_green2_corr_spin_conservation(self):
        """Test that green2_corr preserves spin conservation."""
        gf = GreenFunctionIndices(
            nsite=1, NsiteUC=1, locspinflag=[0],
            is_kondo=False, is_mvmc=False,
        )
        for t in gf.green2_corr():
            # HPhi order: (site1k, sp1, site1k, sp2, site3, sp3, site3, sp4)
            assert t[1] - t[3] + t[5] - t[7] == 0

    def test_green2_corr_mvmc_reorders(self):
        """Test that mVMC mode produces different index ordering."""
        gf_hphi = GreenFunctionIndices(
            nsite=2, NsiteUC=1, locspinflag=[0, 0],
            is_kondo=False, is_mvmc=False,
        )
        gf_mvmc = GreenFunctionIndices(
            nsite=2, NsiteUC=1, locspinflag=[0, 0],
            is_kondo=False, is_mvmc=True,
        )
        assert len(gf_hphi.green2_corr()) == len(gf_mvmc.green2_corr())
        assert set(gf_hphi.green2_corr()) != set(gf_mvmc.green2_corr())

    def test_green2_corr_kondo_doubles_range(self):
        """Test that Kondo mode doubles UC range in 2-body corr."""
        gf_normal = GreenFunctionIndices(
            nsite=4, NsiteUC=1, locspinflag=[0, 0, 0, 0],
            is_kondo=False, is_mvmc=False,
        )
        gf_kondo = GreenFunctionIndices(
            nsite=4, NsiteUC=1, locspinflag=[0, 0, 0, 0],
            is_kondo=True, is_mvmc=False,
        )
        assert len(gf_kondo.green2_corr()) > len(gf_normal.green2_corr())

    def test_green2_raw_1site_itinerant(self):
        """Test 1-site itinerant raw: 16 combinations (2^4)."""
        gf = GreenFunctionIndices(
            nsite=1, NsiteUC=1, locspinflag=[0], is_kondo=False,
        )
        assert len(gf.green2_raw()) == 16

    def test_green2_raw_local_spin_filtering(self):
        """Test that local-spin filtering applies to (site1,site2) and (site3,site4)."""
        gf = GreenFunctionIndices(
            nsite=2, NsiteUC=1, locspinflag=[1, 1], is_kondo=False,
        )
        for s1, sp1, s2, sp2, s3, sp3, s4, sp4 in gf.green2_raw():
            assert s1 == s2
            assert s3 == s4
