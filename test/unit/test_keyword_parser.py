"""Unit tests for keyword_parser module.

Tests for the keyword parsing helpers and solver-specific parsers
extracted from ``stdface_main.py``.
"""
from __future__ import annotations

import math
import cmath

import pytest

from stdface.core.keyword_parser import (
    trim_space_quote,
    store_with_check_dup_s,
    store_with_check_dup_sl,
    store_with_check_dup_i,
    store_with_check_dup_d,
    store_with_check_dup_c,
    parse_common_keyword,
    parse_solver_keyword,
    _COMMON_KEYWORDS,
    _HPHI_KEYWORDS,
    _MVMC_KEYWORDS,
    _UHF_KEYWORDS,
    _HWAVE_KEYWORDS,
    _UHF_BASE_KEYWORDS,
    _BOXSUB_KEYWORDS,
    _SOLVER_KEYWORD_TABLES,
    _apply_keyword_table,
    _j_matrix_keywords,
    _grid3x3_keywords,
    NaN_i,
)
from stdface.core.stdface_vals import StdIntList, SolverType

import numpy as np

NaN_d = float("nan")
NaN_c = complex(float("nan"), 0.0)


def _make_stdi(solver: str = "HPhi") -> StdIntList:
    """Create a StdIntList with sentinel values for testing parsers."""
    StdI = StdIntList()
    StdI.solver = solver
    StdI.pi = math.acos(-1.0)

    # Lattice scalars / vectors
    StdI.a = NaN_d
    StdI.length[:] = NaN_d
    StdI.box[:, :] = NaN_i
    StdI.direct[:, :] = NaN_d
    StdI.Gamma = NaN_d
    StdI.Gamma_y = NaN_d
    StdI.h = NaN_d
    StdI.Height = None
    StdI.L = None
    StdI.W = None

    # Isotropic scalar couplings
    StdI.JAll = NaN_d
    StdI.JpAll = NaN_d
    StdI.JppAll = NaN_d
    StdI.J0All = NaN_d
    StdI.J0pAll = NaN_d
    StdI.J0ppAll = NaN_d
    StdI.J1All = NaN_d
    StdI.J1pAll = NaN_d
    StdI.J1ppAll = NaN_d
    StdI.J2All = NaN_d
    StdI.J2pAll = NaN_d
    StdI.J2ppAll = NaN_d

    # 3x3 coupling matrices
    StdI.J[:, :] = NaN_d
    StdI.Jp[:, :] = NaN_d
    StdI.Jpp[:, :] = NaN_d
    StdI.J0[:, :] = NaN_d
    StdI.J0p[:, :] = NaN_d
    StdI.J0pp[:, :] = NaN_d
    StdI.J1[:, :] = NaN_d
    StdI.J1p[:, :] = NaN_d
    StdI.J1pp[:, :] = NaN_d
    StdI.J2[:, :] = NaN_d
    StdI.J2p[:, :] = NaN_d
    StdI.J2pp[:, :] = NaN_d
    StdI.D[:, :] = 0.0
    StdI.D[2, 2] = NaN_d
    StdI.K = NaN_d

    # Chemical potential / spin
    StdI.mu = NaN_d
    StdI.S2 = None

    # Hopping (complex)
    StdI.t = NaN_c
    StdI.tp = NaN_c
    StdI.tpp = NaN_c
    StdI.t0 = NaN_c
    StdI.t0p = NaN_c
    StdI.t0pp = NaN_c
    StdI.t1 = NaN_c
    StdI.t1p = NaN_c
    StdI.t1pp = NaN_c
    StdI.t2 = NaN_c
    StdI.t2p = NaN_c
    StdI.t2pp = NaN_c

    # Coulomb
    StdI.U = NaN_d
    StdI.V = NaN_d
    StdI.Vp = NaN_d
    StdI.Vpp = NaN_d
    StdI.V0 = NaN_d
    StdI.V0p = NaN_d
    StdI.V0pp = NaN_d
    StdI.V1 = NaN_d
    StdI.V1p = NaN_d
    StdI.V1pp = NaN_d
    StdI.V2 = NaN_d
    StdI.V2p = NaN_d
    StdI.V2pp = NaN_d

    # Phase / boundary
    StdI.phase[:] = NaN_d

    # Strings
    StdI.ncond = None
    StdI.Sz2 = None
    StdI.model = None
    StdI.lattice = None
    StdI.outputmode = None
    StdI.CDataFileHead = None
    StdI.double_counting_mode = None

    # Wannier90 cutoffs
    StdI.cutoff_t = NaN_d
    StdI.cutoff_u = NaN_d
    StdI.cutoff_j = NaN_d
    StdI.cutoff_length_t = NaN_d
    StdI.cutoff_length_U = NaN_d
    StdI.cutoff_length_J = NaN_d
    StdI.lambda_ = NaN_d
    StdI.lambda_U = NaN_d
    StdI.lambda_J = NaN_d
    StdI.alpha = NaN_d
    StdI.cutoff_tR[:] = NaN_i
    StdI.cutoff_UR[:] = NaN_i
    StdI.cutoff_JR[:] = NaN_i
    StdI.cutoff_tVec[:, :] = NaN_d
    StdI.cutoff_UVec[:, :] = NaN_d
    StdI.cutoff_JVec[:, :] = NaN_d

    # HPhi-specific fields
    if solver == "HPhi":
        StdI.LargeValue = NaN_d
        StdI.OmegaMax = NaN_d
        StdI.OmegaMin = NaN_d
        StdI.OmegaOrg = NaN_d
        StdI.OmegaIm = NaN_d
        StdI.Nomega = None
        StdI.SpectrumQ[:] = NaN_d
        StdI.method = None
        StdI.Restart = None
        StdI.EigenVecIO = None
        StdI.InitialVecType = None
        StdI.HamIO = None
        StdI.CalcSpec = None
        StdI.SpectrumType = None
        StdI.OutputExVec = None
        StdI.FlgTemp = 1
        StdI.Lanczos_max = None
        StdI.initial_iv = None
        StdI.nvec = None
        StdI.exct = None
        StdI.LanczosEps = None
        StdI.LanczosTarget = None
        StdI.NumAve = None
        StdI.ExpecInterval = None
        StdI.dt = NaN_d
        StdI.tdump = NaN_d
        StdI.tshift = NaN_d
        StdI.freq = NaN_d
        StdI.Uquench = NaN_d
        StdI.VecPot[:] = NaN_d
        StdI.PumpType = None
        StdI.ExpandCoef = None
        StdI.NGPU = None
        StdI.Scalapack = None

    elif solver == "mVMC":
        StdI.CParaFileHead = None
        StdI.NVMCCalMode = None
        StdI.NLanczosMode = None
        StdI.NDataIdxStart = None
        StdI.NDataQtySmp = None
        StdI.NSPGaussLeg = None
        StdI.NSPStot = None
        StdI.NMPTrans = None
        StdI.NSROptItrStep = None
        StdI.NSROptItrSmp = None
        StdI.DSROptRedCut = NaN_d
        StdI.DSROptStaDel = NaN_d
        StdI.DSROptStepDt = NaN_d
        StdI.NVMCWarmUp = None
        StdI.NVMCInterval = None
        StdI.NVMCSample = None
        StdI.NExUpdatePath = None
        StdI.RndSeed = None
        StdI.NSplitSize = None
        StdI.NStore = None
        StdI.NSRCG = None
        StdI.ComplexType = None
        StdI.boxsub[:, :] = NaN_i
        StdI.Hsub = None
        StdI.Lsub = None
        StdI.Wsub = None

    elif solver == "UHF":
        StdI.NMPTrans = None
        StdI.RndSeed = None
        StdI.mix = NaN_d
        StdI.eps = None
        StdI.eps_slater = None
        StdI.Iteration_max = None
        StdI.boxsub[:, :] = NaN_i
        StdI.Hsub = None
        StdI.Lsub = None
        StdI.Wsub = None

    elif solver == "HWAVE":
        StdI.NMPTrans = None
        StdI.RndSeed = None
        StdI.mix = NaN_d
        StdI.eps = None
        StdI.eps_slater = None
        StdI.Iteration_max = None
        StdI.boxsub[:, :] = NaN_i
        StdI.Hsub = None
        StdI.Lsub = None
        StdI.Wsub = None
        StdI.calcmode = None
        StdI.fileprefix = None
        StdI.export_all = None
        StdI.lattice_gp = None

    return StdI


# =====================================================================
#  Tests for trim_space_quote
# =====================================================================


class TestTrimSpaceQuote:
    """Tests for the trim_space_quote helper."""

    def test_removes_spaces(self):
        """Test that spaces are stripped."""
        assert trim_space_quote("a b c") == "abc"

    def test_removes_colons_semicolons(self):
        """Test that colons and semicolons are stripped."""
        assert trim_space_quote("key:value;end") == "keyvalueend"

    def test_removes_quotes(self):
        """Test that double quotes are stripped."""
        assert trim_space_quote('"hello"') == "hello"

    def test_removes_backslashes(self):
        """Test that backslashes are stripped."""
        assert trim_space_quote("a\\b") == "ab"

    def test_preserves_equals(self):
        """Test that '=' is NOT stripped (needed for key=value parsing)."""
        assert "=" in trim_space_quote("key=value")


# =====================================================================
#  Tests for store_with_check_dup_s
# =====================================================================


class TestStoreWithCheckDupS:
    """Tests for string duplicate-check storage."""

    def test_stores_new_value(self):
        """Test that a new value is stored when current is sentinel."""
        result = store_with_check_dup_s("model", "Hubbard", None)
        assert result == "Hubbard"

    def test_exits_on_duplicate(self):
        """Test that duplicate assignment raises ValueError."""
        with pytest.raises(ValueError):
            store_with_check_dup_s("model", "spin", "hubbard")


# =====================================================================
#  Tests for store_with_check_dup_sl
# =====================================================================


class TestStoreWithCheckDupSl:
    """Tests for string-lowercase duplicate-check storage."""

    def test_stores_and_lowercases(self):
        """Test that value is stored and lowercased."""
        result = store_with_check_dup_sl("lattice", "Chain", None)
        assert result == "chain"

    def test_truncates_long_value(self):
        """Test truncation to maxlen."""
        result = store_with_check_dup_sl("key", "abcdefghij", None, maxlen=5)
        assert result == "abcde"

    def test_exits_on_duplicate(self):
        """Test that duplicate assignment raises ValueError."""
        with pytest.raises(ValueError):
            store_with_check_dup_sl("lattice", "chain", "square")


# =====================================================================
#  Tests for store_with_check_dup_i
# =====================================================================


class TestStoreWithCheckDupI:
    """Tests for integer duplicate-check storage."""

    def test_stores_integer(self):
        """Test that a string is parsed to int."""
        result = store_with_check_dup_i("L", "4", None)
        assert result == 4

    def test_truncates_float_string(self):
        """Test that '2.0' is parsed as 2 (C sscanf behavior)."""
        result = store_with_check_dup_i("W", "2.0", None)
        assert result == 2

    def test_exits_on_duplicate(self):
        """Test that duplicate assignment raises ValueError."""
        with pytest.raises(ValueError):
            store_with_check_dup_i("L", "8", 4)


# =====================================================================
#  Tests for store_with_check_dup_d
# =====================================================================


class TestStoreWithCheckDupD:
    """Tests for float duplicate-check storage."""

    def test_stores_float(self):
        """Test that a string is parsed to float."""
        result = store_with_check_dup_d("U", "4.5", NaN_d)
        assert result == 4.5

    def test_exits_on_duplicate(self):
        """Test that duplicate assignment raises ValueError."""
        with pytest.raises(ValueError):
            store_with_check_dup_d("U", "5.0", 4.5)


# =====================================================================
#  Tests for store_with_check_dup_c
# =====================================================================


class TestStoreWithCheckDupC:
    """Tests for complex duplicate-check storage."""

    def test_stores_real_imag(self):
        """Test 'real,imag' format."""
        result = store_with_check_dup_c("t", "1.5,2.5", NaN_c)
        assert result == complex(1.5, 2.5)

    def test_stores_real_only(self):
        """Test 'real' format (imaginary defaults to 0)."""
        result = store_with_check_dup_c("t", "3.0", NaN_c)
        assert result == complex(3.0, 0.0)

    def test_stores_imag_only(self):
        """Test ',imag' format (real defaults to 0)."""
        result = store_with_check_dup_c("t", ",2.0", NaN_c)
        assert result == complex(0.0, 2.0)

    def test_invalid_real_part_yields_zero_real(self):
        """Non-numeric real token is treated as 0.0 via _safe_float."""
        result = store_with_check_dup_c("t", "bad,1.0", NaN_c)
        assert result == complex(0.0, 1.0)

    def test_invalid_imag_part_yields_zero_imag(self):
        """Non-numeric imaginary token is treated as 0.0 via _safe_float."""
        result = store_with_check_dup_c("t", "1.0,bad", NaN_c)
        assert result == complex(1.0, 0.0)

    def test_exits_on_duplicate(self):
        """Test that duplicate assignment raises ValueError."""
        with pytest.raises(ValueError):
            store_with_check_dup_c("t", "1.0", complex(2.0, 0.0))


# =====================================================================
#  Tests for parse_common_keyword
# =====================================================================


class TestParseCommonKeyword:
    """Tests for the common keyword parser."""

    def test_recognises_model(self):
        """Test that 'model' keyword is recognised."""
        StdI = _make_stdi()
        result = parse_common_keyword("model", "Hubbard", StdI)
        assert result is True
        assert StdI.model == "hubbard"

    def test_recognises_lattice(self):
        """Test that 'lattice' keyword is recognised."""
        StdI = _make_stdi()
        result = parse_common_keyword("lattice", "Chain", StdI)
        assert result is True
        assert StdI.lattice == "chain"

    def test_recognises_l(self):
        """Test that 'l' (length) keyword is parsed as integer."""
        StdI = _make_stdi()
        result = parse_common_keyword("l", "8", StdI)
        assert result is True
        assert StdI.L == 8

    def test_recognises_w(self):
        """Test that 'w' (width) keyword is parsed as integer."""
        StdI = _make_stdi()
        result = parse_common_keyword("w", "4", StdI)
        assert result is True
        assert StdI.W == 4

    def test_recognises_u(self):
        """Test that 'u' (Coulomb U) keyword is parsed as float."""
        StdI = _make_stdi()
        result = parse_common_keyword("u", "4.5", StdI)
        assert result is True
        assert StdI.U == 4.5

    def test_recognises_t(self):
        """Test that 't' (hopping) keyword is parsed as complex."""
        StdI = _make_stdi()
        result = parse_common_keyword("t", "1.0,0.5", StdI)
        assert result is True
        assert StdI.t == complex(1.0, 0.5)

    def test_recognises_j(self):
        """Test that 'j' (exchange coupling) keyword is parsed."""
        StdI = _make_stdi()
        result = parse_common_keyword("j", "1.0", StdI)
        assert result is True
        assert StdI.JAll == 1.0

    def test_recognises_2sz(self):
        """Test that '2sz' keyword is parsed."""
        StdI = _make_stdi()
        result = parse_common_keyword("2sz", "0", StdI)
        assert result is True
        assert StdI.Sz2 == 0

    def test_recognises_phase0(self):
        """Test that 'phase0' keyword is parsed."""
        StdI = _make_stdi()
        result = parse_common_keyword("phase0", "180.0", StdI)
        assert result is True
        assert StdI.phase[0] == 180.0

    def test_unrecognised_returns_false(self):
        """Test that an unknown keyword returns False."""
        StdI = _make_stdi()
        result = parse_common_keyword("nosuchkeyword", "42", StdI)
        assert result is False

    def test_box_a0w(self):
        """Test that 'a0w' supercell keyword is parsed."""
        StdI = _make_stdi()
        result = parse_common_keyword("a0w", "2", StdI)
        assert result is True
        assert int(StdI.box[0, 0]) == 2

    def test_j0_matrix_element(self):
        """Test that 'j0x' sets J0[0,0]."""
        StdI = _make_stdi()
        result = parse_common_keyword("j0x", "0.5", StdI)
        assert result is True
        assert float(StdI.J0[0, 0]) == 0.5


# =====================================================================
#  Tests for parse_solver_keyword
# =====================================================================


class TestParseSolverKeyword:
    """Tests for the solver keyword dispatch function."""

    def test_dispatches_to_hphi(self):
        """Test that HPhi keywords are dispatched correctly."""
        StdI = _make_stdi("HPhi")
        result = parse_solver_keyword("method", "lanczos", StdI, "HPhi")
        assert result is True
        assert StdI.method == "lanczos"

    def test_dispatches_to_mvmc(self):
        """Test that mVMC keywords are dispatched correctly."""
        StdI = _make_stdi("mVMC")
        result = parse_solver_keyword("nvmcsample", "1000", StdI, "mVMC")
        assert result is True
        assert StdI.NVMCSample == 1000

    def test_dispatches_to_uhf(self):
        """Test that UHF keywords are dispatched correctly."""
        StdI = _make_stdi("UHF")
        result = parse_solver_keyword("iteration_max", "100", StdI, "UHF")
        assert result is True
        assert StdI.Iteration_max == 100

    def test_dispatches_to_hwave(self):
        """Test that HWAVE keywords are dispatched correctly."""
        StdI = _make_stdi("HWAVE")
        result = parse_solver_keyword("calcmode", "uhfr", StdI, "HWAVE")
        assert result is True
        assert StdI.calcmode == "uhfr"

    def test_unknown_solver_returns_false(self):
        """Test that an unknown solver returns False."""
        StdI = _make_stdi("HPhi")
        result = parse_solver_keyword("method", "lanczos", StdI, "unknown")
        assert result is False


# =====================================================================
#  Tests for parse_hphi_keyword
# =====================================================================


class TestParseHPhiKeyword:
    """Tests for the HPhi-specific keyword parser via parse_solver_keyword."""

    def test_method(self):
        """Test 'method' keyword."""
        StdI = _make_stdi("HPhi")
        assert parse_solver_keyword("method", "Lanczos", StdI, "HPhi") is True
        assert StdI.method == "lanczos"

    def test_exct(self):
        """Test 'exct' keyword."""
        StdI = _make_stdi("HPhi")
        assert parse_solver_keyword("exct", "2", StdI, "HPhi") is True
        assert StdI.exct == 2

    def test_dt(self):
        """Test 'dt' keyword."""
        StdI = _make_stdi("HPhi")
        assert parse_solver_keyword("dt", "0.01", StdI, "HPhi") is True
        assert StdI.dt == 0.01

    def test_unrecognised(self):
        """Test that unknown HPhi keyword returns False."""
        StdI = _make_stdi("HPhi")
        assert parse_solver_keyword("notahphikey", "1", StdI, "HPhi") is False


# =====================================================================
#  Tests for parse_mvmc_keyword
# =====================================================================


class TestParseMVMCKeyword:
    """Tests for the mVMC-specific keyword parser via parse_solver_keyword."""

    def test_nvmcsample(self):
        """Test 'nvmcsample' keyword."""
        StdI = _make_stdi("mVMC")
        assert parse_solver_keyword("nvmcsample", "5000", StdI, "mVMC") is True
        assert StdI.NVMCSample == 5000

    def test_complextype(self):
        """Test 'complextype' keyword."""
        StdI = _make_stdi("mVMC")
        assert parse_solver_keyword("complextype", "1", StdI, "mVMC") is True
        assert StdI.ComplexType == 1

    def test_boxsub(self):
        """Test 'a0wsub' keyword sets boxsub."""
        StdI = _make_stdi("mVMC")
        assert parse_solver_keyword("a0wsub", "3", StdI, "mVMC") is True
        assert int(StdI.boxsub[0, 0]) == 3

    def test_unrecognised(self):
        """Test that unknown mVMC keyword returns False."""
        StdI = _make_stdi("mVMC")
        assert parse_solver_keyword("notamvmckey", "1", StdI, "mVMC") is False


# =====================================================================
#  Tests for parse_uhf_keyword
# =====================================================================


class TestParseUHFKeyword:
    """Tests for the UHF-specific keyword parser via parse_solver_keyword."""

    def test_iteration_max(self):
        """Test 'iteration_max' keyword."""
        StdI = _make_stdi("UHF")
        assert parse_solver_keyword("iteration_max", "200", StdI, "UHF") is True
        assert StdI.Iteration_max == 200

    def test_mix(self):
        """Test 'mix' keyword."""
        StdI = _make_stdi("UHF")
        assert parse_solver_keyword("mix", "0.5", StdI, "UHF") is True
        assert StdI.mix == 0.5

    def test_unrecognised(self):
        """Test that unknown UHF keyword returns False."""
        StdI = _make_stdi("UHF")
        assert parse_solver_keyword("notauhfkey", "1", StdI, "UHF") is False


# =====================================================================
#  Tests for parse_hwave_keyword
# =====================================================================


class TestParseHWAVEKeyword:
    """Tests for the HWAVE-specific keyword parser via parse_solver_keyword."""

    def test_calcmode(self):
        """Test 'calcmode' keyword."""
        StdI = _make_stdi("HWAVE")
        assert parse_solver_keyword("calcmode", "UHFR", StdI, "HWAVE") is True
        assert StdI.calcmode == "uhfr"

    def test_fileprefix(self):
        """Test 'fileprefix' keyword."""
        StdI = _make_stdi("HWAVE")
        assert parse_solver_keyword("fileprefix", "output", StdI, "HWAVE") is True
        assert StdI.fileprefix == "output"

    def test_export_all(self):
        """Test 'exportall' keyword."""
        StdI = _make_stdi("HWAVE")
        assert parse_solver_keyword("exportall", "1", StdI, "HWAVE") is True
        assert StdI.export_all == 1

    def test_unrecognised(self):
        """Test that unknown HWAVE keyword returns False."""
        StdI = _make_stdi("HWAVE")
        assert parse_solver_keyword("notahwavekey", "1", StdI, "HWAVE") is False


# =====================================================================
#  Tests for keyword dispatch tables
# =====================================================================


class TestKeywordTableStructure:
    """Tests for the solver keyword dispatch table structure."""

    def test_solver_table_has_four_entries(self):
        """Test that _SOLVER_KEYWORD_TABLES has exactly 4 solver entries."""
        assert len(_SOLVER_KEYWORD_TABLES) == 4

    def test_all_solvers_present(self):
        """Test that all solver types are present in the table."""
        for solver in (SolverType.HPhi, SolverType.mVMC,
                       SolverType.UHF, SolverType.HWAVE):
            assert solver in _SOLVER_KEYWORD_TABLES

    def test_hphi_table_is_dict(self):
        """Test that HPhi keyword table is a dict."""
        assert isinstance(_HPHI_KEYWORDS, dict)

    def test_all_entries_are_tuples(self):
        """Test that all table entries are tuples of length 2 or 4."""
        for name, table in _SOLVER_KEYWORD_TABLES.items():
            for kw, entry in table.items():
                assert isinstance(entry, tuple), f"{name}/{kw}: not tuple"
                assert len(entry) in (2, 4), f"{name}/{kw}: len={len(entry)}"

    def test_all_keys_are_lowercase(self):
        """Test that all keyword keys are lowercase strings."""
        for name, table in _SOLVER_KEYWORD_TABLES.items():
            for kw in table:
                assert kw == kw.lower(), f"{name}: key {kw!r} not lowercase"
                assert isinstance(kw, str), f"{name}: key {kw!r} not str"


class TestUHFHWAVEKeywordSharing:
    """Tests for UHF/HWAVE shared keyword base."""

    def test_uhf_is_base(self):
        """Test that UHF keywords is exactly _UHF_BASE_KEYWORDS."""
        assert _UHF_KEYWORDS is _UHF_BASE_KEYWORDS

    def test_hwave_contains_all_uhf_keys(self):
        """Test that HWAVE keywords contain all UHF base keys."""
        for kw in _UHF_BASE_KEYWORDS:
            assert kw in _HWAVE_KEYWORDS, f"Missing in HWAVE: {kw}"

    def test_hwave_has_extra_keys(self):
        """Test that HWAVE has calcmode, fileprefix, exportall, lattice_gp."""
        assert "calcmode" in _HWAVE_KEYWORDS
        assert "fileprefix" in _HWAVE_KEYWORDS
        assert "exportall" in _HWAVE_KEYWORDS
        assert "lattice_gp" in _HWAVE_KEYWORDS

    def test_uhf_lacks_hwave_extras(self):
        """Test that UHF does not have HWAVE-specific keywords."""
        assert "calcmode" not in _UHF_KEYWORDS
        assert "fileprefix" not in _UHF_KEYWORDS
        assert "exportall" not in _UHF_KEYWORDS

    def test_boxsub_shared_in_mvmc(self):
        """Test that mVMC contains all boxsub keywords."""
        for kw in _BOXSUB_KEYWORDS:
            assert kw in _MVMC_KEYWORDS, f"Missing in mVMC: {kw}"

    def test_boxsub_shared_in_uhf(self):
        """Test that UHF contains all boxsub keywords."""
        for kw in _BOXSUB_KEYWORDS:
            assert kw in _UHF_KEYWORDS, f"Missing in UHF: {kw}"


class TestApplyKeywordTable:
    """Tests for the generic _apply_keyword_table function."""

    def test_scalar_int_field(self):
        """Test applying a scalar integer keyword."""
        StdI = _make_stdi("HPhi")
        result = _apply_keyword_table(_HPHI_KEYWORDS, "exct", "5", StdI)
        assert result is True
        assert StdI.exct == 5

    def test_scalar_float_field(self):
        """Test applying a scalar float keyword."""
        StdI = _make_stdi("HPhi")
        result = _apply_keyword_table(_HPHI_KEYWORDS, "dt", "0.05", StdI)
        assert result is True
        assert StdI.dt == 0.05

    def test_scalar_string_field(self):
        """Test applying a scalar string keyword."""
        StdI = _make_stdi("HPhi")
        result = _apply_keyword_table(_HPHI_KEYWORDS, "method", "CG", StdI)
        assert result is True
        assert StdI.method == "cg"

    def test_array_element_field(self):
        """Test applying an array-element keyword (boxsub)."""
        StdI = _make_stdi("mVMC")
        result = _apply_keyword_table(_MVMC_KEYWORDS, "a0wsub", "7", StdI)
        assert result is True
        assert int(StdI.boxsub[0, 0]) == 7

    def test_array_element_spectrumq(self):
        """Test applying SpectrumQ array element keyword."""
        StdI = _make_stdi("HPhi")
        result = _apply_keyword_table(_HPHI_KEYWORDS, "spectrumqw", "1.5", StdI)
        assert result is True
        assert float(StdI.SpectrumQ[0]) == 1.5

    def test_array_element_vecpot(self):
        """Test applying VecPot array element keyword."""
        StdI = _make_stdi("HPhi")
        result = _apply_keyword_table(_HPHI_KEYWORDS, "vecpoth", "0.3", StdI)
        assert result is True
        assert float(StdI.VecPot[2]) == 0.3

    def test_unknown_keyword_returns_false(self):
        """Test that unknown keyword returns False."""
        StdI = _make_stdi("HPhi")
        result = _apply_keyword_table(_HPHI_KEYWORDS, "nonexistent", "1", StdI)
        assert result is False

    def test_empty_table(self):
        """Test that empty table returns False for any keyword."""
        StdI = _make_stdi("HPhi")
        result = _apply_keyword_table({}, "exct", "1", StdI)
        assert result is False


# =====================================================================
#  Tests for _COMMON_KEYWORDS table
# =====================================================================


class TestCommonKeywordsTable:
    """Tests for the _COMMON_KEYWORDS dispatch table structure."""

    def test_is_dict(self):
        """Test that _COMMON_KEYWORDS is a dict."""
        assert isinstance(_COMMON_KEYWORDS, dict)

    def test_all_keys_are_lowercase_strings(self):
        """Test that all keyword keys are lowercase strings."""
        for kw in _COMMON_KEYWORDS:
            assert isinstance(kw, str), f"key {kw!r} not str"
            assert kw == kw.lower(), f"key {kw!r} not lowercase"

    def test_all_entries_are_tuples(self):
        """Test that all entries are tuples of length 2 or 4."""
        for kw, entry in _COMMON_KEYWORDS.items():
            assert isinstance(entry, tuple), f"{kw}: not tuple"
            assert len(entry) in (2, 4), f"{kw}: len={len(entry)}"

    def test_contains_model_and_lattice(self):
        """Test that model and lattice keywords are present."""
        assert "model" in _COMMON_KEYWORDS
        assert "lattice" in _COMMON_KEYWORDS

    def test_contains_nelec_alias(self):
        """Test that 'nelec' keyword maps to 'ncond' field."""
        entry = _COMMON_KEYWORDS["nelec"]
        assert entry[1] == "ncond"

    def test_contains_all_j_families(self):
        """Test that all 12 J-family prefixes are present."""
        prefixes = ["j", "j0", "j0'", "j0''", "j1", "j1'", "j1''",
                    "j2", "j2'", "j2''", "j'", "j''"]
        for prefix in prefixes:
            assert prefix in _COMMON_KEYWORDS, f"Missing J prefix: {prefix}"

    def test_contains_all_t_hopping(self):
        """Test that all 12 hopping keywords are present."""
        keywords = ["t", "t0", "t0'", "t0''", "t1", "t1'", "t1''",
                    "t2", "t2'", "t2''", "t'", "t''"]
        for kw in keywords:
            assert kw in _COMMON_KEYWORDS, f"Missing hopping: {kw}"

    def test_contains_all_v_coulomb(self):
        """Test that all 13 V-Coulomb keywords are present."""
        keywords = ["u", "v", "v0", "v0'", "v0''", "v1", "v1'", "v1''",
                    "v2", "v2'", "v2''", "v'", "v''"]
        for kw in keywords:
            assert kw in _COMMON_KEYWORDS, f"Missing Coulomb: {kw}"

    def test_contains_phase_keywords(self):
        """Test that phase0, phase1, phase2 are present."""
        for i in range(3):
            assert f"phase{i}" in _COMMON_KEYWORDS

    def test_contains_box_keywords(self):
        """Test that all 9 box (supercell) keywords are present."""
        for a in ("a0", "a1", "a2"):
            for c in ("w", "l", "h"):
                assert f"{a}{c}" in _COMMON_KEYWORDS, f"Missing box: {a}{c}"

    def test_contains_cutoff_keywords(self):
        """Test that cutoff_t, cutoff_u, cutoff_j scalar keys are present."""
        assert "cutoff_t" in _COMMON_KEYWORDS
        assert "cutoff_u" in _COMMON_KEYWORDS
        assert "cutoff_j" in _COMMON_KEYWORDS


# =====================================================================
#  Tests for _j_matrix_keywords generator
# =====================================================================


class TestJMatrixKeywords:
    """Tests for the _j_matrix_keywords generator function."""

    def test_returns_10_entries(self):
        """Test that generator produces exactly 10 entries per family."""
        d = _j_matrix_keywords("j0", "J0All", "J0")
        assert len(d) == 10

    def test_has_scalar_entry(self):
        """Test that the isotropic scalar keyword is present."""
        d = _j_matrix_keywords("j0", "J0All", "J0")
        assert "j0" in d
        assert len(d["j0"]) == 2
        assert d["j0"][1] == "J0All"

    def test_has_all_9_components(self):
        """Test that all 9 anisotropic component keywords are present."""
        d = _j_matrix_keywords("j0", "J0All", "J0")
        for suffix in ("x", "xy", "xz", "y", "yx", "yz", "z", "zx", "zy"):
            assert f"j0{suffix}" in d, f"Missing: j0{suffix}"

    def test_component_indices(self):
        """Test that component indices map correctly."""
        d = _j_matrix_keywords("j", "JAll", "J")
        # jx -> (0,0), jy -> (1,1), jz -> (2,2), jxy -> (0,1) etc.
        assert d["jx"][2] == (0, 0)
        assert d["jy"][2] == (1, 1)
        assert d["jz"][2] == (2, 2)
        assert d["jxy"][2] == (0, 1)
        assert d["jxz"][2] == (0, 2)
        assert d["jyx"][2] == (1, 0)
        assert d["jyz"][2] == (1, 2)
        assert d["jzx"][2] == (2, 0)
        assert d["jzy"][2] == (2, 1)

    def test_component_entries_are_4_tuples(self):
        """Test that component entries have length 4."""
        d = _j_matrix_keywords("j1", "J1All", "J1")
        for key in d:
            if key != "j1":
                assert len(d[key]) == 4, f"{key}: len={len(d[key])}"

    def test_prime_prefix_works(self):
        """Test that primed prefixes generate correct keys."""
        d = _j_matrix_keywords("j'", "JpAll", "Jp")
        assert "j'" in d
        assert "j'x" in d
        assert "j'xy" in d
        assert len(d) == 10


# =====================================================================
#  Tests for _cutoff_vec_keywords generator
# =====================================================================


class TestCutoffVecKeywords:
    """Tests for the _grid3x3_keywords generator (cutoff vec variant)."""

    def test_returns_9_entries(self):
        """Test that generator produces exactly 9 entries."""
        d = _grid3x3_keywords("cutoff_j_{a}{c}", "cutoff_JVec", store_with_check_dup_d, float)
        assert len(d) == 9

    def test_has_all_9_keys(self):
        """Test that all a{0,1,2}{w,l,h} combos are present."""
        d = _grid3x3_keywords("cutoff_t_{a}{c}", "cutoff_tVec", store_with_check_dup_d, float)
        for row, aname in enumerate(("a0", "a1", "a2")):
            for col, comp in enumerate(("w", "l", "h")):
                key = f"cutoff_t_{aname}{comp}"
                assert key in d, f"Missing: {key}"

    def test_index_mapping(self):
        """Test that indices map to correct (row, col)."""
        d = _grid3x3_keywords("cutoff_u_{a}{c}", "cutoff_UVec", store_with_check_dup_d, float)
        assert d["cutoff_u_a0w"][2] == (0, 0)
        assert d["cutoff_u_a1l"][2] == (1, 1)
        assert d["cutoff_u_a2h"][2] == (2, 2)
        assert d["cutoff_u_a0l"][2] == (0, 1)
        assert d["cutoff_u_a2w"][2] == (2, 0)

    def test_entries_are_4_tuples(self):
        """Test that all entries have length 4."""
        d = _grid3x3_keywords("cutoff_j_{a}{c}", "cutoff_JVec", store_with_check_dup_d, float)
        for key, entry in d.items():
            assert len(entry) == 4, f"{key}: len={len(entry)}"
            assert entry[1] == "cutoff_JVec"


# =====================================================================
#  Tests for _box_keywords generator
# =====================================================================


class TestBoxKeywords:
    """Tests for the _grid3x3_keywords generator (box variant)."""

    def _box(self):
        return _grid3x3_keywords("{a}{c}", "box", store_with_check_dup_i, int)

    def test_returns_9_entries(self):
        """Test that generator produces exactly 9 entries."""
        assert len(self._box()) == 9

    def test_has_all_9_keys(self):
        """Test that all a{0,1,2}{w,l,h} combos are present."""
        d = self._box()
        for a in ("a0", "a1", "a2"):
            for c in ("w", "l", "h"):
                assert f"{a}{c}" in d

    def test_index_mapping(self):
        """Test that indices map to correct (row, col)."""
        d = self._box()
        assert d["a0w"][2] == (0, 0)
        assert d["a1l"][2] == (1, 1)
        assert d["a2h"][2] == (2, 2)

    def test_entries_use_int_cast(self):
        """Test that box entries use int as cast function."""
        for key, entry in self._box().items():
            assert entry[3] is int, f"{key}: cast is not int"

    def test_entries_target_box_field(self):
        """Test that all entries target the 'box' array."""
        for key, entry in self._box().items():
            assert entry[1] == "box", f"{key}: field is {entry[1]}"
