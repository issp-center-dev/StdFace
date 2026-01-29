"""Unit tests for solver_writer module.

Tests for the SolverWriter class hierarchy and the ``get_solver_writer``
factory function.
"""
from __future__ import annotations

import os
import math
import tempfile

import pytest

from writer.solver_writer import (
    SolverWriter,
    HPhiWriter,
    MVMCWriter,
    UHFWriter,
    HWaveWriter,
    get_solver_writer,
)
from stdface_vals import StdIntList

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
#  Tests for get_solver_writer factory
# =====================================================================


class TestGetSolverWriter:
    """Tests for the get_solver_writer factory function."""

    def test_returns_hphi_writer(self):
        """Test that 'HPhi' returns an HPhiWriter."""
        writer = get_solver_writer("HPhi")
        assert isinstance(writer, HPhiWriter)
        assert writer.name == "HPhi"

    def test_returns_mvmc_writer(self):
        """Test that 'mVMC' returns an MVMCWriter."""
        writer = get_solver_writer("mVMC")
        assert isinstance(writer, MVMCWriter)
        assert writer.name == "mVMC"

    def test_returns_uhf_writer(self):
        """Test that 'UHF' returns a UHFWriter."""
        writer = get_solver_writer("UHF")
        assert isinstance(writer, UHFWriter)
        assert writer.name == "UHF"

    def test_returns_hwave_writer(self):
        """Test that 'HWAVE' returns an HWaveWriter."""
        writer = get_solver_writer("HWAVE")
        assert isinstance(writer, HWaveWriter)
        assert writer.name == "HWAVE"

    def test_raises_on_unknown_solver(self):
        """Test that an unknown solver raises ValueError."""
        with pytest.raises(ValueError, match="Unknown solver"):
            get_solver_writer("unknown")


# =====================================================================
#  Tests for SolverWriter subclasses
# =====================================================================


class TestSolverWriterIsAbstract:
    """Tests for the abstract base class."""

    def test_cannot_instantiate_directly(self):
        """Test that SolverWriter cannot be instantiated directly."""
        with pytest.raises(TypeError):
            SolverWriter("test")


class TestHPhiWriter:
    """Tests for the HPhiWriter class."""

    def test_writes_namelist_def(self):
        """Test that HPhiWriter creates namelist.def."""
        StdI = _make_stdi_for_hphi(nsite=4)
        writer = HPhiWriter()

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                writer.write(StdI)
                assert os.path.exists("namelist.def")
                assert os.path.exists("calcmod.def")
                assert os.path.exists("modpara.def")
                assert os.path.exists("locspn.def")
            finally:
                os.chdir(orig)

    def test_writes_trans_def(self):
        """Test that HPhiWriter creates trans.def."""
        StdI = _make_stdi_for_hphi(nsite=4)
        writer = HPhiWriter()

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                writer.write(StdI)
                assert os.path.exists("trans.def")
            finally:
                os.chdir(orig)

    def test_writes_green_files(self):
        """Test that HPhiWriter creates greenone.def and greentwo.def."""
        StdI = _make_stdi_for_hphi(nsite=4)
        writer = HPhiWriter()

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                writer.write(StdI)
                assert os.path.exists("greenone.def")
                assert os.path.exists("greentwo.def")
            finally:
                os.chdir(orig)


class TestUHFWriter:
    """Tests for the UHFWriter class."""

    def test_writes_expected_files(self):
        """Test that UHFWriter creates the expected set of files."""
        StdI = _make_stdi_for_hphi(nsite=4)
        StdI.solver = "UHF"
        StdI.outputmode = "****"
        writer = UHFWriter()

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                writer.write(StdI)
                assert os.path.exists("locspn.def")
                assert os.path.exists("trans.def")
                assert os.path.exists("modpara.def")
                assert os.path.exists("namelist.def")
                assert os.path.exists("greenone.def")
            finally:
                os.chdir(orig)


class TestHWaveWriter:
    """Tests for the HWaveWriter class."""

    def test_uhfr_mode_writes_trans(self):
        """Test that HWAVE in uhfr mode writes trans.def."""
        StdI = _make_stdi_for_hphi(nsite=4)
        StdI.solver = "HWAVE"
        StdI.calcmode = "uhfr"
        StdI.outputmode = "****"
        writer = HWaveWriter()

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                writer.write(StdI)
                assert os.path.exists("trans.def")
                assert os.path.exists("greenone.def")
            finally:
                os.chdir(orig)
