"""Unit tests for stdface_vals module.

Tests for the Python translation of StdFace_vals.h.
"""
from __future__ import annotations

import numpy as np
import pytest

from stdface_vals import StdIntList


class TestStdIntListDefaults:
    """Tests for default initialization of StdIntList."""

    def test_sentinel_value(self):
        """NaN_i should be initialized to 2147483647."""
        s = StdIntList()
        assert s.NaN_i == 2147483647

    def test_pi_default(self):
        """pi should default to 0.0 (set later by ResetVals)."""
        s = StdIntList()
        assert s.pi == 0.0

    def test_string_defaults(self):
        """String fields should default to empty string."""
        s = StdIntList()
        assert s.lattice == ""
        assert s.model == ""
        assert s.outputmode == ""
        assert s.CDataFileHead == ""
        assert s.solver == ""

    def test_int_defaults(self):
        """Integer fields should default to 0."""
        s = StdIntList()
        assert s.W == 0
        assert s.L == 0
        assert s.Height == 0
        assert s.NCell == 0
        assert s.NsiteUC == 0
        assert s.nsite == 0
        assert s.ntrans == 0
        assert s.nintr == 0

    def test_float_defaults(self):
        """Float fields should default to 0.0."""
        s = StdIntList()
        assert s.a == 0.0
        assert s.mu == 0.0
        assert s.U == 0.0
        assert s.h == 0.0
        assert s.Gamma == 0.0

    def test_complex_defaults(self):
        """Complex hopping fields should default to 0+0j."""
        s = StdIntList()
        assert s.t == 0 + 0j
        assert s.tp == 0 + 0j
        assert s.t0 == 0 + 0j
        assert s.t1 == 0 + 0j
        assert s.t2 == 0 + 0j
        assert s.tpp == 0 + 0j

    def test_dynamic_arrays_none(self):
        """Dynamic (pointer) arrays should be initialized to None."""
        s = StdIntList()
        assert s.Cell is None
        assert s.tau is None
        assert s.locspinflag is None
        assert s.transindx is None
        assert s.trans is None
        assert s.intrindx is None
        assert s.intr is None
        assert s.CintraIndx is None
        assert s.Cintra is None
        assert s.CinterIndx is None
        assert s.Cinter is None
        assert s.HundIndx is None
        assert s.Hund is None
        assert s.ExIndx is None
        assert s.Ex is None
        assert s.PLIndx is None
        assert s.PairLift is None
        assert s.PHIndx is None
        assert s.PairHopp is None


class TestStdIntListArrayShapes:
    """Tests for correct array shapes in StdIntList."""

    def test_length_shape(self):
        """length should be shape (3,) float."""
        s = StdIntList()
        assert s.length.shape == (3,)
        assert s.length.dtype == np.float64

    def test_direct_shape(self):
        """direct should be shape (3,3) float."""
        s = StdIntList()
        assert s.direct.shape == (3, 3)
        assert s.direct.dtype == np.float64

    def test_box_shape_and_dtype(self):
        """box should be shape (3,3) int."""
        s = StdIntList()
        assert s.box.shape == (3, 3)
        assert np.issubdtype(s.box.dtype, np.integer)

    def test_rbox_shape_and_dtype(self):
        """rbox should be shape (3,3) int."""
        s = StdIntList()
        assert s.rbox.shape == (3, 3)
        assert np.issubdtype(s.rbox.dtype, np.integer)

    def test_spin_coupling_matrices(self):
        """All J matrices should be shape (3,3) float."""
        s = StdIntList()
        for name in ["J", "Jp", "J0", "J0p", "J0pp",
                      "J1", "J1p", "J1pp",
                      "J2", "J2p", "J2pp", "Jpp", "D"]:
            arr = getattr(s, name)
            assert arr.shape == (3, 3), f"{name} shape mismatch"
            assert arr.dtype == np.float64, f"{name} dtype mismatch"

    def test_phase_shape(self):
        """phase should be shape (3,) float."""
        s = StdIntList()
        assert s.phase.shape == (3,)
        assert s.phase.dtype == np.float64

    def test_expphase_shape_and_dtype(self):
        """ExpPhase should be shape (3,) complex."""
        s = StdIntList()
        assert s.ExpPhase.shape == (3,)
        assert np.issubdtype(s.ExpPhase.dtype, np.complexfloating)

    def test_antiperiod_shape_and_dtype(self):
        """AntiPeriod should be shape (3,) int."""
        s = StdIntList()
        assert s.AntiPeriod.shape == (3,)
        assert np.issubdtype(s.AntiPeriod.dtype, np.integer)

    def test_cutoff_vectors(self):
        """Cutoff vectors should be shape (3,3) float."""
        s = StdIntList()
        for name in ["cutoff_tVec", "cutoff_UVec", "cutoff_JVec"]:
            arr = getattr(s, name)
            assert arr.shape == (3, 3), f"{name} shape mismatch"

    def test_cutoff_r_vectors(self):
        """Cutoff R-vectors should be shape (3,) int."""
        s = StdIntList()
        for name in ["cutoff_tR", "cutoff_UR", "cutoff_JR"]:
            arr = getattr(s, name)
            assert arr.shape == (3,), f"{name} shape mismatch"
            assert np.issubdtype(arr.dtype, np.integer), f"{name} dtype mismatch"


class TestStdIntListIndependence:
    """Tests that each instance has independent arrays."""

    def test_independent_arrays(self):
        """Two instances should not share the same array objects."""
        s1 = StdIntList()
        s2 = StdIntList()
        s1.J[0, 0] = 99.0
        assert s2.J[0, 0] == 0.0

    def test_independent_length(self):
        """Modifying length in one instance should not affect another."""
        s1 = StdIntList()
        s2 = StdIntList()
        s1.length[0] = 5.0
        assert s2.length[0] == 0.0


class TestStdIntListSolverFields:
    """Tests for solver-specific fields."""

    def test_hphi_fields_exist(self):
        """HPhi-specific fields should exist."""
        s = StdIntList()
        assert hasattr(s, "method")
        assert hasattr(s, "Lanczos_max")
        assert hasattr(s, "CalcSpec")
        assert hasattr(s, "SpectrumQ")
        assert s.SpectrumQ.shape == (3,)

    def test_mvmc_fields_exist(self):
        """mVMC-specific fields should exist."""
        s = StdIntList()
        assert hasattr(s, "CParaFileHead")
        assert hasattr(s, "NVMCCalMode")
        assert hasattr(s, "NSROptItrStep")
        assert hasattr(s, "Orb")
        assert s.Orb is None

    def test_uhf_hwave_fields_exist(self):
        """UHF/HWAVE shared fields should exist."""
        s = StdIntList()
        assert hasattr(s, "mix")
        assert hasattr(s, "eps")
        assert hasattr(s, "Iteration_max")

    def test_hwave_only_fields_exist(self):
        """HWAVE-only fields should exist."""
        s = StdIntList()
        assert hasattr(s, "calcmode")
        assert hasattr(s, "fileprefix")
        assert hasattr(s, "export_all")
        assert hasattr(s, "lattice_gp")

    def test_boxsub_shape(self):
        """boxsub and rboxsub should be (3,3) int arrays."""
        s = StdIntList()
        assert s.boxsub.shape == (3, 3)
        assert np.issubdtype(s.boxsub.dtype, np.integer)
        assert s.rboxsub.shape == (3, 3)
        assert np.issubdtype(s.rboxsub.dtype, np.integer)

    def test_lambda_renamed(self):
        """C 'lambda' field should be 'lambda_' in Python."""
        s = StdIntList()
        assert hasattr(s, "lambda_")
        assert s.lambda_ == 0.0
