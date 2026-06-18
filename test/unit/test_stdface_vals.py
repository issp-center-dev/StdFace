"""Unit tests for stdface_vals module.

Tests for the Python translation of StdFace_vals.h.
"""
from __future__ import annotations

import numpy as np
import pytest

from stdface.core.stdface_vals import (
    StdIntList, ModelType, SolverType, MethodType,
    NaN_i, NaN_d, NaN_c, UNSET_STRING,
    AMPLITUDE_EPS, ZERO_BODY_EPS,
    is_unset_or_trivial_d,
)


class TestIsUnsetOrTrivialD:
    """Tests for is_unset_or_trivial_d (helper for SolverPlugin.validate)."""

    def test_nan_is_unset(self):
        assert is_unset_or_trivial_d(NaN_d) is True

    def test_zero_is_trivial(self):
        assert is_unset_or_trivial_d(0.0) is True

    def test_nonzero_value(self):
        assert is_unset_or_trivial_d(0.5) is False

    def test_custom_trivial(self):
        assert is_unset_or_trivial_d(1.0, trivial=1.0) is True
        assert is_unset_or_trivial_d(2.0, trivial=1.0) is False


class TestStdIntListDefaults:
    """Tests for default initialization of StdIntList."""

    def test_string_defaults(self):
        """String fields should default to UNSET_STRING or empty string."""
        s = StdIntList()
        assert s.lattice == UNSET_STRING
        assert s.model == UNSET_STRING
        assert s.outputmode == UNSET_STRING
        assert s.CDataFileHead == UNSET_STRING
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
        assert len(s.trans_list) == 0
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
        assert s.trans_list == []
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


# ===================================================================
#  ModelType enum
# ===================================================================


class TestModelType:
    """Tests for the ModelType enum."""

    def test_has_three_members(self):
        """ModelType should have exactly 3 members."""
        assert len(ModelType) == 3

    def test_values(self):
        """ModelType members should have the expected string values."""
        assert ModelType.SPIN == "spin"
        assert ModelType.HUBBARD == "hubbard"
        assert ModelType.KONDO == "kondo"

    def test_is_str_subclass(self):
        """ModelType members should be instances of str."""
        assert isinstance(ModelType.SPIN, str)
        assert isinstance(ModelType.HUBBARD, str)
        assert isinstance(ModelType.KONDO, str)

    def test_string_equality(self):
        """ModelType members should compare equal to their string values."""
        assert ModelType.SPIN == "spin"
        assert ModelType.HUBBARD == "hubbard"
        assert ModelType.KONDO == "kondo"
        assert "spin" == ModelType.SPIN
        assert "kondo" == ModelType.KONDO

    def test_string_inequality(self):
        """ModelType members should not equal other strings."""
        assert ModelType.SPIN != "hubbard"
        assert ModelType.HUBBARD != "kondo"
        assert ModelType.KONDO != "spin"

    def test_construction_from_string(self):
        """ModelType can be constructed from its string value."""
        assert ModelType("spin") is ModelType.SPIN
        assert ModelType("hubbard") is ModelType.HUBBARD
        assert ModelType("kondo") is ModelType.KONDO

    def test_invalid_string_raises(self):
        """ModelType raises ValueError for invalid strings."""
        with pytest.raises(ValueError):
            ModelType("invalid")
        with pytest.raises(ValueError):
            ModelType("****")

    def test_can_be_used_as_dict_key(self):
        """ModelType members work as dict keys, interchangeable with strings."""
        d = {ModelType.SPIN: 1, ModelType.HUBBARD: 2}
        assert d["spin"] == 1
        assert d["hubbard"] == 2

    def test_model_field_assignment(self):
        """StdIntList.model can be set to a ModelType value."""
        s = StdIntList()
        s.model = ModelType.HUBBARD
        assert s.model == "hubbard"
        assert s.model == ModelType.HUBBARD


# ===================================================================
#  SolverType enum
# ===================================================================


class TestSolverType:
    """Tests for the SolverType enum."""

    def test_has_six_members(self):
        """SolverType should have exactly 6 members (incl. UHFR/UHFK)."""
        assert len(SolverType) == 6

    def test_values(self):
        """SolverType members should have the expected string values."""
        assert SolverType.HPhi == "HPhi"
        assert SolverType.mVMC == "mVMC"
        assert SolverType.UHF == "UHF"
        assert SolverType.HWAVE == "HWAVE"
        assert SolverType.UHFR == "UHFR"
        assert SolverType.UHFK == "UHFK"

    def test_is_str_subclass(self):
        """SolverType members should be instances of str."""
        assert isinstance(SolverType.HPhi, str)
        assert isinstance(SolverType.mVMC, str)
        assert isinstance(SolverType.UHF, str)
        assert isinstance(SolverType.HWAVE, str)

    def test_string_equality(self):
        """SolverType members should compare equal to their string values."""
        assert SolverType.HPhi == "HPhi"
        assert SolverType.mVMC == "mVMC"
        assert SolverType.UHF == "UHF"
        assert SolverType.HWAVE == "HWAVE"
        assert "HPhi" == SolverType.HPhi
        assert "HWAVE" == SolverType.HWAVE

    def test_string_inequality(self):
        """SolverType members should not equal other strings."""
        assert SolverType.HPhi != "mVMC"
        assert SolverType.mVMC != "HPhi"
        assert SolverType.UHF != "HWAVE"
        assert SolverType.HWAVE != "UHF"

    def test_construction_from_string(self):
        """SolverType can be constructed from its string value."""
        assert SolverType("HPhi") is SolverType.HPhi
        assert SolverType("mVMC") is SolverType.mVMC
        assert SolverType("UHF") is SolverType.UHF
        assert SolverType("HWAVE") is SolverType.HWAVE

    def test_invalid_string_raises(self):
        """SolverType raises ValueError for invalid strings."""
        with pytest.raises(ValueError):
            SolverType("invalid")
        with pytest.raises(ValueError):
            SolverType("hphi")

    def test_can_be_used_as_dict_key(self):
        """SolverType members work as dict keys, interchangeable with strings."""
        d = {SolverType.HPhi: 1, SolverType.mVMC: 2}
        assert d["HPhi"] == 1
        assert d["mVMC"] == 2

    def test_solver_field_assignment(self):
        """StdIntList.solver can be set to a SolverType value."""
        s = StdIntList()
        s.solver = SolverType.HPhi
        assert s.solver == "HPhi"
        assert s.solver == SolverType.HPhi


# ===================================================================
#  MethodType enum
# ===================================================================


class TestMethodType:
    """Tests for the MethodType enum."""

    def test_has_seven_members(self):
        """MethodType should have exactly 7 members."""
        assert len(MethodType) == 7

    def test_values(self):
        """MethodType members should have the expected string values."""
        assert MethodType.LANCZOS == "lanczos"
        assert MethodType.LANCZOS_ENERGY == "lanczosenergy"
        assert MethodType.TPQ == "tpq"
        assert MethodType.FULLDIAG == "fulldiag"
        assert MethodType.CG == "cg"
        assert MethodType.TIME_EVOLUTION == "timeevolution"
        assert MethodType.CTPQ == "ctpq"

    def test_is_str_subclass(self):
        """MethodType members should be instances of str."""
        for member in MethodType:
            assert isinstance(member, str)

    def test_string_equality(self):
        """MethodType members should compare equal to their string values."""
        assert MethodType.LANCZOS == "lanczos"
        assert "lanczos" == MethodType.LANCZOS
        assert MethodType.TIME_EVOLUTION == "timeevolution"
        assert "timeevolution" == MethodType.TIME_EVOLUTION
        assert MethodType.FULLDIAG == "fulldiag"
        assert "fulldiag" == MethodType.FULLDIAG

    def test_string_inequality(self):
        """MethodType members should not equal other strings."""
        assert MethodType.LANCZOS != "tpq"
        assert MethodType.TPQ != "cg"
        assert MethodType.FULLDIAG != "lanczos"
        assert MethodType.TIME_EVOLUTION != "ctpq"

    def test_construction_from_string(self):
        """MethodType can be constructed from its string value."""
        assert MethodType("lanczos") is MethodType.LANCZOS
        assert MethodType("lanczosenergy") is MethodType.LANCZOS_ENERGY
        assert MethodType("tpq") is MethodType.TPQ
        assert MethodType("fulldiag") is MethodType.FULLDIAG
        assert MethodType("cg") is MethodType.CG
        assert MethodType("timeevolution") is MethodType.TIME_EVOLUTION
        assert MethodType("ctpq") is MethodType.CTPQ

    def test_invalid_string_raises(self):
        """MethodType raises ValueError for invalid strings."""
        with pytest.raises(ValueError):
            MethodType("invalid")
        with pytest.raises(ValueError):
            MethodType("Lanczos")  # case-sensitive

    def test_can_be_used_as_dict_key(self):
        """MethodType members work as dict keys, interchangeable with strings."""
        d = {MethodType.LANCZOS: 0, MethodType.TIME_EVOLUTION: 4}
        assert d["lanczos"] == 0
        assert d["timeevolution"] == 4

    def test_method_field_assignment(self):
        """StdIntList.method can be set to a MethodType value."""
        s = StdIntList()
        s.method = MethodType.FULLDIAG
        assert s.method == "fulldiag"
        assert s.method == MethodType.FULLDIAG


# ===================================================================
#  Sentinel constants
# ===================================================================


class TestSentinelConstants:
    """Tests for the module-level sentinel constants."""

    def test_nan_i_value(self):
        """NaN_i should be INT_MAX (2147483647)."""
        assert NaN_i == 2147483647

    def test_nan_i_type(self):
        """NaN_i should be an int."""
        assert isinstance(NaN_i, int)

    def test_nan_d_is_nan(self):
        """NaN_d should be IEEE NaN."""
        import math
        assert math.isnan(NaN_d)

    def test_nan_d_type(self):
        """NaN_d should be a float."""
        assert isinstance(NaN_d, float)

    def test_nan_c_real_is_nan(self):
        """NaN_c should have NaN real part."""
        import math
        assert math.isnan(NaN_c.real)

    def test_nan_c_imag_is_zero(self):
        """NaN_c should have zero imaginary part."""
        assert NaN_c.imag == 0.0

    def test_nan_c_type(self):
        """NaN_c should be a complex."""
        assert isinstance(NaN_c, complex)

    def test_unset_string_value(self):
        """UNSET_STRING should be '****'."""
        assert UNSET_STRING == "****"

    def test_unset_string_type(self):
        """UNSET_STRING should be a str."""
        assert isinstance(UNSET_STRING, str)

    def test_sentinels_are_distinct(self):
        """The three numeric sentinels should have different types."""
        assert type(NaN_i) is int
        assert type(NaN_d) is float
        assert type(NaN_c) is complex


class TestNumericalTolerances:
    """Tests for the AMPLITUDE_EPS and ZERO_BODY_EPS constants."""

    def test_amplitude_eps_value(self):
        """AMPLITUDE_EPS should be 1e-6."""
        assert AMPLITUDE_EPS == 1e-6

    def test_amplitude_eps_type(self):
        """AMPLITUDE_EPS should be a float."""
        assert isinstance(AMPLITUDE_EPS, float)

    def test_zero_body_eps_value(self):
        """ZERO_BODY_EPS should be 1e-12."""
        assert ZERO_BODY_EPS == 1e-12

    def test_zero_body_eps_type(self):
        """ZERO_BODY_EPS should be a float."""
        assert isinstance(ZERO_BODY_EPS, float)

    def test_zero_body_eps_stricter_than_amplitude_eps(self):
        """ZERO_BODY_EPS should be smaller (stricter) than AMPLITUDE_EPS."""
        assert ZERO_BODY_EPS < AMPLITUDE_EPS
