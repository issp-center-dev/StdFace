"""Unit tests for stdface_vals module.

Tests for the Python translation of StdFace_vals.h.
"""
from __future__ import annotations

import numpy as np
import pytest

from stdface.core.stdface_vals import (
    StdIntList, HamiltonianTerms, LatticeGeometry, ModelInput,
    ModelType, SolverType, MethodType,
    NaN_i,
    AMPLITUDE_EPS, ZERO_BODY_EPS,
    is_unset_or_trivial_d,
)


class TestModelInputSplit:
    """C1-3: model Hamiltonian parameters live in a ModelInput sub-object,
    exposed on StdIntList via transparent façade properties."""

    def test_model_subobject_present(self):
        s = StdIntList()
        assert isinstance(s._model, ModelInput)

    def test_facade_reads_subobject(self):
        s = StdIntList()
        assert s.J is s._model.J
        assert s.t is None and s.U is None

    def test_facade_scalar_write(self):
        s = StdIntList()
        s.t = 1.0 + 2.0j
        s.U = 4.0
        s.JAll = 0.3
        assert s._model.t == 1.0 + 2.0j
        assert s._model.U == 4.0
        assert s._model.JAll == 0.3

    def test_facade_inplace_numpy_mutation(self):
        s = StdIntList()
        s.J[0, 0] = 1.5
        s.D[2, 2] = -0.7
        assert s._model.J[0, 0] == 1.5
        assert s._model.D[2, 2] == -0.7

    def test_instances_independent(self):
        s1, s2 = StdIntList(), StdIntList()
        s1.t = 9.0
        s1.J[0, 0] = 5.0
        assert s2.t is None
        assert s2.J[0, 0] == 0.0


class TestLatticeGeometrySplit:
    """C1-2: lattice geometry lives in a LatticeGeometry sub-object,
    exposed on StdIntList via transparent façade properties."""

    def test_lattice_subobject_present(self):
        s = StdIntList()
        assert isinstance(s._lattice, LatticeGeometry)

    def test_facade_reads_subobject(self):
        s = StdIntList()
        assert s.box is s._lattice.box
        assert s.direct is s._lattice.direct
        assert s.W is None and s.nsite == 0

    def test_facade_inplace_numpy_mutation(self):
        s = StdIntList()
        s.box[0, 0] = 5
        s.direct[1, 2] = 3.0
        assert s._lattice.box[0, 0] == 5
        assert s._lattice.direct[1, 2] == 3.0

    def test_facade_scalar_write(self):
        s = StdIntList()
        s.W, s.L, s.Height = 2, 3, 4
        s.nsite = 24
        assert (s._lattice.W, s._lattice.L, s._lattice.Height) == (2, 3, 4)
        assert s._lattice.nsite == 24

    def test_instances_independent(self):
        s1, s2 = StdIntList(), StdIntList()
        s1.box[0, 0] = 9
        s1.NCell = 7
        assert s2.box[0, 0] == 0
        assert s2.NCell == 0


class TestHamiltonianTermsSplit:
    """C1-1: Hamiltonian terms live in a HamiltonianTerms sub-object,
    exposed on StdIntList via transparent façade properties."""

    def test_terms_subobject_present(self):
        s = StdIntList()
        assert isinstance(s._terms, HamiltonianTerms)

    def test_facade_reads_subobject(self):
        s = StdIntList()
        assert s.trans_list is s._terms.trans_list
        assert s.Cintra_list is s._terms.Cintra_list
        assert s.LCintra == s._terms.LCintra == 0

    def test_facade_write_reaches_subobject(self):
        s = StdIntList()
        s.LCinter = 1
        assert s._terms.LCinter == 1
        s.intr_list.append((1 + 0j, 0, 0, 0, 0, 0, 0, 0, 0))
        assert len(s._terms.intr_list) == 1

    def test_subobject_write_visible_via_facade(self):
        s = StdIntList()
        s._terms.Hund_list.append((0.5, 0, 1))
        assert s.Hund_list == [(0.5, 0, 1)]

    def test_instances_independent(self):
        s1, s2 = StdIntList(), StdIntList()
        s1.trans_list.append((1 + 0j, 0, 0, 1, 1))
        assert s2.trans_list == []


class TestIsUnsetOrTrivialD:
    """Tests for is_unset_or_trivial_d (helper for SolverPlugin.validate)."""

    def test_nan_is_unset(self):
        assert is_unset_or_trivial_d(float("nan")) is True

    def test_none_is_unset(self):
        assert is_unset_or_trivial_d(None) is True

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
        """String fields should default to None or empty string."""
        s = StdIntList()
        assert s.lattice is None
        assert s.model is None
        assert s.outputmode is None
        assert s.CDataFileHead is None
        assert s.solver == ""

    def test_int_defaults(self):
        """Sentinel-reset int fields default to None; computed ones to 0."""
        s = StdIntList()
        # A2: lattice dimensions are unset (None) until specified
        assert s.W is None
        assert s.L is None
        assert s.Height is None
        # Computed counters keep their 0 default
        assert s.NCell == 0
        assert s.NsiteUC == 0
        assert s.nsite == 0
        assert len(s.trans_list) == 0
        assert len(s.intr_list) == 0

    def test_float_defaults(self):
        """A2: sentinel-reset float fields default to None (unset)."""
        s = StdIntList()
        assert s.a is None
        assert s.mu is None
        assert s.U is None
        assert s.h is None
        assert s.Gamma is None

    def test_complex_defaults(self):
        """A2: complex hopping fields default to None (unset)."""
        s = StdIntList()
        assert s.t is None
        assert s.tp is None
        assert s.t0 is None
        assert s.t1 is None
        assert s.t2 is None
        assert s.tpp is None

    def test_dynamic_arrays_none(self):
        """Dynamic (pointer) arrays should be initialized to None."""
        s = StdIntList()
        assert s.Cell is None
        assert s.tau is None
        assert s.locspinflag is None
        assert s.trans_list == []
        assert s.Cintra_list == []
        assert s.Cinter_list == []
        assert s.Hund_list == []
        assert s.Ex_list == []
        assert s.PairLift_list == []
        assert s.PairHopp_list == []


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

    def test_hphi_fields_resolve_via_config(self):
        """HPhi-unique fields moved to HPhiConfig (C3-5)."""
        bare = StdIntList()
        assert not hasattr(bare, "method")
        assert not hasattr(bare, "Lanczos_max")

        s = StdIntList()
        s.solver = "HPhi"  # auto-attaches HPhiConfig
        assert s.method is None
        assert s.CalcSpec is None
        assert s.SpectrumQ.shape == (3,)

    def test_mvmc_fields_resolve_via_config(self):
        """mVMC-unique fields moved to MVMCConfig (C3-4).

        They are absent on a bare StdIntList; an mVMC run resolves them
        through the attached config.
        """
        from stdface.core.stdface_main import _reset_vals

        bare = StdIntList()
        assert not hasattr(bare, "CParaFileHead")
        assert not hasattr(bare, "NVMCCalMode")
        assert not hasattr(bare, "Orb")

        s = StdIntList()
        s.solver = "mVMC"
        _reset_vals(s)
        assert s.NVMCCalMode is None
        assert s.Orb is None
        assert s.NSROptFixSmp == 0  # non-None default preserved

    def test_uhf_hwave_fields_resolve_via_config(self):
        """UHF/HWAVE shared fields (mix/eps/...) moved to the solver configs."""
        bare = StdIntList()
        assert not hasattr(bare, "mix")
        assert not hasattr(bare, "eps")

        s = StdIntList()
        s.solver = "UHF"  # auto-attaches UHFConfig
        assert s.mix is None
        assert s.eps is None
        assert s.Iteration_max is None

    def test_calcmode_stays_on_core(self):
        """calcmode is solver-selection metadata and stays on StdIntList."""
        s = StdIntList()
        assert hasattr(s, "calcmode")

    def test_hwave_only_fields_resolve_via_config(self):
        """fileprefix/export_all/lattice_gp moved to HWaveConfig (C3-3).

        They no longer live on a bare StdIntList; for an H-wave run the
        attached config provides them.
        """
        from stdface.core.stdface_main import _reset_vals

        bare = StdIntList()
        assert not hasattr(bare, "fileprefix")
        assert not hasattr(bare, "export_all")
        assert not hasattr(bare, "lattice_gp")

        s = StdIntList()
        s.solver = "HWAVE"
        _reset_vals(s)
        # resolved via the attached HWaveConfig
        assert s.fileprefix is None
        assert s.export_all is None
        assert s.lattice_gp is None

    def test_boxsub_shape(self):
        """boxsub/rboxsub (shared sublattice block) moved to the solver configs."""
        s = StdIntList()
        s.solver = "mVMC"  # auto-attaches MVMCConfig carrying the sublattice copy
        assert s.boxsub.shape == (3, 3)
        assert np.issubdtype(s.boxsub.dtype, np.integer)
        assert s.rboxsub.shape == (3, 3)
        assert np.issubdtype(s.rboxsub.dtype, np.integer)

    def test_lambda_renamed(self):
        """C 'lambda' field should be 'lambda_' in Python."""
        s = StdIntList()
        assert hasattr(s, "lambda_")
        assert s.lambda_ is None

    def test_wannier90_cutoffs_delegate_to_w90(self):
        """W1: cutoff/lambda/alpha fields live in the _w90 sub-object."""
        from stdface.core.stdface_vals import Wannier90Cutoff
        s = StdIntList()
        assert isinstance(s._w90, Wannier90Cutoff)
        # scalar facade: write reaches the sub-object
        s.alpha = 0.5
        assert s._w90.alpha == 0.5
        # in-place numpy mutation through the getter
        s.cutoff_tR[1] = 3
        assert s._w90.cutoff_tR[1] == 3
        assert s.cutoff_tVec.shape == (3, 3)


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
    """Tests for the surviving sentinel constant.

    A2 removed NaN_d / NaN_c / UNSET_STRING; only NaN_i remains, used as
    the unset marker for integer matrices (box / cutoff_*R / boxsub).
    """

    def test_nan_i_value(self):
        """NaN_i should be INT_MAX (2147483647) -- kept for integer matrices."""
        assert NaN_i == 2147483647

    def test_nan_i_type(self):
        """NaN_i should be an int."""
        assert isinstance(NaN_i, int)

    def test_removed_sentinels_are_gone(self):
        """NaN_d / NaN_c / UNSET_STRING were removed in A2."""
        import stdface.core.stdface_vals as vals
        assert not hasattr(vals, "NaN_d")
        assert not hasattr(vals, "NaN_c")
        assert not hasattr(vals, "UNSET_STRING")


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
