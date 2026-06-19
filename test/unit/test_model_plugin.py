"""Unit tests for the ModelPlugin layer (B5) in core/model_plugin.py."""
from __future__ import annotations

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList, SolverType, ModelType
from stdface.core.model_plugin import (
    ModelPlugin, SpinModel, HubbardModel, KondoModel,
    get_model, register_model,
)
from stdface.lattice.interaction_builder import LocalTerms


def _make_stdi(model, **kw):
    s = StdIntList()
    s.solver = SolverType.HPhi
    s.model = model
    s.mu = kw.get("mu", 1.0)
    s.h = kw.get("h", 0.0)
    s.Gamma = kw.get("Gamma", 0.0)
    s.Gamma_y = kw.get("Gamma_y", 0.0)
    s.U = kw.get("U", 4.0)
    s.S2 = kw.get("S2", 1)
    return s


class TestRegistry:
    def test_three_models_registered(self):
        assert isinstance(get_model(ModelType.SPIN), SpinModel)
        assert isinstance(get_model(ModelType.HUBBARD), HubbardModel)
        assert isinstance(get_model(ModelType.KONDO), KondoModel)

    def test_lookup_by_string(self):
        # ModelType is a str-enum, so plain strings resolve too
        assert isinstance(get_model("spin"), SpinModel)

    def test_unknown_model_raises(self):
        with pytest.raises(KeyError):
            get_model("does_not_exist")

    def test_plugins_are_modelplugin(self):
        for m in (ModelType.SPIN, ModelType.HUBBARD, ModelType.KONDO):
            assert isinstance(get_model(m), ModelPlugin)


class TestBuildLocalTerms:
    def test_returns_localterms_without_mutating_stdi(self):
        s = _make_stdi(ModelType.HUBBARD)
        terms = get_model(s.model).build_local_terms(s, 0, 0)
        assert isinstance(terms, LocalTerms)
        # pure: StdI term lists untouched
        assert s.trans_list == [] and s.Cintra_list == []

    def test_hubbard_local_terms(self):
        s = _make_stdi(ModelType.HUBBARD, U=4.0)
        terms = get_model(s.model).build_local_terms(s, 2, 0)
        assert terms.Cintra == [(4.0, 2)]
        # mu contributes the two diagonal transfers
        assert len(terms.trans) == 2

    def test_spin_uses_field_and_D(self):
        s = _make_stdi(ModelType.SPIN, h=1.0, S2=1)
        s.D[2, 2] = 0.5
        terms = get_model(s.model).build_local_terms(s, 0, 0)
        # magnetic field -> transfers; D single-ion -> Hund/Cinter (spin-1/2)
        assert len(terms.trans) > 0
        assert terms.Hund == [(-0.25, 0, 0)]

    def test_kondo_combines_hubbard_J_and_localized_field(self):
        s = _make_stdi(ModelType.KONDO, h=0.2, U=4.0, S2=1)
        s.J[0, 0] = s.J[1, 1] = s.J[2, 2] = 0.5
        terms = get_model(s.model).build_local_terms(s, 0, 3)
        # Cintra from itinerant Hubbard U on isite=0
        assert terms.Cintra == [(4.0, 0)]
        # Kondo J coupling produces Hund between itinerant(0) and localized(3)
        assert terms.Hund and terms.Hund[0][1:] == (0, 3)
        # both itinerant (mu/h) and localized (field) transfers present
        assert len(terms.trans) >= 3

    def test_matches_legacy_add_local_terms(self):
        # End-to-end: ModelPlugin path equals the old direct-call sequence.
        import stdface.lattice.interaction_builder as ib
        s = _make_stdi(ModelType.KONDO, h=0.2, Gamma=0.1, U=4.0, S2=1)
        s.J[0, 0] = s.J[1, 1] = s.J[2, 2] = 0.5
        ib.add_local_terms(s, isite=0, jsite_kondo=3)

        ref = _make_stdi(ModelType.KONDO, h=0.2, Gamma=0.1, U=4.0, S2=1)
        ref.J[0, 0] = ref.J[1, 1] = ref.J[2, 2] = 0.5
        ib.hubbard_local(ref, ref.mu, -ref.h, -ref.Gamma, -ref.Gamma_y, ref.U, 0)
        ib.general_j(ref, ref.J, 1, ref.S2, 0, 3)
        ib.mag_field(ref, ref.S2, -ref.h, -ref.Gamma, -ref.Gamma_y, 3)

        for attr in ("trans_list", "intr_list", "Cintra_list", "Cinter_list",
                     "Hund_list", "Ex_list", "PairLift_list"):
            assert getattr(s, attr) == getattr(ref, attr), attr


class TestRegisterModel:
    def test_register_and_retrieve_custom(self):
        class _Dummy(ModelPlugin):
            name = "dummy_model"

            def build_local_terms(self, StdI, isite, jsite_kondo):
                return LocalTerms()

        register_model(_Dummy())
        try:
            assert isinstance(get_model("dummy_model"), _Dummy)
        finally:
            from stdface.core import model_plugin
            del model_plugin._MODEL_REGISTRY["dummy_model"]
