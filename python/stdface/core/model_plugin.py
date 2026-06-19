"""Model plugins: per-model on-site (local) Hamiltonian term generation.

This is the B5 ``ModelPlugin`` layer.  Each model (SPIN / HUBBARD /
KONDO) knows how to build the on-site local terms for a single site as a
pure :class:`~stdface.lattice.interaction_builder.LocalTerms` value (no
``StdIntList`` mutation).  The lattice builders reach this layer through
``add_local_terms`` (see ``interaction_builder``), which extends the
returned terms into ``StdI``.

License
-------
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from .stdface_vals import StdIntList, ModelType
from ..lattice.interaction_builder import (
    LocalTerms,
    mag_field_terms,
    hubbard_local_terms,
    general_j_terms,
)


class ModelPlugin(ABC):
    """Base class for a model's on-site (local) term generator."""

    #: Canonical model type this plugin handles.
    name: ModelType

    @abstractmethod
    def build_local_terms(
        self, StdI: StdIntList, isite: int, jsite_kondo: int
    ) -> LocalTerms:
        """Return the on-site local terms for site *isite*.

        Parameters
        ----------
        StdI : StdIntList
            Parameter structure (read-only).
        isite : int
            Global site index.  For KONDO this is the itinerant site.
        jsite_kondo : int
            Global index of the localized spin (KONDO only; ignored
            otherwise).
        """
        raise NotImplementedError


class SpinModel(ModelPlugin):
    """SPIN model: magnetic field + single-ion anisotropy ``D``."""

    name = ModelType.SPIN

    def build_local_terms(
        self, StdI: StdIntList, isite: int, jsite_kondo: int
    ) -> LocalTerms:
        terms = LocalTerms()
        terms.trans += mag_field_terms(
            StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, isite)
        terms.merge(general_j_terms(
            StdI.D, StdI.S2, StdI.S2, isite, isite, StdI.solver, StdI.model))
        return terms


class HubbardModel(ModelPlugin):
    """HUBBARD model: Hubbard ``U``, chemical potential and magnetic field."""

    name = ModelType.HUBBARD

    def build_local_terms(
        self, StdI: StdIntList, isite: int, jsite_kondo: int
    ) -> LocalTerms:
        return hubbard_local_terms(
            StdI.mu, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, StdI.U, isite)


class KondoModel(ModelPlugin):
    """KONDO model: itinerant Hubbard terms + Kondo ``J`` + localized field."""

    name = ModelType.KONDO

    def build_local_terms(
        self, StdI: StdIntList, isite: int, jsite_kondo: int
    ) -> LocalTerms:
        terms = hubbard_local_terms(
            StdI.mu, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, StdI.U, isite)
        terms.merge(general_j_terms(
            StdI.J, 1, StdI.S2, isite, jsite_kondo, StdI.solver, StdI.model))
        terms.trans += mag_field_terms(
            StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, jsite_kondo)
        return terms


# ---------------------------------------------------------------------------
#  Registry
# ---------------------------------------------------------------------------

_MODEL_REGISTRY: dict[ModelType, ModelPlugin] = {}


def register_model(plugin: ModelPlugin) -> None:
    """Register *plugin* under its :attr:`ModelPlugin.name`."""
    _MODEL_REGISTRY[plugin.name] = plugin


def get_model(model: ModelType) -> ModelPlugin:
    """Return the registered :class:`ModelPlugin` for *model*.

    Raises
    ------
    KeyError
        If no plugin is registered for *model*.
    """
    return _MODEL_REGISTRY[model]


register_model(SpinModel())
register_model(HubbardModel())
register_model(KondoModel())
