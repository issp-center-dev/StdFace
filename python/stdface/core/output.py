"""Solver output containers (D2).

A :class:`SolverOutput` bundles the per-file ``XxxData`` objects produced
by the ``build_*`` functions into a single value a caller can ``write`` to
a directory or serialise with ``to_dict``.  This is the unit a library
``generate``-style API hands back.

Currently only :class:`WannierModeOutput` (UHFK / RPA) is provided; it is
fully backed by D1 data objects.  ``ExpertModeOutput`` is deferred until
``ModParaData`` / ``NamelistData`` exist.

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
from dataclasses import dataclass
from pathlib import Path

from .stdface_vals import StdIntList
from ..writer.wannier90_writer import (
    WannierGeometryData,
    WannierInteractionData,
    build_wannier_geometry,
    build_wannier_interactions,
    _prefix,
)


class SolverOutput(ABC):
    """Abstract base for all solver output containers.

    Subclasses hold the mode-specific ``XxxData`` objects and implement
    :meth:`write` (emit files) and :meth:`to_dict` (JSON-serialisable
    view).  Gnuplot output is lattice-level and solver-independent, so it
    is **not** part of the solver output; it is handled on its own path.
    """

    @abstractmethod
    def write(self, directory: Path = Path(".")) -> None:
        """Write all output files into *directory*."""
        raise NotImplementedError

    @abstractmethod
    def to_dict(self) -> dict:
        """Return a JSON-serialisable dict of the output contents."""
        raise NotImplementedError


@dataclass
class WannierModeOutput(SolverOutput):
    """Wannier90-mode output (UHFK / RPA): ``geom.dat`` + interaction files."""

    geometry: WannierGeometryData
    geom_fname: str
    interactions: list  # list[WannierInteractionData]

    def write(self, directory: Path = Path(".")) -> None:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.geometry.write(str(directory / self.geom_fname))
        for data in self.interactions:
            data.write(directory)

    def to_dict(self) -> dict:
        return {
            "geometry": self.geometry.to_dict(),
            "geom_fname": self.geom_fname,
            "interactions": [d.to_dict() for d in self.interactions],
        }


def build_wannier_output(StdI: StdIntList) -> WannierModeOutput:
    """Assemble a :class:`WannierModeOutput` from *StdI*."""
    return WannierModeOutput(
        geometry=build_wannier_geometry(StdI),
        geom_fname=_prefix(StdI, "geom.dat"),
        interactions=build_wannier_interactions(StdI),
    )
