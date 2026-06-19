"""Solver output containers (D2).

A :class:`SolverOutput` bundles the per-file ``XxxData`` objects produced
by the ``build_*`` functions into a single value a caller can ``write`` to
a directory or serialise with ``to_dict``.  This is the unit a library
``generate``-style API hands back.

:class:`WannierModeOutput` (UHFK / RPA) and :class:`ExpertModeOutput`
(HPhi / mVMC / UHF) are provided.  For UHF the Expert container is fully
data-backed; HPhi / mVMC additionally emit solver-specific files
(excitation / calcmod / variational) via a plugin hook, which are not yet
represented as data.

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


@dataclass
class ExpertModeOutput(SolverOutput):
    """Expert-mode (.def) output for HPhi / mVMC / UHF.

    Holds the data-backed common files.  Solver-specific files
    (HPhi excitation / calcmod / pump, mVMC variational group) are not
    represented here; they are emitted by the plugin's
    ``write_solver_files`` hook during assembly.
    """

    locspn: object        # LocSpnData
    trans: object         # TransData
    interactions: list    # list[InteractionData | InterAllData]
    modpara: object       # ModParaData
    namelist: object      # NamelistData
    green_one: object | None = None   # GreenOneData
    green_two: object | None = None   # GreenTwoData

    def write(self, directory: Path = Path(".")) -> None:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.locspn.write(directory)
        self.trans.write(directory)
        for data in self.interactions:
            data.write(directory)
        self.modpara.write(directory)
        if self.green_one is not None:
            self.green_one.write(directory)
        if self.green_two is not None:
            self.green_two.write(directory)
        self.namelist.write(directory)

    def to_dict(self) -> dict:
        return {
            "locspn": self.locspn.to_dict(),
            "trans": self.trans.to_dict(),
            "interactions": [d.to_dict() for d in self.interactions],
            "modpara": self.modpara.to_dict(),
            "namelist": self.namelist.to_dict(),
            "green_one": (self.green_one.to_dict()
                          if self.green_one is not None else None),
            "green_two": (self.green_two.to_dict()
                          if self.green_two is not None else None),
        }


# ---------------------------------------------------------------------------
#  Output formats (D3): how a SolverOutput is serialised
# ---------------------------------------------------------------------------


class OutputFormat(ABC):
    """Strategy for rendering a :class:`SolverOutput` to *directory*.

    Separates *what* to output (the :class:`SolverOutput` data) from *how*
    it is serialised (``.def`` files vs a single JSON document).
    """

    @abstractmethod
    def write_output(self, output: SolverOutput,
                     directory: Path = Path(".")) -> None:
        """Serialise *output* into *directory*."""
        raise NotImplementedError


class DefFileFormat(OutputFormat):
    """Legacy ``.def`` / ``.dat`` file format (delegates to ``output.write``)."""

    def write_output(self, output: SolverOutput,
                     directory: Path = Path(".")) -> None:
        output.write(directory)


class JSONFormat(OutputFormat):
    """Bundle the whole output into a single JSON document.

    ``complex`` values are already split into ``[re, im]`` lists inside
    each ``XxxData.to_dict()``, so the standard ``json`` encoder suffices.
    """

    def __init__(self, filename: str = "stdface_output.json") -> None:
        self.filename = filename

    def write_output(self, output: SolverOutput,
                     directory: Path = Path(".")) -> None:
        import json
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        with open(directory / self.filename, "w") as fp:
            json.dump(output.to_dict(), fp, indent=2, ensure_ascii=False)
