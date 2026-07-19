"""Geometry and structure output functions.

This module provides functions for writing lattice geometry files used by
post-processing tools and visualisation programs.

Functions
---------
build_xsf / XsfData
    Build the ``lattice.xsf`` (XCrysDen) data and write it.
build_geometry / GeometryData
    Build the ``geometry.dat`` data and write it.

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

from dataclasses import dataclass
from pathlib import Path


from ..core.stdface_vals import StdIntList, ModelType


def _cell_diff(Cell, iCell: int, jCell: int) -> list[int]:
    """Compute difference between two cell coordinate rows as integers.

    Parameters
    ----------
    Cell : np.ndarray
        ``(NCell, 3)`` array of fractional cell coordinates.
    iCell : int
        Row index of the first (minuend) cell.
    jCell : int
        Row index of the second (subtrahend) cell.

    Returns
    -------
    list of int
        ``[int(Cell[iCell, k] - Cell[jCell, k]) for k in range(3)]``.
    """
    return (Cell[iCell] - Cell[jCell]).astype(int).tolist()


# Lattices whose ``lattice.xsf`` includes a CONVVEC (conventional cell)
# block; others (e.g. wannier90) emit PRIMVEC/PRIMCOORD only.
_XSF_CONVVEC_LATTICES = (
    "orthorhombic", "face-centeredorthorhombic",
    "fcorthorhombic", "fco", "pyrochlore")


@dataclass
class XsfData:
    """XCrysDen structure (``lattice.xsf``) for visualisation.

    Lattice-level, solver-independent output (an independent path, like
    the gnuplot and ``geometry.dat`` files).

    Attributes
    ----------
    primvec : list of (float, float, float)
        Primitive lattice vectors (``box @ direct``), 3 rows.
    convvec : list of (float, float, float) or None
        Conventional cell vectors (orthorhombic / face-centred /
        pyrochlore), or ``None`` when not applicable.
    coords : list of (float, float, float)
        Atomic Cartesian coordinates (all written as element ``H``).
    """

    primvec: list
    convvec: list | None
    coords: list

    def write(self, directory: Path = Path(".")) -> None:
        lines = ["CRYSTAL\n", "PRIMVEC\n"]
        for vec in self.primvec:
            lines.append(f"{vec[0]:15.5f} {vec[1]:15.5f} {vec[2]:15.5f}\n")
        if self.convvec is not None:
            lines.append("CONVVEC\n")
            for vec in self.convvec:
                lines.append(f"{vec[0]:15.5f} {vec[1]:15.5f} {vec[2]:15.5f}\n")
        lines.append("PRIMCOORD\n")
        lines.append(f"{len(self.coords)} 1\n")
        for vec in self.coords:
            lines.append(f"H {vec[0]:15.5f} {vec[1]:15.5f} {vec[2]:15.5f}\n")
        with open(Path(directory) / "lattice.xsf", "w") as fp:
            fp.write("".join(lines))

    def to_dict(self) -> dict:
        return {"primvec": [list(r) for r in self.primvec],
                "convvec": (None if self.convvec is None
                            else [list(r) for r in self.convvec]),
                "coords": [list(c) for c in self.coords]}

    @classmethod
    def from_dict(cls, data: dict) -> "XsfData":
        conv = data["convvec"]
        return cls(primvec=[tuple(r) for r in data["primvec"]],
                   convvec=(None if conv is None else [tuple(r) for r in conv]),
                   coords=[tuple(c) for c in data["coords"]])


def build_xsf(StdI: StdIntList) -> "XsfData":
    """Build :class:`XsfData` from *StdI*.

    Reads ``lattice`` / ``box`` / ``direct`` / ``length`` / ``NCell`` /
    ``NsiteUC`` / ``Cell`` / ``tau``.  Always returns data; the caller
    decides whether the lattice warrants an xsf (3-D lattices only).
    """
    primvec_arr = StdI.box @ StdI.direct
    primvec = [(float(v[0]), float(v[1]), float(v[2])) for v in primvec_arr]

    convvec: "list | None" = None
    if StdI.lattice in _XSF_CONVVEC_LATTICES:
        convvec = []
        for ii, length_val in enumerate(StdI.length):
            row = [0.0, 0.0, 0.0]
            row[ii] = float(length_val)
            convvec.append((row[0], row[1], row[2]))

    coords = []
    for iCell in range(StdI.NCell):
        for isite in range(StdI.NsiteUC):
            frac_coord = StdI.Cell[iCell, :] + StdI.tau[isite, :]
            vec = frac_coord @ StdI.direct
            coords.append((float(vec[0]), float(vec[1]), float(vec[2])))

    return XsfData(primvec=primvec, convvec=convvec, coords=coords)


@dataclass
class GeometryData:
    """Lattice geometry (``geometry.dat``) for correlation post-processing.

    This is a lattice-level, solver-independent output (an independent
    path, like the gnuplot file).

    Attributes
    ----------
    direct : list of (float, float, float)
        Direct lattice vectors (3 rows).
    phase : list of float
        Boundary phase angles (length 3).
    box : list of (int, int, int)
        Supercell box matrix (3 rows).
    sites : list of (int, int, int, int)
        ``(d0, d1, d2, isite)`` per cell/site (Kondo doubling included).
    """

    direct: list
    phase: list
    box: list
    sites: list

    def write(self, directory: Path = Path(".")) -> None:
        lines = []
        for row in self.direct:
            lines.append(f"{row[0]:25.15e} {row[1]:25.15e} {row[2]:25.15e}\n")
        lines.append(f"{self.phase[0]:25.15e} "
                     f"{self.phase[1]:25.15e} "
                     f"{self.phase[2]:25.15e}\n")
        for row in self.box:
            lines.append(f"{row[0]} {row[1]} {row[2]}\n")
        for d0, d1, d2, isite in self.sites:
            lines.append(f"{d0} {d1} {d2} {isite}\n")
        with open(Path(directory) / "geometry.dat", "w") as fp:
            fp.write("".join(lines))

    def to_dict(self) -> dict:
        return {"direct": [list(r) for r in self.direct],
                "phase": list(self.phase),
                "box": [list(r) for r in self.box],
                "sites": [list(s) for s in self.sites]}

    @classmethod
    def from_dict(cls, data: dict) -> "GeometryData":
        return cls(direct=[tuple(r) for r in data["direct"]],
                   phase=list(data["phase"]),
                   box=[tuple(r) for r in data["box"]],
                   sites=[tuple(s) for s in data["sites"]])


def build_geometry(StdI: StdIntList) -> "GeometryData | None":
    """Build :class:`GeometryData`, or ``None`` for explicit Wannier modes.

    ``geometry.dat`` is suppressed for ``calcmode in ("uhfk", "rpa")``;
    an unset ``calcmode`` still produces it.
    """
    if StdI.calcmode in ("uhfk", "rpa"):
        return None

    direct = [(float(r[0]), float(r[1]), float(r[2])) for r in StdI.direct]
    phase = [float(StdI.phase[0]), float(StdI.phase[1]), float(StdI.phase[2])]
    box = [(int(r[0]), int(r[1]), int(r[2])) for r in StdI.box]

    sites: list = []
    for iCell in range(StdI.NCell):
        diff = _cell_diff(StdI.Cell, iCell, 0)
        for isite in range(StdI.NsiteUC):
            sites.append((diff[0], diff[1], diff[2], isite))
    if StdI.model == ModelType.KONDO:
        for iCell in range(StdI.NCell):
            diff = _cell_diff(StdI.Cell, iCell, 0)
            for isite in range(StdI.NsiteUC):
                sites.append((diff[0], diff[1], diff[2], isite + StdI.NsiteUC))

    return GeometryData(direct=direct, phase=phase, box=box, sites=sites)
