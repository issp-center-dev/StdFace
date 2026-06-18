"""Geometry and structure output functions.

This module provides functions for writing lattice geometry files used by
post-processing tools and visualisation programs.

Functions
---------
print_xsf
    Write ``lattice.xsf`` in XCrysDen format.
print_geometry
    Write ``geometry.dat`` for correlation-function post-processing.

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

import numpy as np

from ..core.stdface_vals import StdIntList, ModelType, SolverType


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


def print_xsf(StdI: StdIntList) -> None:
    """Print lattice.xsf file (XCrysDen format).

    Writes a ``lattice.xsf`` file containing primitive vectors,
    optional conventional vectors (for orthorhombic and face-centred
    lattices), and atomic coordinates for visualisation.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.  The following fields are read:

        - ``lattice`` : str -- lattice type name.
        - ``box`` : ndarray (3x3) -- supercell box matrix.
        - ``direct`` : ndarray (3x3) -- direct lattice vectors.
        - ``length`` : ndarray (3,) -- lattice constant lengths.
        - ``NCell`` : int -- number of unit cells.
        - ``NsiteUC`` : int -- sites per unit cell.
        - ``Cell`` : ndarray (NCell, 3) -- cell fractional coordinates.
        - ``tau`` : ndarray (NsiteUC, 3) -- basis positions.
    """
    do_convvec = StdI.lattice in (
        "orthorhombic", "face-centeredorthorhombic",
        "fcorthorhombic", "fco", "pyrochlore")

    with open("lattice.xsf", "w") as fp:
        fp.write("CRYSTAL\n")
        fp.write("PRIMVEC\n")
        # PRIMVEC = box @ direct (each row is a primitive vector)
        primvec = StdI.box @ StdI.direct
        for vec in primvec:
            fp.write(f"{vec[0]:15.5f} {vec[1]:15.5f} {vec[2]:15.5f}\n")

        if do_convvec:
            fp.write("CONVVEC\n")
            for ii, length_val in enumerate(StdI.length):
                row = [0.0, 0.0, 0.0]
                row[ii] = length_val
                fp.write(f"{row[0]:15.5f} {row[1]:15.5f} {row[2]:15.5f}\n")

        fp.write("PRIMCOORD\n")
        fp.write(f"{StdI.NCell * StdI.NsiteUC} 1\n")
        for iCell in range(StdI.NCell):
            for isite in range(StdI.NsiteUC):
                # vec = (Cell[iCell] + tau[isite]) @ direct
                frac_coord = StdI.Cell[iCell, :] + StdI.tau[isite, :]
                vec = frac_coord @ StdI.direct
                fp.write(f"H {vec[0]:15.5f} {vec[1]:15.5f} {vec[2]:15.5f}\n")


def print_geometry(StdI: StdIntList) -> None:
    """Print geometry.dat for post-processing of correlation functions.

    Writes a ``geometry.dat`` file containing direct lattice vectors,
    boundary phases, the supercell box matrix, and site coordinates
    relative to the first cell.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.  The following fields are read:

        - ``solver`` : str -- solver name.
        - ``calcmode`` : str -- calculation mode (HWAVE only).
        - ``direct`` : ndarray (3x3) -- direct lattice vectors.
        - ``phase`` : ndarray (3,) -- boundary phase angles.
        - ``box`` : ndarray (3x3) -- supercell box matrix.
        - ``NCell`` : int -- number of unit cells.
        - ``NsiteUC`` : int -- sites per unit cell.
        - ``Cell`` : ndarray (NCell, 3) -- cell fractional coordinates.
        - ``model`` : str -- model name (``"kondo"`` doubles sites).
    """
    # Suppressed for explicit Wannier modes; an unset calcmode still writes it.
    if StdI.calcmode in ("uhfk", "rpa"):
        return

    with open("geometry.dat", "w") as fp:
        for row in StdI.direct:
            fp.write(f"{row[0]:25.15e} {row[1]:25.15e} {row[2]:25.15e}\n")
        fp.write(f"{StdI.phase[0]:25.15e} "
                 f"{StdI.phase[1]:25.15e} "
                 f"{StdI.phase[2]:25.15e}\n")
        for row in StdI.box:
            fp.write(f"{int(row[0])} {int(row[1])} {int(row[2])}\n")

        for iCell in range(StdI.NCell):
            diff = _cell_diff(StdI.Cell, iCell, 0)
            for isite in range(StdI.NsiteUC):
                fp.write(f"{diff[0]} {diff[1]} {diff[2]} {isite}\n")
        if StdI.model == ModelType.KONDO:
            for iCell in range(StdI.NCell):
                diff = _cell_diff(StdI.Cell, iCell, 0)
                for isite in range(StdI.NsiteUC):
                    fp.write(f"{diff[0]} {diff[1]} {diff[2]} "
                             f"{isite + StdI.NsiteUC}\n")
