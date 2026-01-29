"""Super-cell site initialisation, folding, and labelling utilities.

This module provides functions for setting up the simulation super-cell,
folding site coordinates into the cell, finding site indices, and writing
gnuplot labels for 2-D lattice visualisation.

Functions
---------
init_site
    Initialise the super-cell (box, reciprocal box, cells, gnuplot header).
find_site
    Find global site indices and boundary phase for a pair of sites.
set_label
    Write gnuplot labels and return site indices (2-D lattices).
lattice_gp
    Context manager for ``lattice.gp`` gnuplot output (2-D lattices).
close_lattice_xsf
    Write ``lattice.xsf``, ``geometry.dat``, and finalise 3-D lattice output.

The private helper ``_fold_site`` is also available for use by other
modules that need to fold a coordinate into the original cell.

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

import itertools
from contextlib import contextmanager
from collections.abc import Iterator
from typing import TextIO

import numpy as np

from stdface_vals import StdIntList, ModelType, SolverType, NaN_i, AMPLITUDE_EPS
from param_check import exit_program, print_val_i
from .geometry_output import print_geometry, print_xsf


def _cell_vector(Cell: np.ndarray, idx: int) -> list[int]:
    """Extract a cell row as a list of three ints.

    ``StdI.Cell`` is stored as a ``float`` array, so individual elements
    must be cast to ``int`` whenever they are used in integer arithmetic
    or coordinate comparison.  This helper centralises that conversion.

    Parameters
    ----------
    Cell : np.ndarray
        ``(NCell, 3)`` array of fractional cell coordinates.
    idx : int
        Row index of the cell to extract.

    Returns
    -------
    list of int
        ``[int(Cell[idx, 0]), int(Cell[idx, 1]), int(Cell[idx, 2])]``.
    """
    return [int(Cell[idx, 0]), int(Cell[idx, 1]), int(Cell[idx, 2])]


def _find_cell_index(StdI: StdIntList, cellV: list[int]) -> int:
    """Find the cell index whose coordinates match *cellV*.

    Performs a linear search over ``StdI.Cell`` to find the row whose
    three coordinates match *cellV*.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure containing the ``Cell`` array and
        ``NCell`` count.
    cellV : list of int
        Length-3 fractional coordinate vector to search for.

    Returns
    -------
    int
        Cell index ``k`` such that ``StdI.Cell[k]`` equals *cellV*.
        Returns 0 if no match is found (matching the original C behavior).
    """
    for k in range(StdI.NCell):
        if _cell_vector(StdI.Cell, k) == cellV:
            return k
    return 0


def _fold_to_cell(
    rbox: np.ndarray,
    ncell: int,
    box: np.ndarray,
    iCellV: list[int],
) -> tuple[list[int], list[int]]:
    """Fold a site coordinate into a cell defined by *box*.

    This is the core algorithm shared by :func:`_fold_site` (main
    super-cell) and the sub-lattice folding in ``mvmc_variational``.

    Parameters
    ----------
    rbox : np.ndarray
        3×3 reciprocal (cofactor) matrix of *box*.
    ncell : int
        Determinant of *box* (number of cells).
    box : np.ndarray
        3×3 lattice-vector matrix defining the cell.
    iCellV : list of int
        Fractional coordinate of a site (length 3).

    Returns
    -------
    nBox : list of int
        Index of the periodic image that originally contained this site.
    iCellV_fold : list of int
        Fractional coordinate folded into the cell.
    """
    # (1) Transform to fractional coordinate (times ncell)
    iCellV_arr = np.asarray(iCellV)
    iCellV_frac = rbox @ iCellV_arr

    # (2) Search which periodic image contains this cell
    nBox = (iCellV_frac + ncell * 1000) // ncell - 1000

    # (3) Fractional coordinate in the original cell
    iCellV_frac = iCellV_frac - ncell * nBox

    # (4) Transform back to lattice coordinates and fold
    iCellV_fold = (box.T @ iCellV_frac + ncell * 1000) // ncell - 1000

    return nBox.astype(int).tolist(), iCellV_fold.astype(int).tolist()


def _fold_site(
    StdI: StdIntList,
    iCellV: list[int],
) -> tuple[list[int], list[int]]:
    """Move a site into the original super-cell if it is outside.

    Delegates to :func:`_fold_to_cell` with the main super-cell
    parameters (``rbox``, ``NCell``, ``box``).

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.
    iCellV : list of int
        Fractional coordinate of a site (length 3).

    Returns
    -------
    nBox : list of int
        Index of the super-cell that originally contained this site.
    iCellV_fold : list of int
        Fractional coordinate folded into the original cell.
    """
    return _fold_to_cell(StdI.rbox, StdI.NCell, StdI.box, iCellV)


def _write_gnuplot_header(fp: TextIO, StdI: StdIntList) -> None:
    """Write the gnuplot header for ``lattice.gp`` (2D lattices only).

    Computes the four corner positions of the super-cell from the
    ``direct`` and ``box`` matrices, sets up axis ranges, styles, and
    draws the boundary as four arrows.

    Parameters
    ----------
    fp : file object
        Open file to write gnuplot commands to.
    StdI : StdIntList
        Model parameter structure.  Reads ``direct`` and ``box``.
    """
    pos = np.zeros((4, 2))
    # pos[1:3] = supercell corner positions: box[:2,:2] @ direct[:2,:2]
    pos[1:3, :] = StdI.box[:2, :2] @ StdI.direct[:2, :2]
    pos[3, :] = pos[1, :] + pos[2, :]

    xmin = min(pos[:, 0].min(), pos[:, 1].min()) - 2.0
    xmax = max(pos[:, 0].max(), pos[:, 1].max()) + 2.0

    fp.write("#set terminal pdf color enhanced \\\n")
    fp.write("#dashed dl 1.0 size 20.0cm, 20.0cm \n")
    fp.write('#set output "lattice.pdf"\n')
    fp.write(f"set xrange [{xmin:f}: {xmax:f}]\n")
    fp.write(f"set yrange [{xmin:f}: {xmax:f}]\n")
    fp.write("set size square\n")
    fp.write("unset key\n")
    fp.write("unset tics\n")
    fp.write("unset border\n")
    fp.write("set style line 1 lc 1 lt 1\n")
    fp.write("set style line 2 lc 5 lt 1\n")
    fp.write("set style line 3 lc 0 lt 1\n")

    # Draw boundary arrows: 0→1→3→2→0
    corners = [(0, 1), (1, 3), (3, 2), (2, 0)]
    for src, dst in corners:
        fp.write(f"set arrow from {pos[src][0]:f}, {pos[src][1]:f} "
                 f"to {pos[dst][0]:f}, {pos[dst][1]:f} nohead front ls 3\n")


_LATTICE_GP_FOOTER = "plot '-' w d lc 7\n0.0 0.0\nend\npause -1\n"


@contextmanager
def lattice_gp(StdI: StdIntList) -> Iterator[TextIO | None]:
    """Context manager for ``lattice.gp`` gnuplot output.

    Opens ``lattice.gp`` for writing when the solver is not H-wave, or
    when ``StdI.lattice_gp`` is explicitly set to 1.  Otherwise the
    yielded handle is ``None``.  On exit the gnuplot footer is written,
    the file is closed, and :func:`print_geometry` is called.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.

    Yields
    ------
    fp : TextIO or None
        Open file handle, or ``None`` when output is suppressed.
    """
    fp: TextIO | None = None
    if StdI.solver != SolverType.HWAVE or StdI.lattice_gp == 1:
        fp = open("lattice.gp", "w")
    try:
        yield fp
    finally:
        if fp is not None:
            fp.write(_LATTICE_GP_FOOTER)
            fp.close()
        print_geometry(StdI)


def close_lattice_xsf(StdI: StdIntList) -> None:
    """Write ``lattice.xsf``, ``geometry.dat``, and print geometry.

    This is the 3-D counterpart of :func:`lattice_gp`.  It writes
    the XCrySDen structure file via :func:`print_xsf` and the
    ``geometry.dat`` file via :func:`print_geometry`.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.
    """
    print_xsf(StdI)
    print_geometry(StdI)


def _validate_box_params(
    L: int, W: int, Height: int,
    box: np.ndarray,
    suffix: str = "",
    defaults: np.ndarray | None = None,
) -> tuple[int, int, int]:
    """Validate and normalise L/W/Height vs box specification.

    Exactly one of {L, W, Height} or {box entries} may be specified
    (i.e. differ from ``NaN_i``).  If both are given an error is raised.
    When L/W/Height are specified, the box is filled as a diagonal matrix
    ``diag(W, L, Height)``.  When the box entries are specified (or
    neither is), each entry is validated with ``print_val_i`` using
    *defaults* as the fallback values.

    Parameters
    ----------
    L, W, Height : int
        Scalar lattice dimensions (may be ``NaN_i`` for unset).
    box : np.ndarray
        3x3 box matrix (modified **in-place**).
    suffix : str, optional
        Label suffix appended to parameter names (e.g. ``""`` for main
        lattice, ``"sub"`` for the sub-lattice).  Default ``""``.
    defaults : np.ndarray or None, optional
        3x3 array of default values for the box-entry branch.  If
        ``None``, the identity matrix ``diag(1, 1, 1)`` is used.

    Returns
    -------
    L, W, Height : int
        Validated (and possibly defaulted) values.

    Raises
    ------
    SystemExit
        If both L/W/Height and box entries are specified simultaneously.
    """
    sfx = suffix
    lwh_specified = (L != NaN_i or W != NaN_i or Height != NaN_i)
    box_specified = any(box[i, j] != NaN_i for i in range(3) for j in range(3))

    if lwh_specified and box_specified:
        lwh_names = f"(L{sfx}, W{sfx}, H{sfx + 'eight' if not sfx else sfx})"
        box_names = f"(a0W{sfx}, ..., a2H{sfx})"
        print(f"\nERROR ! {lwh_names} and {box_names} conflict !\n")
        exit_program(-1)
    elif lwh_specified:
        L = print_val_i(f"L{sfx}", L, 1)
        W = print_val_i(f"W{sfx}", W, 1)
        h_label = f"H{sfx + 'eight' if not sfx else sfx}"
        Height = print_val_i(h_label, Height, 1)
        box[:, :] = 0
        box[0, 0] = W
        box[1, 1] = L
        box[2, 2] = Height
    else:
        if defaults is None:
            defaults = np.eye(3, dtype=int)
        _labels = [
            [f"a0W{sfx}", f"a0L{sfx}", f"a0H{sfx}"],
            [f"a1W{sfx}", f"a1L{sfx}", f"a1H{sfx}"],
            [f"a2W{sfx}", f"a2L{sfx}", f"a2H{sfx}"],
        ]
        for i in range(3):
            for j in range(3):
                box[i, j] = print_val_i(
                    _labels[i][j], int(box[i, j]), int(defaults[i, j]))

    return L, W, Height


def _det_and_cofactor(box: np.ndarray) -> tuple[int, np.ndarray]:
    """Compute the determinant and cofactor matrix of a 3x3 integer matrix.

    Uses the Sarrus rule for the determinant and the standard cofactor
    formula for the adjugate (transposed cofactor matrix).  If the
    determinant is negative, both the determinant and cofactor matrix
    are sign-flipped so that the returned determinant is non-negative.

    Parameters
    ----------
    box : np.ndarray
        3x3 integer matrix.

    Returns
    -------
    det : int
        Absolute value of the determinant (always >= 0).
    cofactor : np.ndarray
        3x3 cofactor matrix, sign-adjusted so that ``det >= 0``.
    """
    # Compute determinant using Sarrus rule (equivalent to np.linalg.det for 3x3)
    det = int(round(np.linalg.det(box.astype(float))))

    # Compute cofactor matrix using vectorized indexing
    # cofactor[i,j] = box[(i+1)%3, (j+1)%3] * box[(i+2)%3, (j+2)%3]
    #               - box[(i+1)%3, (j+2)%3] * box[(i+2)%3, (j+1)%3]
    idx = np.array([1, 2, 0])  # (i+1) % 3 for i=0,1,2
    idx2 = np.array([2, 0, 1])  # (i+2) % 3 for i=0,1,2
    cofactor = (box[idx][:, idx] * box[idx2][:, idx2]
                - box[idx][:, idx2] * box[idx2][:, idx])

    if det < 0:
        cofactor = -cofactor
        det = -det

    return det, cofactor


def _compute_reciprocal_box(StdI: StdIntList) -> None:
    """Compute the number of cells and the reciprocal-box matrix.

    Sets ``StdI.NCell`` to the determinant of ``StdI.box`` and
    ``StdI.rbox`` to the cofactor matrix (adjugate transpose) of
    ``StdI.box``.  If the determinant is negative the sign of both
    ``NCell`` and ``rbox`` is flipped so that ``NCell > 0``.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
        Reads ``box``; sets ``NCell`` and ``rbox``.

    Raises
    ------
    SystemExit
        If the determinant is zero (degenerate super-cell).
    """
    StdI.NCell, StdI.rbox = _det_and_cofactor(StdI.box)
    print(f"   Number of Cell = {abs(StdI.NCell)}")
    if StdI.NCell == 0:
        exit_program(-1)


def _enumerate_cells(StdI: StdIntList) -> None:
    """Enumerate all unit cells inside the super-cell.

    Allocates ``StdI.Cell`` of shape ``(NCell, 3)`` and fills it with the
    integer coordinates of each unit cell that lies inside the super-cell.
    The super-cell is defined by the ``box`` matrix; a coordinate belongs
    to the cell if folding it via ``_fold_site`` produces ``nBox == [0,0,0]``.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
        Reads ``NCell`` and ``box``; sets ``Cell``.
    """
    # Find bounding box by checking all 8 cube corners
    box_int = StdI.box.astype(int)
    # All 8 corners of the unit cube: shape (8, 3)
    corners = np.array(list(itertools.product(range(2), repeat=3)))
    # edges[corner, dim] = corner @ box_int[:, dim] => corners @ box_int
    edges = corners @ box_int
    # bound[dim] = [min_edge, max_edge]
    bound = list(zip(edges.min(axis=0).tolist(), edges.max(axis=0).tolist()))

    # Enumerate cells within the bounding box
    # Note: iteration order must be ic2 outermost, ic0 innermost (matching C code)
    StdI.Cell = np.zeros((StdI.NCell, 3), dtype=int)
    jj_idx = 0
    for ic2, ic1, ic0 in itertools.product(
        range(bound[2][0], bound[2][1] + 1),
        range(bound[1][0], bound[1][1] + 1),
        range(bound[0][0], bound[0][1] + 1),
    ):
        iCellV = [ic0, ic1, ic2]
        nBox, iCellV_fold = _fold_site(StdI, iCellV)
        if nBox == [0, 0, 0]:
            StdI.Cell[jj_idx, :] = iCellV
            jj_idx += 1


def init_site(StdI: StdIntList, fp: TextIO | None, dim: int) -> None:
    """Initialize the super-cell where simulation is performed.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    fp : file object or None
        File pointer to ``lattice.gp`` (may be ``None``).
    dim : int
        Dimension of the system.  If 2, the gnuplot header for
        ``lattice.gp`` is written.
    """
    print("\n  @ Super-Lattice setting\n")

    # (1) Check input parameters about the shape of super-cell
    StdI.L, StdI.W, StdI.Height = _validate_box_params(
        StdI.L, StdI.W, StdI.Height, StdI.box)

    if dim == 2:
        StdI.direct[:2, 2] = 0.0   # zero z-component of first two vectors
        StdI.direct[2, :] = [0.0, 0.0, 1.0]  # third vector = unit z

    # (2) Define the phase factor at each boundary
    if dim == 2:
        StdI.phase[2] = 0.0
    StdI.ExpPhase = np.exp(1j * (np.pi / 180.0) * StdI.phase)
    StdI.AntiPeriod = np.where(np.abs(StdI.ExpPhase + 1.0) < AMPLITUDE_EPS, 1, 0)

    # (3) Allocate tau (intrinsic structure of unit-cell)
    StdI.tau = np.zeros((StdI.NsiteUC, 3))

    # (4) Calculate reciprocal lattice vectors and NCell
    _compute_reciprocal_box(StdI)

    # (5) Find cells in the super-cell
    _enumerate_cells(StdI)

    # (6) For 2D, print lattice.gp header
    if dim == 2 and fp is not None:
        _write_gnuplot_header(fp, StdI)


def find_site(
    StdI: StdIntList,
    iW: int, iL: int, iH: int,
    diW: int, diL: int, diH: int,
    isiteUC: int, jsiteUC: int,
) -> tuple[int, int, complex, np.ndarray]:
    """Find the site indices and boundary phase for a pair of sites.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.
    iW, iL, iH : int
        Position of the initial site.
    diW, diL, diH : int
        Translation from the initial site.
    isiteUC : int
        Intrinsic site index of the initial site in the unit cell.
    jsiteUC : int
        Intrinsic site index of the final site in the unit cell.

    Returns
    -------
    isite : int
        Global index of the initial site.
    jsite : int
        Global index of the final site.
    Cphase : complex
        Boundary phase factor.
    dR : numpy.ndarray
        Distance vector R_i - R_j in fractional coordinates (shape ``(3,)``).
    """
    di = np.array([diW, diL, diH], dtype=float)
    dR = -di + StdI.tau[isiteUC, :] - StdI.tau[jsiteUC, :]

    jCellV = [iW + diW, iL + diL, iH + diH]
    nBox, jCellV = _fold_site(StdI, jCellV)
    Cphase = np.prod(StdI.ExpPhase ** np.array(nBox))

    jCell = _find_cell_index(StdI, jCellV)
    iCell = _find_cell_index(StdI, [iW, iL, iH])

    isite = iCell * StdI.NsiteUC + isiteUC
    jsite = jCell * StdI.NsiteUC + jsiteUC
    if StdI.model == ModelType.KONDO:
        isite += StdI.NCell * StdI.NsiteUC
        jsite += StdI.NCell * StdI.NsiteUC

    return isite, jsite, Cphase, dR


def _write_gnuplot_bond(
    fp: TextIO,
    isite: int, jsite: int,
    xi: float, yi: float,
    xj: float, yj: float,
    connect: int,
) -> None:
    """Write gnuplot label and arrow commands for a single bond.

    Writes two ``set label`` commands (one per site endpoint) and,
    when *connect* < 3, a ``set arrow`` command connecting them.

    Parameters
    ----------
    fp : TextIO
        Open gnuplot file handle.
    isite, jsite : int
        Global site indices (used as label text).
    xi, yi : float
        2-D position of site *isite*.
    xj, yj : float
        2-D position of site *jsite*.
    connect : int
        Line-style selector.  The arrow is only drawn when ``connect < 3``.
    """
    if isite < 10:
        fp.write(f'set label "{isite:1d}" at {xi:f}, {yi:f} center front\n')
    else:
        fp.write(f'set label "{isite:2d}" at {xi:f}, {yi:f} center front\n')
    if jsite < 10:
        fp.write(f'set label "{jsite:1d}" at {xj:f}, {yj:f} center front\n')
    else:
        fp.write(f'set label "{jsite:2d}" at {xj:f}, {yj:f} center front\n')
    if connect < 3:
        fp.write(f"set arrow from {xi:f}, {yi:f} to {xj:f}, {yj:f} "
                 f"nohead ls {connect:d}\n")


def set_label(
    StdI: StdIntList,
    fp: TextIO | None,
    iW: int, iL: int,
    diW: int, diL: int,
    isiteUC: int, jsiteUC: int,
    connect: int,
) -> tuple[int, int, complex, np.ndarray]:
    """Set label in the gnuplot display (2D systems only).

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.
    fp : file object or None
        File pointer to ``lattice.gp``.
    iW, iL : int
        Position of the initial site.
    diW, diL : int
        Translation from the initial site.
    isiteUC : int
        Intrinsic site index of the initial site.
    jsiteUC : int
        Intrinsic site index of the final site.
    connect : int
        Connection type (1 for nearest, 2 for 2nd nearest).

    Returns
    -------
    isite : int
        Global index of the initial site.
    jsite : int
        Global index of the final site.
    Cphase : complex
        Boundary phase factor.
    dR : numpy.ndarray
        Distance vector R_i - R_j.
    """
    # First print the reversed one
    isite, jsite, Cphase, dR = find_site(
        StdI, iW, iL, 0, -diW, -diL, 0, jsiteUC, isiteUC)

    # Compute 2D positions via direct[:2,:2].T @ fractional_coords
    D = StdI.direct[:2, :2]
    frac_i = np.array([iW + StdI.tau[jsiteUC, 0], iL + StdI.tau[jsiteUC, 1]])
    frac_j = np.array([iW - diW + StdI.tau[isiteUC, 0], iL - diL + StdI.tau[isiteUC, 1]])
    xi, yi = frac_i @ D
    xj, yj = frac_j @ D

    if fp is not None:
        _write_gnuplot_bond(fp, isite, jsite, xi, yi, xj, yj, connect)

    # Then print the normal one
    isite, jsite, Cphase, dR = find_site(
        StdI, iW, iL, 0, diW, diL, 0, isiteUC, jsiteUC)

    frac_i = np.array([iW + StdI.tau[isiteUC, 0], iL + StdI.tau[isiteUC, 1]])
    frac_j = np.array([iW + diW + StdI.tau[jsiteUC, 0], iL + diL + StdI.tau[jsiteUC, 1]])
    xi, yi = frac_i @ D
    xj, yj = frac_j @ D

    if fp is not None:
        _write_gnuplot_bond(fp, isite, jsite, xi, yi, xj, yj, connect)

    return isite, jsite, Cphase, dR


def set_local_spin_flags(StdI: StdIntList, nsite_base: int) -> None:
    """Set local spin flags and total site count based on the model type.

    This function is shared across all lattice modules. It sets
    ``StdI.nsite`` (doubling for Kondo), allocates ``StdI.locspinflag``,
    and fills the array according to the model type:

    * **SPIN**: all sites get ``S2``
    * **HUBBARD**: all sites get ``0``
    * **KONDO**: first half ``S2``, second half ``0``

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    nsite_base : int
        Number of sites before Kondo doubling (e.g.
        ``NsiteUC * NCell``, ``L``, or ``L * NsiteUC``).
    """
    StdI.nsite = nsite_base
    if StdI.model == ModelType.KONDO:
        StdI.nsite *= 2
    StdI.locspinflag = np.zeros(StdI.nsite, dtype=int)

    if StdI.model == ModelType.SPIN:
        StdI.locspinflag[:] = StdI.S2
    elif StdI.model == ModelType.HUBBARD:
        StdI.locspinflag[:] = 0
    elif StdI.model == ModelType.KONDO:
        half = StdI.nsite // 2
        StdI.locspinflag[:half] = StdI.S2
        StdI.locspinflag[half:] = 0

