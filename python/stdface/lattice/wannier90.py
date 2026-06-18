"""
Standard mode for the Wannier90 interface.

This module sets up the Hamiltonian for the Wannier90 ``*_hr.dat`` format,
supporting hopping, Coulomb, and Hund coupling terms read from
Wannier90/RESPACK output files.

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

import logging
import itertools
import math
from enum import IntEnum
from typing import NamedTuple, TextIO

import numpy as np

from ..core.stdface_vals import StdIntList, ModelType, SolverType, NaN_i, UNSET_STRING, AMPLITUDE_EPS
from ..core.param_check import print_val_d, print_val_i, not_used_d
from .geometry_output import print_geometry, print_xsf
from .interaction_builder import (
    malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb,
)
from .site_util import init_site, find_site, set_local_spin_flags


# ---------------------------------------------------------------------------
#  Internal helpers
# ---------------------------------------------------------------------------


logger = logging.getLogger(__name__)


def _check_in_box(rvec: np.ndarray, inverse_matrix: np.ndarray) -> bool:
    """Check if a lattice vector is inside the unit cell box.

    Parameters
    ----------
    rvec : numpy.ndarray
        Integer lattice vector to check (shape ``(3,)``).
    inverse_matrix : numpy.ndarray
        Inverse of the cutoff lattice vectors (shape ``(3, 3)``).

    Returns
    -------
    bool
        True if inside the box, False otherwise.
    """
    judge_vec = rvec @ inverse_matrix
    return bool(np.all(np.abs(judge_vec) <= 1))


def _skip_degeneracy_weights(fp: TextIO, n_wigner_seitz: int) -> None:
    """Skip the degeneracy-weight lines in a Wannier90 ``*_hr.dat`` file.

    The weights are written as whitespace-separated integers, potentially
    spanning multiple lines.  This helper reads and discards exactly
    *n_wigner_seitz* values.

    Parameters
    ----------
    fp : TextIO
        Open file positioned just after the ``nWSC`` header line.
    n_wigner_seitz : int
        Total number of Wigner-Seitz cells (degeneracy entries to skip).
    """
    count = 0
    while count < n_wigner_seitz:
        count += len(fp.readline().split())


def _geometry_w90(StdI: StdIntList) -> None:
    """Read Wannier90 geometry file.

    Reads lattice vectors and Wannier center positions from the geometry
    file ``<CDataFileHead>_geom.dat``.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters. Modified in-place:
        ``StdI.direct`` (lattice vectors) and ``StdI.tau`` (Wannier centres)
        are populated.
    """
    filename = f"{StdI.CDataFileHead}_geom.dat"
    logger.info(f"    Wannier90 Geometry file = {filename}")

    try:
        fp_geom = open(filename, "r")
    except OSError as exc:
        logger.error("Fail to open the file %s", filename)
        raise FileNotFoundError(filename) from exc

    with fp_geom:
        # Read direct lattice vectors
        for ii in range(3):
            StdI.direct[ii, :] = [float(x) for x in fp_geom.readline().split()[:3]]

        # Read number of correlated sites
        StdI.NsiteUC = int(fp_geom.readline().split()[0])
        logger.info(f"    Number of Correlated Sites = {StdI.NsiteUC}")

        # Allocate and read Wannier centre positions
        StdI.tau = np.zeros((StdI.NsiteUC, 3))
        for isite in range(StdI.NsiteUC):
            StdI.tau[isite, :] = [float(x) for x in fp_geom.readline().split()[:3]]

    logger.info("    Direct lattice vectors:")
    for row in StdI.direct:
        logger.info(f"      {row[0]:10.5f} {row[1]:10.5f} {row[2]:10.5f}")
    logger.info("    Wannier centres:")
    for tau_row in StdI.tau[:StdI.NsiteUC]:
        logger.info(f"      {tau_row[0]:10.5f} {tau_row[1]:10.5f} {tau_row[2]:10.5f}")


def _apply_boundary_weights(
    indx_tot: np.ndarray,
    Weight_tot: np.ndarray,
    nWSC: int,
    StdI: StdIntList,
) -> np.ndarray:
    """Apply boundary-halving weights at model lattice boundaries.

    For periodic models, matrix elements at the boundary of the model
    super-cell are halved to avoid double-counting.

    Parameters
    ----------
    indx_tot : numpy.ndarray
        R-vector indices for each Wigner-Seitz cell, shape ``(nWSC, 3)``.
    Weight_tot : numpy.ndarray
        Weight array for each WSC, shape ``(nWSC,)``.  Modified in-place.
    nWSC : int
        Number of Wigner-Seitz cells.
    StdI : StdIntList
        Structure containing model parameters (``W``, ``L``, ``Height``).

    Returns
    -------
    numpy.ndarray
        Band lattice extent for each dimension, shape ``(3,)``, dtype int.
    """
    # Compute max absolute index per dimension
    Band_lattice = np.max(np.abs(indx_tot[:nWSC]), axis=0).astype(int)

    if StdI.W != NaN_i and StdI.L != NaN_i and StdI.Height != NaN_i:
        dims = np.array([StdI.W, StdI.L, StdI.Height], dtype=int)
        # Model_lattice[i] = dims[i] // 2 if dims[i] is even, else 0
        Model_lattice = np.where(dims % 2 == 0, dims // 2, 0)
        for ii in range(3):
            if Model_lattice[ii] < Band_lattice[ii] and Model_lattice[ii] != 0:
                # Apply 0.5 weight at boundary
                mask = np.abs(indx_tot[:nWSC, ii]) == Model_lattice[ii]
                Weight_tot[:nWSC][mask] *= 0.5

    return Band_lattice


def _count_and_store_terms(
    Mat_tot: np.ndarray,
    indx_tot: np.ndarray,
    Weight_tot: np.ndarray,
    nWSC: int,
    NsiteUC: int,
    cutoff: float,
    itUJ: int,
    NtUJ: list[int],
    tUJindx: list,
    tUJ: list,
) -> None:
    """Apply weights, count effective terms, and store surviving terms.

    Multiplies each matrix element by its weight, counts terms above the
    cutoff threshold, prints them, and packs the surviving terms into
    the output arrays ``tUJ`` and ``tUJindx``.

    Parameters
    ----------
    Mat_tot : numpy.ndarray
        Matrix elements, shape ``(nWSC, nWan, nWan)``, dtype complex.
    indx_tot : numpy.ndarray
        R-vector indices, shape ``(nWSC, 3)``, dtype int.
    Weight_tot : numpy.ndarray
        Per-WSC weights, shape ``(nWSC,)``.
    nWSC : int
        Number of Wigner-Seitz cells.
    NsiteUC : int
        Number of correlated sites in the unit cell.
    cutoff : float
        Threshold for matrix elements.
    itUJ : int
        Interaction type index (0=t, 1=U, 2=J).
    NtUJ : list of int
        Counts per interaction type.  Modified in-place.
    tUJindx : list
        Index arrays.  Modified in-place.
    tUJ : list
        Coefficient arrays.  Modified in-place.
    """
    # Apply weights: broadcast Weight_tot over Wannier indices
    Mat_tot[:nWSC, :, :] *= Weight_tot[:nWSC, np.newaxis, np.newaxis]

    # Print and count effective terms
    logger.info("\n      EFFECTIVE terms:")
    logger.info("           R0   R1   R2 band_i band_f Hamiltonian")
    NtUJ[itUJ] = 0
    for iWSC in range(nWSC):
        for iWan in range(NsiteUC):
            for jWan in range(NsiteUC):
                if cutoff < abs(Mat_tot[iWSC, iWan, jWan]):
                    logger.info(
                        "        %5d%5d%5d%5d%5d%12.6f%12.6f",
                        indx_tot[iWSC, 0], indx_tot[iWSC, 1],
                        indx_tot[iWSC, 2], iWan, jWan,
                        Mat_tot[iWSC, iWan, jWan].real,
                        Mat_tot[iWSC, iWan, jWan].imag,
                    )
                    NtUJ[itUJ] += 1
    logger.info(f"      Total number of EFFECTIVE term = {NtUJ[itUJ]}")

    # Extract surviving terms using numpy masking
    abs_mat = np.abs(Mat_tot[:nWSC, :NsiteUC, :NsiteUC])
    mask = abs_mat > cutoff
    wsc_idx, iwan_idx, jwan_idx = np.nonzero(mask)

    tUJ_arr = Mat_tot[wsc_idx, iwan_idx, jwan_idx].copy()
    tUJindx_arr = np.column_stack([
        indx_tot[wsc_idx, :],
        iwan_idx,
        jwan_idx,
    ])

    # Extend the lists to hold the results
    while len(tUJ) <= itUJ:
        tUJ.append(None)
    while len(tUJindx) <= itUJ:
        tUJindx.append(None)
    tUJ[itUJ] = tUJ_arr
    tUJindx[itUJ] = tUJindx_arr


def _read_w90(
    StdI: StdIntList,
    filename: str,
    cutoff: float,
    cutoff_R: np.ndarray,
    cutoff_Rvec: np.ndarray,
    cutoff_length: float,
    itUJ: int,
    NtUJ: list[int],
    tUJindx: list,
    lam: float,
    tUJ: list,
) -> None:
    """Read Wannier90 hopping/interaction file.

    Reads hopping or interaction matrix elements from Wannier90 files,
    applies cutoffs and stores non-zero terms.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters.
    filename : str
        Input filename (e.g. ``*_hr.dat``, ``*_ur.dat``, ``*_jr.dat``).
    cutoff : float
        Threshold for matrix elements.
    cutoff_R : numpy.ndarray
        Cutoff for R vectors (shape ``(3,)``, int).
    cutoff_Rvec : numpy.ndarray
        Cutoff vectors for unit cell (shape ``(3, 3)``).
    cutoff_length : float
        Real-space cutoff length.
    itUJ : int
        Type of interaction (0: hopping t, 1: Coulomb U, 2: Hund J).
    NtUJ : list of int
        Number of terms for each interaction type. Modified in-place.
    tUJindx : list
        Indices for terms. Modified in-place.
    lam : float
        Scaling factor (lambda).
    tUJ : list
        Matrix elements. Modified in-place.
    """
    flg_vec = int(cutoff_Rvec[0, 0] != NaN_i)

    # Try to open the file
    try:
        fp_hr = open(filename, "r")
    except FileNotFoundError:
        logger.info(f"\n  Skip to read the file {filename}. \n")
        return

    with fp_hr:
        # Header part
        _header_line = fp_hr.readline()  # comment line
        nWan = int(fp_hr.readline().split()[0])
        nWSC = int(fp_hr.readline().split()[0])

        # Skip degeneracy weights
        _skip_degeneracy_weights(fp_hr, nWSC)

        # Allocate arrays
        Weight_tot = np.ones(nWSC)
        Mat_tot = np.zeros((nWSC, nWan, nWan), dtype=complex)
        indx_tot = np.zeros((nWSC, 3), dtype=int)

        if flg_vec:
            inverse_rvec = np.linalg.inv(cutoff_Rvec.astype(float))

        # Read body
        for iWSC in range(nWSC):
            for iWan in range(nWan):
                for jWan in range(nWan):
                    vals = fp_hr.readline().split()
                    indx_tot[iWSC, :] = [int(vals[0]), int(vals[1]), int(vals[2])]
                    iWan0 = int(vals[3])
                    jWan0 = int(vals[4])
                    dtmp_re = float(vals[5])
                    dtmp_im = float(vals[6])
                    # Compute Euclidean length
                    tau_diff = StdI.tau[jWan, :] - StdI.tau[iWan, :] + indx_tot[iWSC, :]
                    dR = StdI.direct.T @ tau_diff
                    length = np.linalg.norm(dR)
                    if length > cutoff_length > 0.0:
                        dtmp_re = 0.0
                        dtmp_im = 0.0

                    if flg_vec:
                        if not _check_in_box(indx_tot[iWSC], inverse_rvec):
                            dtmp_re = 0.0
                            dtmp_im = 0.0
                    else:
                        if np.any(np.abs(indx_tot[iWSC]) > cutoff_R):
                            dtmp_re = 0.0
                            dtmp_im = 0.0

                    if iWan0 <= StdI.NsiteUC and jWan0 <= StdI.NsiteUC:
                        Mat_tot[iWSC, iWan0 - 1, jWan0 - 1] = lam * (dtmp_re + 1j * dtmp_im)

            # Apply inversion symmetry and delete duplication
            if iWSC > 0 and np.any(np.all(indx_tot[iWSC] == -indx_tot[:iWSC], axis=1)):
                Mat_tot[iWSC, :, :] = 0.0

            if np.all(indx_tot[iWSC] == 0):
                for iWan in range(StdI.NsiteUC):
                    Mat_tot[iWSC, iWan, :iWan] = 0.0

    # Apply boundary-halving weights
    _apply_boundary_weights(indx_tot, Weight_tot, nWSC, StdI)

    # Count effective terms, print summary, and store
    _count_and_store_terms(
        Mat_tot, indx_tot, Weight_tot, nWSC,
        StdI.NsiteUC, cutoff, itUJ, NtUJ, tUJindx, tUJ,
    )


def _read_density_matrix(
    StdI: StdIntList,
    filename: str,
) -> dict[tuple[int, int, int], np.ndarray]:
    """Read RESPACK density matrix file.

    Reads density matrix elements from a RESPACK output file.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters.
    filename : str
        Input filename (e.g. ``*_dr.dat``).

    Returns
    -------
    dict of tuple to numpy.ndarray
        Dictionary mapping ``(R0, R1, R2)`` lattice vector tuples to
        2D numpy arrays of shape ``(NsiteUC, NsiteUC)`` containing the
        density matrix elements.
    """


    try:
        fp_dr = open(filename, "r")
    except OSError as exc:
        logger.error("Fail to open the file %s", filename)
        raise FileNotFoundError(filename) from exc

    with fp_dr:
        # Header
        _header_line = fp_dr.readline()
        nWan = int(fp_dr.readline().split()[0])
        nWSC = int(fp_dr.readline().split()[0])

        _skip_degeneracy_weights(fp_dr, nWSC)

        # Allocate
        Mat_tot = np.zeros((nWSC, nWan, nWan), dtype=complex)
        indx_tot = np.zeros((nWSC, 3), dtype=int)

        Rmin = np.zeros(3, dtype=int)
        Rmax = np.zeros(3, dtype=int)

        # Read body
        for iWSC in range(nWSC):
            for iWan in range(nWan):
                for jWan in range(nWan):
                    vals = fp_dr.readline().split()
                    indx_tot[iWSC, :] = [int(vals[0]), int(vals[1]), int(vals[2])]
                    iWan0 = int(vals[3])
                    jWan0 = int(vals[4])
                    dtmp_re = float(vals[5])
                    dtmp_im = float(vals[6])

                    if iWan0 <= StdI.NsiteUC and jWan0 <= StdI.NsiteUC:
                        Mat_tot[iWSC, iWan0 - 1, jWan0 - 1] = dtmp_re + 1j * dtmp_im
                    Rmin = np.minimum(Rmin, indx_tot[iWSC])
                    Rmax = np.maximum(Rmax, indx_tot[iWSC])

    NR = Rmax - Rmin + 1
    logger.info(f"      Minimum R : {Rmin[0]} {Rmin[1]} {Rmin[2]}")
    logger.info(f"      Maximum R : {Rmax[0]} {Rmax[1]} {Rmax[2]}")
    logger.info(f"      Numver of R : {NR[0]} {NR[1]} {NR[2]}")

    # Build dictionary: (R0, R1, R2) -> 2D array
    DenMat: dict[tuple[int, int, int], np.ndarray] = {}
    for i0, i1, i2 in itertools.product(
        range(Rmin[0], Rmax[0] + 1),
        range(Rmin[1], Rmax[1] + 1),
        range(Rmin[2], Rmax[2] + 1),
    ):
        DenMat[(i0, i1, i2)] = np.zeros(
            (StdI.NsiteUC, StdI.NsiteUC), dtype=complex
        )

    for iWSC in range(nWSC):
        key = tuple(indx_tot[iWSC].astype(int))
        DenMat[key][:, :] = Mat_tot[iWSC, :nWan, :nWan]

    return DenMat


def _print_uhf_initial(
    StdI: StdIntList,
    NtUJ: list[int],
    tUJ: list[np.ndarray],
    DenMat: dict[tuple[int, int, int], np.ndarray],
    tUJindx: list[np.ndarray],
) -> None:
    """Print initial UHF guess to ``initial.def``.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters.
    NtUJ : list of int
        Number of terms for each interaction type.
    tUJ : list of numpy.ndarray
        Matrix elements.
    DenMat : dict
        Density matrix elements keyed by R-vector tuples.
    tUJindx : list of numpy.ndarray
        Indices for interaction terms.
    """
    IniGuess = np.zeros((StdI.nsite, StdI.nsite), dtype=complex)

    for kCell in range(StdI.NCell):
        cell_w = StdI.Cell[kCell, 0]
        cell_l = StdI.Cell[kCell, 1]
        iH = StdI.Cell[kCell, 2]

        # Diagonal term
        for isite_uc in range(StdI.NsiteUC):
            jsite = isite_uc + StdI.NsiteUC * kCell
            IniGuess[jsite, jsite] = DenMat[(0, 0, 0)][isite_uc, isite_uc]

        # Coulomb integral (U) and Exchange integral (J)
        for idx in (1, 2):
            for it in range(NtUJ[idx]):
                row = tUJindx[idx][it]
                isite, jsite, Cphase, dR = find_site(
                    StdI, cell_w, cell_l, iH,
                    int(row[0]), int(row[1]), int(row[2]),
                    int(row[3]), int(row[4]),
                )
                key = (int(row[0]), int(row[1]), int(row[2]))
                dm_val = DenMat[key][int(row[3]), int(row[4])]
                IniGuess[isite, jsite] = dm_val
                IniGuess[jsite, isite] = np.conj(dm_val)

    mask = np.abs(IniGuess) > AMPLITUDE_EPS
    NIniGuess = int(np.count_nonzero(mask))

    with open("initial.def", "w") as fp:
        fp.write("======================== \n")
        fp.write(f"NInitialGuess {NIniGuess * 2:7d}  \n")
        fp.write("======================== \n")
        fp.write("========i_j_s_tijs====== \n")
        fp.write("======================== \n")

        rows, cols = np.nonzero(mask)
        for isite, jsite in zip(rows, cols):
            val = 0.5 * IniGuess[isite, jsite]
            for ispin in range(2):
                fp.write(
                    f"{jsite:5d} {ispin:5d} {isite:5d} {ispin:5d} "
                    f"{val.real:25.15f} "
                    f"{val.imag:25.15f}\n"
                )

    logger.info("      initial.def is written.")


# ---------------------------------------------------------------------------
#  Double-counting mode enum
# ---------------------------------------------------------------------------


class _DCMode(IntEnum):
    """Double-counting correction mode."""
    NOTCORRECT = 0
    HARTREE = 1
    HARTREE_U = 2
    FULL = 3


_DC_MODE_MAP: dict[str, _DCMode] = {
    "none":      _DCMode.NOTCORRECT,
    UNSET_STRING: _DCMode.NOTCORRECT,
    "hartree":   _DCMode.HARTREE,
    "hartree_u": _DCMode.HARTREE_U,
    "full":      _DCMode.FULL,
}
"""Maps double-counting mode strings to :class:`_DCMode` enum members.

The sentinel :data:`UNSET_STRING` (``"****"``) is treated the same as
``"none"`` (no correction).
"""


def _parse_double_counting_mode(mode_str: str) -> _DCMode:
    """Convert a double-counting mode string to the corresponding enum.

    Uses the :data:`_DC_MODE_MAP` dispatch table for lookup.

    Parameters
    ----------
    mode_str : str
        Mode specification from the input file.  Recognised values (case
        sensitive) are ``"none"``, ``"hartree"``, ``"hartree_u"`` and
        ``"full"``.  The sentinel :data:`UNSET_STRING` is treated the same
        as ``"none"``.

    Returns
    -------
    _DCMode
        The matching enum member.

    Raises
    ------
    ValueError
        If *mode_str* is not one of the recognised values.
    """
    result = _DC_MODE_MAP.get(mode_str)
    if result is not None:
        return result


    msg = (
        "the word of doublecounting is not correct "
        "(select from none, hartree, hartree_u, full)."
    )
    logger.error(msg)
    raise ValueError(msg)


# ---------------------------------------------------------------------------
#  Cutoff-parameter setup + file read helper
# ---------------------------------------------------------------------------


def _read_w90_with_cutoff(
    StdI: StdIntList,
    label: str,
    label_suffix: str,
    file_suffix: str,
    cutoff_val: float,
    cutoff_length: float,
    cutoff_R: np.ndarray,
    cutoff_Vec: np.ndarray,
    cutoff_length_default: float,
    cutoff_R_defaults: tuple[int | None, int | None, int | None],
    itUJ: int,
    NtUJ: list[int],
    tUJindx: list,
    lam: float,
    tUJ: list,
) -> tuple[float, float]:
    """Set cutoff parameters and read a Wannier90 interaction file.

    Prints parameter values, sets cutoff thresholds for R-vectors and
    real-space length, then calls :func:`_read_w90` to read the file.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters.
    label : str
        Lowercase label for the cutoff threshold name
        (e.g. ``"t"``, ``"u"``, ``"j"``).
    label_suffix : str
        Suffix for length/R/Vec parameter names (e.g. ``"t"``, ``"U"``,
        ``"J"``).  May differ in case from *label*.
    file_suffix : str
        File suffix including underscore (e.g. ``"_hr.dat"``).
    cutoff_val : float
        Current cutoff threshold value (may be NaN if unset).
    cutoff_length : float
        Current cutoff length value (may be NaN if unset).
    cutoff_R : numpy.ndarray
        Cutoff R-vector array (shape ``(3,)``, int).  Modified in-place.
    cutoff_Vec : numpy.ndarray
        Cutoff vector matrix (shape ``(3, 3)``).  Modified in-place.
    cutoff_length_default : float
        Default value for cutoff length if unset.
    cutoff_R_defaults : tuple of (int or None)
        Default values for cutoff_R[0], cutoff_R[1], cutoff_R[2].
        Use ``None`` to skip a dimension (leave unchanged).
    itUJ : int
        Interaction type index (0=hopping, 1=Coulomb, 2=Hund).
    NtUJ : list of int
        Number of terms per interaction type.  Modified in-place.
    tUJindx : list
        Indices per interaction type.  Modified in-place.
    lam : float
        Scaling factor (lambda).
    tUJ : list
        Matrix elements per interaction type.  Modified in-place.

    Returns
    -------
    tuple of (float, float)
        Updated ``(cutoff_val, cutoff_length)`` after applying defaults.
    """
    cutoff_name = f"cutoff_{label}"
    cutoff_length_name = f"cutoff_length_{label_suffix}"
    cutoff_R_name = f"cutoff_{label_suffix}R"
    cutoff_Vec_name = f"cutoff_{label_suffix}Vec"

    cutoff_val = print_val_d(cutoff_name, cutoff_val, 1.0e-8)
    cutoff_length = print_val_d(cutoff_length_name, cutoff_length, cutoff_length_default)

    for dim in range(3):
        if cutoff_R_defaults[dim] is not None:
            cutoff_R[dim] = print_val_i(
                f"{cutoff_R_name}[{dim}]", int(cutoff_R[dim]), cutoff_R_defaults[dim]
            )

    for i in range(3):
        for j in range(3):
            if StdI.box[i, j] != NaN_i:
                param_name = f"{cutoff_Vec_name}[{i}][{j}]"
            cutoff_Vec[i, j] = print_val_d(
                param_name, cutoff_Vec[i, j], float(StdI.box[i, j]) * 0.5
            )

    filename = f"{StdI.CDataFileHead}{file_suffix}"
    _read_w90(
        StdI, filename,
        cutoff_val, cutoff_R, cutoff_Vec, cutoff_length,
        itUJ, NtUJ, tUJindx, lam, tUJ,
    )

    return cutoff_val, cutoff_length


# ---------------------------------------------------------------------------
#  Per-cell interaction helpers
# ---------------------------------------------------------------------------


def _apply_hopping_terms(
    StdI: StdIntList,
    kCell: int,
    cell_w: int,
    cell_l: int,
    iH: int,
    NtUJ: list[int],
    tUJ: list,
    tUJindx: list,
    Uspin: np.ndarray | None,
) -> None:
    """Apply hopping transfer terms for one unit cell.

    Processes all hopping (t) terms for the unit cell at position
    ``(cell_w, cell_l, iH)``.  Local terms contribute on-site energies (Hubbard)
    and non-local terms contribute either super-exchange (spin) or
    hopping integrals (Hubbard).

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters.  Modified in-place.
    kCell : int
        Linear index of the current unit cell.
    cell_w, cell_l, iH : int
        Unit-cell coordinates.
    NtUJ : list of int
        Number of terms per interaction type.
    tUJ : list
        Matrix elements per interaction type.
    tUJindx : list
        Indices per interaction type.
    Uspin : numpy.ndarray or None
        On-site Coulomb values per orbital (spin model only).
    """
    if tUJindx[0] is None:
        return
    for it in range(NtUJ[0]):
        # Local term
        if (tUJindx[0][it, 0] == 0 and tUJindx[0][it, 1] == 0
                and tUJindx[0][it, 2] == 0
                and tUJindx[0][it, 3] == tUJindx[0][it, 4]):
            if StdI.model == ModelType.HUBBARD:
                isite = StdI.NsiteUC * kCell + int(tUJindx[0][it, 3])
                for ispin in range(2):
                    StdI.trans[StdI.ntrans] = -tUJ[0][it]
                    StdI.transindx[StdI.ntrans, 0] = isite
                    StdI.transindx[StdI.ntrans, 1] = ispin
                    StdI.transindx[StdI.ntrans, 2] = isite
                    StdI.transindx[StdI.ntrans, 3] = ispin
                    StdI.ntrans += 1
        else:
            # Non-local term
            isite, jsite, Cphase, dR = find_site(
                StdI, cell_w, cell_l, iH,
                int(tUJindx[0][it, 0]), int(tUJindx[0][it, 1]),
                int(tUJindx[0][it, 2]),
                int(tUJindx[0][it, 3]), int(tUJindx[0][it, 4]),
            )
            if StdI.model == ModelType.SPIN:
                diag_val = (
                    2.0 * tUJ[0][it] * np.conj(tUJ[0][it])
                    * (1.0 / Uspin[int(tUJindx[0][it, 3])]
                       + 1.0 / Uspin[int(tUJindx[0][it, 4])])
                ).real
                Jtmp = np.diag([diag_val, diag_val, diag_val])
                general_j(StdI, Jtmp, StdI.S2, StdI.S2, isite, jsite)
            else:
                hopping(StdI, -Cphase * tUJ[0][it], jsite, isite, dR)


def _apply_coulomb_terms(
    StdI: StdIntList,
    kCell: int,
    cell_w: int,
    cell_l: int,
    iH: int,
    NtUJ: list[int],
    tUJ: list,
    tUJindx: list,
    idcmode: _DCMode,
    DenMat: dict[tuple[int, int, int], np.ndarray] | None,
) -> None:
    """Apply Coulomb (U) interaction terms for one unit cell.

    Processes all Coulomb terms for the unit cell at position
    ``(cell_w, cell_l, iH)``, including local intra-site Coulomb, non-local
    inter-site Coulomb, and double-counting corrections when enabled.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters.  Modified in-place.
    kCell : int
        Linear index of the current unit cell.
    cell_w, cell_l, iH : int
        Unit-cell coordinates.
    NtUJ : list of int
        Number of terms per interaction type.
    tUJ : list
        Matrix elements per interaction type.
    tUJindx : list
        Indices per interaction type.
    idcmode : _DCMode
        Double-counting correction mode.
    DenMat : dict or None
        Density matrix elements keyed by R-vector tuples.
    """
    if tUJindx[1] is None:
        return
    for it in range(NtUJ[1]):
        # Local term
        if (tUJindx[1][it, 0] == 0 and tUJindx[1][it, 1] == 0
                and tUJindx[1][it, 2] == 0
                and tUJindx[1][it, 3] == tUJindx[1][it, 4]):
            StdI.Cintra[StdI.NCintra] = tUJ[1][it].real
            StdI.CintraIndx[StdI.NCintra, 0] = (
                StdI.NsiteUC * kCell + int(tUJindx[1][it, 3])
            )
            StdI.NCintra += 1

            # Double-counting correction
            if idcmode != _DCMode.NOTCORRECT:
                isite = StdI.NsiteUC * kCell + int(tUJindx[1][it, 3])
                for ispin in range(2):
                    DenMat0 = DenMat[(0, 0, 0)][
                        int(tUJindx[1][it, 3]), int(tUJindx[1][it, 3])
                    ]
                    StdI.trans[StdI.ntrans] = (
                        StdI.alpha * tUJ[1][it].real * DenMat0
                    )
                    StdI.transindx[StdI.ntrans, 0] = isite
                    StdI.transindx[StdI.ntrans, 1] = ispin
                    StdI.transindx[StdI.ntrans, 2] = isite
                    StdI.transindx[StdI.ntrans, 3] = ispin
                    StdI.ntrans += 1
        else:
            # Non-local term
            isite, jsite, Cphase, dR = find_site(
                StdI, cell_w, cell_l, iH,
                int(tUJindx[1][it, 0]), int(tUJindx[1][it, 1]),
                int(tUJindx[1][it, 2]),
                int(tUJindx[1][it, 3]), int(tUJindx[1][it, 4]),
            )
            coulomb(StdI, tUJ[1][it].real, isite, jsite)

            # Double-counting correction
            if idcmode != _DCMode.NOTCORRECT:
                for ispin in range(2):
                    # U_{Rij} D_{0jj} (Local)
                    DenMat0 = DenMat[(0, 0, 0)][
                        int(tUJindx[1][it, 4]), int(tUJindx[1][it, 4])
                    ]
                    StdI.trans[StdI.ntrans] = tUJ[1][it].real * DenMat0
                    StdI.transindx[StdI.ntrans, 0] = isite
                    StdI.transindx[StdI.ntrans, 1] = ispin
                    StdI.transindx[StdI.ntrans, 2] = isite
                    StdI.transindx[StdI.ntrans, 3] = ispin
                    StdI.ntrans += 1

                    # U_{Rij} D_{0ii} (Local)
                    DenMat0 = DenMat[(0, 0, 0)][
                        int(tUJindx[1][it, 3]), int(tUJindx[1][it, 3])
                    ]
                    StdI.trans[StdI.ntrans] = tUJ[1][it].real * DenMat0
                    StdI.transindx[StdI.ntrans, 0] = jsite
                    StdI.transindx[StdI.ntrans, 1] = ispin
                    StdI.transindx[StdI.ntrans, 2] = jsite
                    StdI.transindx[StdI.ntrans, 3] = ispin
                    StdI.ntrans += 1

                # Hartree-Fock correction
                if idcmode == _DCMode.FULL:
                    key = (
                        int(tUJindx[1][it, 0]),
                        int(tUJindx[1][it, 1]),
                        int(tUJindx[1][it, 2]),
                    )
                    DenMat0 = DenMat[key][
                        int(tUJindx[1][it, 3]), int(tUJindx[1][it, 4])
                    ]
                    hopping(
                        StdI,
                        -0.5 * Cphase * tUJ[1][it].real * DenMat0,
                        jsite, isite, dR,
                    )


def _apply_hund_terms(
    StdI: StdIntList,
    kCell: int,
    cell_w: int,
    cell_l: int,
    iH: int,
    NtUJ: list[int],
    tUJ: list,
    tUJindx: list,
    idcmode: _DCMode,
    DenMat: dict[tuple[int, int, int], np.ndarray] | None,
) -> None:
    """Apply Hund (J) coupling terms for one unit cell.

    Processes all Hund coupling terms for the unit cell at position
    ``(cell_w, cell_l, iH)``, including exchange, pair-hopping (Hubbard), and
    double-counting corrections when enabled.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters.  Modified in-place.
    kCell : int
        Linear index of the current unit cell.
    cell_w, cell_l, iH : int
        Unit-cell coordinates.
    NtUJ : list of int
        Number of terms per interaction type.
    tUJ : list
        Matrix elements per interaction type.
    tUJindx : list
        Indices per interaction type.
    idcmode : _DCMode
        Double-counting correction mode.
    DenMat : dict or None
        Density matrix elements keyed by R-vector tuples.
    """
    if tUJindx[2] is None:
        return
    for it in range(NtUJ[2]):
        # Local term should not be computed
        if (tUJindx[2][it, 0] != 0 or tUJindx[2][it, 1] != 0
                or tUJindx[2][it, 2] != 0
                or tUJindx[2][it, 3] != tUJindx[2][it, 4]):
            isite, jsite, Cphase, dR = find_site(
                StdI, cell_w, cell_l, iH,
                int(tUJindx[2][it, 0]), int(tUJindx[2][it, 1]),
                int(tUJindx[2][it, 2]),
                int(tUJindx[2][it, 3]), int(tUJindx[2][it, 4]),
            )

            StdI.Hund[StdI.NHund] = tUJ[2][it].real
            StdI.HundIndx[StdI.NHund, 0] = isite
            StdI.HundIndx[StdI.NHund, 1] = jsite
            StdI.NHund += 1

            if StdI.model == ModelType.HUBBARD:
                StdI.Ex[StdI.NEx] = tUJ[2][it].real
                StdI.ExIndx[StdI.NEx, 0] = isite
                StdI.ExIndx[StdI.NEx, 1] = jsite
                StdI.NEx += 1

                StdI.PairHopp[StdI.NPairHopp] = tUJ[2][it].real
                StdI.PHIndx[StdI.NPairHopp, 0] = isite
                StdI.PHIndx[StdI.NPairHopp, 1] = jsite
                StdI.NPairHopp += 1

                # Double-counting correction
                if idcmode != _DCMode.NOTCORRECT and idcmode != _DCMode.HARTREE_U:
                    for ispin in range(2):
                        # -0.5 J_{Rij} D_{0jj}
                        DenMat0 = DenMat[(0, 0, 0)][
                            int(tUJindx[2][it, 4]), int(tUJindx[2][it, 4])
                        ]
                        StdI.trans[StdI.ntrans] = (
                            -(1.0 - StdI.alpha) * tUJ[2][it].real * DenMat0
                        )
                        StdI.transindx[StdI.ntrans, 0] = isite
                        StdI.transindx[StdI.ntrans, 1] = ispin
                        StdI.transindx[StdI.ntrans, 2] = isite
                        StdI.transindx[StdI.ntrans, 3] = ispin
                        StdI.ntrans += 1

                        # -0.5 J_{Rij} D_{0ii}
                        DenMat0 = DenMat[(0, 0, 0)][
                            int(tUJindx[2][it, 3]), int(tUJindx[2][it, 3])
                        ]
                        StdI.trans[StdI.ntrans] = (
                            -(1.0 - StdI.alpha) * tUJ[2][it].real * DenMat0
                        )
                        StdI.transindx[StdI.ntrans, 0] = jsite
                        StdI.transindx[StdI.ntrans, 1] = ispin
                        StdI.transindx[StdI.ntrans, 2] = jsite
                        StdI.transindx[StdI.ntrans, 3] = ispin
                        StdI.ntrans += 1

                    # Hartree-Fock correction
                    if idcmode == _DCMode.FULL:
                        key = (
                            int(tUJindx[2][it, 0]),
                            int(tUJindx[2][it, 1]),
                            int(tUJindx[2][it, 2]),
                        )
                        DenMat0 = DenMat[key][
                            int(tUJindx[2][it, 3]), int(tUJindx[2][it, 4])
                        ]
                        hopping(
                            StdI,
                            0.5 * Cphase * tUJ[2][it].real
                            * (DenMat0 + 2.0 * DenMat0.real),
                            jsite, isite, dR,
                        )
            else:
                # spin model
                if StdI.solver == SolverType.mVMC:
                    StdI.Ex[StdI.NEx] = tUJ[2][it].real
                else:
                    StdI.Ex[StdI.NEx] = -tUJ[2][it].real
                StdI.ExIndx[StdI.NEx, 0] = isite
                StdI.ExIndx[StdI.NEx, 1] = jsite
                StdI.NEx += 1


# ---------------------------------------------------------------------------
#  High-level helpers (called by wannier90)
# ---------------------------------------------------------------------------


def _validate_wannier_params(StdI: StdIntList) -> None:
    """Check and store Hamiltonian parameters for the Wannier90 lattice.

    Validates model-specific parameters (``S2`` for spin, ``mu`` for Hubbard)
    and reports unused parameters.  Exits on unsupported Kondo model.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters.  Modified in-place.
    """
    not_used_d("K", StdI.K)
    StdI.h = print_val_d("h", StdI.h, 0.0)
    StdI.Gamma = print_val_d("Gamma", StdI.Gamma, 0.0)
    StdI.Gamma_y = print_val_d("Gamma_y", StdI.Gamma_y, 0.0)
    not_used_d("U", StdI.U)

    if StdI.model == ModelType.SPIN:
        StdI.S2 = print_val_i("2S", StdI.S2, 1)
    elif StdI.model == ModelType.HUBBARD:
        StdI.mu = print_val_d("mu", StdI.mu, 0.0)
    else:
        msg = "wannier + Kondo is not available !"
        logger.error(msg)
        raise ValueError(msg)


def _build_wannier_interactions(
    StdI: StdIntList,
    NtUJ: list[int],
    tUJ: list,
    tUJindx: list,
    idcmode: _DCMode,
    DenMat: dict | None,
) -> None:
    """Allocate interaction arrays and populate transfer / interaction terms.

    Computes upper bounds for transfer and interaction arrays, allocates
    memory via :func:`malloc_interactions`, then loops over super-cells to
    apply hopping, Coulomb, and Hund terms.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice.  Modified in-place.
    NtUJ : list of int
        Number of hopping / Coulomb / Hund terms (length 3).
    tUJ : list
        Complex coefficient arrays for hopping / Coulomb / Hund (length 3).
    tUJindx : list
        Index arrays for hopping / Coulomb / Hund (length 3).
    idcmode : _DCMode
        Double-counting mode enum.
    DenMat : dict or None
        Density matrix (keyed by ``(R0, R1, R2)``), or ``None``.
    """
    # Compute upper limits for Transfer & Interaction arrays
    if StdI.model == ModelType.SPIN:
        ntransMax = StdI.nsite * (StdI.S2 + 1 + 2 * StdI.S2)
        nintrMax = StdI.NCell * (
            StdI.NsiteUC + NtUJ[0] + NtUJ[1] + NtUJ[2]
        ) * (3 * StdI.S2 + 1) * (3 * StdI.S2 + StdI.NsiteUC)
    elif StdI.model == ModelType.HUBBARD:
        ntransMax = StdI.NCell * 2 * (
            2 * StdI.NsiteUC + NtUJ[0] * 2
            + NtUJ[1] * 2 * 3 + NtUJ[2] * 2 * 2
        )
        nintrMax = StdI.NCell * (NtUJ[1] + NtUJ[2] + StdI.NsiteUC)

    malloc_interactions(StdI, ntransMax, nintrMax)

    # For spin systems, compute super-exchange interaction on-site U
    Uspin = None
    if StdI.model == ModelType.SPIN:
        Uspin = np.zeros(StdI.NsiteUC)
        if tUJindx[1] is not None:
            for it in range(NtUJ[1]):
                if (tUJindx[1][it, 0] == 0 and tUJindx[1][it, 1] == 0
                        and tUJindx[1][it, 2] == 0
                        and tUJindx[1][it, 3] == tUJindx[1][it, 4]):
                    Uspin[int(tUJindx[1][it, 3])] = tUJ[1][it].real

    # Main cell loop — apply all interaction terms
    for kCell in range(StdI.NCell):
        cell_w = StdI.Cell[kCell, 0]
        cell_l = StdI.Cell[kCell, 1]
        iH = StdI.Cell[kCell, 2]

        # Local term
        if StdI.model == ModelType.SPIN:
            for isite in range(StdI.NsiteUC * kCell, StdI.NsiteUC * (kCell + 1)):
                mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, isite)
        else:
            for isite in range(StdI.NsiteUC * kCell, StdI.NsiteUC * (kCell + 1)):
                hubbard_local(StdI, StdI.mu, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, 0.0, isite)

        # Hopping
        _apply_hopping_terms(StdI, kCell, cell_w, cell_l, iH, NtUJ, tUJ, tUJindx, Uspin)

        # Coulomb integral (U)
        _apply_coulomb_terms(StdI, kCell, cell_w, cell_l, iH, NtUJ, tUJ, tUJindx, idcmode, DenMat)

        # Hund coupling (J)
        _apply_hund_terms(StdI, kCell, cell_w, cell_l, iH, NtUJ, tUJ, tUJindx, idcmode, DenMat)


def _write_wan2site(StdI: StdIntList) -> None:
    """Write ``wan2site.dat`` mapping Wannier orbitals to super-cell sites.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing lattice and cell information.
    """
    with open("wan2site.dat", "w") as fp:
        fp.write("======================== \n")
        fp.write(f"Total site number {StdI.NCell * StdI.NsiteUC:7d}  \n")
        fp.write("======================== \n")
        fp.write("========site nx ny nz norb====== \n")
        fp.write("======================== \n")

        for kCell in range(StdI.NCell):
            nx = StdI.Cell[kCell, 0]
            ny = StdI.Cell[kCell, 1]
            nz = StdI.Cell[kCell, 2]
            for it in range(StdI.NsiteUC):
                isite = StdI.NsiteUC * kCell + it
                fp.write(f"{isite:5d}{nx:5d}{ny:5d}{nz:5d}{it:5d}\n")


def _validate_interaction_params(StdI: StdIntList) -> None:
    """Validate and default lambda/alpha interaction-strength parameters.

    Sets ``lambda_U``, ``lambda_J``, and ``alpha`` on *StdI*, using
    ``lambda_`` as a shared default when it is not NaN.  Exits if
    any value is out of range.

    Parameters
    ----------
    StdI : StdIntList
        Parameter structure (modified in place).

    Raises
    ------
    ValueError
        If ``lambda_U`` or ``lambda_J`` is negative, or ``alpha`` is
        outside [0, 1].
    """
    if math.isnan(StdI.lambda_):
        StdI.lambda_U = print_val_d("lambda_U", StdI.lambda_U, 1.0)
        StdI.lambda_J = print_val_d("lambda_J", StdI.lambda_J, 1.0)
    else:
        StdI.lambda_U = print_val_d("lambda_U", StdI.lambda_U, StdI.lambda_)
        StdI.lambda_J = print_val_d("lambda_J", StdI.lambda_J, StdI.lambda_)

    if StdI.lambda_U < 0.0 or StdI.lambda_J < 0.0:
        msg = (
            "the value of lambda_U / lambda_J must be "
            "greater than or equal to 0."
        )
        logger.error(msg)
        raise ValueError(msg)

    StdI.alpha = print_val_d("alpha", StdI.alpha, 0.5)
    if StdI.alpha > 1.0 or StdI.alpha < 0.0:
        msg = "the value of alpha must be in the range 0 <= alpha <= 1."
        logger.error(msg)
        raise ValueError(msg)


# ---------------------------------------------------------------------------
#  Main entry point
# ---------------------------------------------------------------------------


class _W90Channel(NamedTuple):
    """Configuration for one Wannier90 interaction channel (hopping/Coulomb/Hund)."""

    label: str
    key_lower: str
    key_upper: str
    file_suffix: str
    cutoff_length_default: float
    cutoff_R_defaults: tuple[int | None, ...]
    itUJ: int
    lam: float


# Mapping from channel key_lower to (cutoff_attr, cutoff_length_attr, cutoffR_attr, cutoffVec_attr)
_W90_FIELD_MAP: dict[str, tuple[str, str, str, str]] = {
    "t": ("cutoff_t", "cutoff_length_t", "cutoff_tR", "cutoff_tVec"),
    "u": ("cutoff_u", "cutoff_length_U", "cutoff_UR", "cutoff_UVec"),
    "j": ("cutoff_j", "cutoff_length_J", "cutoff_JR", "cutoff_JVec"),
}


def _read_w90_channels(
    StdI: StdIntList,
    channels: tuple[_W90Channel, ...],
    NtUJ: int,
    tUJindx: np.ndarray,
    tUJ: np.ndarray,
) -> None:
    """Read Wannier90 interaction files for all channels.

    Parameters
    ----------
    StdI : StdIntList
        Standard interface data; cutoff fields are updated in place.
    channels : tuple of _W90Channel
        Channel configurations (hopping, Coulomb, Hund).
    NtUJ : int
        Number of interaction entries.
    tUJindx : numpy.ndarray
        Interaction index array.
    tUJ : numpy.ndarray
        Interaction value array.
    """
    for ch in channels:
        logger.info(f"\n  @ Wannier90 {ch.label} \n")
        co_attr, cl_attr, cr_attr, cv_attr = _W90_FIELD_MAP[ch.key_lower]
        cutoff, cutoff_length = _read_w90_with_cutoff(
            StdI, ch.key_lower, ch.key_upper, ch.file_suffix,
            getattr(StdI, co_attr), getattr(StdI, cl_attr),
            getattr(StdI, cr_attr), getattr(StdI, cv_attr),
            cutoff_length_default=ch.cutoff_length_default,
            cutoff_R_defaults=ch.cutoff_R_defaults,
            itUJ=ch.itUJ, NtUJ=NtUJ, tUJindx=tUJindx, lam=ch.lam, tUJ=tUJ,
        )
        setattr(StdI, co_attr, cutoff)
        setattr(StdI, cl_attr, cutoff_length)


def wannier90(StdI: StdIntList) -> None:
    """Set up a Hamiltonian for the Wannier90 ``*_hr.dat``.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.

    Notes
    -----
    This function performs the following steps:

    1. Compute the shape of the super-cell and sites in the super-cell.
    2. Read Wannier90 geometry, hopping, Coulomb and Hund files.
    3. Validate and store Hamiltonian parameters.
    4. Set local spin flags and number of sites.
    5. Allocate memory for interactions.
    6. Set up transfers and interactions between sites.
    7. Write ``lattice.xsf``, ``geometry.dat`` and ``wan2site.dat``.
    """
    NtUJ = [0, 0, 0]
    tUJ: list = [None, None, None]
    tUJindx: list = [None, None, None]

    # (1) Compute the shape of the super-cell and sites in the super-cell
    with open("lattice.xsf", "w") as fp_xsf:
        StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
        StdI.phase[1] = print_val_d("phase1", StdI.phase[1], 0.0)
        StdI.phase[2] = print_val_d("phase2", StdI.phase[2], 0.0)
        StdI.NsiteUC = 1
        init_site(StdI, fp_xsf, 3)
    logger.info("\n  @ Wannier90 Geometry \n")
    _geometry_w90(StdI)

    _validate_interaction_params(StdI)
    idcmode = _parse_double_counting_mode(StdI.double_counting_mode)

    # Read hopping, Coulomb, and Hund interaction files
    hopping_R_defaults = (
        (StdI.W - 1) // 2 if StdI.W != NaN_i else None,
        (StdI.L - 1) // 2 if StdI.L != NaN_i else None,
        (StdI.Height - 1) // 2 if StdI.Height != NaN_i else None,
    )
    _W90_CHANNELS = (
        _W90Channel("hopping", "t", "t", "_hr.dat", -1.0, hopping_R_defaults, 0, 1.0),
        _W90Channel("Coulomb", "u", "U", "_ur.dat", 0.3, (0, 0, 0), 1, StdI.lambda_U),
        _W90Channel("Hund", "j", "J", "_jr.dat", 0.3, (0, 0, 0), 2, StdI.lambda_J),
    )
    _read_w90_channels(StdI, _W90_CHANNELS, NtUJ, tUJindx, tUJ)

    # Read Density matrix
    DenMat = None
    if idcmode != _DCMode.NOTCORRECT:
        logger.info("\n  @ Wannier90 Density-matrix \n")
        filename = f"{StdI.CDataFileHead}_dr.dat"
        DenMat = _read_density_matrix(StdI, filename)

    # (2) Check and store parameters of Hamiltonian
    logger.info("\n  @ Hamiltonian \n")
    _validate_wannier_params(StdI)

    logger.info("\n  @ Numerical conditions\n")

    # (3) Set local spin flag and number of sites
    set_local_spin_flags(StdI, StdI.NsiteUC * StdI.NCell)

    # (4)-(5) Allocate arrays and populate transfer / interaction terms
    _build_wannier_interactions(StdI, NtUJ, tUJ, tUJindx, idcmode, DenMat)

    if idcmode != _DCMode.NOTCORRECT:
        _print_uhf_initial(StdI, NtUJ, tUJ, DenMat, tUJindx)
    print_xsf(StdI)
    print_geometry(StdI)

    # Write wan2site.dat
    _write_wan2site(StdI)


# ---------------------------------------------------------------------------
#  Lattice plugin registration
# ---------------------------------------------------------------------------

from . import LatticePlugin, register_lattice

_wannier90_setup = wannier90


class Wannier90Plugin(LatticePlugin):
    """Plugin for the Wannier90 interface."""

    @property
    def name(self) -> str:
        return "wannier90"

    @property
    def aliases(self) -> list[str]:
        return ["wannier90"]

    @property
    def ndim(self) -> int:
        return 3

    def setup(self, StdI: StdIntList) -> None:
        """Delegate to the wannier90() function."""
        _wannier90_setup(StdI)


register_lattice(Wannier90Plugin())
