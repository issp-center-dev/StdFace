"""Input-file readers for the wannier90 lattice.

Pure reading layer for the Wannier90 / RESPACK data files consumed by
:mod:`stdface.lattice.wannier90`:

- ``_geometry_w90``       -- ``*_geom.dat`` (lattice vectors, Wannier centres)
- ``_read_w90``           -- ``*_hr.dat`` / ``*_ur.dat`` / ``*_jr.dat``
- ``_read_density_matrix``-- ``*_dr.dat``

File names are resolved against ``StdI.input_dir`` (see ``_input_path``).
Split out of ``wannier90.py`` (F-1).

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
import logging
import os
from typing import Iterator, TextIO

import numpy as np

from ..core.stdface_vals import StdIntList, NaN_i

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


def _records(fp: TextIO) -> Iterator[list[str]]:
    """Yield the whitespace-separated fields of each non-blank line.

    Blank and whitespace-only lines are not records (C reader, #72).
    CRLF line ends are handled by text-mode universal newlines.
    """
    for line in fp:
        fields = line.split()
        if fields:
            yield fields


def _input_error(filename: str, reason: str) -> ValueError:
    """Log and build the error for a malformed input file."""
    msg = f"Wannier90 input {filename}: {reason}"
    logger.error(msg)
    return ValueError(msg)


def _next_record(records: Iterator[list[str]], filename: str) -> list[str]:
    """Return the next record, or raise ``ValueError`` at end of file."""
    fields = next(records, None)
    if fields is None:
        raise _input_error(filename, "unexpected end of file")
    return fields


def _read_matrix_file(
    fp: TextIO,
    filename: str,
    nsite: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Read a Wannier90-format matrix file (``*_hr/_ur/_jr/_dr.dat``).

    The layout is a comment line, ``num_wann``, ``nrpts``, the degeneracy
    list (any number of lines), then ``nrpts`` blocks of
    ``num_wann * num_wann`` records ``R0 R1 R2 i j re im``.  Blocks are
    contiguous, but the records inside a block may come in any order
    (Wannier90 writes ``i`` fastest, RESPACK ``j``).  Only the leading
    *nsite* orbitals (the geometry's count) are kept.

    Parameters
    ----------
    fp : TextIO
        Open file positioned at the comment line.
    filename : str
        File name used in diagnostics.
    nsite : int
        Number of orbitals selected by the geometry file.

    Returns
    -------
    indx_tot : numpy.ndarray
        R vector of each block, shape ``(nWSC, 3)``, dtype int.
    Mat_tot : numpy.ndarray
        Matrix elements indexed by the records' own orbital numbers,
        shape ``(nWSC, nsite, nsite)``, dtype complex.
    """
    fp.readline()  # comment line
    records = _records(fp)
    nWan = int(_next_record(records, filename)[0])
    if nWan < nsite:
        raise _input_error(
            filename,
            f"file has fewer orbitals ({nWan}) than geometry ({nsite})")
    nWSC = int(_next_record(records, filename)[0])

    # Degeneracy weights are not used (boundary weights are recomputed).
    count = 0
    while count < nWSC:
        count += len(_next_record(records, filename))

    Mat_tot = np.zeros((nWSC, nsite, nsite), dtype=complex)
    indx_tot = np.zeros((nWSC, 3), dtype=int)
    for iWSC in range(nWSC):
        for _ in range(nWan * nWan):
            vals = _next_record(records, filename)
            if len(vals) < 7:
                raise _input_error(filename, "incomplete matrix record")
            indx_tot[iWSC, :] = [int(vals[0]), int(vals[1]), int(vals[2])]
            iWan0 = int(vals[3])
            jWan0 = int(vals[4])
            if iWan0 <= nsite and jWan0 <= nsite:
                Mat_tot[iWSC, iWan0 - 1, jWan0 - 1] = (
                    float(vals[5]) + 1j * float(vals[6])
                )
    return indx_tot, Mat_tot


def _input_path(StdI: StdIntList, fname: str) -> str:
    """Resolve an auxiliary input file against ``StdI.input_dir``.

    ``input_dir`` is set by ``generate()`` to the caller's cwd before it
    chdirs into the output directory; ``None`` (the CLI flow) keeps the
    plain relative name.
    """
    return os.path.join(StdI.input_dir, fname) if StdI.input_dir else fname


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
    filename = _input_path(StdI, f"{StdI.CDataFileHead}_geom.dat")
    logger.info(f"    Wannier90 Geometry file = {filename}")

    try:
        fp_geom = open(filename, "r")
    except OSError as exc:
        logger.error("Fail to open the file %s", filename)
        raise FileNotFoundError(filename) from exc

    with fp_geom:
        records = _records(fp_geom)
        # Read direct lattice vectors
        for ii in range(3):
            StdI.direct[ii, :] = [
                float(x) for x in _next_record(records, filename)[:3]
            ]

        # Read number of correlated sites
        StdI.NsiteUC = int(_next_record(records, filename)[0])
        logger.info(f"    Number of Correlated Sites = {StdI.NsiteUC}")

        # A count below the matrix files' num_wann selects leading orbitals
        StdI.tau = np.zeros((StdI.NsiteUC, 3))
        for isite in range(StdI.NsiteUC):
            StdI.tau[isite, :] = [
                float(x) for x in _next_record(records, filename)[:3]
            ]

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

    if StdI.W is not None and StdI.L is not None and StdI.Height is not None:
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
        Matrix elements, shape ``(nWSC, NsiteUC, NsiteUC)``, dtype complex.
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
    filename = _input_path(StdI, filename)
    try:
        fp_hr = open(filename, "r")
    except FileNotFoundError:
        logger.info(f"\n  Skip to read the file {filename}. \n")
        return

    with fp_hr:
        indx_tot, Mat_tot = _read_matrix_file(fp_hr, filename, StdI.NsiteUC)
    nWSC = indx_tot.shape[0]
    Weight_tot = np.ones(nWSC)

    if flg_vec:
        inverse_rvec = np.linalg.inv(cutoff_Rvec.astype(float))

    for iWSC in range(nWSC):
        for iWan in range(StdI.NsiteUC):
            for jWan in range(StdI.NsiteUC):
                # Python complex, so lam * dtmp keeps the pre-#72 arithmetic
                dtmp = complex(Mat_tot[iWSC, iWan, jWan])
                # Distance between the record's own orbitals (C: #72)
                tau_diff = StdI.tau[jWan, :] - StdI.tau[iWan, :] + indx_tot[iWSC, :]
                dR = StdI.direct.T @ tau_diff
                length = np.linalg.norm(dR)
                if length > cutoff_length > 0.0:
                    dtmp = 0j

                if flg_vec:
                    if not _check_in_box(indx_tot[iWSC], inverse_rvec):
                        dtmp = 0j
                else:
                    if np.any(np.abs(indx_tot[iWSC]) > cutoff_R):
                        dtmp = 0j

                Mat_tot[iWSC, iWan, jWan] = lam * dtmp

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


class _DensityMatrix(dict):
    """Density matrix per R vector.

    An R vector that the density file does not list has zero density,
    whether or not it lies inside the R range of the file.
    """

    def __init__(self, nsite: int) -> None:
        super().__init__()
        self._zero = np.zeros((nsite, nsite), dtype=complex)
        self._zero.setflags(write=False)

    def __missing__(self, key: tuple[int, int, int]) -> np.ndarray:
        return self._zero


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
        density matrix elements.  Looking up an R vector the file does
        not list gives a zero matrix.
    """


    filename = _input_path(StdI, filename)
    try:
        fp_dr = open(filename, "r")
    except OSError as exc:
        logger.error("Fail to open the file %s", filename)
        raise FileNotFoundError(filename) from exc

    with fp_dr:
        indx_tot, Mat_tot = _read_matrix_file(fp_dr, filename, StdI.NsiteUC)
    nWSC = indx_tot.shape[0]

    Rmin = indx_tot.min(axis=0, initial=0)
    Rmax = indx_tot.max(axis=0, initial=0)
    NR = Rmax - Rmin + 1
    logger.info(f"      Minimum R : {Rmin[0]} {Rmin[1]} {Rmin[2]}")
    logger.info(f"      Maximum R : {Rmax[0]} {Rmax[1]} {Rmax[2]}")
    logger.info(f"      Numver of R : {NR[0]} {NR[1]} {NR[2]}")

    # Build dictionary: (R0, R1, R2) -> 2D array
    DenMat = _DensityMatrix(StdI.NsiteUC)
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
        DenMat[key][:, :] = Mat_tot[iWSC]

    return DenMat
