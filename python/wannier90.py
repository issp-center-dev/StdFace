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

import math
from enum import IntEnum
from typing import TextIO

import numpy as np

from stdface_vals import StdIntList
from stdface_model_util import (
    exit_program,
    print_val_d,
    print_val_i,
    not_used_d,
    init_site,
    find_site,
    malloc_interactions,
    mag_field,
    general_j,
    hubbard_local,
    hopping,
    coulomb,
    print_geometry,
    print_xsf,
)


# ---------------------------------------------------------------------------
#  Internal helpers
# ---------------------------------------------------------------------------


def _calc_inverse_matrix(cutoff_rvec: np.ndarray) -> np.ndarray:
    """Calculate inverse of a 3x3 matrix.

    Parameters
    ----------
    cutoff_rvec : numpy.ndarray
        Input 3x3 matrix (shape ``(3, 3)``).

    Returns
    -------
    numpy.ndarray
        Inverse matrix (shape ``(3, 3)``).
    """
    N = cutoff_rvec.copy().astype(float)

    det = (N[0, 0] * N[1, 1] * N[2, 2]
           + N[1, 0] * N[2, 1] * N[0, 2]
           + N[2, 0] * N[0, 1] * N[1, 2]
           - N[2, 0] * N[1, 1] * N[0, 2]
           - N[1, 0] * N[0, 1] * N[2, 2]
           - N[0, 0] * N[2, 1] * N[1, 2])

    inv = np.zeros((3, 3))
    inv[0, 0] = N[1, 1] * N[2, 2] - N[1, 2] * N[2, 1]
    inv[0, 1] = -(N[0, 1] * N[2, 2] - N[0, 2] * N[2, 1])
    inv[0, 2] = N[0, 1] * N[1, 2] - N[0, 2] * N[1, 1]

    inv[1, 0] = -(N[1, 0] * N[2, 2] - N[2, 0] * N[1, 2])
    inv[1, 1] = N[0, 0] * N[2, 2] - N[0, 2] * N[2, 0]
    inv[1, 2] = -(N[0, 0] * N[1, 2] - N[0, 2] * N[1, 0])

    inv[2, 0] = N[1, 0] * N[2, 1] - N[2, 0] * N[1, 1]
    inv[2, 1] = -(N[0, 0] * N[2, 1] - N[2, 0] * N[0, 1])
    inv[2, 2] = N[0, 0] * N[1, 1] - N[0, 1] * N[1, 0]

    inv /= det
    return inv


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
    judge_vec = np.zeros(3)
    for i in range(3):
        for j in range(3):
            judge_vec[i] += rvec[j] * inverse_matrix[j, i]
    return bool(
        abs(judge_vec[0]) <= 1
        and abs(judge_vec[1]) <= 1
        and abs(judge_vec[2]) <= 1
    )


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
    print(f"    Wannier90 Geometry file = {filename}")

    try:
        fp = open(filename, "r")
    except FileNotFoundError:
        import sys
        print(f"\n  Error: Fail to open the file {filename}. \n", file=sys.stderr)
        exit_program(-1)

    # Read direct lattice vectors
    for ii in range(3):
        vals = fp.readline().split()
        StdI.direct[ii, 0] = float(vals[0])
        StdI.direct[ii, 1] = float(vals[1])
        StdI.direct[ii, 2] = float(vals[2])

    # Read number of correlated sites
    StdI.NsiteUC = int(fp.readline().split()[0])
    print(f"    Number of Correlated Sites = {StdI.NsiteUC}")

    # Allocate and read Wannier centre positions
    StdI.tau = np.zeros((StdI.NsiteUC, 3))
    for isite in range(StdI.NsiteUC):
        vals = fp.readline().split()
        StdI.tau[isite, 0] = float(vals[0])
        StdI.tau[isite, 1] = float(vals[1])
        StdI.tau[isite, 2] = float(vals[2])

    fp.close()

    print("    Direct lattice vectors:")
    for ii in range(3):
        print(f"      {StdI.direct[ii, 0]:10.5f} {StdI.direct[ii, 1]:10.5f} {StdI.direct[ii, 2]:10.5f}")
    print("    Wannier centres:")
    for isite in range(StdI.NsiteUC):
        print(f"      {StdI.tau[isite, 0]:10.5f} {StdI.tau[isite, 1]:10.5f} {StdI.tau[isite, 2]:10.5f}")


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
    NaN_i = StdI.NaN_i
    flg_vec = int(cutoff_Rvec[0, 0] != NaN_i)

    # Try to open the file
    try:
        fp = open(filename, "r")
    except FileNotFoundError:
        print(f"\n  Skip to read the file {filename}. \n")
        return

    # Header part
    _header_line = fp.readline()  # comment line
    nWan = int(fp.readline().split()[0])
    nWSC = int(fp.readline().split()[0])

    # Read degeneracy weights (skip them, only needed for count)
    count = 0
    while count < nWSC:
        line = fp.readline().split()
        count += len(line)

    # Allocate arrays
    Weight_tot = np.ones(nWSC)
    Band_lattice = np.zeros(3, dtype=int)
    Model_lattice = np.zeros(3, dtype=int)

    Mat_tot = np.zeros((nWSC, nWan, nWan), dtype=complex)
    indx_tot = np.zeros((nWSC, 3), dtype=int)

    if flg_vec:
        inverse_rvec = _calc_inverse_matrix(cutoff_Rvec)

    # Read body
    for iWSC in range(nWSC):
        for iWan in range(nWan):
            for jWan in range(nWan):
                vals = fp.readline().split()
                indx_tot[iWSC, 0] = int(vals[0])
                indx_tot[iWSC, 1] = int(vals[1])
                indx_tot[iWSC, 2] = int(vals[2])
                iWan0 = int(vals[3])
                jWan0 = int(vals[4])
                dtmp_re = float(vals[5])
                dtmp_im = float(vals[6])

                # Compute Euclidean length
                dR = np.zeros(3)
                for ii in range(3):
                    for jj in range(3):
                        dR[ii] += StdI.direct[jj, ii] * (
                            StdI.tau[jWan, jj] - StdI.tau[iWan, jj] + indx_tot[iWSC, jj]
                        )
                length = math.sqrt(dR[0] ** 2 + dR[1] ** 2 + dR[2] ** 2)
                if length > cutoff_length > 0.0:
                    dtmp_re = 0.0
                    dtmp_im = 0.0

                if flg_vec:
                    if not _check_in_box(indx_tot[iWSC], inverse_rvec):
                        dtmp_re = 0.0
                        dtmp_im = 0.0
                else:
                    if (abs(indx_tot[iWSC, 0]) > cutoff_R[0]
                            or abs(indx_tot[iWSC, 1]) > cutoff_R[1]
                            or abs(indx_tot[iWSC, 2]) > cutoff_R[2]):
                        dtmp_re = 0.0
                        dtmp_im = 0.0

                if iWan0 <= StdI.NsiteUC and jWan0 <= StdI.NsiteUC:
                    Mat_tot[iWSC, iWan0 - 1, jWan0 - 1] = lam * (dtmp_re + 1j * dtmp_im)

        # Apply inversion symmetry and delete duplication
        for jWSC in range(iWSC):
            if (indx_tot[iWSC, 0] == -indx_tot[jWSC, 0]
                    and indx_tot[iWSC, 1] == -indx_tot[jWSC, 1]
                    and indx_tot[iWSC, 2] == -indx_tot[jWSC, 2]):
                for iWan in range(StdI.NsiteUC):
                    for jWan in range(StdI.NsiteUC):
                        Mat_tot[iWSC, iWan, jWan] = 0.0

        if (indx_tot[iWSC, 0] == 0
                and indx_tot[iWSC, 1] == 0
                and indx_tot[iWSC, 2] == 0):
            for iWan in range(StdI.NsiteUC):
                for jWan in range(iWan):
                    Mat_tot[iWSC, iWan, jWan] = 0.0

    fp.close()

    # Apply weight - get lattice length
    for iWSC in range(nWSC):
        for ii in range(3):
            if abs(indx_tot[iWSC, ii]) > Band_lattice[ii]:
                Band_lattice[ii] = abs(indx_tot[iWSC, ii])

    if StdI.W != NaN_i and StdI.L != NaN_i and StdI.Height != NaN_i:
        Model_lattice[0] = StdI.W // 2 if StdI.W % 2 == 0 else 0
        Model_lattice[1] = StdI.L // 2 if StdI.L % 2 == 0 else 0
        Model_lattice[2] = StdI.Height // 2 if StdI.Height % 2 == 0 else 0
        for ii in range(3):
            if Model_lattice[ii] < Band_lattice[ii] and Model_lattice[ii] != 0:
                for iWSC in range(nWSC):
                    if abs(indx_tot[iWSC, ii]) == Model_lattice[ii]:
                        Weight_tot[iWSC] *= 0.5

    # Count effective terms
    print("\n      EFFECTIVE terms:")
    print("           R0   R1   R2 band_i band_f Hamiltonian")
    NtUJ[itUJ] = 0
    for iWSC in range(nWSC):
        for iWan in range(StdI.NsiteUC):
            for jWan in range(StdI.NsiteUC):
                Mat_tot[iWSC, iWan, jWan] *= Weight_tot[iWSC]
                if cutoff < abs(Mat_tot[iWSC, iWan, jWan]):
                    print(
                        f"        {indx_tot[iWSC, 0]:5d}{indx_tot[iWSC, 1]:5d}"
                        f"{indx_tot[iWSC, 2]:5d}{iWan:5d}{jWan:5d}"
                        f"{Mat_tot[iWSC, iWan, jWan].real:12.6f}"
                        f"{Mat_tot[iWSC, iWan, jWan].imag:12.6f}"
                    )
                    NtUJ[itUJ] += 1
    print(f"      Total number of EFFECTIVE term = {NtUJ[itUJ]}")

    # Store terms
    tUJ_arr = np.zeros(NtUJ[itUJ], dtype=complex)
    tUJindx_arr = np.zeros((NtUJ[itUJ], 5), dtype=int)

    count = 0
    for iWSC in range(nWSC):
        for iWan in range(StdI.NsiteUC):
            for jWan in range(StdI.NsiteUC):
                if cutoff < abs(Mat_tot[iWSC, iWan, jWan]):
                    tUJindx_arr[count, 0] = indx_tot[iWSC, 0]
                    tUJindx_arr[count, 1] = indx_tot[iWSC, 1]
                    tUJindx_arr[count, 2] = indx_tot[iWSC, 2]
                    tUJindx_arr[count, 3] = iWan
                    tUJindx_arr[count, 4] = jWan
                    tUJ_arr[count] = Mat_tot[iWSC, iWan, jWan]
                    count += 1

    # Extend the lists to hold the results
    while len(tUJ) <= itUJ:
        tUJ.append(None)
    while len(tUJindx) <= itUJ:
        tUJindx.append(None)
    tUJ[itUJ] = tUJ_arr
    tUJindx[itUJ] = tUJindx_arr


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
    import sys

    try:
        fp = open(filename, "r")
    except FileNotFoundError:
        print(f"\n  Error: Fail to open the file {filename}. \n", file=sys.stderr)
        exit_program(-1)

    # Header
    _header_line = fp.readline()
    nWan = int(fp.readline().split()[0])
    nWSC = int(fp.readline().split()[0])

    count = 0
    while count < nWSC:
        line = fp.readline().split()
        count += len(line)

    # Allocate
    Mat_tot = np.zeros((nWSC, nWan, nWan), dtype=complex)
    indx_tot = np.zeros((nWSC, 3), dtype=int)

    Rmin = np.zeros(3, dtype=int)
    Rmax = np.zeros(3, dtype=int)

    # Read body
    for iWSC in range(nWSC):
        for iWan in range(nWan):
            for jWan in range(nWan):
                vals = fp.readline().split()
                indx_tot[iWSC, 0] = int(vals[0])
                indx_tot[iWSC, 1] = int(vals[1])
                indx_tot[iWSC, 2] = int(vals[2])
                iWan0 = int(vals[3])
                jWan0 = int(vals[4])
                dtmp_re = float(vals[5])
                dtmp_im = float(vals[6])

                if iWan0 <= StdI.NsiteUC and jWan0 <= StdI.NsiteUC:
                    Mat_tot[iWSC, iWan0 - 1, jWan0 - 1] = dtmp_re + 1j * dtmp_im
                for ii in range(3):
                    if indx_tot[iWSC, ii] < Rmin[ii]:
                        Rmin[ii] = indx_tot[iWSC, ii]
                    if indx_tot[iWSC, ii] > Rmax[ii]:
                        Rmax[ii] = indx_tot[iWSC, ii]
    fp.close()

    NR = Rmax - Rmin + 1
    print(f"      Minimum R : {Rmin[0]} {Rmin[1]} {Rmin[2]}")
    print(f"      Maximum R : {Rmax[0]} {Rmax[1]} {Rmax[2]}")
    print(f"      Numver of R : {NR[0]} {NR[1]} {NR[2]}")

    # Build dictionary: (R0, R1, R2) -> 2D array
    DenMat: dict[tuple[int, int, int], np.ndarray] = {}
    for i0 in range(Rmin[0], Rmax[0] + 1):
        for i1 in range(Rmin[1], Rmax[1] + 1):
            for i2 in range(Rmin[2], Rmax[2] + 1):
                DenMat[(i0, i1, i2)] = np.zeros(
                    (StdI.NsiteUC, StdI.NsiteUC), dtype=complex
                )

    for iWSC in range(nWSC):
        key = (int(indx_tot[iWSC, 0]), int(indx_tot[iWSC, 1]), int(indx_tot[iWSC, 2]))
        for iWan in range(nWan):
            for jWan in range(nWan):
                DenMat[key][iWan, jWan] = Mat_tot[iWSC, iWan, jWan]

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
        iW = StdI.Cell[kCell, 0]
        iL = StdI.Cell[kCell, 1]
        iH = StdI.Cell[kCell, 2]

        # Diagonal term
        for isite_uc in range(StdI.NsiteUC):
            jsite = isite_uc + StdI.NsiteUC * kCell
            IniGuess[jsite, jsite] = DenMat[(0, 0, 0)][isite_uc, isite_uc]

        # Coulomb integral (U)
        for it in range(NtUJ[1]):
            isite, jsite, Cphase, dR = find_site(
                StdI, iW, iL, iH,
                int(tUJindx[1][it, 0]), int(tUJindx[1][it, 1]), int(tUJindx[1][it, 2]),
                int(tUJindx[1][it, 3]), int(tUJindx[1][it, 4]),
            )
            key = (int(tUJindx[1][it, 0]), int(tUJindx[1][it, 1]), int(tUJindx[1][it, 2]))
            IniGuess[isite, jsite] = DenMat[key][int(tUJindx[1][it, 3]), int(tUJindx[1][it, 4])]
            IniGuess[jsite, isite] = np.conj(
                DenMat[key][int(tUJindx[1][it, 3]), int(tUJindx[1][it, 4])]
            )

        # Exchange integral (J)
        for it in range(NtUJ[2]):
            isite, jsite, Cphase, dR = find_site(
                StdI, iW, iL, iH,
                int(tUJindx[2][it, 0]), int(tUJindx[2][it, 1]), int(tUJindx[2][it, 2]),
                int(tUJindx[2][it, 3]), int(tUJindx[2][it, 4]),
            )
            key = (int(tUJindx[2][it, 0]), int(tUJindx[2][it, 1]), int(tUJindx[2][it, 2]))
            IniGuess[isite, jsite] = DenMat[key][int(tUJindx[2][it, 3]), int(tUJindx[2][it, 4])]
            IniGuess[jsite, isite] = np.conj(
                DenMat[key][int(tUJindx[2][it, 3]), int(tUJindx[2][it, 4])]
            )

    NIniGuess = 0
    for isite in range(StdI.nsite):
        for jsite in range(StdI.nsite):
            if abs(IniGuess[isite, jsite]) > 1.0e-6:
                NIniGuess += 1

    with open("initial.def", "w") as fp:
        fp.write("======================== \n")
        fp.write(f"NInitialGuess {NIniGuess * 2:7d}  \n")
        fp.write("======================== \n")
        fp.write("========i_j_s_tijs====== \n")
        fp.write("======================== \n")

        for isite in range(StdI.nsite):
            for jsite in range(StdI.nsite):
                if abs(IniGuess[isite, jsite]) > 1.0e-6:
                    for ispin in range(2):
                        fp.write(
                            f"{jsite:5d} {ispin:5d} {isite:5d} {ispin:5d} "
                            f"{0.5 * IniGuess[isite, jsite].real:25.15f} "
                            f"{0.5 * IniGuess[isite, jsite].imag:25.15f}\n"
                        )

    print("      initial.def is written.")


# ---------------------------------------------------------------------------
#  Double-counting mode enum
# ---------------------------------------------------------------------------


class _DCMode(IntEnum):
    """Double-counting correction mode."""
    NOTCORRECT = 0
    HARTREE = 1
    HARTREE_U = 2
    FULL = 3


# ---------------------------------------------------------------------------
#  Main entry point
# ---------------------------------------------------------------------------


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
    NaN_i = StdI.NaN_i
    NtUJ = [0, 0, 0]
    tUJ: list = [None, None, None]
    tUJindx: list = [None, None, None]

    # (1) Compute the shape of the super-cell and sites in the super-cell
    fp_xsf = open("lattice.xsf", "w")

    StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
    StdI.phase[1] = print_val_d("phase1", StdI.phase[1], 0.0)
    StdI.phase[2] = print_val_d("phase2", StdI.phase[2], 0.0)
    StdI.NsiteUC = 1
    init_site(StdI, fp_xsf, 3)
    print("\n  @ Wannier90 Geometry \n")
    _geometry_w90(StdI)

    # Set parameters to tune the strength of interactions
    if math.isnan(StdI.lambda_):
        # Lambda is not defined
        StdI.lambda_U = print_val_d("lambda_U", StdI.lambda_U, 1.0)
        StdI.lambda_J = print_val_d("lambda_J", StdI.lambda_J, 1.0)
    else:
        StdI.lambda_U = print_val_d("lambda_U", StdI.lambda_U, StdI.lambda_)
        StdI.lambda_J = print_val_d("lambda_J", StdI.lambda_J, StdI.lambda_)

    if StdI.lambda_U < 0.0 or StdI.lambda_J < 0.0:
        import sys
        print(
            "\n  Error: the value of lambda_U / lambda_J must be "
            "greater than or equal to 0. \n",
            file=sys.stderr,
        )
        exit_program(-1)

    # Determine double-counting mode
    dcm = StdI.double_counting_mode
    if dcm in ("none", "****"):
        idcmode = _DCMode.NOTCORRECT
    elif dcm == "hartree":
        idcmode = _DCMode.HARTREE
    elif dcm == "hartree_u":
        idcmode = _DCMode.HARTREE_U
    elif dcm == "full":
        idcmode = _DCMode.FULL
    else:
        import sys
        print(
            "\n  Error: the word of doublecounting is not correct "
            "(select from none, hartree, hartree_u, full). \n",
            file=sys.stderr,
        )
        exit_program(-1)

    StdI.alpha = print_val_d("alpha", StdI.alpha, 0.5)
    if StdI.alpha > 1.0 or StdI.alpha < 0.0:
        import sys
        print(
            "\n  Error: the value of alpha must be in the range 0<= alpha <= 1. \n",
            file=sys.stderr,
        )
        exit_program(-1)

    # Read Hopping
    print("\n  @ Wannier90 hopping \n")
    StdI.cutoff_t = print_val_d("cutoff_t", StdI.cutoff_t, 1.0e-8)
    StdI.cutoff_length_t = print_val_d("cutoff_length_t", StdI.cutoff_length_t, -1.0)
    if StdI.W != NaN_i:
        StdI.cutoff_tR[0] = print_val_i("cutoff_tR[0]", int(StdI.cutoff_tR[0]), (StdI.W - 1) // 2)
    if StdI.L != NaN_i:
        StdI.cutoff_tR[1] = print_val_i("cutoff_tR[1]", int(StdI.cutoff_tR[1]), (StdI.L - 1) // 2)
    if StdI.Height != NaN_i:
        StdI.cutoff_tR[2] = print_val_i("cutoff_tR[2]", int(StdI.cutoff_tR[2]), (StdI.Height - 1) // 2)

    for i in range(3):
        for j in range(3):
            if StdI.box[i, j] != NaN_i:
                tempwords = f"cutoff_tVec[{i}][{j}]"
                StdI.cutoff_tVec[i, j] = print_val_d(
                    tempwords, StdI.cutoff_tVec[i, j], float(StdI.box[i, j]) * 0.5
                )

    filename = f"{StdI.CDataFileHead}_hr.dat"
    _read_w90(
        StdI, filename,
        StdI.cutoff_t, StdI.cutoff_tR, StdI.cutoff_tVec, StdI.cutoff_length_t,
        0, NtUJ, tUJindx, 1.0, tUJ,
    )

    # Read Coulomb
    print("\n  @ Wannier90 Coulomb \n")
    StdI.cutoff_u = print_val_d("cutoff_u", StdI.cutoff_u, 1.0e-8)
    StdI.cutoff_length_U = print_val_d("cutoff_length_U", StdI.cutoff_length_U, 0.3)
    StdI.cutoff_UR[0] = print_val_i("cutoff_UR[0]", int(StdI.cutoff_UR[0]), 0)
    StdI.cutoff_UR[1] = print_val_i("cutoff_UR[1]", int(StdI.cutoff_UR[1]), 0)
    StdI.cutoff_UR[2] = print_val_i("cutoff_UR[2]", int(StdI.cutoff_UR[2]), 0)
    for i in range(3):
        for j in range(3):
            if StdI.box[i, j] != NaN_i:
                tempwords = f"cutoff_UVec[{i}][{j}]"
            StdI.cutoff_UVec[i, j] = print_val_d(
                tempwords, StdI.cutoff_UVec[i, j], float(StdI.box[i, j]) * 0.5
            )

    filename = f"{StdI.CDataFileHead}_ur.dat"
    _read_w90(
        StdI, filename,
        StdI.cutoff_u, StdI.cutoff_UR, StdI.cutoff_UVec, StdI.cutoff_length_U,
        1, NtUJ, tUJindx, StdI.lambda_U, tUJ,
    )

    # Read Hund
    print("\n  @ Wannier90 Hund \n")
    StdI.cutoff_j = print_val_d("cutoff_j", StdI.cutoff_j, 1.0e-8)
    StdI.cutoff_length_J = print_val_d("cutoff_length_J", StdI.cutoff_length_J, 0.3)
    StdI.cutoff_JR[0] = print_val_i("cutoff_JR[0]", int(StdI.cutoff_JR[0]), 0)
    StdI.cutoff_JR[1] = print_val_i("cutoff_JR[1]", int(StdI.cutoff_JR[1]), 0)
    StdI.cutoff_JR[2] = print_val_i("cutoff_JR[2]", int(StdI.cutoff_JR[2]), 0)
    for i in range(3):
        for j in range(3):
            if StdI.box[i, j] != NaN_i:
                tempwords = f"cutoff_JVec[{i}][{j}]"
            StdI.cutoff_JVec[i, j] = print_val_d(
                tempwords, StdI.cutoff_JVec[i, j], float(StdI.box[i, j]) * 0.5
            )

    filename = f"{StdI.CDataFileHead}_jr.dat"
    _read_w90(
        StdI, filename,
        StdI.cutoff_j, StdI.cutoff_JR, StdI.cutoff_JVec, StdI.cutoff_length_J,
        2, NtUJ, tUJindx, StdI.lambda_J, tUJ,
    )

    # Read Density matrix
    DenMat = None
    if idcmode != _DCMode.NOTCORRECT:
        print("\n  @ Wannier90 Density-matrix \n")
        filename = f"{StdI.CDataFileHead}_dr.dat"
        DenMat = _read_density_matrix(StdI, filename)

    # (2) Check and store parameters of Hamiltonian
    print("\n  @ Hamiltonian \n")
    not_used_d("K", StdI.K)
    StdI.h = print_val_d("h", StdI.h, 0.0)
    StdI.Gamma = print_val_d("Gamma", StdI.Gamma, 0.0)
    StdI.Gamma_y = print_val_d("Gamma_y", StdI.Gamma_y, 0.0)
    not_used_d("U", StdI.U)

    if StdI.model == "spin":
        StdI.S2 = print_val_i("2S", StdI.S2, 1)
    elif StdI.model == "hubbard":
        StdI.mu = print_val_d("mu", StdI.mu, 0.0)
    else:
        print("wannier + Kondo is not available !")
        exit_program(-1)

    print("\n  @ Numerical conditions\n")

    # (3) Set local spin flag and number of sites
    StdI.nsite = StdI.NsiteUC * StdI.NCell
    StdI.locspinflag = np.zeros(StdI.nsite, dtype=int)

    if StdI.model == "spin":
        for isite in range(StdI.nsite):
            StdI.locspinflag[isite] = StdI.S2
    elif StdI.model == "hubbard":
        for isite in range(StdI.nsite):
            StdI.locspinflag[isite] = 0

    # (4) Compute the upper limit of the number of Transfer & Interaction
    if StdI.model == "spin":
        ntransMax = StdI.nsite * (StdI.S2 + 1 + 2 * StdI.S2)
        nintrMax = StdI.NCell * (
            StdI.NsiteUC + NtUJ[0] + NtUJ[1] + NtUJ[2]
        ) * (3 * StdI.S2 + 1) * (3 * StdI.S2 + StdI.NsiteUC)
    elif StdI.model == "hubbard":
        ntransMax = StdI.NCell * 2 * (
            2 * StdI.NsiteUC + NtUJ[0] * 2
            + NtUJ[1] * 2 * 3 + NtUJ[2] * 2 * 2
        )
        nintrMax = StdI.NCell * (NtUJ[1] + NtUJ[2] + StdI.NsiteUC)

    malloc_interactions(StdI, ntransMax, nintrMax)

    # (4.5) For spin system, compute super exchange interaction
    Uspin = None
    if StdI.model == "spin":
        Uspin = np.zeros(StdI.NsiteUC)
        if tUJindx[1] is not None:
            for it in range(NtUJ[1]):
                if (tUJindx[1][it, 0] == 0 and tUJindx[1][it, 1] == 0
                        and tUJindx[1][it, 2] == 0
                        and tUJindx[1][it, 3] == tUJindx[1][it, 4]):
                    Uspin[int(tUJindx[1][it, 3])] = tUJ[1][it].real

    # (5) Set Transfer & Interaction
    for kCell in range(StdI.NCell):
        iW = StdI.Cell[kCell, 0]
        iL = StdI.Cell[kCell, 1]
        iH = StdI.Cell[kCell, 2]

        # Local term 1
        if StdI.model == "spin":
            for isite in range(StdI.NsiteUC * kCell, StdI.NsiteUC * (kCell + 1)):
                mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, isite)
        else:
            for isite in range(StdI.NsiteUC * kCell, StdI.NsiteUC * (kCell + 1)):
                hubbard_local(StdI, StdI.mu, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, 0.0, isite)

        # Hopping
        if tUJindx[0] is not None:
            for it in range(NtUJ[0]):
                # Local term
                if (tUJindx[0][it, 0] == 0 and tUJindx[0][it, 1] == 0
                        and tUJindx[0][it, 2] == 0
                        and tUJindx[0][it, 3] == tUJindx[0][it, 4]):
                    if StdI.model == "hubbard":
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
                        StdI, iW, iL, iH,
                        int(tUJindx[0][it, 0]), int(tUJindx[0][it, 1]),
                        int(tUJindx[0][it, 2]),
                        int(tUJindx[0][it, 3]), int(tUJindx[0][it, 4]),
                    )
                    if StdI.model == "spin":
                        Jtmp = np.zeros((3, 3))
                        for ii in range(3):
                            Jtmp[ii, ii] = (
                                2.0 * tUJ[0][it] * np.conj(tUJ[0][it])
                                * (1.0 / Uspin[int(tUJindx[0][it, 3])]
                                   + 1.0 / Uspin[int(tUJindx[0][it, 4])])
                            ).real
                        general_j(StdI, Jtmp, StdI.S2, StdI.S2, isite, jsite)
                    else:
                        hopping(StdI, -Cphase * tUJ[0][it], jsite, isite, dR)

        # Coulomb integral (U)
        if tUJindx[1] is not None:
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
                        StdI, iW, iL, iH,
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

        # Hund coupling (J)
        if tUJindx[2] is not None:
            for it in range(NtUJ[2]):
                # Local term should not be computed
                if (tUJindx[2][it, 0] != 0 or tUJindx[2][it, 1] != 0
                        or tUJindx[2][it, 2] != 0
                        or tUJindx[2][it, 3] != tUJindx[2][it, 4]):
                    isite, jsite, Cphase, dR = find_site(
                        StdI, iW, iL, iH,
                        int(tUJindx[2][it, 0]), int(tUJindx[2][it, 1]),
                        int(tUJindx[2][it, 2]),
                        int(tUJindx[2][it, 3]), int(tUJindx[2][it, 4]),
                    )

                    StdI.Hund[StdI.NHund] = tUJ[2][it].real
                    StdI.HundIndx[StdI.NHund, 0] = isite
                    StdI.HundIndx[StdI.NHund, 1] = jsite
                    StdI.NHund += 1

                    if StdI.model == "hubbard":
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
                        if StdI.solver == "mVMC":
                            StdI.Ex[StdI.NEx] = tUJ[2][it].real
                        else:
                            StdI.Ex[StdI.NEx] = -tUJ[2][it].real
                        StdI.ExIndx[StdI.NEx, 0] = isite
                        StdI.ExIndx[StdI.NEx, 1] = jsite
                        StdI.NEx += 1

    fp_xsf.close()
    if idcmode != _DCMode.NOTCORRECT:
        _print_uhf_initial(StdI, NtUJ, tUJ, DenMat, tUJindx)
    print_xsf(StdI)
    print_geometry(StdI)

    # Write wan2site.dat
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
