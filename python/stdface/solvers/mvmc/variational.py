"""mVMC variational parameter generation functions.

This module contains functions for generating variational wave-function
parameter files used by the mVMC solver:

- Quantum number projection (``qptransidx.def``)
- Orbital indices (``Orb`` / ``AntiOrb`` arrays)
- Jastrow factor indices (``jastrowidx.def``)

These functions were extracted from ``stdface_model_util.py`` to improve
module cohesion.

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

from ...core.stdface_vals import StdIntList, ModelType, NaN_i
from ...core.param_check import exit_program
from ...lattice.site_util import (
    _cell_vector, _fold_to_cell, _fold_site, _find_cell_index,
    _validate_box_params, _det_and_cofactor, find_site,
)


def _anti_period_dot(AntiPeriod: np.ndarray, nBox: list[int]) -> int:
    """Compute the dot product of anti-period flags and box indices.

    Parameters
    ----------
    AntiPeriod : np.ndarray
        Length-3 array of anti-periodic boundary flags (0 or 1).
    nBox : list of int
        Length-3 super-cell image index.

    Returns
    -------
    int
        Sum ``AntiPeriod[0]*nBox[0] + AntiPeriod[1]*nBox[1] +
        AntiPeriod[2]*nBox[2]``.
    """
    return (int(AntiPeriod[0]) * nBox[0]
            + int(AntiPeriod[1]) * nBox[1]
            + int(AntiPeriod[2]) * nBox[2])


def _parity_sign(value: int) -> int:
    """Convert an integer to +1 (even) or -1 (odd).

    Used to translate an anti-periodic boundary count into a phase sign.

    Parameters
    ----------
    value : int
        Integer whose parity determines the sign.

    Returns
    -------
    int
        ``+1`` if *value* is even, ``-1`` if odd.
    """
    return 1 if value % 2 == 0 else -1


def _check_commensurate(rbox_sub: np.ndarray, box: np.ndarray, ncell_sub: int) -> bool:
    """Check whether a sublattice is commensurate with the main lattice.

    The sublattice is commensurate if for all i, j in {0,1,2}:
    ``(rbox_sub[i, :] · box[j, :]) % ncell_sub == 0``

    Parameters
    ----------
    rbox_sub : np.ndarray
        Reciprocal box matrix of the sublattice (3x3).
    box : np.ndarray
        Box matrix of the main lattice (3x3).
    ncell_sub : int
        Determinant of the sublattice box (number of cells).

    Returns
    -------
    bool
        ``True`` if the sublattice is commensurate, ``False`` otherwise.
    """
    # Compute all dot products: prod[i,j] = rbox_sub[i,:] · box[j,:]
    prod = rbox_sub.astype(int) @ box.astype(int).T
    return bool(np.all(prod % ncell_sub == 0))


def _fold_site_sub(
    StdI: StdIntList,
    iCellV: list[int],
) -> tuple[list[int], list[int]]:
    """Fold site into the sub-lattice cell (mVMC only).

    Delegates to :func:`~lattice.site_util._fold_to_cell` with the
    sub-lattice parameters (``rboxsub``, ``NCellsub``, ``boxsub``).

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.
    iCellV : list of int
        Fractional coordinate of a site (length 3).

    Returns
    -------
    nBox : list of int
        Super-cell index.
    iCellV_fold : list of int
        Folded fractional coordinate.
    """
    return _fold_to_cell(StdI.rboxsub, StdI.NCellsub, StdI.boxsub, iCellV)


def proj(StdI: StdIntList) -> None:
    """Print quantum number projection file ``qptransidx.def`` (mVMC only).

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    """
    Sym = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
    Anti = np.zeros((StdI.nsite, StdI.nsite), dtype=int)

    StdI.NSym = 0
    for iCell in range(StdI.NCell):
        iCV = _cell_vector(StdI.Cell, iCell)
        nBox, iCellV = _fold_site_sub(StdI, iCV)
        nBox, iCellV = _fold_site(StdI, iCellV)

        if iCellV == iCV:
            for jCell in range(StdI.NCell):
                jCellV = [c + v for c, v in zip(
                    _cell_vector(StdI.Cell, jCell), iCellV)]
                nBox, jCellV = _fold_site(StdI, jCellV)

                kCell = _find_cell_index(StdI, jCellV)
                ap_dot = _anti_period_dot(StdI.AntiPeriod, nBox)
                for jsite in range(StdI.NsiteUC):
                    Sym[StdI.NSym][jCell * StdI.NsiteUC + jsite] = (
                        kCell * StdI.NsiteUC + jsite)
                    Anti[StdI.NSym][jCell * StdI.NsiteUC + jsite] = ap_dot

                    if StdI.model == ModelType.KONDO:
                        half = StdI.nsite // 2
                        Sym[StdI.NSym][half + jCell * StdI.NsiteUC + jsite] = (
                            half + kCell * StdI.NsiteUC + jsite)
                        Anti[StdI.NSym][half + jCell * StdI.NsiteUC + jsite] = ap_dot
            StdI.NSym += 1

    with open("qptransidx.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NQPTrans {StdI.NSym:10d}\n")
        fp.write("=============================================\n")
        fp.write("======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n")
        fp.write("=============================================\n")
        for iSym in range(StdI.NSym):
            fp.write(f"{iSym} {1.0:10.5f}\n")
        for iSym in range(StdI.NSym):
            for jsite in range(StdI.nsite):
                a = _parity_sign(Anti[iSym][jsite])
                fp.write(f"{iSym:5d}  {jsite:5d}  {Sym[iSym][jsite]:5d}  {a:5d}\n")
    print("    qptransidx.def is written.")


def _init_site_sub(StdI: StdIntList) -> None:
    """Initialize sub-cell for mVMC/UHF/HWAVE.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    """
    StdI.Lsub, StdI.Wsub, StdI.Hsub = _validate_box_params(
        StdI.Lsub, StdI.Wsub, StdI.Hsub, StdI.boxsub,
        suffix="sub", defaults=StdI.box)

    # Calculate reciprocal lattice vectors
    StdI.NCellsub, StdI.rboxsub = _det_and_cofactor(StdI.boxsub)
    print(f"         Number of Cell in the sublattice: {abs(StdI.NCellsub)}")
    if StdI.NCellsub == 0:
        exit_program(-1)

    # Check commensurate
    if not _check_commensurate(StdI.rboxsub, StdI.box, StdI.NCellsub):
        print("\n ERROR ! Sublattice is INCOMMENSURATE !\n")
        exit_program(-1)


def _assign_orb_sector(
    StdI: StdIntList,
    iOrb: int,
    anti_val: int,
    iCell: int, jCell: int,
    iCell2: int, jCell2: int,
    isite: int, jsite: int,
    i_off: int, j_off: int,
    is_new: bool,
) -> int:
    """Assign orbital and anti-orbital indices for one sector of a cell pair.

    For each ``(i_off, j_off)`` sector (itinerant/local-spin combination
    in the Kondo model, or ``(0, 0)`` for the base sector), this function:

    1.  If *is_new*, assigns a fresh ``iOrb`` to the reference cell pair
        ``(iCell2, jCell2)`` and records the anti-periodic sign.
    2.  Always copies the reference orbital index to the current cell pair
        ``(iCell, jCell)`` and records the anti-periodic sign.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    iOrb : int
        Next available orbital index.
    anti_val : int
        Anti-periodic sign (+1 or -1).
    iCell, jCell : int
        Current cell indices.
    iCell2, jCell2 : int
        Reference (reduced) cell indices.
    isite, jsite : int
        Intra-unit-cell site indices.
    i_off, j_off : int
        Site-index offsets for this sector (0 for itinerant sites,
        ``nsite // 2`` for local-spin sites in the Kondo model).
    is_new : bool
        Whether this ``(iCell2, jCell2)`` pair is being visited for the
        first time.

    Returns
    -------
    int
        Updated ``iOrb`` (incremented by 1 if *is_new*).
    """
    nu = StdI.NsiteUC
    row2 = i_off + iCell2 * nu + isite
    col2 = j_off + jCell2 * nu + jsite

    if is_new:
        StdI.Orb[row2, col2] = iOrb
        StdI.AntiOrb[row2, col2] = anti_val
        iOrb += 1

    row = i_off + iCell * nu + isite
    col = j_off + jCell * nu + jsite
    StdI.Orb[row, col] = StdI.Orb[row2, col2]
    StdI.AntiOrb[row, col] = anti_val

    return iOrb


def generate_orb(StdI: StdIntList) -> None:
    """Generate orbital index for mVMC.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    """
    _init_site_sub(StdI)

    StdI.Orb = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
    StdI.AntiOrb = np.zeros((StdI.nsite, StdI.nsite), dtype=int)
    CellDone = np.zeros((StdI.NCell, StdI.NCell), dtype=int)

    iOrb = 0
    for iCell in range(StdI.NCell):
        iCV = _cell_vector(StdI.Cell, iCell)
        nBox, iCellV = _fold_site_sub(StdI, iCV)
        nBox, iCellV = _fold_site(StdI, iCellV)

        iCell2 = _find_cell_index(StdI, iCellV)

        for jCell in range(StdI.NCell):
            jCV = _cell_vector(StdI.Cell, jCell)
            jCellV = [jc + v - ic for jc, v, ic in zip(jCV, iCellV, iCV)]
            nBox, jCellV = _fold_site(StdI, jCellV)

            jCell2 = _find_cell_index(StdI, jCellV)

            # AntiPeriodic factor
            dCellV = [jc - ic for jc, ic in zip(jCV, iCV)]
            nBox_d, _ = _fold_site(StdI, dCellV)
            anti_val = _parity_sign(
                _anti_period_dot(StdI.AntiPeriod, nBox_d))

            # Build list of (i_offset, j_offset) sector pairs
            sectors = [(0, 0)]
            if StdI.model == ModelType.KONDO:
                half = StdI.nsite // 2
                sectors += [(half, 0), (0, half), (half, half)]

            for isite in range(StdI.NsiteUC):
                for jsite in range(StdI.NsiteUC):
                    for i_off, j_off in sectors:
                        iOrb = _assign_orb_sector(
                            StdI, iOrb, anti_val,
                            iCell, jCell, iCell2, jCell2,
                            isite, jsite, i_off, j_off,
                            is_new=(CellDone[iCell2, jCell2] == 0))

            CellDone[iCell2, jCell2] = 1

    StdI.NOrb = iOrb


def _jastrow_momentum_projected(
    StdI: StdIntList,
    Jastrow: np.ndarray,
) -> tuple[int, np.ndarray]:
    """Compute Jastrow indices using momentum-projected (symmetrised) orbital.

    Steps:

    1. Copy the orbital index matrix into *Jastrow*.
    2. Symmetrise: for each orbital index, set ``J[j,i] = J[i,j]``.
    3. Exclude local-spin sites (set their rows/columns to -1).
    4. Renumber: walk the strict lower triangle, assign negative
       temporaries, then invert so indices become non-negative.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (read-only).
    Jastrow : numpy.ndarray
        ``(nsite, nsite)`` integer array, modified **in place** during
        steps 1-3.  Replaced by ``-1 - Jastrow`` in step 4 (the
        returned array may be a new object).

    Returns
    -------
    NJastrow : int
        Number of unique Jastrow indices.
    Jastrow : numpy.ndarray
        Final renumbered Jastrow index matrix.
    """
    # (1) Copy Orbital index
    for isite in range(StdI.nsite):
        for jsite in range(StdI.nsite):
            Jastrow[isite, jsite] = StdI.Orb[isite, jsite]

    # (2) Symmetrize
    for iorb in range(StdI.NOrb):
        for isite in range(StdI.nsite):
            for jsite in range(StdI.nsite):
                if Jastrow[isite, jsite] == iorb:
                    Jastrow[jsite, isite] = Jastrow[isite, jsite]

    # (3) Exclude local-spin sites and renumber
    NJastrow = 0 if StdI.model == ModelType.HUBBARD else -1
    for isite in range(StdI.nsite):
        if StdI.locspinflag[isite] != 0:
            Jastrow[isite, :] = -1
            Jastrow[:, isite] = -1
            continue

        for jsite in range(isite):
            if Jastrow[isite, jsite] >= 0:
                iJastrow = Jastrow[isite, jsite]
                NJastrow -= 1
                mask = (Jastrow == iJastrow)
                Jastrow[mask] = NJastrow

    NJastrow = -NJastrow
    Jastrow = -1 - Jastrow
    return NJastrow, Jastrow


def _jastrow_global_optimization(
    StdI: StdIntList,
    Jastrow: np.ndarray,
) -> int:
    """Compute Jastrow indices using global (cell-based) optimisation.

    For the Spin model, all pairs share a single Jastrow index.
    For Hubbard/Kondo, unique indices are assigned per cell-displacement
    pair, respecting reversal symmetry and excluding the on-site term.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (read-only except for lattice lookups).
    Jastrow : numpy.ndarray
        ``(nsite, nsite)`` integer array, modified **in place**.

    Returns
    -------
    int
        Number of unique Jastrow indices (``NJastrow``).
    """
    if StdI.model == ModelType.SPIN:
        Jastrow[:, :] = 0
        return 1

    NJastrow = 0
    if StdI.model == ModelType.KONDO:
        half = StdI.nsite // 2
        Jastrow[:, :half] = 0
        Jastrow[:half, :] = 0
        NJastrow += 1

    for dCell in range(StdI.NCell):
        dCV = _cell_vector(StdI.Cell, dCell)
        isite, jsite, Cphase, dR_arr = find_site(
            StdI, 0, 0, 0, -dCV[0], -dCV[1], -dCV[2], 0, 0)
        if StdI.model == ModelType.KONDO:
            jsite -= StdI.NCell * StdI.NsiteUC
        iCell_j = jsite // StdI.NsiteUC
        if iCell_j < dCell:
            continue
        reversal = 1 if iCell_j == dCell else 0

        for isiteUC in range(StdI.NsiteUC):
            for jsiteUC in range(StdI.NsiteUC):
                if reversal == 1 and jsiteUC > isiteUC:
                    continue
                if isiteUC == jsiteUC and dCV == [0, 0, 0]:
                    continue

                for iCell_idx in range(StdI.NCell):
                    iCV = _cell_vector(StdI.Cell, iCell_idx)
                    i_s, j_s, _, _ = find_site(
                        StdI,
                        iCV[0], iCV[1], iCV[2],
                        dCV[0], dCV[1], dCV[2],
                        isiteUC, jsiteUC)
                    Jastrow[i_s, j_s] = NJastrow
                    Jastrow[j_s, i_s] = NJastrow

                NJastrow += 1

    return NJastrow


def print_jastrow(StdI: StdIntList) -> None:
    """Output Jastrow factor index file ``jastrowidx.def`` (mVMC only).

    Delegates computation to :func:`_jastrow_momentum_projected` or
    :func:`_jastrow_global_optimization` based on ``NMPTrans``, then
    writes the results to ``jastrowidx.def``.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.
    """
    Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)

    if abs(StdI.NMPTrans) == 1 or StdI.NMPTrans == NaN_i:
        NJastrow, Jastrow = _jastrow_momentum_projected(StdI, Jastrow)
    else:
        NJastrow = _jastrow_global_optimization(StdI, Jastrow)

    with open("jastrowidx.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NJastrowIdx {NJastrow:10d}\n")
        fp.write(f"ComplexType {0:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")

        for isite in range(StdI.nsite):
            for jsite in range(StdI.nsite):
                if isite == jsite:
                    continue
                fp.write(f"{isite:5d}  {jsite:5d}  {Jastrow[isite, jsite]:5d}\n")

        for iJastrow in range(NJastrow):
            if StdI.model == ModelType.HUBBARD or iJastrow > 0:
                fp.write(f"{iJastrow:5d}  {1:5d}\n")
            else:
                fp.write(f"{iJastrow:5d}  {0:5d}\n")
    print("    jastrowidx.def is written.")
