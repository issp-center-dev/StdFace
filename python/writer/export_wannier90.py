"""
Export functions for the Wannier90/HWAVE format.

This module contains functions to export lattice geometry and interaction
parameters in a format compatible with Wannier90. The main public functions
are:

- ``export_geometry()`` -- Exports lattice vectors and orbital positions
- ``export_interaction()`` -- Exports hopping and interaction parameters

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

from dataclasses import dataclass, field

import numpy as np

from stdface_vals import StdIntList, NaN_i, UNSET_STRING
from param_check import exit_program
from lattice.site_util import _cell_vector

# -----------------------------------------------------------------------
#  Module-level constants
# -----------------------------------------------------------------------

_EPS = 1.0e-8
"""Small number for floating-point comparisons."""


# -----------------------------------------------------------------------
#  Control flags
# -----------------------------------------------------------------------

_is_export_all = 1
"""Default for the exportall parameter (1 = export zero elements too)."""


# -----------------------------------------------------------------------
#  Internal data structures
# -----------------------------------------------------------------------

@dataclass
class _IntrItem:
    """Store an interaction parameter between orbitals.

    Attributes
    ----------
    r : list of int
        Relative coordinate between orbitals (length 3).
    a : int
        First orbital index.
    b : int
        Second orbital index.
    s : int
        First spin index.
    t : int
        Second spin index.
    v : complex
        Interaction strength.
    """
    r: list[int] = field(default_factory=lambda: [0, 0, 0])
    a: int = 0
    b: int = 0
    s: int = 0
    t: int = 0
    v: complex = 0.0 + 0.0j


# -----------------------------------------------------------------------
#  Helper: fatal error
# -----------------------------------------------------------------------

def _fatal(msg: str) -> None:
    """Print an error message and exit.

    Parameters
    ----------
    msg : str
        Error description.
    """
    import sys
    print(f"ERROR: {msg}", file=sys.stderr)
    exit_program(-1)


# -----------------------------------------------------------------------
#  Write geometry file
# -----------------------------------------------------------------------

def _write_geometry(StdI: StdIntList, fname: str) -> None:
    """Write geometry data to file in Wannier90 format.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters containing geometry info.
    fname : str
        Output filename.

    Notes
    -----
    Writes:
    - Primitive lattice vectors
    - Number of orbitals per unit cell
    - Orbital positions in fractional coordinates
    """
    try:
        fp_out = open(fname, "w")
    except OSError:
        _fatal(f"cannot open file for output: {fname}")
        return  # unreachable, but keeps type checker happy

    with fp_out:
        # Print primitive vectors
        for row in StdI.direct:
            fp_out.write(f"{row[0]:16.12f} {row[1]:16.12f} {row[2]:16.12f}\n")

        # Print number of orbits
        fp_out.write(f"{StdI.NsiteUC}\n")

        # Print centre of orbits
        for tau_row in StdI.tau[:StdI.NsiteUC]:
            fp_out.write(f"{tau_row[0]:25.15e} "
                         f"{tau_row[1]:25.15e} "
                         f"{tau_row[2]:25.15e}\n")
    print(f"{fname:>24s} is written.")


# -----------------------------------------------------------------------
#  Index macro (replaces C's _index)
# -----------------------------------------------------------------------

def _compute_index(rx: int, ry: int, rz: int,
                   a: int, b: int, s: int, t: int,
                   rr: list[int], nsiteuc: int, nspin: int) -> int:
    """Compute flat index into the interaction matrix.

    Parameters
    ----------
    rx, ry, rz : int
        Relative coordinates.
    a, b : int
        Orbital indices.
    s, t : int
        Spin indices.
    rr : list of int
        Half-ranges for each direction.
    nsiteuc : int
        Number of sites per unit cell.
    nspin : int
        Number of spin states.

    Returns
    -------
    int
        Flat index into the matrix array.
    """
    return (t + nspin * (
        s + nspin * (
            b + nsiteuc * (
                a + nsiteuc * (
                    (rz + rr[2]) + (rr[2] * 2 + 1) * (
                        (ry + rr[1]) + (rr[1] * 2 + 1) * (
                            rx + rr[0]
                        )))))))


# -----------------------------------------------------------------------
#  Write Wannier90-format interaction file
# -----------------------------------------------------------------------

def _build_wannier_matrix(
    nintr_table: int,
    intr_table: list[_IntrItem],
    nsiteuc: int,
    nspin: int,
) -> tuple[list[int], int, np.ndarray]:
    """Build the flat interaction matrix from a list of interaction items.

    Scans the items to determine the coordinate half-ranges, allocates
    a flat complex matrix, and populates it — including Hermitian-conjugate
    entries for any reverse-direction slot that is still empty.

    Parameters
    ----------
    nintr_table : int
        Number of interaction terms.
    intr_table : list of _IntrItem
        Array of interaction parameters.
    nsiteuc : int
        Number of sites per unit cell.
    nspin : int
        Number of spin states.

    Returns
    -------
    rr : list of int
        Half-ranges ``[rr0, rr1, rr2]`` for each lattice direction.
    nvol : int
        Total number of real-space unit cells in the range.
    matrix : numpy.ndarray
        Flat complex array of length ``nvol * nsiteuc**2 * nspin**2``.
    """
    rmin = list(intr_table[0].r)
    rmax = list(intr_table[0].r)

    for item in intr_table[1:]:
        for i, r in enumerate(item.r):
            if r < rmin[i]:
                rmin[i] = r
            if r > rmax[i]:
                rmax[i] = r

    rr = [max(abs(lo), abs(hi)) for lo, hi in zip(rmin, rmax)]

    nvol = (rr[0] * 2 + 1) * (rr[1] * 2 + 1) * (rr[2] * 2 + 1)
    matrix_size = nvol * nsiteuc * nsiteuc * nspin * nspin
    matrix = np.zeros(matrix_size, dtype=complex)

    for k in range(nintr_table):
        idx = _compute_index(
            intr_table[k].r[0], intr_table[k].r[1], intr_table[k].r[2],
            intr_table[k].a, intr_table[k].b,
            intr_table[k].s, intr_table[k].t,
            rr, nsiteuc, nspin)
        matrix[idx] = intr_table[k].v

        ridx = _compute_index(
            -intr_table[k].r[0], -intr_table[k].r[1], -intr_table[k].r[2],
            intr_table[k].b, intr_table[k].a,
            intr_table[k].t, intr_table[k].s,
            rr, nsiteuc, nspin)

        if abs(matrix[ridx]) < _EPS:
            matrix[ridx] = np.conj(intr_table[k].v)

    return rr, nvol, matrix


def _write_wannier_body(
    fp: 'TextIO',
    rr: list[int],
    nvol: int,
    nsiteuc: int,
    nspin: int,
    matrix: np.ndarray,
) -> None:
    """Write the matrix body of a Wannier90-format interaction file.

    Iterates over all real-space cells and orbital/spin pairs, writing
    one line per entry.  When ``nspin > 1`` the extended format with
    explicit spin columns ``s`` and ``t`` is used; otherwise the compact
    format without spin columns is written.

    Parameters
    ----------
    fp : file object
        Open file to write matrix entries to.
    rr : list of int
        Half-ranges for each lattice direction.
    nvol : int
        Total number of real-space unit cells.
    nsiteuc : int
        Number of sites per unit cell.
    nspin : int
        Number of spin states.
    matrix : numpy.ndarray
        Flat complex interaction matrix.
    """
    for r in range(nvol):
        rz = r % (rr[2] * 2 + 1) - rr[2]
        ry = (r // (rr[2] * 2 + 1)) % (rr[1] * 2 + 1) - rr[1]
        rx = (r // ((rr[2] * 2 + 1) * (rr[1] * 2 + 1))) % (rr[0] * 2 + 1) - rr[0]

        for a in range(nsiteuc):
            for b in range(nsiteuc):

                if nspin > 1:
                    # Extended format
                    for s in range(nspin):
                        for t in range(nspin):
                            idx = _compute_index(rx, ry, rz, a, b, s, t,
                                                 rr, nsiteuc, nspin)
                            if _is_export_all or abs(matrix[idx]) > _EPS:
                                fp.write(
                                    f"{rx:4d} {ry:4d} {rz:4d} "
                                    f"{a + 1:4d} {b + 1:4d} "
                                    f"{s:4d} {t:4d} "
                                    f"{matrix[idx].real:16.12f} "
                                    f"{matrix[idx].imag:16.12f}\n")
                else:
                    s = 0
                    t = 0
                    idx = _compute_index(rx, ry, rz, a, b, s, t,
                                         rr, nsiteuc, nspin)
                    if _is_export_all or abs(matrix[idx]) > _EPS:
                        fp.write(
                            f"{rx:4d} {ry:4d} {rz:4d} "
                            f"{a + 1:4d} {b + 1:4d} "
                            f"{matrix[idx].real:16.12f} "
                            f"{matrix[idx].imag:16.12f}\n")


def _write_wannier90(nintr_table: int, intr_table: list[_IntrItem],
                     nsiteuc: int, nspin: int,
                     fname: str, tagname: str) -> None:
    """Write interaction parameters to file in Wannier90 format.

    Builds the interaction matrix from the item list, then writes the
    Wannier90-format file with header and body.

    Parameters
    ----------
    nintr_table : int
        Number of interaction terms.
    intr_table : list of _IntrItem
        Array of interaction parameters.
    nsiteuc : int
        Number of sites per unit cell.
    nspin : int
        Number of spin states.
    fname : str
        Output filename.
    tagname : str
        Tag identifying interaction type.
    """
    rr, nvol, matrix = _build_wannier_matrix(
        nintr_table, intr_table, nsiteuc, nspin)

    try:
        fp_out = open(fname, "w")
    except OSError:
        _fatal(f"cannot open file: {fname}")
        return

    with fp_out:
        # Write header
        fp_out.write(f"{tagname} in wannier90-like format for uhfk\n")
        fp_out.write(f"{nsiteuc}\n{nvol}\n")
        for i in range(nvol):
            if i > 0 and i % 15 == 0:
                fp_out.write("\n")
            fp_out.write(f" {1}")
        fp_out.write("\n")

        _write_wannier_body(fp_out, rr, nvol, nsiteuc, nspin, matrix)
    print(f"{fname:>24s} is written.")


# -----------------------------------------------------------------------
#  Unfold site coordinates
# -----------------------------------------------------------------------

def _unfold_site(StdI: StdIntList, v_in: list[int]) -> list[int]:
    """Convert site coordinates from [0, N] to [-N/2, N/2] range.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters.
    v_in : list of int
        Input coordinates (length 3).

    Returns
    -------
    list of int
        Output coordinates in the unfolded range (length 3).
    """
    # v = rbox @ v_in (matrix-vector product)
    v = StdI.rbox.astype(int) @ np.array(v_in)

    # Fold to [-N/2, N/2] range using vectorized operations
    vv = v / StdI.NCell
    v = np.where(vv > 0.5, v - StdI.NCell, v)
    v = np.where(vv <= -0.5, v + StdI.NCell, v)

    # w = (v @ box) // NCell
    w = (v @ StdI.box.astype(int)) // StdI.NCell

    return w.tolist()


# -----------------------------------------------------------------------
#  Key generation / comparison helpers
# -----------------------------------------------------------------------

def _generate_key(keylen: int, index: list[int], ordered: int) -> list[int]:
    """Generate key for interaction table entry.

    Parameters
    ----------
    keylen : int
        Length of key (1, 2, or 4).
    index : list of int
        Input indices.
    ordered : int
        Whether to order indices (1 = yes).

    Returns
    -------
    list of int
        Generated key.
    """
    if keylen == 1:
        return [index[0]]
    elif keylen == 2:
        i, j = index[0], index[1]
        if ordered == 1 and i > j:
            return [j, i]
        else:
            return [i, j]
    elif keylen == 4:
        i, s, j, t = index[0], index[1], index[2], index[3]
        if ordered == 1 and i > j:
            return [j, t, i, s]
        else:
            return [i, s, j, t]
    else:
        _fatal(f"unsupported keylen: {keylen}")
        return []  # unreachable


# -----------------------------------------------------------------------
#  Accumulate entries
# -----------------------------------------------------------------------

def _accumulate_list(keylen: int,
                     ntbl: int,
                     tbl_index: np.ndarray,
                     tbl_value: np.ndarray,
                     ordered: int
                     ) -> tuple:
    """Accumulate interaction table entries with same key.

    Parameters
    ----------
    keylen : int
        Length of key (1, 2, or 4).
    ntbl : int
        Number of input entries.
    tbl_index : np.ndarray
        Input indices, shape ``(ntbl, keylen)``.
    tbl_value : np.ndarray
        Input values, shape ``(ntbl,)``.
    ordered : int
        Whether to order indices (1 = yes).

    Returns
    -------
    nintr : int
        Number of output entries.
    intr_index : list of list of int
        Output indices (each entry is a list of length *keylen*).
    intr_value : list of complex
        Output values.
    """
    # Use a dict keyed by tuple for O(1) lookup instead of O(n) linear scan
    accum: dict[tuple[int, ...], complex] = {}
    key_order: list[tuple[int, ...]] = []

    for k in range(ntbl):
        idx = tuple(_generate_key(keylen, list(tbl_index[k, :keylen]), ordered))
        val = complex(tbl_value[k])

        if idx in accum:
            accum[idx] += val
        else:
            accum[idx] = val
            key_order.append(idx)

    # Filter out near-zero entries
    intr_index: list[list[int]] = []
    intr_value: list[complex] = []
    for key in key_order:
        if abs(accum[key]) >= _EPS:
            intr_index.append(list(key))
            intr_value.append(accum[key])

    return len(intr_index), intr_index, intr_value


# -----------------------------------------------------------------------
#  Export: inter-site interaction (complex)
# -----------------------------------------------------------------------

def _build_inter_table(
    StdI: StdIntList,
    nintr: int,
    intr_index: list[list[int]],
    intr_value: np.ndarray,
) -> list[_IntrItem]:
    """Build a deduplicated interaction table in relative coordinates.

    Converts accumulated inter-site interaction entries from absolute
    site indices to relative-coordinate ``_IntrItem`` entries,
    deduplicating by the key ``(rr, isite, jsite)``.  Spin indices
    are fixed to ``(0, 0)``.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters (reads ``NsiteUC``, ``Cell``).
    nintr : int
        Number of accumulated entries.
    intr_index : list of list of int
        Accumulated indices, each row is ``[idx_i, idx_j]``.
    intr_value : np.ndarray
        Accumulated complex values, shape ``(nintr,)``.

    Returns
    -------
    list of _IntrItem
        Deduplicated interaction entries in relative coordinates.
    """
    intr_table: list[_IntrItem] = []
    seen: dict[tuple, int] = {}  # key -> index in intr_table

    for k in range(nintr):
        idx_i = intr_index[k][0]
        idx_j = intr_index[k][1]

        icell = idx_i // StdI.NsiteUC
        isite = idx_i % StdI.NsiteUC
        jcell = idx_j // StdI.NsiteUC
        jsite = idx_j % StdI.NsiteUC

        jCV = _cell_vector(StdI.Cell, jcell)
        iCV = _cell_vector(StdI.Cell, icell)
        rr = [j - i for j, i in zip(jCV, iCV)]
        rr = _unfold_site(StdI, rr)

        lookup_key = (tuple(rr), isite, jsite)
        if lookup_key in seen:
            existing = intr_table[seen[lookup_key]]
            if abs(existing.v - intr_value[k]) > _EPS:
                print(f"WARNING: not uniform. "
                      f"expected=({existing.v.real},{existing.v.imag}), "
                      f"found=({intr_value[k].real},{intr_value[k].imag}) "
                      f"for index {idx_i},{idx_j}")
        else:
            seen[lookup_key] = len(intr_table)
            intr_table.append(_IntrItem(
                r=list(rr), a=isite, b=jsite,
                s=0, t=0, v=intr_value[k]))

    return intr_table


def _export_inter(StdI: StdIntList,
                  ntbl: int,
                  tbl_index: np.ndarray,
                  tbl_value: np.ndarray,
                  fname: str, tagname: str) -> None:
    """Export interaction coefficients from complex array.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters.
    ntbl : int
        Number of interaction terms.
    tbl_index : np.ndarray
        Interaction indices, shape ``(ntbl, 2)``.
    tbl_value : np.ndarray
        Interaction values, shape ``(ntbl,)`` (complex).
    fname : str
        Output filename.
    tagname : str
        Tag identifying interaction type.
    """
    if ntbl == 0:
        print(f"{fname:>24s} is skipped.")
        return

    nintr, intr_index, intr_value = _accumulate_list(
        2, ntbl, tbl_index, tbl_value, 1)

    if nintr > 0:
        intr_table = _build_inter_table(StdI, nintr, intr_index, intr_value)

        if len(intr_table) > 0:
            _write_wannier90(len(intr_table), intr_table, StdI.NsiteUC, 1,
                             fname, tagname)
        else:
            print(f"{fname:>24s} is skipped.")
    else:
        print(f"{fname:>24s} is skipped.")


# -----------------------------------------------------------------------
#  Export: inter-site interaction (real -> complex wrapper)
# -----------------------------------------------------------------------

def _export_inter_real(StdI: StdIntList,
                       ntbl: int,
                       tbl_index: np.ndarray,
                       tbl_value: np.ndarray,
                       fname: str, tagname: str) -> None:
    """Export interaction coefficients from real array.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters.
    ntbl : int
        Number of interaction terms.
    tbl_index : np.ndarray
        Interaction indices, shape ``(ntbl, 2)``.
    tbl_value : np.ndarray
        Interaction values, shape ``(ntbl,)`` (real).
    fname : str
        Output filename.
    tagname : str
        Tag identifying interaction type.

    Notes
    -----
    Converts real array to complex and calls ``_export_inter()``.
    """
    if tbl_value is not None:
        buf = tbl_value.astype(complex)
    else:
        buf = np.array([], dtype=complex)
    _export_inter(StdI, ntbl, tbl_index, buf, fname, tagname)


# -----------------------------------------------------------------------
#  Export: transfer (hopping) coefficients
# -----------------------------------------------------------------------

def _build_transfer_table(
    StdI: StdIntList,
    nintr: int,
    intr_index: list[list[int]],
    intr_value: np.ndarray,
    spin_dep: int,
) -> list[_IntrItem]:
    """Build a deduplicated transfer table in relative coordinates.

    Converts accumulated transfer entries from absolute site indices
    to relative-coordinate ``_IntrItem`` entries, deduplicating by
    the key ``(rr, isite, jsite, ispin, jspin)``.  Values are
    sign-flipped (multiplied by −1) per the Wannier90 convention.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters (reads ``NsiteUC``, ``Cell``).
    nintr : int
        Number of accumulated entries.
    intr_index : list of list of int
        Accumulated indices, each row is ``[idx_i, ispin, idx_j, jspin]``.
    intr_value : np.ndarray
        Accumulated complex values, shape ``(nintr,)``.
        **Modified in-place** (sign-flipped).
    spin_dep : int
        Whether transfer is spin-dependent (1 = yes, 0 = no).
        When 0, entries where ``(ispin, jspin) != (0, 0)`` are skipped.

    Returns
    -------
    list of _IntrItem
        Deduplicated transfer entries in relative coordinates.
    """
    intr_table: list[_IntrItem] = []
    seen: dict[tuple, int] = {}  # key -> index in intr_table

    for k in range(nintr):
        idx_i = intr_index[k][0]
        ispin = intr_index[k][1]
        idx_j = intr_index[k][2]
        jspin = intr_index[k][3]

        icell = idx_i // StdI.NsiteUC
        isite = idx_i % StdI.NsiteUC
        jcell = idx_j // StdI.NsiteUC
        jsite = idx_j % StdI.NsiteUC

        jCV = _cell_vector(StdI.Cell, jcell)
        iCV = _cell_vector(StdI.Cell, icell)
        rr = [j - i for j, i in zip(jCV, iCV)]
        rr = _unfold_site(StdI, rr)

        intr_value[k] *= -1  # by convention

        lookup_key = (tuple(rr), isite, jsite, ispin, jspin)
        if lookup_key in seen:
            existing = intr_table[seen[lookup_key]]
            if abs(existing.v - intr_value[k]) > _EPS:
                print(f"WARNING: not uniform. "
                      f"expected=({existing.v.real},{existing.v.imag}), "
                      f"found=({intr_value[k].real},{intr_value[k].imag}) "
                      f"for index {idx_i},{idx_j}")
        else:
            if spin_dep == 0 and not (ispin == 0 and jspin == 0):
                continue  # skip

            seen[lookup_key] = len(intr_table)
            intr_table.append(_IntrItem(
                r=list(rr), a=isite, b=jsite,
                s=ispin, t=jspin, v=intr_value[k]))

    return intr_table


def _export_transfer(StdI: StdIntList,
                     ntbl: int,
                     tbl_index: np.ndarray,
                     tbl_value: np.ndarray,
                     fname: str, tagname: str,
                     spin_dep: int) -> None:
    """Export transfer (hopping) coefficients.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters.
    ntbl : int
        Number of transfer terms.
    tbl_index : np.ndarray
        Transfer indices, shape ``(ntbl, 4)``.
    tbl_value : np.ndarray
        Transfer values, shape ``(ntbl,)`` (complex).
    fname : str
        Output filename.
    tagname : str
        Tag identifying transfer type.
    spin_dep : int
        Whether transfer is spin-dependent (1 = yes, 0 = no).
    """
    if ntbl == 0:
        print(f"{fname:>24s} is skipped.")
        return

    # Accumulate entries of the same index pair
    nintr, intr_index, intr_value = _accumulate_list(
        4, ntbl, tbl_index, tbl_value, 0)

    if nintr > 0:
        intr_table = _build_transfer_table(
            StdI, nintr, intr_index, intr_value, spin_dep)

        if len(intr_table) > 0:
            _write_wannier90(len(intr_table), intr_table, StdI.NsiteUC,
                             2 if spin_dep == 1 else 1,
                             fname, tagname)
        else:
            print(f"{fname:>24s} is skipped.")
    else:
        print(f"{fname:>24s} is skipped.")


# -----------------------------------------------------------------------
#  Export: on-site Coulomb
# -----------------------------------------------------------------------

def _build_coulomb_intra_table(
    StdI: StdIntList,
    nintr: int,
    intr_index: list[list[int]],
    intr_value: np.ndarray,
) -> list[_IntrItem]:
    """Build a deduplicated on-site Coulomb table.

    Converts accumulated on-site Coulomb entries from absolute site
    indices to ``_IntrItem`` entries, deduplicating by the unit-cell
    site index ``isite``.  Each entry has ``rr = [0,0,0]`` and
    ``a == b == isite``.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters (reads ``NsiteUC``).
    nintr : int
        Number of accumulated entries.
    intr_index : list of list of int
        Accumulated indices, each row is ``[idx_i]``.
    intr_value : np.ndarray
        Accumulated complex values, shape ``(nintr,)``.

    Returns
    -------
    list of _IntrItem
        Deduplicated on-site Coulomb entries.
    """
    intr_table: list[_IntrItem] = []

    for k in range(nintr):
        idx_i = intr_index[k][0]
        isite = idx_i % StdI.NsiteUC

        is_found = False
        for item in intr_table:
            if item.a == isite:
                is_found = True
                if abs(item.v - intr_value[k]) > _EPS:
                    print(f"WARNING: not uniform. "
                          f"expected=({item.v.real},{item.v.imag}), "
                          f"found=({intr_value[k].real},{intr_value[k].imag}) "
                          f"for index {idx_i}")
                break

        if not is_found:
            item = _IntrItem(
                r=[0, 0, 0], a=isite, b=isite,
                s=0, t=0, v=intr_value[k])
            intr_table.append(item)

    return intr_table


def _export_coulomb_intra(StdI: StdIntList,
                          ntbl: int,
                          tbl_index: np.ndarray,
                          tbl_value: np.ndarray,
                          fname: str, tagname: str) -> None:
    """Export on-site Coulomb term coefficients.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters.
    ntbl : int
        Number of Coulomb terms.
    tbl_index : np.ndarray
        Site indices, shape ``(ntbl, 1)``.
    tbl_value : np.ndarray
        Coulomb coefficients, shape ``(ntbl,)`` (real).
    fname : str
        Output filename.
    tagname : str
        Tag identifying the term type.
    """
    if ntbl == 0:
        print(f"{fname:>24s} is skipped.")
        return

    tbl_value_c = tbl_value.astype(complex)

    nintr, intr_index, intr_value = _accumulate_list(
        1, ntbl, tbl_index, tbl_value_c, 1)

    if nintr > 0:
        intr_table = _build_coulomb_intra_table(
            StdI, nintr, intr_index, intr_value)

        if len(intr_table) > 0:
            _write_wannier90(len(intr_table), intr_table, StdI.NsiteUC, 1,
                             fname, tagname)
        else:
            print(f"{fname:>24s} is skipped.")
    else:
        print(f"{fname:>24s} is skipped.")


# -----------------------------------------------------------------------
#  Filename prefix helper
# -----------------------------------------------------------------------

def _prefix(StdI: StdIntList, fname: str) -> str:
    """Add prefix to filename if ``StdI.fileprefix`` is set.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters.
    fname : str
        Base filename.

    Returns
    -------
    str
        Filename with optional prefix prepended.

    Notes
    -----
    In the C code, the sentinel for "no prefix" is ``"****"``.
    In Python, the sentinel is the empty string ``""``.
    """
    if StdI.fileprefix == "" or StdI.fileprefix == UNSET_STRING:
        return fname
    else:
        return f"{StdI.fileprefix}_{fname}"


# =======================================================================
#  Public API
# =======================================================================

def export_geometry(StdI: StdIntList) -> None:
    """Export geometry information to file.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters containing lattice geometry.

    Notes
    -----
    This function is only meaningful when ``StdI.solver == "HWAVE"``.
    It writes the file ``geom.dat`` (optionally with prefix).
    """
    _write_geometry(StdI, _prefix(StdI, "geom.dat"))


def export_interaction(StdI: StdIntList) -> None:
    """Export interaction term coefficients to files.

    Parameters
    ----------
    StdI : StdIntList
        Standard input parameters containing interaction terms.

    Notes
    -----
    This function is only meaningful when ``StdI.solver == "HWAVE"``.
    It writes the following files (optionally with prefix):

    - ``transfer.dat``
    - ``coulombintra.dat``
    - ``coulombinter.dat``
    - ``hund.dat``
    - ``exchange.dat``
    - ``pairlift.dat``
    - ``pairhopp.dat``
    """
    global _is_export_all

    if StdI.export_all != NaN_i:
        _is_export_all = StdI.export_all

    _export_transfer(
        StdI,
        StdI.ntrans, StdI.transindx, StdI.trans,
        _prefix(StdI, "transfer.dat"), "Transfer",
        0)

    _export_coulomb_intra(
        StdI,
        StdI.NCintra, StdI.CintraIndx, StdI.Cintra,
        _prefix(StdI, "coulombintra.dat"), "CoulombIntra")

    _export_inter_real(
        StdI,
        StdI.NCinter, StdI.CinterIndx, StdI.Cinter,
        _prefix(StdI, "coulombinter.dat"), "CoulombInter")

    _export_inter_real(
        StdI,
        StdI.NHund, StdI.HundIndx, StdI.Hund,
        _prefix(StdI, "hund.dat"), "Hund")

    _export_inter_real(
        StdI,
        StdI.NEx, StdI.ExIndx, StdI.Ex,
        _prefix(StdI, "exchange.dat"), "Exchange")

    _export_inter_real(
        StdI,
        StdI.NPairLift, StdI.PLIndx, StdI.PairLift,
        _prefix(StdI, "pairlift.dat"), "PairLift")

    _export_inter_real(
        StdI,
        StdI.NPairHopp, StdI.PHIndx, StdI.PairHopp,
        _prefix(StdI, "pairhopp.dat"), "PairHopp")
