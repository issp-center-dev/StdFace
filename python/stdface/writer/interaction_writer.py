"""Interaction-term merge and file-output functions.

This module contains :func:`build_interactions`, which processes all
interaction types (CoulombIntra, CoulombInter, Hund, Exchange, PairLift,
PairHopp, InterAll) -- merging duplicate terms and counting non-zero
terms -- and returns the corresponding ``XxxData`` objects.

Functions
---------
build_interactions
    Merge duplicate interaction terms and return the interaction data.

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
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import numpy as np

from ..core.stdface_vals import StdIntList, AMPLITUDE_EPS

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
#  Helpers for the repeated merge → count → write pattern
# ---------------------------------------------------------------------------

def _merge_1idx(nterms: int, indx: list, coeff: list) -> None:
    """Merge duplicate terms for a 1-index interaction (CoulombIntra).

    Two terms are duplicates if they share the same single site index.
    The coefficient of the first occurrence is updated and the duplicate
    is zeroed.  Uses dict-based O(n) lookup instead of O(n²) pairwise scan.

    Parameters
    ----------
    nterms : int
        Number of interaction terms.
    indx : list
        Index array; each element is a list whose first element is the
        site index.
    coeff : list
        Coefficient array, modified in place.
    """
    seen: dict[int, int] = {}  # site index -> first occurrence position
    for k in range(nterms):
        key = indx[k][0]
        if key in seen:
            coeff[seen[key]] += coeff[k]
            coeff[k] = 0.0
        else:
            seen[key] = k


def _merge_2idx(nterms: int, indx: list, coeff: list) -> None:
    """Merge duplicate terms for a 2-index interaction.

    Two terms are duplicates if they share the same pair of site indices
    (in either order, i.e. symmetric).  Uses dict-based O(n) lookup
    instead of O(n²) pairwise scan.

    Parameters
    ----------
    nterms : int
        Number of interaction terms.
    indx : list
        Index array; each element is a list of two site indices.
    coeff : list
        Coefficient array, modified in place.
    """
    seen: dict[tuple[int, int], int] = {}  # canonical pair -> first occurrence
    for k in range(nterms):
        i0, i1 = indx[k][0], indx[k][1]
        key = (min(i0, i1), max(i0, i1))
        if key in seen:
            coeff[seen[key]] += coeff[k]
            coeff[k] = 0.0
        else:
            seen[key] = k


def _count_nonzero(nterms: int, coeff: list) -> int:
    """Count the number of terms with ``|coeff| > AMPLITUDE_EPS``.

    Parameters
    ----------
    nterms : int
        Total number of terms.
    coeff : list
        Coefficient array.

    Returns
    -------
    int
        Number of non-zero terms.
    """
    return sum(1 for k in range(nterms) if abs(coeff[k]) > AMPLITUDE_EPS)


@dataclass
class InteractionData:
    """One real-valued interaction ``.def`` file (1- or 2-index family).

    Attributes
    ----------
    filename : str
        Output file name (e.g. ``"coulombintra.def"``).
    count_label : str
        Count header label (e.g. ``"NCoulombIntra"``).
    banner : str
        Description banner line.
    rows : list of tuple
        ``(*site_indices, coeff)`` per non-zero term (``coeff`` real).
    """

    filename: str
    count_label: str
    banner: str
    rows: list

    def write(self, directory: Path = Path(".")) -> None:
        lines = ["=============================================\n",
                 f"{self.count_label} {len(self.rows):10d}\n",
                 "=============================================\n",
                 f"{self.banner}\n",
                 "=============================================\n"]
        for row in self.rows:
            *idx, coeff = row
            idx_str = " ".join(f"{i:5d}" for i in idx)
            lines.append(f"{idx_str} {coeff:25.15f}\n")
        with open(Path(directory) / self.filename, "w") as fp:
            fp.write("".join(lines))
        logger.info("    %s is written.", self.filename)

    def to_dict(self) -> dict:
        return {"filename": self.filename, "count_label": self.count_label,
                "banner": self.banner, "rows": [list(r) for r in self.rows]}

    @classmethod
    def from_dict(cls, data: dict) -> "InteractionData":
        return cls(data["filename"], data["count_label"], data["banner"],
                   [tuple(r) for r in data["rows"]])


def _build_interaction(
    StdI: StdIntList, meta: "_InteractionMeta"
) -> "InteractionData | None":
    """Merge one interaction family, set its ``L*`` flag, and return data.

    Returns ``None`` (and sets the flag to 0) when the family is empty or
    suppressed by boost mode.  Otherwise sets the flag to 1 and returns an
    :class:`InteractionData` with the merged non-zero rows.
    """
    terms = getattr(StdI, meta.list_attr)
    nterms = len(terms)
    coeff = [t[0] for t in terms]
    indx = [list(t[1:]) for t in terms]

    if meta.n_indices == 1:
        _merge_1idx(nterms, indx, coeff)
    else:
        _merge_2idx(nterms, indx, coeff)

    nintr0 = _count_nonzero(nterms, coeff)
    if nintr0 == 0 or StdI.lBoost == 1:
        setattr(StdI, meta.flag_attr, 0)
        return None
    setattr(StdI, meta.flag_attr, 1)

    rows = [tuple(int(x) for x in indx[k]) + (float(coeff[k]),)
            for k in range(nterms) if abs(coeff[k]) > AMPLITUDE_EPS]
    return InteractionData(meta.filename, meta.count_label, meta.banner, rows)


# ---------------------------------------------------------------------------
#  Interaction-type metadata table
# ---------------------------------------------------------------------------


class _InteractionMeta(NamedTuple):
    """Metadata for one interaction type (attribute names, file info).

    Attributes
    ----------
    list_attr : str
        Name of the term-list attribute (e.g. ``"Cintra_list"``).
    flag_attr : str
        Name of the output-flag attribute.
    filename : str
        Output ``.def`` file name.
    count_label : str
        Label for the count header line.
    banner : str
        Description banner in the header.
    n_indices : int
        Number of site indices per term (1 or 2).
    """

    list_attr: str
    flag_attr: str
    filename: str
    count_label: str
    banner: str
    n_indices: int


_INTERACTION_TYPES: list[_InteractionMeta] = [
    _InteractionMeta("Cintra_list", "LCintra",
                     "coulombintra.def", "NCoulombIntra",
                     "================== CoulombIntra ================", 1),
    _InteractionMeta("Cinter_list", "LCinter",
                     "coulombinter.def", "NCoulombInter",
                     "================== CoulombInter ================", 2),
    _InteractionMeta("Hund_list", "LHund",
                     "hund.def", "NHund",
                     "=============== Hund coupling ===============", 2),
    _InteractionMeta("Ex_list", "LEx",
                     "exchange.def", "NExchange",
                     "====== ExchangeCoupling coupling ============", 2),
    _InteractionMeta("PairLift_list", "LPairLift",
                     "pairlift.def", "NPairLift",
                     "====== Pair-Lift term ============", 2),
    _InteractionMeta("PairHopp_list", "LPairHopp",
                     "pairhopp.def", "NPairHopp",
                     "====== Pair-Hopping term ============", 2),
]
"""Metadata for the 6 standard interaction types processed by
:func:`build_interactions`.  Each entry maps attribute names, file
names, and header strings so that :func:`_process_interaction` can
handle all types generically.
"""


# ---------------------------------------------------------------------------
#  InterAll helpers
# ---------------------------------------------------------------------------

def _merge_interall_equivalent(
    nintr: int, intrindx: list, intr: list
) -> None:
    """Merge equivalent InterAll terms (Pass 1).

    For each pair of terms, check three patterns:

    - **Case A**: All 8 indices match exactly → **add** coefficients.
    - **Case B**: Hermitian conjugate match (indices 0–3 ↔ 4–7) with
      two exclusion conditions → **add** coefficients.
    - **Case C**: Partial swap patterns → **subtract** coefficients.

    In all cases the duplicate term is zeroed.

    Parameters
    ----------
    nintr : int
        Number of InterAll terms.
    intrindx : list
        Index array; each element is a list of 8 integers.
    intr : list
        Coefficient array (complex), modified in place.
    """
    for jintr in range(nintr):
        j = intrindx[jintr]
        for kintr in range(jintr + 1, nintr):
            k = intrindx[kintr]

            # Case A: all 8 indices match exactly
            # Case B: Hermitian conjugate match (indices 0-3 <-> 4-7)
            if (
                (j[0] == k[0] and j[1] == k[1]
                 and j[2] == k[2] and j[3] == k[3]
                 and j[4] == k[4] and j[5] == k[5]
                 and j[6] == k[6] and j[7] == k[7])
                or
                (j[0] == k[4] and j[1] == k[5]
                 and j[2] == k[6] and j[3] == k[7]
                 and j[4] == k[0] and j[5] == k[1]
                 and j[6] == k[2] and j[7] == k[3]
                 and not (j[0] == j[6] and j[1] == j[7])
                 and not (j[2] == j[4] and j[3] == j[5]))
            ):
                intr[jintr] = intr[jintr] + intr[kintr]
                intr[kintr] = 0.0

            # Case C: Partial swap patterns -> SUBTRACT
            elif (
                (j[0] == k[4] and j[1] == k[5]
                 and j[2] == k[2] and j[3] == k[3]
                 and j[4] == k[0] and j[5] == k[1]
                 and j[6] == k[6] and j[7] == k[7]
                 and not (j[2] == j[0] and j[3] == j[1])
                 and not (j[2] == j[4] and j[3] == j[5]))
                or
                (j[0] == k[0] and j[1] == k[1]
                 and j[2] == k[6] and j[3] == k[7]
                 and j[4] == k[4] and j[5] == k[5]
                 and j[6] == k[2] and j[7] == k[3]
                 and not (j[4] == j[2] and j[5] == j[3])
                 and not (j[4] == j[6] and j[5] == j[7]))
            ):
                intr[jintr] = intr[jintr] - intr[kintr]
                intr[kintr] = 0.0


def _reorder_interall_hermitian(
    nintr: int, intrindx: list, intr: list
) -> None:
    """Force Hermitian ordering on InterAll terms (Pass 2).

    For the two-body operator ``(c1† c2 c3† c4)†  = c4† c3 c2† c1``,
    certain index-pair patterns require reordering the 8 indices of
    *kintr* to match the Hermitian ordering of *jintr*.

    Two reorder patterns are detected:

    - **Pattern 1**: Direct Hermitian conjugate → reorder indices.
    - **Pattern 2**: Two sub-cases with sign flip → reorder indices
      and negate coefficient.

    Parameters
    ----------
    nintr : int
        Number of InterAll terms.
    intrindx : list
        Index array (8 integers per term), modified in place.
    intr : list
        Coefficient array (complex), modified in place.
    """
    for jintr in range(nintr):
        j = intrindx[jintr]
        for kintr in range(jintr + 1, nintr):
            k = intrindx[kintr]

            # Pattern 1: direct Hermitian conjugate match
            if (j[6] == k[4] and j[7] == k[5]
                    and j[4] == k[6] and j[5] == k[7]
                    and j[2] == k[0] and j[3] == k[1]
                    and j[0] == k[2] and j[1] == k[3]
                    and not (k[0] == k[6] and k[1] == k[7])
                    and not (k[2] == k[4] and k[3] == k[5])):
                intrindx[kintr][0] = j[6]
                intrindx[kintr][1] = j[7]
                intrindx[kintr][2] = j[4]
                intrindx[kintr][3] = j[5]
                intrindx[kintr][4] = j[2]
                intrindx[kintr][5] = j[3]
                intrindx[kintr][6] = j[0]
                intrindx[kintr][7] = j[1]

            # Pattern 2: two sub-cases, with sign flip
            elif (
                (j[6] == k[4] and j[7] == k[5]
                 and j[4] == k[2] and j[5] == k[3]
                 and j[2] == k[0] and j[3] == k[1]
                 and j[0] == k[6] and j[1] == k[7]
                 and not (k[2] == k[0] and k[3] == k[1])
                 and not (k[2] == k[4] and k[3] == k[5]))
                or
                (j[6] == k[0] and j[7] == k[1]
                 and j[4] == k[6] and j[5] == k[7]
                 and j[2] == k[4] and j[3] == k[5]
                 and j[0] == k[2] and j[1] == k[3]
                 and not (k[4] == k[2] and k[5] == k[3])
                 and not (k[4] == k[6] and k[5] == k[7]))
            ):
                intrindx[kintr][0] = j[6]
                intrindx[kintr][1] = j[7]
                intrindx[kintr][2] = j[4]
                intrindx[kintr][3] = j[5]
                intrindx[kintr][4] = j[2]
                intrindx[kintr][5] = j[3]
                intrindx[kintr][6] = j[0]
                intrindx[kintr][7] = j[1]

                intr[kintr] = -intr[kintr]


def _remove_interall_diagonal(
    nintr: int, intrindx: list, intr: list
) -> None:
    """Remove spurious diagonal InterAll terms (Pass 3).

    A term is zeroed if it has a diagonal pair — i.e.
    ``(site0, spin0) == (site4, spin4)`` or
    ``(site2, spin2) == (site6, spin6)`` — **unless** any of the four
    possible diagonal-matching conditions holds.

    Parameters
    ----------
    nintr : int
        Number of InterAll terms.
    intrindx : list
        Index array (8 integers per term).
    intr : list
        Coefficient array, modified in place.
    """
    for jintr in range(nintr):
        idx = intrindx[jintr]

        has_diagonal_pair = (
            (idx[0] == idx[4] and idx[1] == idx[5])
            or (idx[2] == idx[6] and idx[3] == idx[7])
        )

        if has_diagonal_pair:
            any_matching_diagonal = (
                (idx[0] == idx[2] and idx[1] == idx[3])
                or (idx[0] == idx[6] and idx[1] == idx[7])
                or (idx[4] == idx[2] and idx[5] == idx[3])
                or (idx[4] == idx[6] and idx[5] == idx[7])
            )
            if not any_matching_diagonal:
                intr[jintr] = 0.0


@dataclass
class InterAllData:
    """Complex two-body ``interall.def`` file.

    Attributes
    ----------
    rows : list of tuple
        ``(i0, s0, i1, s1, i2, s2, i3, s3, re, im)`` per non-zero term.
    """

    rows: list

    def write(self, directory: Path = Path(".")) -> None:
        lines = ["====================== \n",
                 f"NInterAll {len(self.rows):7d}  \n",
                 "====================== \n",
                 "========zInterAll===== \n",
                 "====================== \n"]
        for i0, s0, i1, s1, i2, s2, i3, s3, re, im in self.rows:
            lines.append(
                f"{i0:5d} {s0:5d} "
                f"{i1:5d} {s1:5d} "
                f"{i2:5d} {s2:5d} "
                f"{i3:5d} {s3:5d} "
                f"{re:25.15f}  {im:25.15f}\n"
            )
        with open(Path(directory) / "interall.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("    interall.def is written.")

    def to_dict(self) -> dict:
        return {"rows": [list(r) for r in self.rows]}

    @classmethod
    def from_dict(cls, data: dict) -> "InterAllData":
        return cls(rows=[tuple(r) for r in data["rows"]])


def _build_interall(StdI: StdIntList) -> "InterAllData | None":
    """Run the three InterAll merge passes, set ``Lintr``, and return data.

    Returns ``None`` (and sets ``Lintr`` to 0) when there are no non-zero
    terms or boost mode suppresses output.
    """
    nintr = len(StdI.intr_list)
    intr_indx = np.array([t[1:9] for t in StdI.intr_list], dtype=int).reshape(nintr, 8)
    intr_val = np.array([t[0] for t in StdI.intr_list], dtype=complex)
    _merge_interall_equivalent(nintr, intr_indx, intr_val)
    _reorder_interall_hermitian(nintr, intr_indx, intr_val)
    _remove_interall_diagonal(nintr, intr_indx, intr_val)

    nintr0 = _count_nonzero(nintr, intr_val)
    if nintr0 == 0 or StdI.lBoost == 1:
        StdI.Lintr = 0
        return None
    StdI.Lintr = 1

    rows = []
    for kintr in range(nintr):
        val = intr_val[kintr]
        if abs(val) > AMPLITUDE_EPS:
            i0, s0, i1, s1, i2, s2, i3, s3 = intr_indx[kintr]
            rows.append((int(i0), int(s0), int(i1), int(s1),
                         int(i2), int(s2), int(i3), int(s3),
                         float(val.real), float(val.imag)))
    return InterAllData(rows=rows)


def build_interactions(StdI: StdIntList):
    """Merge all interaction types, set the ``L*`` flags, and return data.

    Returns a list of :class:`InteractionData` (for each enabled 1-/2-index
    family) and an :class:`InterAllData` (when enabled).  The ``L*`` /
    ``Lintr`` flags are set on *StdI* as a side effect because the namelist
    writer reads them; only file writing is deferred to the returned
    objects' :meth:`write`.
    """
    out: list = []
    # Standard 1-/2-index families (data-driven)
    for spec in _INTERACTION_TYPES:
        data = _build_interaction(StdI, spec)
        if data is not None:
            out.append(data)
    # InterAll (three merge passes inside _build_interall)
    interall = _build_interall(StdI)
    if interall is not None:
        out.append(interall)
    return out
