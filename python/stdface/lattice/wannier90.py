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
import math
from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path
from typing import NamedTuple

import numpy as np

from ..core.stdface_vals import StdIntList, ModelType, SolverType, NaN_i, AMPLITUDE_EPS
from ..core.param_check import print_val_d, print_val_i, not_used_d
from .interaction_builder import (
    malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb,
)
from .site_util import init_site, find_site, set_local_spin_flags
from .wannier90_io import _geometry_w90, _read_w90, _read_density_matrix


# ---------------------------------------------------------------------------
#  Internal helpers
# ---------------------------------------------------------------------------


logger = logging.getLogger(__name__)


@dataclass
class UHFInitialData:
    """Initial UHF guess (``initial.def``).

    Attributes
    ----------
    rows : list of tuple
        ``(jsite, isite, re, im)`` per non-negligible entry; each row is
        written once per spin with the value halved at build time.
    """

    rows: list

    def write(self, directory: Path = Path(".")) -> None:
        lines = ["======================== \n",
                 f"NInitialGuess {len(self.rows) * 2:7d}  \n",
                 "======================== \n",
                 "========i_j_s_tijs====== \n",
                 "======================== \n"]
        for jsite, isite, re, im in self.rows:
            for ispin in range(2):
                lines.append(
                    f"{jsite:5d} {ispin:5d} {isite:5d} {ispin:5d} "
                    f"{re:25.15f} "
                    f"{im:25.15f}\n"
                )
        with open(Path(directory) / "initial.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("      initial.def is written.")

    def to_dict(self) -> dict:
        return {"rows": [list(r) for r in self.rows]}

    @classmethod
    def from_dict(cls, data: dict) -> "UHFInitialData":
        return cls(rows=[tuple(r) for r in data["rows"]])


@dataclass
class Wan2SiteData:
    """Wannier-orbital to super-cell-site mapping (``wan2site.dat``).

    Attributes
    ----------
    rows : list of tuple
        ``(isite, nx, ny, nz, iorb)`` per site.
    """

    rows: list

    def write(self, directory: Path = Path(".")) -> None:
        lines = ["======================== \n",
                 f"Total site number {len(self.rows):7d}  \n",
                 "======================== \n",
                 "========site nx ny nz norb====== \n",
                 "======================== \n"]
        for isite, nx, ny, nz, iorb in self.rows:
            lines.append(f"{isite:5d}{nx:5d}{ny:5d}{nz:5d}{iorb:5d}\n")
        with open(Path(directory) / "wan2site.dat", "w") as fp:
            fp.write("".join(lines))

    def to_dict(self) -> dict:
        return {"rows": [list(r) for r in self.rows]}

    @classmethod
    def from_dict(cls, data: dict) -> "Wan2SiteData":
        return cls(rows=[tuple(r) for r in data["rows"]])


def _build_uhf_initial(
    StdI: StdIntList,
    NtUJ: list[int],
    tUJ: list[np.ndarray],
    DenMat: dict[tuple[int, int, int], np.ndarray],
    tUJindx: list[np.ndarray],
) -> UHFInitialData:
    """Build the initial UHF guess (``initial.def``) as data.

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
    out_rows = []
    rows, cols = np.nonzero(mask)
    for isite, jsite in zip(rows, cols):
        val = 0.5 * IniGuess[isite, jsite]
        out_rows.append((int(jsite), int(isite),
                         float(val.real), float(val.imag)))
    return UHFInitialData(rows=out_rows)


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
    None: _DCMode.NOTCORRECT,
    "hartree":   _DCMode.HARTREE,
    "hartree_u": _DCMode.HARTREE_U,
    "full":      _DCMode.FULL,
}
"""Maps double-counting mode strings to :class:`_DCMode` enum members.

The sentinel ``None`` (unset) is treated the same as
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
        ``"full"``.  The sentinel ``None`` is treated the same
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
                cutoff_Vec[i, j] = print_val_d(
                    f"{cutoff_Vec_name}[{i}][{j}]",
                    cutoff_Vec[i, j], float(StdI.box[i, j]) * 0.5
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
                    StdI.trans_list.append(
                        (-tUJ[0][it], isite, ispin, isite, ispin))
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
            StdI.Cintra_list.append((
                tUJ[1][it].real,
                StdI.NsiteUC * kCell + int(tUJindx[1][it, 3]),
            ))

            # Double-counting correction
            if idcmode != _DCMode.NOTCORRECT:
                isite = StdI.NsiteUC * kCell + int(tUJindx[1][it, 3])
                for ispin in range(2):
                    DenMat0 = DenMat[(0, 0, 0)][
                        int(tUJindx[1][it, 3]), int(tUJindx[1][it, 3])
                    ]
                    StdI.trans_list.append(
                        (StdI.alpha * tUJ[1][it].real * DenMat0,
                         isite, ispin, isite, ispin))
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
                    StdI.trans_list.append(
                        (tUJ[1][it].real * DenMat0, isite, ispin, isite, ispin))

                    # U_{Rij} D_{0ii} (Local)
                    DenMat0 = DenMat[(0, 0, 0)][
                        int(tUJindx[1][it, 3]), int(tUJindx[1][it, 3])
                    ]
                    StdI.trans_list.append(
                        (tUJ[1][it].real * DenMat0, jsite, ispin, jsite, ispin))

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

            StdI.Hund_list.append((tUJ[2][it].real, isite, jsite))

            if StdI.model == ModelType.HUBBARD:
                StdI.Ex_list.append((tUJ[2][it].real, isite, jsite))
                StdI.PairHopp_list.append((tUJ[2][it].real, isite, jsite))

                # Double-counting correction
                if idcmode != _DCMode.NOTCORRECT and idcmode != _DCMode.HARTREE_U:
                    for ispin in range(2):
                        # -0.5 J_{Rij} D_{0jj}
                        DenMat0 = DenMat[(0, 0, 0)][
                            int(tUJindx[2][it, 4]), int(tUJindx[2][it, 4])
                        ]
                        StdI.trans_list.append(
                            (-(1.0 - StdI.alpha) * tUJ[2][it].real * DenMat0,
                             isite, ispin, isite, ispin))

                        # -0.5 J_{Rij} D_{0ii}
                        DenMat0 = DenMat[(0, 0, 0)][
                            int(tUJindx[2][it, 3]), int(tUJindx[2][it, 3])
                        ]
                        StdI.trans_list.append(
                            (-(1.0 - StdI.alpha) * tUJ[2][it].real * DenMat0,
                             jsite, ispin, jsite, ispin))

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
                    ex_val = tUJ[2][it].real
                else:
                    ex_val = -tUJ[2][it].real
                StdI.Ex_list.append((ex_val, isite, jsite))


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
    # Compute upper limit on the number of transfer terms (for pump arrays)
    if StdI.model == ModelType.SPIN:
        ntransMax = StdI.nsite * (StdI.S2 + 1 + 2 * StdI.S2)
    elif StdI.model == ModelType.HUBBARD:
        ntransMax = StdI.NCell * 2 * (
            2 * StdI.NsiteUC + NtUJ[0] * 2
            + NtUJ[1] * 2 * 3 + NtUJ[2] * 2 * 2
        )

    malloc_interactions(StdI, ntransMax)

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


def _build_wan2site(StdI: StdIntList) -> Wan2SiteData:
    """Build the Wannier-orbital to super-cell-site mapping as data.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing lattice and cell information.
    """
    rows = []
    for kCell in range(StdI.NCell):
        nx = int(StdI.Cell[kCell, 0])
        ny = int(StdI.Cell[kCell, 1])
        nz = int(StdI.Cell[kCell, 2])
        for it in range(StdI.NsiteUC):
            isite = StdI.NsiteUC * kCell + it
            rows.append((isite, nx, ny, nz, it))
    return Wan2SiteData(rows=rows)


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
    if StdI.lambda_ is None or math.isnan(StdI.lambda_):
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
    7. Build the auxiliary outputs (``initial.def`` when double-counting
       correction is on, and ``wan2site.dat``) into ``StdI._aux_outputs``;
       the main flow writes them alongside gnuplot / geometry / xsf.
    """
    NtUJ = [0, 0, 0]
    tUJ: list = [None, None, None]
    tUJindx: list = [None, None, None]

    # (1) Compute the shape of the super-cell and sites in the super-cell
    StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
    StdI.phase[1] = print_val_d("phase1", StdI.phase[1], 0.0)
    StdI.phase[2] = print_val_d("phase2", StdI.phase[2], 0.0)
    StdI.NsiteUC = 1
    init_site(StdI, 3)
    logger.info("\n  @ Wannier90 Geometry \n")
    _geometry_w90(StdI)

    _validate_interaction_params(StdI)
    idcmode = _parse_double_counting_mode(StdI.double_counting_mode)

    # Read hopping, Coulomb, and Hund interaction files
    hopping_R_defaults = (
        (StdI.W - 1) // 2 if StdI.W is not None else None,
        (StdI.L - 1) // 2 if StdI.L is not None else None,
        (StdI.Height - 1) // 2 if StdI.Height is not None else None,
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

    # (7) Auxiliary outputs are built here but written by the main flow
    # (StdI._aux_outputs, alongside gnuplot / geometry / xsf); lattice.xsf
    # itself is emitted on the lattice-level independent path (build_xsf).
    aux: list = []
    if idcmode != _DCMode.NOTCORRECT:
        aux.append(_build_uhf_initial(StdI, NtUJ, tUJ, DenMat, tUJindx))
    aux.append(_build_wan2site(StdI))
    # Extend (not assign): parse-stage builders (HPhi potential.dat) may
    # already have queued outputs here.
    StdI._aux_outputs = (StdI._aux_outputs or []) + aux


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
