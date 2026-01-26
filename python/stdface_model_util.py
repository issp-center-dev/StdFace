"""
Utility functions for constructing lattice models in Standard mode.

This module provides various utility functions for building interactions,
handling site labelling, input validation, and memory allocation that are
shared across all lattice implementations.

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
import sys
from typing import TextIO

import numpy as np

from stdface_vals import StdIntList


# ---------------------------------------------------------------------------
#  Exit wrapper
# ---------------------------------------------------------------------------


def exit_program(errorcode: int) -> None:
    """Terminate the program with the given error code.

    Parameters
    ----------
    errorcode : int
        Exit code passed to ``sys.exit``.

    Notes
    -----
    In the original C code this was ``StdFace_exit`` which wrapped
    ``MPI_Abort`` / ``MPI_Finalize``.  The Python version simply calls
    ``sys.exit``.
    """
    sys.stdout.flush()
    sys.stderr.flush()
    sys.exit(errorcode)


# ---------------------------------------------------------------------------
#  Transfer / one-body terms
# ---------------------------------------------------------------------------


def trans(
    StdI: StdIntList,
    trans0: complex,
    isite: int,
    ispin: int,
    jsite: int,
    jspin: int,
) -> None:
    """Add a transfer (one-body) term to the list.

    Appends to ``StdI.trans`` and ``StdI.transindx`` and increments
    ``StdI.ntrans``.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    trans0 : complex
        Hopping integral t, mu, etc.
    isite : int
        Site index i for c†_{i σ}.
    ispin : int
        Spin index σ for c†_{i σ}.
    jsite : int
        Site index j for c_{j σ'}.
    jspin : int
        Spin index σ' for c_{j σ'}.
    """
    if abs(trans0) < 1.0e-12:
        return
    StdI.trans[StdI.ntrans] = trans0
    StdI.transindx[StdI.ntrans][0] = isite
    StdI.transindx[StdI.ntrans][1] = ispin
    StdI.transindx[StdI.ntrans][2] = jsite
    StdI.transindx[StdI.ntrans][3] = jspin
    StdI.ntrans += 1


def hopping(
    StdI: StdIntList,
    trans0: complex,
    isite: int,
    jsite: int,
    dR: np.ndarray,
) -> None:
    """Add hopping for both spin channels.

    Both c†_{i σ} c_{j σ} and c†_{j σ} c_{i σ} for every spin channel
    (σ) are added.  For HPhi time-evolution pump mode, the pump arrays
    are populated instead.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    trans0 : complex
        Hopping integral t.
    isite : int
        Site index i.
    jsite : int
        Site index j.
    dR : numpy.ndarray
        Distance vector R_i - R_j (shape ``(3,)``).
    """
    if (StdI.solver == "HPhi"
            and StdI.method == "timeevolution"
            and StdI.PumpBody == 1):
        for it in range(StdI.Lanczos_max):
            Cphase = 0.0
            for ii in range(3):
                Cphase += StdI.At[it][ii] * dR[ii]
            coef = math.cos(Cphase) + 1j * math.sin(-Cphase)
            for ispin in range(2):
                StdI.pump[it][StdI.npump[it]] = coef * trans0
                StdI.pumpindx[it][StdI.npump[it]][0] = isite
                StdI.pumpindx[it][StdI.npump[it]][1] = ispin
                StdI.pumpindx[it][StdI.npump[it]][2] = jsite
                StdI.pumpindx[it][StdI.npump[it]][3] = ispin
                StdI.npump[it] += 1

                StdI.pump[it][StdI.npump[it]] = np.conj(coef * trans0)
                StdI.pumpindx[it][StdI.npump[it]][0] = jsite
                StdI.pumpindx[it][StdI.npump[it]][1] = ispin
                StdI.pumpindx[it][StdI.npump[it]][2] = isite
                StdI.pumpindx[it][StdI.npump[it]][3] = ispin
                StdI.npump[it] += 1
    else:
        for ispin in range(2):
            trans(StdI, trans0, jsite, ispin, isite, ispin)
            trans(StdI, np.conj(trans0), isite, ispin, jsite, ispin)


def hubbard_local(
    StdI: StdIntList,
    mu0: float,
    h0: float,
    Gamma0: float,
    Gamma0_y: float,
    U0: float,
    isite: int,
) -> None:
    """Add intra-Coulomb, magnetic field and chemical potential for itinerant electrons.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    mu0 : float
        Chemical potential.
    h0 : float
        Longitudinal magnetic field.
    Gamma0 : float
        Transverse magnetic field (x).
    Gamma0_y : float
        Transverse magnetic field (y).
    U0 : float
        Intra-site Coulomb potential.
    isite : int
        Site index.
    """
    trans(StdI, mu0 - 0.5 * h0, isite, 0, isite, 0)
    trans(StdI, mu0 + 0.5 * h0, isite, 1, isite, 1)
    trans(StdI, -0.5 * Gamma0, isite, 1, isite, 0)
    trans(StdI, -0.5 * Gamma0, isite, 0, isite, 1)
    trans(StdI, -0.5 * 1j * Gamma0_y, isite, 1, isite, 0)
    trans(StdI, 0.5 * 1j * Gamma0_y, isite, 0, isite, 1)

    StdI.Cintra[StdI.NCintra] = U0
    StdI.CintraIndx[StdI.NCintra][0] = isite
    StdI.NCintra += 1


def mag_field(
    StdI: StdIntList,
    S2: int,
    h: float,
    Gamma: float,
    Gamma_y: float,
    isite: int,
) -> None:
    """Add longitudinal and transverse magnetic field terms.

    Uses the Bogoliubov representation for arbitrary spin S.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    S2 : int
        Twice the spin moment (2S) at site *isite*.
    h : float
        Longitudinal magnetic field.
    Gamma : float
        Transverse magnetic field (x).
    Gamma_y : float
        Transverse magnetic field (y).
    isite : int
        Site index.
    """
    S = S2 * 0.5
    for ispin in range(S2 + 1):
        Sz = S - float(ispin)
        # Longitudinal part: -h σ c†_{i σ} c_{i σ}
        trans(StdI, -h * Sz, isite, ispin, isite, ispin)
        # Transverse part
        if ispin > 0:
            factor = math.sqrt(S * (S + 1.0) - Sz * (Sz + 1.0))
            trans(StdI, -0.5 * Gamma * factor - 0.5 * 1j * Gamma_y * factor,
                  isite, ispin, isite, ispin - 1)
            trans(StdI, -0.5 * Gamma * factor + 0.5 * 1j * Gamma_y * factor,
                  isite, ispin - 1, isite, ispin)


# ---------------------------------------------------------------------------
#  Interaction / two-body terms
# ---------------------------------------------------------------------------


def intr(
    StdI: StdIntList,
    intr0: complex,
    site1: int, spin1: int,
    site2: int, spin2: int,
    site3: int, spin3: int,
    site4: int, spin4: int,
) -> None:
    """Add a general two-body (InterAll) interaction term to the list.

    Appends to ``StdI.intr`` and ``StdI.intrindx`` and increments
    ``StdI.nintr``.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    intr0 : complex
        Interaction coefficient U, V, J, etc.
    site1, spin1 : int
        Indices for c†_{i₁ σ₁}.
    site2, spin2 : int
        Indices for c_{i₂ σ₂}.
    site3, spin3 : int
        Indices for c†_{i₃ σ₃}.
    site4, spin4 : int
        Indices for c_{i₄ σ₄}.
    """
    if abs(intr0) < 1.0e-12:
        return
    StdI.intr[StdI.nintr] = intr0
    idx = StdI.intrindx[StdI.nintr]
    idx[0] = site1; idx[1] = spin1
    idx[2] = site2; idx[3] = spin2
    idx[4] = site3; idx[5] = spin3
    idx[6] = site4; idx[7] = spin4
    StdI.nintr += 1


def general_j(
    StdI: StdIntList,
    J: np.ndarray,
    Si2: int,
    Sj2: int,
    isite: int,
    jsite: int,
) -> None:
    """Treat J as a 3×3 matrix for general spin interactions.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    J : numpy.ndarray
        3×3 spin interaction matrix.
    Si2 : int
        Twice the spin moment (2S) at site *isite*.
    Sj2 : int
        Twice the spin moment (2S) at site *jsite*.
    isite : int
        First site index.
    jsite : int
        Second site index.
    """
    ZGeneral = 1
    ExGeneral = 1

    if Si2 == 1 or Sj2 == 1:
        ZGeneral = 0

        StdI.Hund[StdI.NHund] = -0.5 * J[2, 2]
        StdI.HundIndx[StdI.NHund][0] = isite
        StdI.HundIndx[StdI.NHund][1] = jsite
        StdI.NHund += 1

        StdI.Cinter[StdI.NCinter] = -0.25 * J[2, 2]
        StdI.CinterIndx[StdI.NCinter][0] = isite
        StdI.CinterIndx[StdI.NCinter][1] = jsite
        StdI.NCinter += 1

        cond_offdiag = (abs(J[0, 1]) < 1e-6 and abs(J[1, 0]) < 1e-6)
        if StdI.solver == "mVMC":
            cond_offdiag = cond_offdiag and (abs(J[0, 0] - J[1, 1]) < 1e-6)

        if cond_offdiag:
            ExGeneral = 0

            if StdI.solver == "mVMC":
                StdI.Ex[StdI.NEx] = -0.25 * (J[0, 0] + J[1, 1])
            else:
                if StdI.model == "kondo":
                    StdI.Ex[StdI.NEx] = -0.25 * (J[0, 0] + J[1, 1])
                else:
                    StdI.Ex[StdI.NEx] = 0.25 * (J[0, 0] + J[1, 1])
            StdI.ExIndx[StdI.NEx][0] = isite
            StdI.ExIndx[StdI.NEx][1] = jsite
            StdI.NEx += 1

            StdI.PairLift[StdI.NPairLift] = 0.25 * (J[0, 0] - J[1, 1])
            StdI.PLIndx[StdI.NPairLift][0] = isite
            StdI.PLIndx[StdI.NPairLift][1] = jsite
            StdI.NPairLift += 1

    Si = 0.5 * Si2
    Sj = 0.5 * Sj2

    for ispin in range(Si2 + 1):
        Siz = Si - float(ispin)
        for jspin in range(Sj2 + 1):
            Sjz = Sj - float(jspin)

            # (1) J_z S_{iz} S_{jz}
            if ZGeneral == 1:
                intr0 = J[2, 2] * Siz * Sjz
                intr(StdI, intr0,
                     isite, ispin, isite, ispin,
                     jsite, jspin, jsite, jspin)

            if ispin > 0 and jspin > 0 and ExGeneral == 1:
                # (2) S_i^+ S_j^- + h.c.
                intr0 = (0.25 * (J[0, 0] + J[1, 1] + 1j * (J[0, 1] - J[1, 0]))
                         * math.sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0))
                         * math.sqrt(Sj * (Sj + 1.0) - Sjz * (Sjz + 1.0)))
                intr(StdI, intr0,
                     isite, ispin - 1, isite, ispin,
                     jsite, jspin, jsite, jspin - 1)
                intr(StdI, np.conj(intr0),
                     isite, ispin, isite, ispin - 1,
                     jsite, jspin - 1, jsite, jspin)

                # (3) S_i^+ S_j^+ + h.c.
                intr0 = (0.25 * (J[0, 0] - J[1, 1] - 1j * (J[0, 1] + J[1, 0]))
                         * math.sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0))
                         * math.sqrt(Sj * (Sj + 1.0) - Sjz * (Sjz + 1.0)))
                intr(StdI, intr0,
                     isite, ispin - 1, isite, ispin,
                     jsite, jspin - 1, jsite, jspin)
                intr(StdI, np.conj(intr0),
                     isite, ispin, isite, ispin - 1,
                     jsite, jspin, jsite, jspin - 1)

            # (4) S_i^+ S_{jz} + h.c.
            if ispin > 0:
                intr0 = (0.5 * (J[0, 2] - 1j * J[1, 2])
                         * math.sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0)) * Sjz)
                intr(StdI, intr0,
                     isite, ispin - 1, isite, ispin,
                     jsite, jspin, jsite, jspin)
                intr(StdI, np.conj(intr0),
                     jsite, jspin, jsite, jspin,
                     isite, ispin, isite, ispin - 1)

            # (5) S_{iz} S_j^+ + h.c.
            if jspin > 0:
                intr0 = (0.5 * (J[2, 0] - 1j * J[2, 1])
                         * Siz
                         * math.sqrt(Sj * (Sj + 1.0) - Sjz * (Sjz + 1.0)))
                intr(StdI, intr0,
                     isite, ispin, isite, ispin,
                     jsite, jspin - 1, jsite, jspin)
                intr(StdI, np.conj(intr0),
                     jsite, jspin, jsite, jspin - 1,
                     isite, ispin, isite, ispin)


def coulomb(StdI: StdIntList, V: float, isite: int, jsite: int) -> None:
    """Add an onsite/offsite Coulomb interaction term.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    V : float
        Coulomb integral U, V, etc.
    isite : int
        First site index.
    jsite : int
        Second site index.
    """
    StdI.Cinter[StdI.NCinter] = V
    StdI.CinterIndx[StdI.NCinter][0] = isite
    StdI.CinterIndx[StdI.NCinter][1] = jsite
    StdI.NCinter += 1


# ---------------------------------------------------------------------------
#  PrintVal / NotUsed / Required helpers
# ---------------------------------------------------------------------------


def print_val_d(valname: str, val: float, val0: float) -> float:
    """Print and optionally set a real-valued parameter.

    If *val* is NaN, set it to the default *val0* and print with a
    ``DEFAULT VALUE`` tag.

    Parameters
    ----------
    valname : str
        Name of the variable (for display).
    val : float
        Current value (may be NaN if not specified).
    val0 : float
        Default value to use when *val* is NaN.

    Returns
    -------
    float
        The (possibly updated) value.
    """
    if math.isnan(val):
        val = val0
        print(f"  {valname:>15s} = {val:<10.5f}  ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"  {valname:>15s} = {val:<10.5f}")
    return val


def print_val_dd(valname: str, val: float, val0: float, val1: float) -> float:
    """Print and optionally set a real-valued parameter with two defaults.

    If *val* is NaN, use the primary default *val0* if specified,
    otherwise fall back to the secondary default *val1*.

    Parameters
    ----------
    valname : str
        Name of the variable (for display).
    val : float
        Current value (may be NaN).
    val0 : float
        Primary default (may itself be NaN).
    val1 : float
        Secondary default.

    Returns
    -------
    float
        The (possibly updated) value.
    """
    if math.isnan(val):
        if math.isnan(val0):
            val = val1
        else:
            val = val0
        print(f"  {valname:>15s} = {val:<10.5f}  ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"  {valname:>15s} = {val:<10.5f}")
    return val


def print_val_c(valname: str, val: complex, val0: complex) -> complex:
    """Print and optionally set a complex-valued parameter.

    Parameters
    ----------
    valname : str
        Name of the variable (for display).
    val : complex
        Current value (real part NaN if not specified).
    val0 : complex
        Default value.

    Returns
    -------
    complex
        The (possibly updated) value.
    """
    if math.isnan(val.real):
        val = val0
        print(f"  {valname:>15s} = {val.real:<10.5f} {val.imag:<10.5f}"
              f"  ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"  {valname:>15s} = {val.real:<10.5f} {val.imag:<10.5f}")
    return val


def print_val_i(valname: str, val: int, val0: int) -> int:
    """Print and optionally set an integer parameter.

    If *val* equals the sentinel ``2147483647`` (NaN_i), set it to *val0*.

    Parameters
    ----------
    valname : str
        Name of the variable (for display).
    val : int
        Current value (sentinel if not specified).
    val0 : int
        Default value.

    Returns
    -------
    int
        The (possibly updated) value.
    """
    NaN_i = 2147483647
    if val == NaN_i:
        val = val0
        print(f"  {valname:>15s} = {val:<10d}  ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"  {valname:>15s} = {val:<10d}")
    return val


def not_used_d(valname: str, val: float | complex) -> None:
    """Abort if a real parameter is specified but will not be used.

    Parameters
    ----------
    valname : str
        Name of the variable.
    val : float or complex
        Value to check (abort if not NaN).  If complex, the real part
        is tested (matching the C behaviour of implicitly casting
        ``double complex`` to ``double``).
    """
    check = val.real if isinstance(val, complex) else val
    if not math.isnan(check):
        print(f"\n Check !  {valname} is SPECIFIED but will NOT be USED. ")
        print("            Please COMMENT-OUT this line ")
        print("            or check this input is REALLY APPROPRIATE for your purpose !\n")
        exit_program(-1)


def not_used_c(valname: str, val: complex) -> None:
    """Abort if a complex parameter is specified but will not be used.

    Parameters
    ----------
    valname : str
        Name of the variable.
    val : complex
        Value to check (abort if real part is not NaN).
    """
    if not math.isnan(val.real):
        print(f"\n Check !  {valname} is SPECIFIED but will NOT be USED. ")
        print("            Please COMMENT-OUT this line ")
        print("            or check this input is REALLY APPROPRIATE for your purpose !\n")
        exit_program(-1)


def not_used_j(valname: str, JAll: float, J: np.ndarray) -> None:
    """Abort if any component of a J-type interaction is specified but unused.

    Parameters
    ----------
    valname : str
        Base name of the variable (e.g. ``"J0"``).
    JAll : float
        Scalar (isotropic) part.
    J : numpy.ndarray
        3×3 matrix of anisotropic components.
    """
    suffixes = [["x", "xy", "xz"],
                ["yx", "y", "yz"],
                ["zx", "zy", "z"]]
    not_used_d(valname, JAll)
    for i1 in range(3):
        for i2 in range(3):
            not_used_d(f"{valname}{suffixes[i1][i2]}", J[i1, i2])


def not_used_i(valname: str, val: int) -> None:
    """Abort if an integer parameter is specified but will not be used.

    Parameters
    ----------
    valname : str
        Name of the variable.
    val : int
        Value to check (abort if not the sentinel 2147483647).
    """
    NaN_i = 2147483647
    if val != NaN_i:
        print(f"\n Check !  {valname} is SPECIFIED but will NOT be USED. ")
        print("            Please COMMENT-OUT this line ")
        print("            or check this input is REALLY APPROPRIATE for your purpose !\n")
        exit_program(-1)


def required_val_i(valname: str, val: int) -> None:
    """Abort if a required integer parameter is missing.

    Parameters
    ----------
    valname : str
        Name of the variable.
    val : int
        Value to check (abort if equals the sentinel 2147483647).
    """
    NaN_i = 2147483647
    if val == NaN_i:
        print(f"ERROR ! {valname} is NOT specified !")
        exit_program(-1)
    else:
        print(f"  {valname:>15s} = {val:<3d}")


# ---------------------------------------------------------------------------
#  Fold / Site utilities
# ---------------------------------------------------------------------------


def _fold_site(
    StdI: StdIntList,
    iCellV: list[int],
) -> tuple[list[int], list[int]]:
    """Move a site into the original super-cell if it is outside.

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
    # (1) Transform to fractional coordinate (times NCell)
    iCellV_frac = [0, 0, 0]
    for ii in range(3):
        for jj in range(3):
            iCellV_frac[ii] += StdI.rbox[ii, jj] * iCellV[jj]

    # (2) Search which super-cell contains this cell
    nBox = [0, 0, 0]
    for ii in range(3):
        nBox[ii] = (iCellV_frac[ii] + StdI.NCell * 1000) // StdI.NCell - 1000

    # (3) Fractional coordinate in the original super-cell
    for ii in range(3):
        iCellV_frac[ii] -= StdI.NCell * nBox[ii]

    iCellV_fold = [0, 0, 0]
    for ii in range(3):
        for jj in range(3):
            iCellV_fold[ii] += StdI.box[jj, ii] * iCellV_frac[jj]
        iCellV_fold[ii] = (iCellV_fold[ii] + StdI.NCell * 1000) // StdI.NCell - 1000

    return nBox, iCellV_fold


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
    NaN_i = StdI.NaN_i

    print("\n  @ Super-Lattice setting\n")

    # (1) Check input parameters about the shape of super-cell
    lwh_specified = (StdI.L != NaN_i or StdI.W != NaN_i or StdI.Height != NaN_i)
    box_specified = any(StdI.box[i, j] != NaN_i for i in range(3) for j in range(3))

    if lwh_specified and box_specified:
        print("\nERROR ! (L, W, Height) and (a0W, ..., a2H) conflict !\n")
        exit_program(-1)
    elif lwh_specified:
        StdI.L = print_val_i("L", StdI.L, 1)
        StdI.W = print_val_i("W", StdI.W, 1)
        StdI.Height = print_val_i("Height", StdI.Height, 1)
        StdI.box[:, :] = 0
        StdI.box[0, 0] = StdI.W
        StdI.box[1, 1] = StdI.L
        StdI.box[2, 2] = StdI.Height
    else:
        StdI.box[0, 0] = print_val_i("a0W", int(StdI.box[0, 0]), 1)
        StdI.box[0, 1] = print_val_i("a0L", int(StdI.box[0, 1]), 0)
        StdI.box[0, 2] = print_val_i("a0H", int(StdI.box[0, 2]), 0)
        StdI.box[1, 0] = print_val_i("a1W", int(StdI.box[1, 0]), 0)
        StdI.box[1, 1] = print_val_i("a1L", int(StdI.box[1, 1]), 1)
        StdI.box[1, 2] = print_val_i("a1H", int(StdI.box[1, 2]), 0)
        StdI.box[2, 0] = print_val_i("a2W", int(StdI.box[2, 0]), 0)
        StdI.box[2, 1] = print_val_i("a2L", int(StdI.box[2, 1]), 0)
        StdI.box[2, 2] = print_val_i("a2H", int(StdI.box[2, 2]), 1)

    if dim == 2:
        StdI.direct[0, 2] = 0.0
        StdI.direct[1, 2] = 0.0
        StdI.direct[2, 0] = 0.0
        StdI.direct[2, 1] = 0.0
        StdI.direct[2, 2] = 1.0

    # (2) Define the phase factor at each boundary
    if dim == 2:
        StdI.phase[2] = 0.0
    for ii in range(3):
        StdI.ExpPhase[ii] = (math.cos(StdI.pi180 * StdI.phase[ii])
                             + 1j * math.sin(StdI.pi180 * StdI.phase[ii]))
        if abs(StdI.ExpPhase[ii] + 1.0) < 1e-6:
            StdI.AntiPeriod[ii] = 1
        else:
            StdI.AntiPeriod[ii] = 0

    # (3) Allocate tau (intrinsic structure of unit-cell)
    StdI.tau = np.zeros((StdI.NsiteUC, 3))

    # (4) Calculate reciprocal lattice vectors and NCell
    StdI.NCell = 0
    for ii in range(3):
        StdI.NCell += (int(StdI.box[0, ii])
                       * int(StdI.box[1, (ii + 1) % 3])
                       * int(StdI.box[2, (ii + 2) % 3])
                       - int(StdI.box[0, ii])
                       * int(StdI.box[1, (ii + 2) % 3])
                       * int(StdI.box[2, (ii + 1) % 3]))
    print(f"   Number of Cell = {abs(StdI.NCell)}")
    if StdI.NCell == 0:
        exit_program(-1)

    for ii in range(3):
        for jj in range(3):
            StdI.rbox[ii, jj] = (int(StdI.box[(ii + 1) % 3, (jj + 1) % 3])
                                 * int(StdI.box[(ii + 2) % 3, (jj + 2) % 3])
                                 - int(StdI.box[(ii + 1) % 3, (jj + 2) % 3])
                                 * int(StdI.box[(ii + 2) % 3, (jj + 1) % 3]))
    if StdI.NCell < 0:
        StdI.rbox *= -1
        StdI.NCell *= -1

    # (5) Find cells in the super-cell
    # (5-1) Find bounds
    bound = [[0, 0], [0, 0], [0, 0]]
    for ii in range(3):
        for n2 in range(2):
            for n1 in range(2):
                for n0 in range(2):
                    nBox = [n0, n1, n2]
                    edge = sum(nBox[jj] * int(StdI.box[jj, ii]) for jj in range(3))
                    if edge < bound[ii][0]:
                        bound[ii][0] = edge
                    if edge > bound[ii][1]:
                        bound[ii][1] = edge

    # (5-2) Find cells in the super-cell
    StdI.Cell = np.zeros((StdI.NCell, 3), dtype=int)
    jj_idx = 0
    for ic2 in range(bound[2][0], bound[2][1] + 1):
        for ic1 in range(bound[1][0], bound[1][1] + 1):
            for ic0 in range(bound[0][0], bound[0][1] + 1):
                iCellV = [ic0, ic1, ic2]
                nBox, iCellV_fold = _fold_site(StdI, iCellV)
                if nBox[0] == 0 and nBox[1] == 0 and nBox[2] == 0:
                    StdI.Cell[jj_idx, :] = iCellV
                    jj_idx += 1

    # (6) For 2D, print lattice.gp header
    if dim == 2 and fp is not None:
        pos = np.zeros((4, 2))
        pos[1, 0] = StdI.direct[0, 0] * StdI.box[0, 0] + StdI.direct[1, 0] * StdI.box[0, 1]
        pos[1, 1] = StdI.direct[0, 1] * StdI.box[0, 0] + StdI.direct[1, 1] * StdI.box[0, 1]
        pos[2, 0] = StdI.direct[0, 0] * StdI.box[1, 0] + StdI.direct[1, 0] * StdI.box[1, 1]
        pos[2, 1] = StdI.direct[0, 1] * StdI.box[1, 0] + StdI.direct[1, 1] * StdI.box[1, 1]
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
        for i in range(4):
            j = (i + 1) if i < 3 else 2
            if i == 0:
                fp.write(f"set arrow from {pos[0][0]:f}, {pos[0][1]:f} "
                         f"to {pos[1][0]:f}, {pos[1][1]:f} nohead front ls 3\n")
            elif i == 1:
                fp.write(f"set arrow from {pos[1][0]:f}, {pos[1][1]:f} "
                         f"to {pos[3][0]:f}, {pos[3][1]:f} nohead front ls 3\n")
            elif i == 2:
                fp.write(f"set arrow from {pos[3][0]:f}, {pos[3][1]:f} "
                         f"to {pos[2][0]:f}, {pos[2][1]:f} nohead front ls 3\n")
            elif i == 3:
                fp.write(f"set arrow from {pos[2][0]:f}, {pos[2][1]:f} "
                         f"to {pos[0][0]:f}, {pos[0][1]:f} nohead front ls 3\n")


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
    dR = np.zeros(3)
    dR[0] = -float(diW) + StdI.tau[isiteUC, 0] - StdI.tau[jsiteUC, 0]
    dR[1] = -float(diL) + StdI.tau[isiteUC, 1] - StdI.tau[jsiteUC, 1]
    dR[2] = -float(diH) + StdI.tau[isiteUC, 2] - StdI.tau[jsiteUC, 2]

    jCellV = [iW + diW, iL + diL, iH + diH]
    nBox, jCellV = _fold_site(StdI, jCellV)
    Cphase = 1.0 + 0j
    for ii in range(3):
        Cphase *= StdI.ExpPhase[ii] ** nBox[ii]

    iCell = 0
    jCell = 0
    for kCell in range(StdI.NCell):
        if (jCellV[0] == StdI.Cell[kCell, 0]
                and jCellV[1] == StdI.Cell[kCell, 1]
                and jCellV[2] == StdI.Cell[kCell, 2]):
            jCell = kCell
        if (iW == StdI.Cell[kCell, 0]
                and iL == StdI.Cell[kCell, 1]
                and iH == StdI.Cell[kCell, 2]):
            iCell = kCell

    isite = iCell * StdI.NsiteUC + isiteUC
    jsite = jCell * StdI.NsiteUC + jsiteUC
    if StdI.model == "kondo":
        isite += StdI.NCell * StdI.NsiteUC
        jsite += StdI.NCell * StdI.NsiteUC

    return isite, jsite, Cphase, dR


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

    xi = (StdI.direct[0, 0] * (iW + StdI.tau[jsiteUC, 0])
          + StdI.direct[1, 0] * (iL + StdI.tau[jsiteUC, 1]))
    yi = (StdI.direct[0, 1] * (iW + StdI.tau[jsiteUC, 0])
          + StdI.direct[1, 1] * (iL + StdI.tau[jsiteUC, 1]))
    xj = (StdI.direct[0, 0] * (iW - diW + StdI.tau[isiteUC, 0])
          + StdI.direct[1, 0] * (iL - diL + StdI.tau[isiteUC, 1]))
    yj = (StdI.direct[0, 1] * (iW - diW + StdI.tau[isiteUC, 0])
          + StdI.direct[1, 1] * (iL - diL + StdI.tau[isiteUC, 1]))

    if fp is not None:
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

    # Then print the normal one
    isite, jsite, Cphase, dR = find_site(
        StdI, iW, iL, 0, diW, diL, 0, isiteUC, jsiteUC)

    xi = (StdI.direct[1, 0] * (iL + StdI.tau[isiteUC, 1])
          + StdI.direct[0, 0] * (iW + StdI.tau[isiteUC, 0]))
    yi = (StdI.direct[1, 1] * (iL + StdI.tau[isiteUC, 1])
          + StdI.direct[0, 1] * (iW + StdI.tau[isiteUC, 0]))
    xj = (StdI.direct[0, 0] * (iW + diW + StdI.tau[jsiteUC, 0])
          + StdI.direct[1, 0] * (iL + diL + StdI.tau[jsiteUC, 1]))
    yj = (StdI.direct[0, 1] * (iW + diW + StdI.tau[jsiteUC, 0])
          + StdI.direct[1, 1] * (iL + diL + StdI.tau[jsiteUC, 1]))

    if fp is not None:
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

    return isite, jsite, Cphase, dR


# ---------------------------------------------------------------------------
#  XSF / geometry output
# ---------------------------------------------------------------------------


def print_xsf(StdI: StdIntList) -> None:
    """Print lattice.xsf file (XCrysDen format).

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.
    """
    do_convvec = StdI.lattice in (
        "orthorhombic", "face-centeredorthorhombic",
        "fcorthorhombic", "fco", "pyrochlore")

    with open("lattice.xsf", "w") as fp:
        fp.write("CRYSTAL\n")
        fp.write("PRIMVEC\n")
        for ii in range(3):
            vec = [0.0, 0.0, 0.0]
            for jj in range(3):
                for kk in range(3):
                    vec[jj] += float(StdI.box[ii, kk]) * StdI.direct[kk, jj]
            fp.write(f"{vec[0]:15.5f} {vec[1]:15.5f} {vec[2]:15.5f}\n")

        if do_convvec:
            fp.write("CONVVEC\n")
            for ii in range(3):
                row = [0.0, 0.0, 0.0]
                for jj in range(3):
                    if ii == jj:
                        row[jj] = StdI.length[ii]
                fp.write(f"{row[0]:15.5f} {row[1]:15.5f} {row[2]:15.5f}\n")

        fp.write("PRIMCOORD\n")
        fp.write(f"{StdI.NCell * StdI.NsiteUC} 1\n")
        for iCell in range(StdI.NCell):
            for isite in range(StdI.NsiteUC):
                vec = [0.0, 0.0, 0.0]
                for jj in range(3):
                    for kk in range(3):
                        vec[jj] += ((float(StdI.Cell[iCell, kk])
                                     + StdI.tau[isite, kk])
                                    * StdI.direct[kk, jj])
                fp.write(f"H {vec[0]:15.5f} {vec[1]:15.5f} {vec[2]:15.5f}\n")


# ---------------------------------------------------------------------------
#  Input spin / Coulomb / hopping helpers
# ---------------------------------------------------------------------------


def input_spin_nn(
    J: np.ndarray,
    JAll: float,
    J0: np.ndarray,
    J0All: float,
    J0name: str,
) -> None:
    """Input nearest-neighbour spin-spin interaction.

    Resolves conflicts between isotropic and anisotropic specifications
    and fills in the output matrix *J0* in-place.

    Parameters
    ----------
    J : numpy.ndarray
        Anisotropic spin interaction (3×3).
    JAll : float
        Isotropic interaction.
    J0 : numpy.ndarray
        Output anisotropic spin interaction (3×3, modified in-place).
    J0All : float
        Isotropic interaction for this bond.
    J0name : str
        Name of this spin interaction (e.g. ``"J1"``).
    """
    suffixes = [["x", "xy", "xz"],
                ["yx", "y", "yz"],
                ["zx", "zy", "z"]]

    # Conflict checks
    if not math.isnan(JAll) and not math.isnan(J0All):
        print(f"\n ERROR! {J0name} conflict !\n")
        exit_program(-1)

    for i1 in range(3):
        for i2 in range(3):
            if not math.isnan(JAll) and not math.isnan(J[i1, i2]):
                print(f"\n ERROR! J{suffixes[i1][i2]} conflict !\n")
                exit_program(-1)
            elif not math.isnan(J0All) and not math.isnan(J[i1, i2]):
                print(f"\n ERROR! {J0name} and J{suffixes[i1][i2]} conflict !\n")
                exit_program(-1)
            elif not math.isnan(J0All) and not math.isnan(J0[i1, i2]):
                print(f"\n ERROR! {J0name} and {J0name}{suffixes[i1][i2]} conflict !\n")
                exit_program(-1)
            elif not math.isnan(J0[i1, i2]) and not math.isnan(JAll):
                print(f"\n ERROR! {J0name}{suffixes[i1][i2]} conflict !\n")
                exit_program(-1)

    for i1 in range(3):
        for i2 in range(3):
            for i3 in range(3):
                for i4 in range(3):
                    if not math.isnan(J0[i1, i2]) and not math.isnan(J[i3, i4]):
                        print(f"\n ERROR! {J0name}{suffixes[i1][i2]} "
                              f"and J{suffixes[i3][i4]} conflict !\n")
                        exit_program(-1)

    # Set values
    for i1 in range(3):
        for i2 in range(3):
            if not math.isnan(J0[i1, i2]):
                print(f"  {J0name + suffixes[i1][i2]:>14s} = {J0[i1, i2]:<10.5f}")
            elif not math.isnan(J[i1, i2]):
                J0[i1, i2] = J[i1, i2]
                print(f"  {J0name + suffixes[i1][i2]:>14s} = {J0[i1, i2]:<10.5f}")
            elif i1 == i2 and not math.isnan(J0All):
                J0[i1, i2] = J0All
                print(f"  {J0name + suffixes[i1][i2]:>14s} = {J0[i1, i2]:<10.5f}")
            elif i1 == i2 and not math.isnan(JAll):
                J0[i1, i2] = JAll
                print(f"  {J0name + suffixes[i1][i2]:>14s} = {J0[i1, i2]:<10.5f}")
            else:
                J0[i1, i2] = 0.0


def input_spin(Jp: np.ndarray, JpAll: float, Jpname: str) -> None:
    """Input spin-spin interaction other than nearest-neighbour.

    Parameters
    ----------
    Jp : numpy.ndarray
        Fully anisotropic spin interaction (3×3, modified in-place).
    JpAll : float
        Isotropic interaction value.
    Jpname : str
        Name of this spin interaction (e.g. ``"J'"``).
    """
    suffixes = [["x", "xy", "xz"],
                ["yx", "y", "yz"],
                ["zx", "zy", "z"]]

    for i1 in range(3):
        for i2 in range(3):
            if not math.isnan(JpAll) and not math.isnan(Jp[i1, i2]):
                print(f"\n ERROR! {Jpname} and {Jpname}{suffixes[i1][i2]} conflict !\n")
                exit_program(-1)

    for i1 in range(3):
        for i2 in range(3):
            if not math.isnan(Jp[i1, i2]):
                print(f"  {Jpname + suffixes[i1][i2]:>14s} = {Jp[i1, i2]:<10.5f}")
            elif i1 == i2 and not math.isnan(JpAll):
                Jp[i1, i2] = JpAll
                print(f"  {Jpname + suffixes[i1][i2]:>14s} = {Jp[i1, i2]:<10.5f}")
            else:
                Jp[i1, i2] = 0.0


def input_coulomb_v(V: float, V0: float, V0name: str) -> float:
    """Input off-site Coulomb interaction from the input file.

    Parameters
    ----------
    V : float
        Isotropic Coulomb interaction (may be NaN).
    V0 : float
        Specific Coulomb parameter (may be NaN).
    V0name : str
        Name of the parameter (e.g. ``"V1"``).

    Returns
    -------
    float
        The resolved value of *V0*.
    """
    if not math.isnan(V) and not math.isnan(V0):
        print(f"\n ERROR! {V0name} conflicts !\n")
        exit_program(-1)
    elif not math.isnan(V0):
        print(f"  {V0name:>15s} = {V0:<10.5f}")
    elif not math.isnan(V):
        V0 = V
        print(f"  {V0name:>15s} = {V0:<10.5f}")
    else:
        V0 = 0.0
    return V0


def input_hopp(t: complex, t0: complex, t0name: str) -> complex:
    """Input hopping integral from the input file.

    Parameters
    ----------
    t : complex
        Isotropic hopping (may have NaN real part).
    t0 : complex
        Specific hopping parameter (may have NaN real part).
    t0name : str
        Name of the parameter (e.g. ``"t1"``).

    Returns
    -------
    complex
        The resolved value of *t0*.
    """
    if not math.isnan(t.real) and not math.isnan(t0.real):
        print(f"\n ERROR! {t0name} conflicts !\n")
        exit_program(-1)
    elif not math.isnan(t0.real):
        print(f"  {t0name:>15s} = {t0.real:<10.5f} {t0.imag:<10.5f}")
    elif not math.isnan(t.real):
        t0 = t
        print(f"  {t0name:>15s} = {t0.real:<10.5f} {t0.imag:<10.5f}")
    else:
        t0 = 0.0 + 0j
    return t0


# ---------------------------------------------------------------------------
#  Geometry output
# ---------------------------------------------------------------------------


def print_geometry(StdI: StdIntList) -> None:
    """Print geometry.dat for post-processing of correlation functions.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.
    """
    if (StdI.solver == "HWAVE"
            and StdI.calcmode in ("uhfk", "rpa")):
        return

    with open("geometry.dat", "w") as fp:
        for ii in range(3):
            fp.write(f"{StdI.direct[ii, 0]:25.15e} "
                     f"{StdI.direct[ii, 1]:25.15e} "
                     f"{StdI.direct[ii, 2]:25.15e}\n")
        fp.write(f"{StdI.phase[0]:25.15e} "
                 f"{StdI.phase[1]:25.15e} "
                 f"{StdI.phase[2]:25.15e}\n")
        for ii in range(3):
            fp.write(f"{int(StdI.box[ii, 0])} "
                     f"{int(StdI.box[ii, 1])} "
                     f"{int(StdI.box[ii, 2])}\n")

        for iCell in range(StdI.NCell):
            for isite in range(StdI.NsiteUC):
                fp.write(f"{StdI.Cell[iCell, 0] - StdI.Cell[0, 0]} "
                         f"{StdI.Cell[iCell, 1] - StdI.Cell[0, 1]} "
                         f"{StdI.Cell[iCell, 2] - StdI.Cell[0, 2]} "
                         f"{isite}\n")
        if StdI.model == "kondo":
            for iCell in range(StdI.NCell):
                for isite in range(StdI.NsiteUC):
                    fp.write(f"{StdI.Cell[iCell, 0] - StdI.Cell[0, 0]} "
                             f"{StdI.Cell[iCell, 1] - StdI.Cell[0, 1]} "
                             f"{StdI.Cell[iCell, 2] - StdI.Cell[0, 2]} "
                             f"{isite + StdI.NsiteUC}\n")


# ---------------------------------------------------------------------------
#  Memory allocation for interactions
# ---------------------------------------------------------------------------


def malloc_interactions(StdI: StdIntList, ntransMax: int, nintrMax: int) -> None:
    """Allocate arrays for interactions.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    ntransMax : int
        Upper limit of the number of transfer terms.
    nintrMax : int
        Upper limit of the number of interaction terms.
    """
    # (1) Transfer
    StdI.transindx = np.zeros((ntransMax, 4), dtype=int)
    StdI.trans = np.zeros(ntransMax, dtype=complex)
    StdI.ntrans = 0

    # HPhi pump arrays
    if (StdI.solver == "HPhi"
            and StdI.method == "timeevolution"
            and StdI.PumpBody == 1):
        StdI.npump = np.zeros(StdI.Lanczos_max, dtype=int)
        StdI.pumpindx = np.zeros((StdI.Lanczos_max, ntransMax, 4), dtype=int)
        StdI.pump = np.zeros((StdI.Lanczos_max, ntransMax), dtype=complex)

    # (2) InterAll
    StdI.intrindx = np.zeros((nintrMax, 8), dtype=int)
    StdI.intr = np.zeros(nintrMax, dtype=complex)
    StdI.nintr = 0

    # (3) Coulomb intra
    StdI.CintraIndx = np.zeros((nintrMax, 1), dtype=int)
    StdI.Cintra = np.zeros(nintrMax)
    StdI.NCintra = 0

    # (4) Coulomb inter
    StdI.CinterIndx = np.zeros((nintrMax, 2), dtype=int)
    StdI.Cinter = np.zeros(nintrMax)
    StdI.NCinter = 0

    # (5) Hund
    StdI.HundIndx = np.zeros((nintrMax, 2), dtype=int)
    StdI.Hund = np.zeros(nintrMax)
    StdI.NHund = 0

    # (6) Exchange
    StdI.ExIndx = np.zeros((nintrMax, 2), dtype=int)
    StdI.Ex = np.zeros(nintrMax)
    StdI.NEx = 0

    # (7) PairLift
    StdI.PLIndx = np.zeros((nintrMax, 2), dtype=int)
    StdI.PairLift = np.zeros(nintrMax)
    StdI.NPairLift = 0

    # (8) PairHopp
    StdI.PHIndx = np.zeros((nintrMax, 2), dtype=int)
    StdI.PairHopp = np.zeros(nintrMax)
    StdI.NPairHopp = 0


# ---------------------------------------------------------------------------
#  mVMC-specific functions
# ---------------------------------------------------------------------------


def _fold_site_sub(
    StdI: StdIntList,
    iCellV: list[int],
) -> tuple[list[int], list[int]]:
    """Fold site into the sub-lattice cell (mVMC only).

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
    iCellV_frac = [0, 0, 0]
    for ii in range(3):
        for jj in range(3):
            iCellV_frac[ii] += StdI.rboxsub[ii, jj] * iCellV[jj]

    nBox = [0, 0, 0]
    for ii in range(3):
        nBox[ii] = (iCellV_frac[ii] + StdI.NCellsub * 1000) // StdI.NCellsub - 1000

    for ii in range(3):
        iCellV_frac[ii] -= StdI.NCellsub * nBox[ii]

    iCellV_fold = [0, 0, 0]
    for ii in range(3):
        for jj in range(3):
            iCellV_fold[ii] += StdI.boxsub[jj, ii] * iCellV_frac[jj]
        iCellV_fold[ii] = (iCellV_fold[ii] + StdI.NCellsub * 1000) // StdI.NCellsub - 1000

    return nBox, iCellV_fold


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
        nBox, iCellV = _fold_site_sub(
            StdI, [int(StdI.Cell[iCell, k]) for k in range(3)])
        nBox, iCellV = _fold_site(StdI, iCellV)

        if (iCellV[0] == StdI.Cell[iCell, 0]
                and iCellV[1] == StdI.Cell[iCell, 1]
                and iCellV[2] == StdI.Cell[iCell, 2]):
            for jCell in range(StdI.NCell):
                jCellV = [int(StdI.Cell[jCell, k]) + iCellV[k] for k in range(3)]
                nBox, jCellV = _fold_site(StdI, jCellV)

                for kCell in range(StdI.NCell):
                    if (jCellV[0] == StdI.Cell[kCell, 0]
                            and jCellV[1] == StdI.Cell[kCell, 1]
                            and jCellV[2] == StdI.Cell[kCell, 2]):
                        for jsite in range(StdI.NsiteUC):
                            Sym[StdI.NSym][jCell * StdI.NsiteUC + jsite] = (
                                kCell * StdI.NsiteUC + jsite)
                            Anti[StdI.NSym][jCell * StdI.NsiteUC + jsite] = (
                                StdI.AntiPeriod[0] * nBox[0]
                                + StdI.AntiPeriod[1] * nBox[1]
                                + StdI.AntiPeriod[2] * nBox[2])

                            if StdI.model == "kondo":
                                half = StdI.nsite // 2
                                Sym[StdI.NSym][half + jCell * StdI.NsiteUC + jsite] = (
                                    half + kCell * StdI.NsiteUC + jsite)
                                Anti[StdI.NSym][half + jCell * StdI.NsiteUC + jsite] = (
                                    StdI.AntiPeriod[0] * nBox[0]
                                    + StdI.AntiPeriod[1] * nBox[1]
                                    + StdI.AntiPeriod[2] * nBox[2])
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
                a = Anti[iSym][jsite]
                a = 1 if a % 2 == 0 else -1
                fp.write(f"{iSym:5d}  {jsite:5d}  {Sym[iSym][jsite]:5d}  {a:5d}\n")
    print("    qptransidx.def is written.")


def _init_site_sub(StdI: StdIntList) -> None:
    """Initialize sub-cell for mVMC/UHF/HWAVE.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    """
    NaN_i = StdI.NaN_i

    lwh_sub = (StdI.Lsub != NaN_i or StdI.Wsub != NaN_i or StdI.Hsub != NaN_i)
    box_sub = any(StdI.boxsub[i, j] != NaN_i for i in range(3) for j in range(3))

    if lwh_sub and box_sub:
        print("\nERROR ! (Lsub, Wsub, Hsub) and (a0Wsub, ..., a2Hsub) conflict !\n")
        exit_program(-1)
    elif lwh_sub:
        StdI.Lsub = print_val_i("Lsub", StdI.Lsub, 1)
        StdI.Wsub = print_val_i("Wsub", StdI.Wsub, 1)
        StdI.Hsub = print_val_i("Hsub", StdI.Hsub, 1)
        StdI.boxsub[:, :] = 0
        StdI.boxsub[0, 0] = StdI.Wsub
        StdI.boxsub[1, 1] = StdI.Lsub
        StdI.boxsub[2, 2] = StdI.Hsub
    else:
        StdI.boxsub[0, 0] = print_val_i("a0Wsub", int(StdI.boxsub[0, 0]), int(StdI.box[0, 0]))
        StdI.boxsub[0, 1] = print_val_i("a0Lsub", int(StdI.boxsub[0, 1]), int(StdI.box[0, 1]))
        StdI.boxsub[0, 2] = print_val_i("a0Hsub", int(StdI.boxsub[0, 2]), int(StdI.box[0, 2]))
        StdI.boxsub[1, 0] = print_val_i("a1Wsub", int(StdI.boxsub[1, 0]), int(StdI.box[1, 0]))
        StdI.boxsub[1, 1] = print_val_i("a1Lsub", int(StdI.boxsub[1, 1]), int(StdI.box[1, 1]))
        StdI.boxsub[1, 2] = print_val_i("a1Hsub", int(StdI.boxsub[1, 2]), int(StdI.box[1, 2]))
        StdI.boxsub[2, 0] = print_val_i("a2Wsub", int(StdI.boxsub[2, 0]), int(StdI.box[2, 0]))
        StdI.boxsub[2, 1] = print_val_i("a2Lsub", int(StdI.boxsub[2, 1]), int(StdI.box[2, 1]))
        StdI.boxsub[2, 2] = print_val_i("a2Hsub", int(StdI.boxsub[2, 2]), int(StdI.box[2, 2]))

    # Calculate reciprocal lattice vectors
    StdI.NCellsub = 0
    for ii in range(3):
        StdI.NCellsub += (int(StdI.boxsub[0, ii])
                          * int(StdI.boxsub[1, (ii + 1) % 3])
                          * int(StdI.boxsub[2, (ii + 2) % 3])
                          - int(StdI.boxsub[0, ii])
                          * int(StdI.boxsub[1, (ii + 2) % 3])
                          * int(StdI.boxsub[2, (ii + 1) % 3]))
    print(f"         Number of Cell in the sublattice: {abs(StdI.NCellsub)}")
    if StdI.NCellsub == 0:
        exit_program(-1)

    for ii in range(3):
        for jj in range(3):
            StdI.rboxsub[ii, jj] = (
                int(StdI.boxsub[(ii + 1) % 3, (jj + 1) % 3])
                * int(StdI.boxsub[(ii + 2) % 3, (jj + 2) % 3])
                - int(StdI.boxsub[(ii + 1) % 3, (jj + 2) % 3])
                * int(StdI.boxsub[(ii + 2) % 3, (jj + 1) % 3]))
    if StdI.NCellsub < 0:
        StdI.rboxsub *= -1
        StdI.NCellsub *= -1

    # Check commensurate
    for ii in range(3):
        for jj in range(3):
            prod = 0
            for kk in range(3):
                prod += int(StdI.rboxsub[ii, kk]) * int(StdI.box[jj, kk])
            if prod % StdI.NCellsub != 0:
                print("\n ERROR ! Sublattice is INCOMMENSURATE !\n")
                exit_program(-1)


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
        nBox, iCellV = _fold_site_sub(
            StdI, [int(StdI.Cell[iCell, k]) for k in range(3)])
        nBox, iCellV = _fold_site(StdI, iCellV)

        iCell2 = 0
        for kCell in range(StdI.NCell):
            if (iCellV[0] == StdI.Cell[kCell, 0]
                    and iCellV[1] == StdI.Cell[kCell, 1]
                    and iCellV[2] == StdI.Cell[kCell, 2]):
                iCell2 = kCell

        for jCell in range(StdI.NCell):
            jCellV = [int(StdI.Cell[jCell, k]) + iCellV[k]
                      - int(StdI.Cell[iCell, k]) for k in range(3)]
            nBox, jCellV = _fold_site(StdI, jCellV)

            jCell2 = 0
            for kCell in range(StdI.NCell):
                if (jCellV[0] == StdI.Cell[kCell, 0]
                        and jCellV[1] == StdI.Cell[kCell, 1]
                        and jCellV[2] == StdI.Cell[kCell, 2]):
                    jCell2 = kCell

            # AntiPeriodic factor
            dCellV = [int(StdI.Cell[jCell, k]) - int(StdI.Cell[iCell, k])
                      for k in range(3)]
            nBox_d, _ = _fold_site(StdI, dCellV)
            anti_val = sum(int(StdI.AntiPeriod[k]) * nBox_d[k] for k in range(3))
            anti_val = 1 if anti_val % 2 == 0 else -1

            for isite in range(StdI.NsiteUC):
                for jsite in range(StdI.NsiteUC):
                    if CellDone[iCell2, jCell2] == 0:
                        StdI.Orb[iCell2 * StdI.NsiteUC + isite,
                                 jCell2 * StdI.NsiteUC + jsite] = iOrb
                        StdI.AntiOrb[iCell2 * StdI.NsiteUC + isite,
                                     jCell2 * StdI.NsiteUC + jsite] = anti_val
                        iOrb += 1

                    StdI.Orb[iCell * StdI.NsiteUC + isite,
                             jCell * StdI.NsiteUC + jsite] = (
                        StdI.Orb[iCell2 * StdI.NsiteUC + isite,
                                 jCell2 * StdI.NsiteUC + jsite])
                    StdI.AntiOrb[iCell * StdI.NsiteUC + isite,
                                 jCell * StdI.NsiteUC + jsite] = anti_val

                    if StdI.model == "kondo":
                        half = StdI.nsite // 2
                        if CellDone[iCell2, jCell2] == 0:
                            StdI.Orb[half + iCell2 * StdI.NsiteUC + isite,
                                     jCell2 * StdI.NsiteUC + jsite] = iOrb
                            StdI.AntiOrb[half + iCell2 * StdI.NsiteUC + isite,
                                         jCell2 * StdI.NsiteUC + jsite] = anti_val
                            iOrb += 1
                            StdI.Orb[iCell2 * StdI.NsiteUC + isite,
                                     half + jCell2 * StdI.NsiteUC + jsite] = iOrb
                            StdI.AntiOrb[iCell2 * StdI.NsiteUC + isite,
                                         half + jCell2 * StdI.NsiteUC + jsite] = anti_val
                            iOrb += 1
                            StdI.Orb[half + iCell2 * StdI.NsiteUC + isite,
                                     half + jCell2 * StdI.NsiteUC + jsite] = iOrb
                            StdI.AntiOrb[half + iCell2 * StdI.NsiteUC + isite,
                                         half + jCell2 * StdI.NsiteUC + jsite] = anti_val
                            iOrb += 1

                        StdI.Orb[half + iCell * StdI.NsiteUC + isite,
                                 jCell * StdI.NsiteUC + jsite] = (
                            StdI.Orb[half + iCell2 * StdI.NsiteUC + isite,
                                     jCell2 * StdI.NsiteUC + jsite])
                        StdI.AntiOrb[half + iCell * StdI.NsiteUC + isite,
                                     jCell * StdI.NsiteUC + jsite] = anti_val
                        StdI.Orb[iCell * StdI.NsiteUC + isite,
                                 half + jCell * StdI.NsiteUC + jsite] = (
                            StdI.Orb[iCell2 * StdI.NsiteUC + isite,
                                     half + jCell2 * StdI.NsiteUC + jsite])
                        StdI.AntiOrb[iCell * StdI.NsiteUC + isite,
                                     half + jCell * StdI.NsiteUC + jsite] = anti_val
                        StdI.Orb[half + iCell * StdI.NsiteUC + isite,
                                 half + jCell * StdI.NsiteUC + jsite] = (
                            StdI.Orb[half + iCell2 * StdI.NsiteUC + isite,
                                     half + jCell2 * StdI.NsiteUC + jsite])
                        StdI.AntiOrb[half + iCell * StdI.NsiteUC + isite,
                                     half + jCell * StdI.NsiteUC + jsite] = anti_val

            CellDone[iCell2, jCell2] = 1

    StdI.NOrb = iOrb


def print_jastrow(StdI: StdIntList) -> None:
    """Output Jastrow factor index file ``jastrowidx.def`` (mVMC only).

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.
    """
    Jastrow = np.zeros((StdI.nsite, StdI.nsite), dtype=int)

    if abs(StdI.NMPTrans) == 1 or StdI.NMPTrans == StdI.NaN_i:
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

        NJastrow = 0 if StdI.model == "hubbard" else -1
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
    else:
        if StdI.model == "spin":
            NJastrow = 1
            Jastrow[:, :] = 0
        else:
            NJastrow = 0
            if StdI.model == "kondo":
                half = StdI.nsite // 2
                Jastrow[:, :half] = 0
                Jastrow[:half, :] = 0
                NJastrow += 1

            for dCell in range(StdI.NCell):
                isite, jsite, Cphase, dR_arr = find_site(
                    StdI, 0, 0, 0,
                    -int(StdI.Cell[dCell, 0]),
                    -int(StdI.Cell[dCell, 1]),
                    -int(StdI.Cell[dCell, 2]),
                    0, 0)
                if StdI.model == "kondo":
                    jsite -= StdI.NCell * StdI.NsiteUC
                iCell_j = jsite // StdI.NsiteUC
                if iCell_j < dCell:
                    continue
                reversal = 1 if iCell_j == dCell else 0

                for isiteUC in range(StdI.NsiteUC):
                    for jsiteUC in range(StdI.NsiteUC):
                        if reversal == 1 and jsiteUC > isiteUC:
                            continue
                        if (isiteUC == jsiteUC
                                and StdI.Cell[dCell, 0] == 0
                                and StdI.Cell[dCell, 1] == 0
                                and StdI.Cell[dCell, 2] == 0):
                            continue

                        for iCell_idx in range(StdI.NCell):
                            i_s, j_s, _, _ = find_site(
                                StdI,
                                int(StdI.Cell[iCell_idx, 0]),
                                int(StdI.Cell[iCell_idx, 1]),
                                int(StdI.Cell[iCell_idx, 2]),
                                int(StdI.Cell[dCell, 0]),
                                int(StdI.Cell[dCell, 1]),
                                int(StdI.Cell[dCell, 2]),
                                isiteUC, jsiteUC)
                            Jastrow[i_s, j_s] = NJastrow
                            Jastrow[j_s, i_s] = NJastrow

                        NJastrow += 1

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
            if StdI.model == "hubbard" or iJastrow > 0:
                fp.write(f"{iJastrow:5d}  {1:5d}\n")
            else:
                fp.write(f"{iJastrow:5d}  {0:5d}\n")
    print("    jastrowidx.def is written.")
