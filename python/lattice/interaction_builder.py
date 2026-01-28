"""Interaction term builder functions.

This module provides functions for constructing one-body (transfer) and
two-body (interaction) Hamiltonian terms, as well as the memory allocation
helper that prepares the arrays these builders populate.

Functions
---------
trans
    Add a single transfer (one-body) term.
hopping
    Add hopping for both spin channels.
hubbard_local
    Add intra-Coulomb, magnetic field, and chemical potential.
mag_field
    Add longitudinal and transverse magnetic field terms.
intr
    Add a general two-body (InterAll) interaction term.
general_j
    Add a general 3x3 spin-spin interaction.
coulomb
    Add an off-site Coulomb interaction term.
compute_max_interactions
    Compute upper limits for transfer and interaction arrays.
malloc_interactions
    Allocate arrays for transfer and interaction terms.
add_neighbor_interaction
    Label a neighbor bond and add spin or electron interaction terms (2D).
add_neighbor_interaction_3d
    Find a neighbor site and add spin or electron interaction terms (3D).
add_local_terms
    Add on-site (magnetic field, anisotropy, Hubbard U, Kondo) terms for one site.

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

import numpy as np

from typing import TextIO

from stdface_vals import StdIntList, ModelType, SolverType, MethodType, ZERO_BODY_EPS, AMPLITUDE_EPS


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
    if abs(trans0) < ZERO_BODY_EPS:
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
    if (StdI.solver == SolverType.HPhi
            and StdI.method == MethodType.TIME_EVOLUTION
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
    if abs(intr0) < ZERO_BODY_EPS:
        return
    StdI.intr[StdI.nintr] = intr0
    idx = StdI.intrindx[StdI.nintr]
    idx[0] = site1; idx[1] = spin1
    idx[2] = site2; idx[3] = spin2
    idx[4] = site3; idx[5] = spin3
    idx[6] = site4; idx[7] = spin4
    StdI.nintr += 1


def _spin_ladder_factor(S: float, Sz: float) -> float:
    """Compute the spin-ladder matrix element sqrt(S(S+1) - Sz(Sz+1)).

    This is the factor ``<S, Sz+1 | S^+ | S, Sz>`` appearing in the
    spin-raising and lowering operator matrix elements.

    Parameters
    ----------
    S : float
        Total spin quantum number.
    Sz : float
        z-component of spin.

    Returns
    -------
    float
        ``sqrt(S(S+1) - Sz(Sz+1))``
    """
    return math.sqrt(S * (S + 1.0) - Sz * (Sz + 1.0))


def _add_spin_half_terms(
    StdI: StdIntList,
    J: np.ndarray,
    isite: int,
    jsite: int,
) -> tuple[bool, bool]:
    """Add shortcut Hund/Cinter/Ex/PairLift terms for spin-1/2 sites.

    When at least one site is spin-1/2, the Ising (Jzz) part is handled
    via dedicated Hund and Cinter arrays, and the transverse part may be
    handled via Ex and PairLift arrays if the off-diagonal J components
    are negligible.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    J : numpy.ndarray
        3x3 spin interaction matrix.
    isite : int
        First site index.
    jsite : int
        Second site index.

    Returns
    -------
    tuple of (bool, bool)
        ``(use_z_general, use_ex_general)`` — flags indicating whether the
        general intr() path should still handle the Jzz and exchange terms,
        respectively.
    """
    # Hund: -0.5 * Jzz
    StdI.Hund[StdI.NHund] = -0.5 * J[2, 2]
    StdI.HundIndx[StdI.NHund][0] = isite
    StdI.HundIndx[StdI.NHund][1] = jsite
    StdI.NHund += 1

    # Cinter: -0.25 * Jzz
    StdI.Cinter[StdI.NCinter] = -0.25 * J[2, 2]
    StdI.CinterIndx[StdI.NCinter][0] = isite
    StdI.CinterIndx[StdI.NCinter][1] = jsite
    StdI.NCinter += 1

    # Check whether off-diagonal J elements allow Ex/PairLift shortcut
    cond_offdiag = (abs(J[0, 1]) < AMPLITUDE_EPS and abs(J[1, 0]) < AMPLITUDE_EPS)
    if StdI.solver == SolverType.mVMC:
        cond_offdiag = cond_offdiag and (abs(J[0, 0] - J[1, 1]) < AMPLITUDE_EPS)

    if not cond_offdiag:
        return False, True  # z handled by shortcut, ex goes through general

    # Exchange
    if StdI.solver == SolverType.mVMC or StdI.model == ModelType.KONDO:
        StdI.Ex[StdI.NEx] = -0.25 * (J[0, 0] + J[1, 1])
    else:
        StdI.Ex[StdI.NEx] = 0.25 * (J[0, 0] + J[1, 1])
    StdI.ExIndx[StdI.NEx][0] = isite
    StdI.ExIndx[StdI.NEx][1] = jsite
    StdI.NEx += 1

    # Pair lift
    StdI.PairLift[StdI.NPairLift] = 0.25 * (J[0, 0] - J[1, 1])
    StdI.PLIndx[StdI.NPairLift][0] = isite
    StdI.PLIndx[StdI.NPairLift][1] = jsite
    StdI.NPairLift += 1

    return False, False  # both handled by shortcut


def general_j(
    StdI: StdIntList,
    J: np.ndarray,
    Si2: int,
    Sj2: int,
    isite: int,
    jsite: int,
) -> None:
    """Treat J as a 3x3 matrix for general spin interactions.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    J : numpy.ndarray
        3x3 spin interaction matrix.
    Si2 : int
        Twice the spin moment (2S) at site *isite*.
    Sj2 : int
        Twice the spin moment (2S) at site *jsite*.
    isite : int
        First site index.
    jsite : int
        Second site index.
    """
    if Si2 == 1 or Sj2 == 1:
        use_z, use_ex = _add_spin_half_terms(StdI, J, isite, jsite)
    else:
        use_z, use_ex = True, True

    Si = 0.5 * Si2
    Sj = 0.5 * Sj2

    for ispin in range(Si2 + 1):
        Siz = Si - float(ispin)
        for jspin in range(Sj2 + 1):
            Sjz = Sj - float(jspin)

            # (1) J_z S_{iz} S_{jz}
            if use_z:
                intr0 = J[2, 2] * Siz * Sjz
                intr(StdI, intr0,
                     isite, ispin, isite, ispin,
                     jsite, jspin, jsite, jspin)

            if ispin > 0 and jspin > 0 and use_ex:
                fi = _spin_ladder_factor(Si, Siz)
                fj = _spin_ladder_factor(Sj, Sjz)

                # (2) S_i^+ S_j^- + h.c.
                intr0 = 0.25 * (J[0, 0] + J[1, 1] + 1j * (J[0, 1] - J[1, 0])) * fi * fj
                intr(StdI, intr0,
                     isite, ispin - 1, isite, ispin,
                     jsite, jspin, jsite, jspin - 1)
                intr(StdI, np.conj(intr0),
                     isite, ispin, isite, ispin - 1,
                     jsite, jspin - 1, jsite, jspin)

                # (3) S_i^+ S_j^+ + h.c.
                intr0 = 0.25 * (J[0, 0] - J[1, 1] - 1j * (J[0, 1] + J[1, 0])) * fi * fj
                intr(StdI, intr0,
                     isite, ispin - 1, isite, ispin,
                     jsite, jspin - 1, jsite, jspin)
                intr(StdI, np.conj(intr0),
                     isite, ispin, isite, ispin - 1,
                     jsite, jspin, jsite, jspin - 1)

            # (4) S_i^+ S_{jz} + h.c.
            if ispin > 0:
                fi = _spin_ladder_factor(Si, Siz)
                intr0 = 0.5 * (J[0, 2] - 1j * J[1, 2]) * fi * Sjz
                intr(StdI, intr0,
                     isite, ispin - 1, isite, ispin,
                     jsite, jspin, jsite, jspin)
                intr(StdI, np.conj(intr0),
                     jsite, jspin, jsite, jspin,
                     isite, ispin, isite, ispin - 1)

            # (5) S_{iz} S_j^+ + h.c.
            if jspin > 0:
                fj = _spin_ladder_factor(Sj, Sjz)
                intr0 = 0.5 * (J[2, 0] - 1j * J[2, 1]) * Siz * fj
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


def compute_max_interactions(
    StdI: StdIntList,
    n_bonds: int,
) -> tuple[int, int]:
    """Compute upper limits for the transfer and interaction arrays.

    Uses the standard formula shared by most lattice builders.  The
    lattice-specific information is captured by *n_bonds*, the total
    number of distinct neighbor bond types per unit cell (sum of
    nearest, next-nearest, and third-nearest neighbor bonds).

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (read-only; uses ``model``, ``nsite``,
        ``NCell``, ``NsiteUC``, and ``S2``).
    n_bonds : int
        Total number of neighbor bond types per unit cell
        (e.g., ``3 + 6 + 4 = 13`` for orthorhombic).

    Returns
    -------
    tuple of (int, int)
        ``(ntransMax, nintrMax)`` — upper limits for allocation.
    """
    if StdI.model == ModelType.SPIN:
        ntransMax = StdI.nsite * (StdI.S2 + 1 + 2 * StdI.S2)
        nintrMax = (StdI.NCell * (StdI.NsiteUC + n_bonds)
                    * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1))
    else:
        ntransMax = StdI.NCell * 2 * (2 * StdI.NsiteUC + 2 * n_bonds)
        nintrMax = StdI.NCell * (StdI.NsiteUC + 4 * n_bonds)
        if StdI.model == ModelType.KONDO:
            ntransMax += StdI.nsite // 2 * (StdI.S2 + 1 + 2 * StdI.S2)
            nintrMax += (StdI.nsite // 2
                         * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1))
    return ntransMax, nintrMax


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
    if (StdI.solver == SolverType.HPhi
            and StdI.method == MethodType.TIME_EVOLUTION
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


def add_neighbor_interaction(
    StdI: StdIntList,
    fp: TextIO | None,
    iW: int, iL: int,
    diW: int, diL: int,
    isiteUC: int, jsiteUC: int,
    connect: int,
    J: np.ndarray,
    t: complex,
    V: float,
) -> tuple[int, int, complex, np.ndarray]:
    """Label a neighbor bond and add the appropriate model interaction.

    Combines ``set_label`` with the model-dependent dispatch that appears
    in every lattice builder: for spin models calls ``general_j``; for
    electron models calls ``hopping`` + ``coulomb``.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    fp : TextIO or None
        Gnuplot file handle (may be ``None``).
    iW, iL : int
        Cell position of the initial site.
    diW, diL : int
        Translation to the neighbor.
    isiteUC, jsiteUC : int
        Unit-cell site indices for initial and final sites.
    connect : int
        Connection type for gnuplot (1 = nearest, 2 = 2nd, 3+ = no arrow).
    J : numpy.ndarray
        3×3 spin-coupling matrix (used when model is SPIN).
    t : complex
        Hopping amplitude *before* phase multiplication (used otherwise).
    V : float
        Coulomb repulsion (used otherwise).

    Returns
    -------
    isite : int
        Global index of the initial site.
    jsite : int
        Global index of the final site.
    Cphase : complex
        Boundary phase factor.
    dR : numpy.ndarray
        Distance vector R_i − R_j.
    """
    from .site_util import set_label  # local import to avoid circular dependency

    isite, jsite, Cphase, dR = set_label(
        StdI, fp, iW, iL, diW, diL, isiteUC, jsiteUC, connect)
    if StdI.model == ModelType.SPIN:
        general_j(StdI, J, StdI.S2, StdI.S2, isite, jsite)
    else:
        hopping(StdI, Cphase * t, isite, jsite, dR)
        coulomb(StdI, V, isite, jsite)
    return isite, jsite, Cphase, dR


def add_neighbor_interaction_3d(
    StdI: StdIntList,
    iW: int, iL: int, iH: int,
    diW: int, diL: int, diH: int,
    isiteUC: int, jsiteUC: int,
    J: np.ndarray,
    t: complex,
    V: float,
) -> tuple[int, int, complex, np.ndarray]:
    """Find a 3D neighbor site and add the appropriate model interaction.

    Combines ``find_site`` with the model-dependent dispatch: for spin
    models calls ``general_j``; for electron models calls ``hopping`` +
    ``coulomb``.  This is the 3D counterpart of
    :func:`add_neighbor_interaction`.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    iW, iL, iH : int
        Cell position of the initial site.
    diW, diL, diH : int
        Translation to the neighbor.
    isiteUC, jsiteUC : int
        Unit-cell site indices for initial and final sites.
    J : numpy.ndarray
        3×3 spin-coupling matrix (used when model is SPIN).
    t : complex
        Hopping amplitude *before* phase multiplication (used otherwise).
    V : float
        Coulomb repulsion (used otherwise).

    Returns
    -------
    isite : int
        Global index of the initial site.
    jsite : int
        Global index of the final site.
    Cphase : complex
        Boundary phase factor.
    dR : numpy.ndarray
        Distance vector R_i − R_j.
    """
    from .site_util import find_site  # local import to avoid circular dependency

    isite, jsite, Cphase, dR = find_site(
        StdI, iW, iL, iH, diW, diL, diH, isiteUC, jsiteUC)
    if StdI.model == ModelType.SPIN:
        general_j(StdI, J, StdI.S2, StdI.S2, isite, jsite)
    else:
        hopping(StdI, Cphase * t, isite, jsite, dR)
        coulomb(StdI, V, isite, jsite)
    return isite, jsite, Cphase, dR


def add_local_terms(
    StdI: StdIntList,
    isite: int,
    jsite_kondo: int,
) -> None:
    """Add on-site interaction terms for a single site.

    Handles the model-dependent local terms that every lattice builder
    applies to each site:

    * **SPIN**: magnetic field + single-ion anisotropy ``D``.
    * **HUBBARD**: Hubbard ``U``, chemical potential, and magnetic field.
    * **KONDO**: same as HUBBARD for the itinerant site, plus a Kondo
      coupling ``J`` and magnetic field on the localized spin site
      *jsite_kondo*.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    isite : int
        Global site index.  For KONDO models this is the itinerant
        (electron) site; for SPIN this is the spin site.
    jsite_kondo : int
        Global site index of the localized spin in the Kondo model.
        Ignored when the model is not KONDO.
    """
    if StdI.model == ModelType.SPIN:
        mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, isite)
        general_j(StdI, StdI.D, StdI.S2, StdI.S2, isite, isite)
    else:
        hubbard_local(StdI, StdI.mu, -StdI.h, -StdI.Gamma, -StdI.Gamma_y,
                      StdI.U, isite)
        if StdI.model == ModelType.KONDO:
            general_j(StdI, StdI.J, 1, StdI.S2, isite, jsite_kondo)
            mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y,
                      jsite_kondo)
