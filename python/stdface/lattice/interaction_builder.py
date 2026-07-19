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
from dataclasses import dataclass, field

import numpy as np

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .site_util import GnuplotBuffer

from ..core.stdface_vals import StdIntList, ModelType, SolverType, MethodType, ZERO_BODY_EPS, AMPLITUDE_EPS


@dataclass
class LocalTerms:
    """Pure-data container for the on-site (local) Hamiltonian terms.

    Produced by the pure ``*_terms`` builders and by
    :meth:`ModelPlugin.build_local_terms`.  Each field mirrors the
    corresponding ``HamiltonianTerms`` list; :meth:`extend_into` appends
    them to a :class:`StdIntList`.
    """

    trans: list = field(default_factory=list)
    intr: list = field(default_factory=list)
    Cintra: list = field(default_factory=list)
    Cinter: list = field(default_factory=list)
    Hund: list = field(default_factory=list)
    Ex: list = field(default_factory=list)
    PairLift: list = field(default_factory=list)

    def extend_into(self, StdI: StdIntList) -> None:
        """Append every term list to the matching ``StdI`` list."""
        StdI.trans_list.extend(self.trans)
        StdI.intr_list.extend(self.intr)
        StdI.Cintra_list.extend(self.Cintra)
        StdI.Cinter_list.extend(self.Cinter)
        StdI.Hund_list.extend(self.Hund)
        StdI.Ex_list.extend(self.Ex)
        StdI.PairLift_list.extend(self.PairLift)

    def merge(self, other: "LocalTerms") -> None:
        """Append every term list from *other* into this container."""
        self.trans += other.trans
        self.intr += other.intr
        self.Cintra += other.Cintra
        self.Cinter += other.Cinter
        self.Hund += other.Hund
        self.Ex += other.Ex
        self.PairLift += other.PairLift


def _trans_term(amp: complex, isite: int, ispin: int,
                jsite: int, jspin: int) -> list:
    """Return ``[(amp, i, si, j, sj)]`` or ``[]`` if below ``ZERO_BODY_EPS``."""
    if abs(amp) < ZERO_BODY_EPS:
        return []
    return [(amp, isite, ispin, jsite, jspin)]


def _intr_term(amp: complex, i1: int, s1: int, i2: int, s2: int,
               i3: int, s3: int, i4: int, s4: int) -> list:
    """Return ``[(amp, …8 indices)]`` or ``[]`` if below ``ZERO_BODY_EPS``."""
    if abs(amp) < ZERO_BODY_EPS:
        return []
    return [(amp, i1, s1, i2, s2, i3, s3, i4, s4)]


def trans(
    StdI: StdIntList,
    trans0: complex,
    isite: int,
    ispin: int,
    jsite: int,
    jspin: int,
) -> None:
    """Add a transfer (one-body) term to the list.

    Appends ``(trans0, isite, ispin, jsite, jspin)`` to ``StdI.trans_list``.

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
    StdI.trans_list.append((trans0, isite, ispin, jsite, jspin))


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
            Cphase = np.dot(StdI.At[it], dR)
            coef = np.exp(-1j * Cphase)
            for ispin in range(2):
                n = StdI.npump[it]
                StdI.pump[it][n] = coef * trans0
                StdI.pumpindx[it][n] = [isite, ispin, jsite, ispin]
                n += 1

                StdI.pump[it][n] = np.conj(coef * trans0)
                StdI.pumpindx[it][n] = [jsite, ispin, isite, ispin]
                StdI.npump[it] = n + 1
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
    hubbard_local_terms(mu0, h0, Gamma0, Gamma0_y, U0, isite).extend_into(StdI)


def hubbard_local_terms(
    mu0: float,
    h0: float,
    Gamma0: float,
    Gamma0_y: float,
    U0: float,
    isite: int,
) -> LocalTerms:
    """Pure builder for :func:`hubbard_local` (returns :class:`LocalTerms`)."""
    terms = LocalTerms()
    terms.trans += _trans_term(mu0 - 0.5 * h0, isite, 0, isite, 0)
    terms.trans += _trans_term(mu0 + 0.5 * h0, isite, 1, isite, 1)
    terms.trans += _trans_term(-0.5 * Gamma0, isite, 1, isite, 0)
    terms.trans += _trans_term(-0.5 * Gamma0, isite, 0, isite, 1)
    terms.trans += _trans_term(-0.5 * 1j * Gamma0_y, isite, 1, isite, 0)
    terms.trans += _trans_term(0.5 * 1j * Gamma0_y, isite, 0, isite, 1)
    terms.Cintra.append((U0, isite))
    return terms


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
    StdI.trans_list.extend(mag_field_terms(S2, h, Gamma, Gamma_y, isite))


def mag_field_terms(
    S2: int,
    h: float,
    Gamma: float,
    Gamma_y: float,
    isite: int,
) -> list:
    """Pure builder for :func:`mag_field` (returns a list of transfer terms)."""
    out: list = []
    S = S2 * 0.5
    for ispin in range(S2 + 1):
        Sz = S - float(ispin)
        # Longitudinal part: -h σ c†_{i σ} c_{i σ}
        out += _trans_term(-h * Sz, isite, ispin, isite, ispin)
        # Transverse part
        if ispin > 0:
            factor = math.sqrt(S * (S + 1.0) - Sz * (Sz + 1.0))
            out += _trans_term(-0.5 * Gamma * factor - 0.5 * 1j * Gamma_y * factor,
                               isite, ispin, isite, ispin - 1)
            out += _trans_term(-0.5 * Gamma * factor + 0.5 * 1j * Gamma_y * factor,
                               isite, ispin - 1, isite, ispin)
    return out


def intr(
    StdI: StdIntList,
    intr0: complex,
    site1: int, spin1: int,
    site2: int, spin2: int,
    site3: int, spin3: int,
    site4: int, spin4: int,
) -> None:
    """Add a general two-body (InterAll) interaction term to the list.

    Appends ``(intr0, i1, s1, i2, s2, i3, s3, i4, s4)`` to ``StdI.intr_list``.

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
    StdI.intr_list.append(
        (intr0, site1, spin1, site2, spin2, site3, spin3, site4, spin4))


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
    terms, use_z, use_ex = _spin_half_terms(J, isite, jsite, StdI.solver, StdI.model)
    terms.extend_into(StdI)
    return use_z, use_ex


def _spin_half_terms(
    J: np.ndarray,
    isite: int,
    jsite: int,
    solver: SolverType,
    model: ModelType,
) -> tuple[LocalTerms, bool, bool]:
    """Pure builder for :func:`_add_spin_half_terms`.

    Returns ``(terms, use_z_general, use_ex_general)``.
    """
    terms = LocalTerms()
    # Hund: -0.5 * Jzz
    terms.Hund.append((-0.5 * J[2, 2], isite, jsite))
    # Cinter: -0.25 * Jzz
    terms.Cinter.append((-0.25 * J[2, 2], isite, jsite))

    # Check whether off-diagonal J elements allow Ex/PairLift shortcut
    cond_offdiag = (abs(J[0, 1]) < AMPLITUDE_EPS and abs(J[1, 0]) < AMPLITUDE_EPS)
    if solver == SolverType.mVMC:
        cond_offdiag = cond_offdiag and (abs(J[0, 0] - J[1, 1]) < AMPLITUDE_EPS)

    if not cond_offdiag:
        return terms, False, True  # z via shortcut, ex via general

    # Exchange
    if solver == SolverType.mVMC or model == ModelType.KONDO:
        ex_val = -0.25 * (J[0, 0] + J[1, 1])
    else:
        ex_val = 0.25 * (J[0, 0] + J[1, 1])
    terms.Ex.append((ex_val, isite, jsite))

    # Pair lift
    terms.PairLift.append((0.25 * (J[0, 0] - J[1, 1]), isite, jsite))

    return terms, False, False  # both handled by shortcut


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
    general_j_terms(J, Si2, Sj2, isite, jsite,
                    StdI.solver, StdI.model).extend_into(StdI)


def general_j_terms(
    J: np.ndarray,
    Si2: int,
    Sj2: int,
    isite: int,
    jsite: int,
    solver: SolverType,
    model: ModelType,
) -> LocalTerms:
    """Pure builder for :func:`general_j` (returns :class:`LocalTerms`).

    *solver* and *model* select the spin-1/2 shortcut branch (see
    :func:`_spin_half_terms`).
    """
    if Si2 == 1 or Sj2 == 1:
        terms, use_z, use_ex = _spin_half_terms(J, isite, jsite, solver, model)
    else:
        terms, use_z, use_ex = LocalTerms(), True, True

    Si = 0.5 * Si2
    Sj = 0.5 * Sj2

    # Precompute spin ladder factors for all needed (S, Sz) pairs
    ladder_i: dict[int, float] = {}
    for ispin in range(1, Si2 + 1):
        Siz = Si - float(ispin)
        ladder_i[ispin] = _spin_ladder_factor(Si, Siz)
    ladder_j: dict[int, float] = {}
    for jspin in range(1, Sj2 + 1):
        Sjz = Sj - float(jspin)
        ladder_j[jspin] = _spin_ladder_factor(Sj, Sjz)

    for ispin in range(Si2 + 1):
        Siz = Si - float(ispin)
        for jspin in range(Sj2 + 1):
            Sjz = Sj - float(jspin)

            # (1) J_z S_{iz} S_{jz}
            if use_z:
                intr0 = J[2, 2] * Siz * Sjz
                terms.intr += _intr_term(intr0,
                                         isite, ispin, isite, ispin,
                                         jsite, jspin, jsite, jspin)

            if ispin > 0 and jspin > 0 and use_ex:
                fi = ladder_i[ispin]
                fj = ladder_j[jspin]

                # (2) S_i^+ S_j^- + h.c.
                intr0 = 0.25 * (J[0, 0] + J[1, 1] + 1j * (J[0, 1] - J[1, 0])) * fi * fj
                terms.intr += _intr_term(intr0,
                                         isite, ispin - 1, isite, ispin,
                                         jsite, jspin, jsite, jspin - 1)
                terms.intr += _intr_term(np.conj(intr0),
                                         isite, ispin, isite, ispin - 1,
                                         jsite, jspin - 1, jsite, jspin)

                # (3) S_i^+ S_j^+ + h.c.
                intr0 = 0.25 * (J[0, 0] - J[1, 1] - 1j * (J[0, 1] + J[1, 0])) * fi * fj
                terms.intr += _intr_term(intr0,
                                         isite, ispin - 1, isite, ispin,
                                         jsite, jspin - 1, jsite, jspin)
                terms.intr += _intr_term(np.conj(intr0),
                                         isite, ispin, isite, ispin - 1,
                                         jsite, jspin, jsite, jspin - 1)

            # (4) S_i^+ S_{jz} + h.c.
            if ispin > 0:
                fi = ladder_i[ispin]
                intr0 = 0.5 * (J[0, 2] - 1j * J[1, 2]) * fi * Sjz
                terms.intr += _intr_term(intr0,
                                         isite, ispin - 1, isite, ispin,
                                         jsite, jspin, jsite, jspin)
                terms.intr += _intr_term(np.conj(intr0),
                                         jsite, jspin, jsite, jspin,
                                         isite, ispin, isite, ispin - 1)

            # (5) S_{iz} S_j^+ + h.c.
            if jspin > 0:
                fj = ladder_j[jspin]
                intr0 = 0.5 * (J[2, 0] - 1j * J[2, 1]) * Siz * fj
                terms.intr += _intr_term(intr0,
                                         isite, ispin, isite, ispin,
                                         jsite, jspin - 1, jsite, jspin)
                terms.intr += _intr_term(np.conj(intr0),
                                         jsite, jspin, jsite, jspin - 1,
                                         isite, ispin, isite, ispin)
    return terms


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
    StdI.Cinter_list.append((V, isite, jsite))


def compute_max_interactions(
    StdI: StdIntList,
    n_bonds: int,
) -> int:
    """Compute the upper limit on the number of transfer terms.

    Uses the standard formula shared by most lattice builders.  The
    lattice-specific information is captured by *n_bonds*, the total
    number of distinct neighbor bond types per unit cell (sum of
    nearest, next-nearest, and third-nearest neighbor bonds).

    The result is only used to size the HPhi time-evolution pump arrays;
    the transfer and interaction terms themselves are list-based.

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
    int
        ``ntransMax`` — upper limit on the number of transfer terms.
    """
    if StdI.model == ModelType.SPIN:
        ntransMax = StdI.nsite * (StdI.S2 + 1 + 2 * StdI.S2)
    else:
        ntransMax = StdI.NCell * 2 * (2 * StdI.NsiteUC + 2 * n_bonds)
        if StdI.model == ModelType.KONDO:
            ntransMax += StdI.nsite // 2 * (StdI.S2 + 1 + 2 * StdI.S2)
    return ntransMax


def malloc_interactions(StdI: StdIntList, ntransMax: int) -> None:
    """Allocate the HPhi time-evolution pump arrays, if needed.

    Transfer and interaction terms are list-based (see ``trans_list`` /
    ``intr_list``), so no allocation is required for them.  ``ntransMax``
    is only used to size the pump arrays.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    ntransMax : int
        Upper limit on the number of transfer terms.
    """
    if (StdI.solver == SolverType.HPhi
            and StdI.method == MethodType.TIME_EVOLUTION
            and StdI.PumpBody == 1):
        StdI.npump = np.zeros(StdI.Lanczos_max, dtype=int)
        StdI.pumpindx = np.zeros((StdI.Lanczos_max, ntransMax, 4), dtype=int)
        StdI.pump = np.zeros((StdI.Lanczos_max, ntransMax), dtype=complex)

    # (3)-(8) Two-body shortcut term lists (A1: list-based)
    StdI.Cintra_list = []
    StdI.Cinter_list = []
    StdI.Hund_list = []
    StdI.Ex_list = []
    StdI.PairLift_list = []
    StdI.PairHopp_list = []


def _dispatch_bond_interaction(
    StdI: StdIntList,
    isite: int, jsite: int,
    Cphase: complex, dR: np.ndarray,
    J: np.ndarray, t: complex, V: float,
) -> None:
    """Apply the model-dependent interaction for a single bond.

    For spin models calls :func:`general_j`; for electron models calls
    :func:`hopping` and :func:`coulomb`.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in-place).
    isite, jsite : int
        Global site indices of the two bond endpoints.
    Cphase : complex
        Boundary phase factor.
    dR : numpy.ndarray
        Distance vector R_i − R_j (shape ``(3,)``).
    J : numpy.ndarray
        3×3 spin-coupling matrix (used when model is SPIN).
    t : complex
        Hopping amplitude *before* phase multiplication (used otherwise).
    V : float
        Coulomb repulsion (used otherwise).
    """
    if StdI.model == ModelType.SPIN:
        general_j(StdI, J, StdI.S2, StdI.S2, isite, jsite)
    else:
        hopping(StdI, Cphase * t, isite, jsite, dR)
        coulomb(StdI, V, isite, jsite)


def add_neighbor_interaction(
    StdI: StdIntList,
    buf: "GnuplotBuffer | None",
    cell_w: int, cell_l: int,
    delta_w: int, delta_l: int,
    uc_i: int, uc_j: int,
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
    buf : GnuplotBuffer or None
        Gnuplot bond buffer (may be ``None`` to suppress gnuplot output).
    cell_w, cell_l : int
        Cell position of the initial site.
    delta_w, delta_l : int
        Translation to the neighbor.
    uc_i, uc_j : int
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
        StdI, buf, cell_w, cell_l, delta_w, delta_l, uc_i, uc_j, connect)
    _dispatch_bond_interaction(StdI, isite, jsite, Cphase, dR, J, t, V)
    return isite, jsite, Cphase, dR


def add_neighbor_interaction_3d(
    StdI: StdIntList,
    cell_w: int, cell_l: int, iH: int,
    delta_w: int, delta_l: int, diH: int,
    uc_i: int, uc_j: int,
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
    cell_w, cell_l, iH : int
        Cell position of the initial site.
    delta_w, delta_l, diH : int
        Translation to the neighbor.
    uc_i, uc_j : int
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
        StdI, cell_w, cell_l, iH, delta_w, delta_l, diH, uc_i, uc_j)
    _dispatch_bond_interaction(StdI, isite, jsite, Cphase, dR, J, t, V)
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
    # Local import avoids the interaction_builder <-> model_plugin cycle.
    from ..core.model_plugin import get_model
    terms = get_model(StdI.model).build_local_terms(StdI, isite, jsite_kondo)
    terms.extend_into(StdI)


def expand_bonds_2d(
    StdI: StdIntList,
    buf: "GnuplotBuffer | None",
    bonds,
) -> None:
    """Expand a 2-D lattice's relative bond table over the whole super-cell.

    Generic Layer-2 expander shared by the 2-D lattice builders (L1).  For
    each cell it adds the model-dependent on-site (local) terms for every
    unit-cell sublattice, then walks the relative *bonds* table, deferring
    the absolute-index / boundary-phase work to
    :func:`add_neighbor_interaction` (i.e. :func:`find_site`).

    The traversal order — cell-outer, local-terms-then-bonds, bond-inner —
    matches the hand-written per-lattice loops it replaces, so the produced
    ``trans_list`` / interaction-list ordering is byte-identical.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in place).
    buf : GnuplotBuffer or None
        Gnuplot bond buffer (``None`` suppresses gnuplot output).
    bonds : iterable of tuple
        Relative bond table; each entry is
        ``(delta_w, delta_l, uc_i, uc_j, connect, J, t, V)`` — the same shape
        as the ``_BONDS`` tables in the lattice modules.
    """
    StdI._rel_bonds = bonds
    StdI._rel_dim = 2
    StdI._rel_local_fn = None
    kondo_off = (StdI.NsiteUC * StdI.NCell
                 if StdI.model == ModelType.KONDO else 0)
    for kCell in range(StdI.NCell):
        cell_w = StdI.Cell[kCell, 0]
        cell_l = StdI.Cell[kCell, 1]
        base = StdI.NsiteUC * kCell
        for uc in range(StdI.NsiteUC):
            add_local_terms(StdI, base + uc + kondo_off, base + uc)
        for dW, dL, si, sj, nn, J, t, V in bonds:
            add_neighbor_interaction(
                StdI, buf, cell_w, cell_l, dW, dL, si, sj, nn, J, t, V)


def expand_bonds_3d(
    StdI: StdIntList,
    bonds,
    local_fn,
) -> None:
    """Expand a 3-D lattice's relative bond table over the whole super-cell.

    Generic Layer-2 expander shared by the 3-D lattice builders (L1).  For
    each cell it runs the lattice-supplied *local_fn* (which adds the
    model-dependent on-site terms — these differ between 3-D lattices, e.g.
    pyrochlore's Kondo coupling), then walks the relative *bonds* table via
    :func:`add_neighbor_interaction_3d` (i.e. :func:`find_site`).

    The traversal order — cell-outer, local-then-bonds, bond-inner — matches
    the hand-written per-lattice loops it replaces, so the produced
    ``trans_list`` / interaction-list ordering is byte-identical.  3-D
    lattices have no gnuplot output, so there is no buffer argument.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure (modified in place).
    bonds : iterable of tuple
        Relative bond table; each entry is
        ``(delta_w, delta_l, delta_h, uc_i, uc_j, J, t, V)`` — the same shape
        as the ``_BONDS`` tables in the 3-D lattice modules.
    local_fn : callable
        ``local_fn(StdI, kCell)`` adds the on-site (local) terms for cell
        *kCell*.  Supplied by each lattice so its exact local-term handling
        is preserved.
    """
    StdI._rel_bonds = bonds
    StdI._rel_dim = 3
    StdI._rel_local_fn = local_fn
    for kCell in range(StdI.NCell):
        cell_w = StdI.Cell[kCell, 0]
        cell_l = StdI.Cell[kCell, 1]
        iH = StdI.Cell[kCell, 2]
        local_fn(StdI, kCell)
        for dW, dL, dH, si, sj, J, t, V in bonds:
            add_neighbor_interaction_3d(
                StdI, cell_w, cell_l, iH, dW, dL, dH, si, sj, J, t, V)


# ---------------------------------------------------------------------------
#  UHFk supercell normalization (L1)
# ---------------------------------------------------------------------------

_TERM_LIST_ATTRS = (
    "trans_list", "intr_list", "Cintra_list", "Cinter_list",
    "Hund_list", "Ex_list", "PairLift_list", "PairHopp_list",
)


def normalize_supercell_for_wannier(StdI: StdIntList) -> None:
    """Re-expand on a supercell large enough that no bond wraps (L1, UHFk).

    UHFk emits a Wannier90 unit-cell Hamiltonian.  On a small supercell,
    several lattice bonds fold onto the same ``(R, a, b)`` and accumulate
    (e.g. 2x2 square: ``+W`` and ``-W`` both connect the same pair -> 2t),
    making the output depend on the supercell size.  The correct,
    size-independent result needs each cell dimension larger than twice the
    longest bond range in that direction.

    This rebuilds the interaction terms on such a normalized (diagonal)
    supercell, derived from the relative bond model recorded by the expander.
    Dimensions with no bond extent (e.g. the chain's W) are kept at size 1.
    The unit cell (``NsiteUC`` / ``tau`` / ``direct``) and therefore
    ``geom.dat`` are untouched.

    No-op when no relative model was recorded (the ``wannier90`` lattice,
    which reads ``*_hr.dat``; the hand-written ladder).  A boundary phase is
    already rejected upstream (``HWavePlugin.validate``), so ``Cphase == 1``.
    """
    bonds = StdI._rel_bonds
    if bonds is None:
        return
    from .site_util import _compute_reciprocal_box, _enumerate_cells

    dim = StdI._rel_dim
    ndelta = 3 if dim == 3 else 2
    # Index of (delta_w, delta_l[, delta_h]) within each bond tuple.
    maxd = [0, 0, 0]
    for b in bonds:
        for i in range(ndelta):
            maxd[i] = max(maxd[i], abs(int(b[i])))

    box = np.zeros((3, 3), dtype=int)
    for i in range(3):
        if i < ndelta and maxd[i] > 0:
            box[i, i] = 2 * maxd[i] + 1   # no two deltas coincide mod size
        else:
            box[i, i] = max(1, int(StdI.box[i, i]))
    StdI.box = box
    _compute_reciprocal_box(StdI)
    _enumerate_cells(StdI)

    for attr in _TERM_LIST_ATTRS:
        setattr(StdI, attr, [])
    if dim == 2:
        expand_bonds_2d(StdI, None, bonds)
    else:
        expand_bonds_3d(StdI, bonds, StdI._rel_local_fn)
