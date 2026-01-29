"""mVMC solver-specific output functions.

This module contains the mVMC (many-variable Variational Monte Carlo)
solver-specific output functions extracted from ``stdface_main.py``.
These functions write the orbital and Gutzwiller variational-parameter
definition files required by mVMC.

Functions
---------
print_orb
    Write the anti-parallel orbital index file ``orbitalidx.def``.
print_orb_para
    Write parallel orbital index files ``orbitalidxpara.def`` and
    ``orbitalidxgen.def``.
print_gutzwiller
    Write the Gutzwiller variational-parameter file ``gutzwilleridx.def``.

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

from ...core.stdface_vals import StdIntList, ModelType, NaN_i


def _has_anti_period(StdI: StdIntList) -> bool:
    """Return whether any lattice direction has anti-periodic boundary.

    Parameters
    ----------
    StdI : StdIntList
        Model parameter structure.

    Returns
    -------
    bool
        ``True`` if any of the three ``AntiPeriod`` flags equals 1.
    """
    return any(ap == 1 for ap in StdI.AntiPeriod)


def print_orb(StdI: StdIntList) -> None:
    """Write the anti-parallel orbital index file ``orbitalidx.def``.

    The file records the orbital pairing indices used by mVMC for the
    anti-parallel-spin part of the variational wave function.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``nsite`` -- total number of sites
        - ``NOrb`` -- number of orbital indices
        - ``ComplexType`` -- 0 for real, 1 for complex
        - ``Orb`` -- ``nsite x nsite`` orbital index matrix
        - ``AntiOrb`` -- ``nsite x nsite`` anti-periodic sign matrix
        - ``AntiPeriod`` -- length-3 array of anti-periodic boundary flags

    Notes
    -----
    Translated from the C function ``PrintOrb()`` in ``StdFace_main.c``.
    """
    with open("orbitalidx.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NOrbitalIdx {StdI.NOrb:10d}\n")
        fp.write(f"ComplexType {StdI.ComplexType:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")

        has_anti = _has_anti_period(StdI)

        for isite in range(StdI.nsite):
            for jsite in range(StdI.nsite):
                if has_anti:
                    fp.write(f"{isite:5d}  {jsite:5d}  "
                             f"{StdI.Orb[isite][jsite]:5d}  "
                             f"{StdI.AntiOrb[isite][jsite]:5d}\n")
                else:
                    fp.write(f"{isite:5d}  {jsite:5d}  "
                             f"{StdI.Orb[isite][jsite]:5d}\n")

        for iOrb in range(StdI.NOrb):
            fp.write(f"{iOrb:5d}  {1:5d}\n")

    print("    orbitalidx.def is written.")


def _compute_parallel_orbitals(
    nsite: int,
    NOrb: int,
    Orb,
    AntiOrb,
) -> tuple[list[list[int]], list[list[int]], int]:
    """Compute parallel orbital indices from anti-parallel orbital data.

    Performs three steps:

    1. **Copy**: copy the anti-parallel orbital matrix (``Orb``) into a
       local ``OrbGC`` matrix and the sign matrix (``AntiOrb``) into
       ``reverse``.
    2. **Symmetrise**: for each known orbital index, set
       ``OrbGC[j][i] = OrbGC[i][j]`` and ``reverse[j][i] = -reverse[i][j]``.
    3. **Renumber**: walk the strict lower triangle (``isite > jsite``),
       assign negative temporaries to each newly seen orbital, then invert
       all indices so they become non-negative.

    Parameters
    ----------
    nsite : int
        Total number of sites.
    NOrb : int
        Number of anti-parallel orbital indices.
    Orb : array-like
        ``(nsite, nsite)`` orbital index matrix (read-only).
    AntiOrb : array-like
        ``(nsite, nsite)`` anti-periodic sign matrix (read-only).

    Returns
    -------
    OrbGC : list of list of int
        ``(nsite, nsite)`` renumbered parallel orbital indices.
    reverse : list of list of int
        ``(nsite, nsite)`` sign-reversal matrix.
    NOrbGC : int
        Number of unique parallel orbital indices.
    """
    import numpy as np

    # (1) Copy into numpy arrays
    OrbGC = np.asarray(Orb, dtype=int).copy()
    reverse = np.asarray(AntiOrb, dtype=int).copy()

    # (2) Symmetrise — process each orbital in ascending order.
    #     For each iorb, find all (i,j) with OrbGC[i,j]==iorb and set
    #     OrbGC[j,i]=iorb, reverse[j,i]=-reverse[i,j].
    #     Order matters: later iorb can overwrite earlier iorb's writes.
    #     Process sequentially to preserve the original semantics where
    #     each (i,j) match writes to (j,i) immediately.
    for iorb in range(NOrb):
        rows, cols = np.where(OrbGC == iorb)
        for r, c in zip(rows, cols):
            OrbGC[c, r] = iorb
            reverse[c, r] = -reverse[r, c]

    # (3) Renumber — lower triangle (isite > jsite).
    #     Replace each newly-seen positive orbital value with a negative
    #     temporary across the entire matrix (vectorised).
    NOrbGC = 0
    for isite in range(nsite):
        for jsite in range(isite):
            if OrbGC[isite, jsite] >= 0:
                iOrbGC = OrbGC[isite, jsite]
                NOrbGC -= 1
                OrbGC[OrbGC == iOrbGC] = NOrbGC

    NOrbGC = -NOrbGC
    OrbGC = -1 - OrbGC

    # Convert back to list-of-lists for downstream compatibility
    OrbGC_list = OrbGC.tolist()
    reverse_list = reverse.tolist()

    return OrbGC_list, reverse_list, NOrbGC


def _write_orbitalidxpara(
    nsite: int,
    ComplexType: int,
    OrbGC: list,
    reverse: list,
    NOrbGC: int,
) -> None:
    """Write ``orbitalidxpara.def`` from pre-computed parallel orbital data.

    Parameters
    ----------
    nsite : int
        Total number of sites.
    ComplexType : int
        0 for real, 1 for complex variational parameters.
    OrbGC : array-like
        ``nsite x nsite`` parallel orbital index matrix.
    reverse : array-like
        ``nsite x nsite`` sign-reversal matrix.
    NOrbGC : int
        Number of parallel orbital indices.
    """
    with open("orbitalidxpara.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NOrbitalIdx {NOrbGC:10d}\n")
        fp.write(f"ComplexType {ComplexType:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")

        for isite in range(nsite):
            for jsite in range(nsite):
                if isite >= jsite:
                    continue
                fp.write(f"{isite:5d}  {jsite:5d}  "
                         f"{OrbGC[isite][jsite]:5d}  "
                         f"{reverse[isite][jsite]:5d}\n")

        for iOrbGC in range(NOrbGC):
            fp.write(f"{iOrbGC:5d}  {1:5d}\n")


def _write_orbitalidxgen(
    StdI: StdIntList,
    OrbGC: list,
    reverse: list,
    NOrbGC: int,
) -> None:
    """Write ``orbitalidxgen.def`` combining anti-parallel and parallel orbitals.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Reads ``nsite``, ``NOrb``,
        ``ComplexType``, ``Orb``, ``AntiOrb``, and ``AntiPeriod``.
    OrbGC : array-like
        ``nsite x nsite`` parallel orbital index matrix.
    reverse : array-like
        ``nsite x nsite`` sign-reversal matrix.
    NOrbGC : int
        Number of parallel orbital indices.
    """
    nsite = StdI.nsite
    has_anti = _has_anti_period(StdI)

    with open("orbitalidxgen.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NOrbitalIdx {StdI.NOrb + 2 * NOrbGC:10d}\n")
        fp.write(f"ComplexType {StdI.ComplexType:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")

        # -- anti-parallel section --
        for isite in range(nsite):
            for jsite in range(nsite):
                if has_anti:
                    fp.write(f"{isite:5d}  0  {jsite:5d}  1  "
                             f"{StdI.Orb[isite][jsite]:5d}  "
                             f"{StdI.AntiOrb[isite][jsite]:5d}\n")
                else:
                    fp.write(f"{isite:5d}  0  {jsite:5d}  1  "
                             f"{StdI.Orb[isite][jsite]:5d}  {1:5d}\n")

        # -- parallel section (upper triangle) --
        for isite in range(nsite):
            for jsite in range(nsite):
                if isite >= jsite:
                    continue
                fp.write(f"{isite:5d}  0  {jsite:5d}  0  "
                         f"{OrbGC[isite][jsite] + StdI.NOrb:5d}  "
                         f"{reverse[isite][jsite]:5d}\n")
                fp.write(f"{isite:5d}  1  {jsite:5d}  1  "
                         f"{OrbGC[isite][jsite] + StdI.NOrb + NOrbGC:5d}  "
                         f"{reverse[isite][jsite]:5d}\n")

        for iOrbGC in range(StdI.NOrb):
            fp.write(f"{iOrbGC:5d}  {1:5d}\n")

        for iOrbGC in range(NOrbGC * 2):
            fp.write(f"{iOrbGC + StdI.NOrb:5d}  {1:5d}\n")


def print_orb_para(StdI: StdIntList) -> None:
    """Write parallel orbital index files.

    Writes ``orbitalidxpara.def`` via :func:`_write_orbitalidxpara` and
    ``orbitalidxgen.def`` via :func:`_write_orbitalidxgen`.  Computation
    of parallel orbital indices is delegated to
    :func:`_compute_parallel_orbitals`.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``nsite`` -- total number of sites
        - ``NOrb`` -- number of anti-parallel orbital indices
        - ``ComplexType`` -- 0 for real, 1 for complex
        - ``Orb`` -- ``nsite x nsite`` orbital index matrix
        - ``AntiOrb`` -- ``nsite x nsite`` anti-periodic sign matrix
        - ``AntiPeriod`` -- length-3 array of anti-periodic boundary flags

    Notes
    -----
    Translated from the C function ``PrintOrbPara()`` in ``StdFace_main.c``.
    """
    OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(
        StdI.nsite, StdI.NOrb, StdI.Orb, StdI.AntiOrb)

    _write_orbitalidxpara(
        StdI.nsite, StdI.ComplexType, OrbGC, reverse, NOrbGC)
    print("    orbitalidxpara.def is written.")

    _write_orbitalidxgen(StdI, OrbGC, reverse, NOrbGC)
    print("    orbitalidxgen.def is written.")


def _gutzwiller_momentum_projected(
    StdI: StdIntList,
    Gutz: list[int],
) -> int:
    """Compute Gutzwiller indices in momentum-projected mode.

    For the Hubbard model, ``NGutzwiller`` starts at 0; for other models
    it starts at -1.  Diagonal orbital indices (``Orb[i][i]``) are used;
    local-spin sites are excluded (set to -1).  Unique Gutzwiller indices
    are renumbered with negative temporaries and then inverted.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.
    Gutz : list of int
        Mutable per-site Gutzwiller index array (modified in place).

    Returns
    -------
    int
        Number of unique Gutzwiller parameters (``NGutzwiller``).
    """
    nsite = StdI.nsite

    if StdI.model == ModelType.HUBBARD:
        NGutzwiller = 0
    else:
        NGutzwiller = -1

    for isite in range(nsite):
        Gutz[isite] = int(StdI.Orb[isite][isite])

    for isite in range(nsite):
        if StdI.locspinflag[isite] != 0:
            Gutz[isite] = -1
            continue
        if Gutz[isite] >= 0:
            iGutz = Gutz[isite]
            NGutzwiller -= 1
            for jsite in range(nsite):
                if Gutz[jsite] == iGutz:
                    Gutz[jsite] = NGutzwiller

    NGutzwiller = -NGutzwiller
    for isite in range(nsite):
        Gutz[isite] = -1 - Gutz[isite]

    return NGutzwiller


def _gutzwiller_global_optimization(
    StdI: StdIntList,
    Gutz: list[int],
) -> int:
    """Compute Gutzwiller indices in global-optimisation mode.

    - Hubbard: ``NGutzwiller = NsiteUC``, site index modulo ``NsiteUC``.
    - Spin: ``NGutzwiller = 1``, all sites map to index 0.
    - Kondo: ``NGutzwiller = NsiteUC + 1``, conduction sites map to 0,
      localised sites map to ``isite + 1``.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.
    Gutz : list of int
        Mutable per-site Gutzwiller index array (modified in place).

    Returns
    -------
    int
        Number of unique Gutzwiller parameters (``NGutzwiller``).
    """
    if StdI.model == ModelType.HUBBARD:
        NGutzwiller = StdI.NsiteUC
    elif StdI.model == ModelType.SPIN:
        NGutzwiller = 1
    else:
        NGutzwiller = StdI.NsiteUC + 1

    for iCell in range(StdI.NCell):
        for isite in range(StdI.NsiteUC):
            if StdI.model == ModelType.HUBBARD:
                Gutz[isite + StdI.NsiteUC * iCell] = isite
            elif StdI.model == ModelType.SPIN:
                Gutz[isite + StdI.NsiteUC * iCell] = 0
            else:
                Gutz[isite + StdI.NsiteUC * iCell] = 0
                Gutz[isite + StdI.NsiteUC * (iCell + StdI.NCell)] = isite + 1

    return NGutzwiller


def _write_gutzwiller_file(
    StdI: StdIntList,
    NGutzwiller: int,
    Gutz: list[int],
) -> None:
    """Write ``gutzwilleridx.def`` from pre-computed Gutzwiller indices.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Reads ``nsite`` and ``model``.
    NGutzwiller : int
        Number of unique Gutzwiller parameters.
    Gutz : list of int
        Per-site Gutzwiller index array.
    """
    with open("gutzwilleridx.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NGutzwillerIdx {NGutzwiller:10d}\n")
        fp.write(f"ComplexType {0:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")

        for isite in range(StdI.nsite):
            fp.write(f"{isite:5d}  {Gutz[isite]:5d}\n")

        for iGutz in range(NGutzwiller):
            flag = int(StdI.model == ModelType.HUBBARD or iGutz > 0)
            fp.write(f"{iGutz:5d}  {flag:5d}\n")


def print_gutzwiller(StdI: StdIntList) -> None:
    """Write the Gutzwiller variational-parameter file ``gutzwilleridx.def``.

    Delegates computation to :func:`_gutzwiller_momentum_projected` or
    :func:`_gutzwiller_global_optimization` based on ``NMPTrans``, then
    writes the results via :func:`_write_gutzwiller_file`.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``nsite`` -- total number of sites
        - ``NMPTrans`` -- momentum-projection control
        - ``model`` -- model type
        - ``Orb`` -- ``nsite x nsite`` orbital index matrix
        - ``locspinflag`` -- per-site local-spin flag array
        - ``NsiteUC`` -- number of sites per unit cell
        - ``NCell`` -- number of unit cells

    Notes
    -----
    Translated from the C function ``PrintGutzwiller()`` in
    ``StdFace_main.c``.
    """
    Gutz = [0] * StdI.nsite

    if abs(StdI.NMPTrans) == 1 or StdI.NMPTrans == NaN_i:
        NGutzwiller = _gutzwiller_momentum_projected(StdI, Gutz)
    else:
        NGutzwiller = _gutzwiller_global_optimization(StdI, Gutz)

    _write_gutzwiller_file(StdI, NGutzwiller, Gutz)
    print("    gutzwilleridx.def is written.")
