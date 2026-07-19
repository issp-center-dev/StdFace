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

import logging

from ...core.stdface_vals import StdIntList, ModelType
from ...core.param_check import print_val_d, print_val_i, not_used_i


logger = logging.getLogger(__name__)


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
        lines = [
            "=============================================\n",
            f"NOrbitalIdx {StdI.NOrb:10d}\n",
            f"ComplexType {StdI.ComplexType:10d}\n",
            "=============================================\n",
            "=============================================\n",
        ]

        has_anti = _has_anti_period(StdI)

        for isite in range(StdI.nsite):
            for jsite in range(StdI.nsite):
                if has_anti:
                    lines.append(f"{isite:5d}  {jsite:5d}  "
                                 f"{StdI.Orb[isite][jsite]:5d}  "
                                 f"{StdI.AntiOrb[isite][jsite]:5d}\n")
                else:
                    lines.append(f"{isite:5d}  {jsite:5d}  "
                                 f"{StdI.Orb[isite][jsite]:5d}\n")

        for iOrb in range(StdI.NOrb):
            lines.append(f"{iOrb:5d}  {1:5d}\n")
        fp.write("".join(lines))

    logger.info("    orbitalidx.def is written.")


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
    #     Build a reverse map (value -> list of positions) in a single
    #     pass, then iterate orbitals in order.  This replaces NOrb
    #     full-matrix scans with one scan + targeted updates.
    from collections import defaultdict
    orb_positions: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for r in range(nsite):
        for c in range(nsite):
            v = int(OrbGC[r, c])
            if 0 <= v < NOrb:
                orb_positions[v].append((r, c))

    for iorb in range(NOrb):
        for r, c in orb_positions.get(iorb, ()):
            if OrbGC[r, c] == iorb:
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
        lines = [
            "=============================================\n",
            f"NOrbitalIdx {NOrbGC:10d}\n",
            f"ComplexType {ComplexType:10d}\n",
            "=============================================\n",
            "=============================================\n",
        ]

        for isite in range(nsite):
            for jsite in range(isite + 1, nsite):
                lines.append(f"{isite:5d}  {jsite:5d}  "
                             f"{OrbGC[isite][jsite]:5d}  "
                             f"{reverse[isite][jsite]:5d}\n")

        for iOrbGC in range(NOrbGC):
            lines.append(f"{iOrbGC:5d}  {1:5d}\n")
        fp.write("".join(lines))


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
        lines = [
            "=============================================\n",
            f"NOrbitalIdx {StdI.NOrb + 2 * NOrbGC:10d}\n",
            f"ComplexType {StdI.ComplexType:10d}\n",
            "=============================================\n",
            "=============================================\n",
        ]

        # -- anti-parallel section --
        for isite in range(nsite):
            for jsite in range(nsite):
                if has_anti:
                    lines.append(f"{isite:5d}  0  {jsite:5d}  1  "
                                 f"{StdI.Orb[isite][jsite]:5d}  "
                                 f"{StdI.AntiOrb[isite][jsite]:5d}\n")
                else:
                    lines.append(f"{isite:5d}  0  {jsite:5d}  1  "
                                 f"{StdI.Orb[isite][jsite]:5d}  {1:5d}\n")

        # -- parallel section (upper triangle) --
        for isite in range(nsite):
            for jsite in range(isite + 1, nsite):
                lines.append(f"{isite:5d}  0  {jsite:5d}  0  "
                             f"{OrbGC[isite][jsite] + StdI.NOrb:5d}  "
                             f"{reverse[isite][jsite]:5d}\n")
                lines.append(f"{isite:5d}  1  {jsite:5d}  1  "
                             f"{OrbGC[isite][jsite] + StdI.NOrb + NOrbGC:5d}  "
                             f"{reverse[isite][jsite]:5d}\n")

        for iOrbGC in range(StdI.NOrb):
            lines.append(f"{iOrbGC:5d}  {1:5d}\n")

        for iOrbGC in range(NOrbGC * 2):
            lines.append(f"{iOrbGC + StdI.NOrb:5d}  {1:5d}\n")
        fp.write("".join(lines))


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
    logger.info("    orbitalidxpara.def is written.")

    _write_orbitalidxgen(StdI, OrbGC, reverse, NOrbGC)
    logger.info("    orbitalidxgen.def is written.")


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
        lines = [
            "=============================================\n",
            f"NGutzwillerIdx {NGutzwiller:10d}\n",
            f"ComplexType {0:10d}\n",
            "=============================================\n",
            "=============================================\n",
        ]

        for isite in range(StdI.nsite):
            lines.append(f"{isite:5d}  {Gutz[isite]:5d}\n")

        for iGutz in range(NGutzwiller):
            flag = int(StdI.model == ModelType.HUBBARD or iGutz > 0)
            lines.append(f"{iGutz:5d}  {flag:5d}\n")
        fp.write("".join(lines))


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

    if StdI.NMPTrans is None or abs(StdI.NMPTrans) == 1:
        NGutzwiller = _gutzwiller_momentum_projected(StdI, Gutz)
    else:
        NGutzwiller = _gutzwiller_global_optimization(StdI, Gutz)

    _write_gutzwiller_file(StdI, NGutzwiller, Gutz)
    logger.info("    gutzwilleridx.def is written.")


# -----------------------------------------------------------------------
#  modpara / namelist bodies (moved from writer/common_writer.py)
# -----------------------------------------------------------------------

MODEL_GC_TO_EX_UPDATE_PATH: dict[tuple, int] = {
    (ModelType.HUBBARD, 0): 0,
    (ModelType.HUBBARD, 1): 0,
    (ModelType.SPIN, 0):    2,
    (ModelType.SPIN, 1):    2,
    (ModelType.KONDO, 0):   1,
    (ModelType.KONDO, 1):   3,
}
"""Maps ``(ModelType, lGC)`` to the ``NExUpdatePath`` integer for mVMC.

Hubbard always gets 0, Spin always gets 2, and Kondo depends on whether
the grand-canonical flag (``lGC``) is set (1 → 3, 0 → 1).
"""


def namelist_entries(StdI: StdIntList) -> list:
    """Return mVMC-specific ``(keyword, filename)`` namelist entries."""
    entries: list = [
        ("Gutzwiller", "gutzwilleridx.def"),
        ("Jastrow", "jastrowidx.def"),
        ("Orbital", "orbitalidx.def"),
    ]
    if StdI.lGC == 1 or (StdI.Sz2 != 0 and StdI.Sz2 is not None):
        entries.append(("OrbitalParallel", "orbitalidxpara.def"))
        entries.append(("# OrbitalGeneral", "orbitalidxgen.def"))
    entries.append(("TransSym", "qptransidx.def"))
    return entries


def modpara_lines(StdI: StdIntList) -> list:
    """Return the mVMC-specific body line descriptors of ``modpara.def``."""
    lines: list = [
        ("raw", "VMC_Cal_Parameters"),
        ("sep",),
        ("kv", "CDataFileHead", StdI.CDataFileHead, ""),
        ("kv", "CParaFileHead", StdI.CParaFileHead, ""),
        ("sep",),
        ("kv", "NVMCCalMode", StdI.NVMCCalMode, ""),
        ("kv", "NLanczosMode", StdI.NLanczosMode, ""),
        ("sep",),
        ("kv", "NDataIdxStart", StdI.NDataIdxStart, ""),
        ("kv", "NDataQtySmp", StdI.NDataQtySmp, ""),
        ("sep",),
        ("kv", "Nsite", StdI.nsite, ""),
        ("kv", "Ncond", StdI.ncond, "<5d"),
    ]
    if StdI.Sz2 is not None:
        lines.append(("kv", "2Sz", StdI.Sz2, ""))
    if StdI.NSPGaussLeg is not None:
        lines.append(("kv", "NSPGaussLeg", StdI.NSPGaussLeg, ""))
    if StdI.NSPStot is not None:
        lines.append(("kv", "NSPStot", StdI.NSPStot, ""))
    lines += [
        ("kv", "NMPTrans", StdI.NMPTrans, ""),
        ("kv", "NSROptItrStep", StdI.NSROptItrStep, ""),
        ("kv", "NSROptItrSmp", StdI.NSROptItrSmp, ""),
        ("kv", "DSROptRedCut", StdI.DSROptRedCut, ".10f"),
        ("kv", "DSROptStaDel", StdI.DSROptStaDel, ".10f"),
        ("kv", "DSROptStepDt", StdI.DSROptStepDt, ".10f"),
        ("kv", "NVMCWarmUp", StdI.NVMCWarmUp, ""),
        ("kv", "NVMCInterval", StdI.NVMCInterval, ""),
        ("kv", "NVMCSample", StdI.NVMCSample, ""),
        ("kv", "NExUpdatePath", StdI.NExUpdatePath, ""),
        ("kv", "RndSeed", StdI.RndSeed, ""),
        ("kv", "NSplitSize", StdI.NSplitSize, ""),
        ("kv", "NStore", StdI.NStore, ""),
        ("kv", "NSRCG", StdI.NSRCG, ""),
    ]
    return lines


def set_modpara_defaults(StdI: StdIntList) -> None:
    """Set mVMC-specific default model parameters.

    Handles VMC sampling parameters, exchange-update path count,
    and SR-optimisation settings.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure, modified in place.
    """
    if StdI.CParaFileHead is None:
        StdI.CParaFileHead = "zqp"
        logger.info(
            "    CParaFileHead = %-12s######  DEFAULT VALUE IS USED  ######",
            StdI.CParaFileHead,
        )
    else:
        logger.info("    CParaFileHead = %s", StdI.CParaFileHead)

    StdI.NVMCCalMode = print_val_i("NVMCCalMode", StdI.NVMCCalMode, 0)
    StdI.NLanczosMode = print_val_i("NLanczosMode", StdI.NLanczosMode, 0)
    StdI.NDataIdxStart = print_val_i("NDataIdxStart", StdI.NDataIdxStart, 1)

    if StdI.NVMCCalMode == 0:
        not_used_i("NDataQtySmp", StdI.NDataQtySmp)
    StdI.NDataQtySmp = print_val_i("NDataQtySmp", StdI.NDataQtySmp, 1)

    if StdI.lGC == 0 and (StdI.Sz2 == 0 or StdI.Sz2 is None):
        StdI.NSPGaussLeg = print_val_i("NSPGaussLeg", StdI.NSPGaussLeg, 8)
        StdI.NSPStot = print_val_i("NSPStot", StdI.NSPStot, 0)
    else:
        not_used_i("NSPGaussLeg", StdI.NSPGaussLeg)
        not_used_i("NSPStot", StdI.NSPStot)

    StdI.NMPTrans = print_val_i("NMPTrans", StdI.NMPTrans, -1)

    StdI.NSROptItrStep = print_val_i("NSROptItrStep", StdI.NSROptItrStep, 1000)

    if StdI.NVMCCalMode == 1:
        not_used_i("NSROptItrSmp", StdI.NSROptItrSmp)
    StdI.NSROptItrSmp = print_val_i(
        "NSROptItrSmp", StdI.NSROptItrSmp, StdI.NSROptItrStep // 10
    )

    StdI.NVMCWarmUp = print_val_i("NVMCWarmUp", StdI.NVMCWarmUp, 10)
    StdI.NVMCInterval = print_val_i("NVMCInterval", StdI.NVMCInterval, 1)
    StdI.NVMCSample = print_val_i("NVMCSample", StdI.NVMCSample, 1000)

    key = (StdI.model, StdI.lGC)
    ex_path = MODEL_GC_TO_EX_UPDATE_PATH.get(key)
    if ex_path is not None:
        StdI.NExUpdatePath = ex_path
    logger.info("  %15s = %-10d", "NExUpdatePath", StdI.NExUpdatePath)

    StdI.RndSeed = print_val_i("RndSeed", StdI.RndSeed, 123456789)
    StdI.NSplitSize = print_val_i("NSplitSize", StdI.NSplitSize, 1)
    StdI.NStore = print_val_i("NStore", StdI.NStore, 1)
    StdI.NSRCG = print_val_i("NSRCG", StdI.NSRCG, 0)

    StdI.DSROptRedCut = print_val_d("DSROptRedCut", StdI.DSROptRedCut, 0.001)
    StdI.DSROptStaDel = print_val_d("DSROptStaDel", StdI.DSROptStaDel, 0.02)
    StdI.DSROptStepDt = print_val_d("DSROptStepDt", StdI.DSROptStepDt, 0.02)
