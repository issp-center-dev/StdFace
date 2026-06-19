"""Common output writer functions shared across all solvers.

This module contains functions that generate definition files common to all
supported solvers (HPhi, mVMC, UHF, H-wave).  These were extracted from
``stdface_main.py`` during refactoring.

Classes
-------
GreenFunctionIndices
    Generator for one- and two-body Green-function index lists.

Functions
---------
print_loc_spin
    Write ``locspn.def`` listing local-spin flags for every site.
print_trans
    Write ``trans.def`` listing one-body transfer integrals.
print_namelist
    Write ``namelist.def`` that lists all definition files for the solver.
print_mod_para
    Write ``modpara.def`` containing model / calculation parameters.
print_1_green
    Write ``greenone.def`` listing one-body Green-function indices.
print_2_green
    Write ``greentwo.def`` listing two-body Green-function indices.
unsupported_system
    Print an error message and abort for an unsupported model/lattice pair.
check_output_mode
    Verify and set the integer output-mode flag from the string keyword.
check_mod_para
    Validate and set default values for solver-specific model parameters.

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
from collections.abc import Callable
from itertools import product
from typing import NamedTuple

import numpy as np

from ..core.stdface_vals import (
    StdIntList, ModelType, SolverType, MethodType,
    NaN_i, AMPLITUDE_EPS,
)
from ..core.param_check import print_val_i, print_val_d, required_val_i, not_used_i

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
#  String → integer dispatch tables
# ---------------------------------------------------------------------------

OUTPUT_MODE_TO_INT: dict[str, int] = {
    "non":         0,
    "none":        0,
    "off":         0,
    "cor":         1,
    "corr":        1,
    "correlation": 1,
    "raw":         2,
    "all":         2,
    "full":        2,
}
"""Maps ``outputmode`` keyword strings to the ``ioutputmode`` integer code.

- 0 = no correlation output
- 1 = correlation functions only
- 2 = raw (all) output
"""

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


def _merge_duplicate_terms(indx, vals, n: int) -> int:
    """Merge duplicate index quadruples and count non-negligible entries.

    For each pair of terms with identical 4-element index quadruples,
    the first term's amplitude absorbs the second's, and the second is
    zeroed.  After merging, the number of entries whose absolute value
    exceeds ``AMPLITUDE_EPS`` is returned.

    Parameters
    ----------
    indx : array-like
        2-D index array, shape ``(n, 4)``.  Each row is a 4-element
        index quadruple ``(i, s_i, j, s_j)``.
    vals : array-like
        1-D amplitude array, length ``n``.  Modified **in place** during
        merging (duplicates are zeroed).
    n : int
        Number of entries to process.

    Returns
    -------
    int
        Count of entries with ``abs(val) > AMPLITUDE_EPS`` after merging.
    """
    # Merge duplicates: first occurrence absorbs all later ones with same key
    seen: dict[tuple, int] = {}
    for k in range(n):
        key = tuple(indx[k])
        if key in seen:
            vals[seen[key]] += vals[k]
            vals[k] = 0.0
        else:
            seen[key] = k

    return sum(1 for k in range(n) if abs(vals[k]) > AMPLITUDE_EPS)


def print_loc_spin(StdI: StdIntList) -> None:
    """Write ``locspn.def`` listing local-spin flags for every site.

    This is the Python translation of the C function ``PrintLocSpin()``.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``nsite`` : int -- total number of sites.
        - ``locspinflag`` : list of int -- per-site flag (0 = itinerant
          electron, nonzero = local spin with :math:`S` given by the value).
    """
    nlocspin = int(np.count_nonzero(StdI.locspinflag[:StdI.nsite]))

    lines = ["================================ \n",
             f"NlocalSpin {nlocspin:5d}  \n",
             "================================ \n",
             "========i_1LocSpn_0IteElc ====== \n",
             "================================ \n"]
    for isite in range(StdI.nsite):
        lines.append(f"{isite:5d}  {StdI.locspinflag[isite]:5d}\n")
    with open("locspn.def", "w") as fp:
        fp.write("".join(lines))

    logger.info("    locspn.def is written.")


def print_trans(StdI: StdIntList) -> None:
    """Write ``trans.def`` listing one-body transfer integrals.

    This is the Python translation of the C function ``PrintTrans()``.
    Duplicate index quadruples are merged (their amplitudes summed) and
    entries whose absolute value is below ``AMPLITUDE_EPS`` are suppressed.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read
        and (for merging) modified **in place**:

        - ``trans_list`` : list of ``(amp, i, s_i, j, s_j)`` tuples.
    """
    ntrans = len(StdI.trans_list)
    vals = np.array([t[0] for t in StdI.trans_list], dtype=complex)
    indx = np.array([t[1:5] for t in StdI.trans_list], dtype=int).reshape(ntrans, 4)
    ntrans0 = _merge_duplicate_terms(indx, vals, ntrans)

    # --- write file ---
    lines = ["======================== \n",
             f"NTransfer {ntrans0:7d}  \n",
             "======================== \n",
             "========i_j_s_tijs====== \n",
             "======================== \n"]
    for ktrans in range(ntrans):
        val = vals[ktrans]
        if abs(val) > AMPLITUDE_EPS:
            i0, s0, i1, s1 = indx[ktrans]
            lines.append(
                f"{i0:5d} {s0:5d} {i1:5d} {s1:5d} "
                f"{val.real:25.15f} {val.imag:25.15f}\n"
            )
    with open("trans.def", "w") as fp:
        fp.write("".join(lines))

    logger.info("      trans.def is written.")


def _write_namelist_hphi(fp, StdI: StdIntList) -> None:
    """Write HPhi-specific entries in ``namelist.def``.

    Parameters
    ----------
    fp : file object
        Open file handle for ``namelist.def``.
    StdI : StdIntList
        The global parameter structure.
    """
    fp.write("         CalcMod  calcmod.def\n")
    if StdI.SpectrumBody == 1:
        fp.write("SingleExcitation  single.def\n")
    else:
        fp.write("  PairExcitation  pair.def\n")
    if StdI.method == MethodType.TIME_EVOLUTION:
        if StdI.PumpBody == 1:
            fp.write("       TEOneBody  teone.def\n")
        elif StdI.PumpBody == 2:
            fp.write("       TETwoBody  tetwo.def\n")
    fp.write(f"     SpectrumVec  {StdI.CDataFileHead}_eigenvec_0\n")
    if StdI.lBoost == 1:
        fp.write("           Boost  boost.def\n")


def _write_namelist_mvmc(fp, StdI: StdIntList) -> None:
    """Write mVMC-specific entries in ``namelist.def``.

    Parameters
    ----------
    fp : file object
        Open file handle for ``namelist.def``.
    StdI : StdIntList
        The global parameter structure.
    """
    fp.write("      Gutzwiller  gutzwilleridx.def\n")
    fp.write("         Jastrow  jastrowidx.def\n")
    fp.write("         Orbital  orbitalidx.def\n")
    if StdI.lGC == 1 or (StdI.Sz2 != 0 and StdI.Sz2 is not None):
        fp.write(" OrbitalParallel  orbitalidxpara.def\n")
        fp.write("# OrbitalGeneral  orbitalidxgen.def\n")
    fp.write("        TransSym  qptransidx.def\n")


# Interaction file flags and their namelist lines
_INTERACTION_FLAGS: list[tuple[str, str]] = [
    ("LCintra",   "    CoulombIntra  coulombintra.def\n"),
    ("LCinter",   "    CoulombInter  coulombinter.def\n"),
    ("LHund",     "            Hund  hund.def\n"),
    ("LEx",       "        Exchange  exchange.def\n"),
    ("LPairLift", "        PairLift  pairlift.def\n"),
    ("LPairHopp", "         PairHop  pairhopp.def\n"),
    ("Lintr",     "        InterAll  interall.def\n"),
]
"""Maps ``StdIntList`` flag attribute names to their namelist entry lines.

Each flag, when equal to 1, causes the corresponding definition file to be
listed in ``namelist.def``.
"""


def print_namelist(StdI: StdIntList) -> None:
    """Write ``namelist.def`` that lists all definition files for the solver.

    This is the Python translation of the C function ``PrintNamelist()``.
    The content depends on which solver is active (``StdI.solver``).

    The common prefix (ModPara, LocSpin, Trans, conditional interaction
    files, Green-function files) is written for all solvers.  Solver-specific
    entries are delegated to
    :meth:`ExpertModeSolverPlugin.write_namelist_body`.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.
    """
    with open("namelist.def", "w") as fp:
        fp.write("         ModPara  modpara.def\n")
        fp.write("         LocSpin  locspn.def\n")
        fp.write("           Trans  trans.def\n")

        for flag_attr, line in _INTERACTION_FLAGS:
            if getattr(StdI, flag_attr) == 1:
                fp.write(line)

        from ..plugin import get_plugin, ExpertModeSolverPlugin
        plugin = get_plugin(StdI.solver)
        assert isinstance(plugin, ExpertModeSolverPlugin)

        if StdI.ioutputmode != 0:
            fp.write("        OneBodyG  greenone.def\n")
            if plugin.has_two_body_green(StdI):
                fp.write("        TwoBodyG  greentwo.def\n")

        plugin.write_namelist_body(fp, StdI)

    logger.info("    namelist.def is written.")


def _write_modpara_hphi(fp, StdI: StdIntList) -> None:
    """Write the HPhi-specific body of ``modpara.def``.

    Parameters
    ----------
    fp : file object
        Open file handle for ``modpara.def``.
    StdI : StdIntList
        The global parameter structure.
    """
    fp.write("HPhi_Cal_Parameters\n")
    fp.write("--------------------\n")
    fp.write(f"CDataFileHead  {StdI.CDataFileHead}\n")
    fp.write("CParaFileHead  zqp\n")
    fp.write("--------------------\n")
    fp.write(f"Nsite          {StdI.nsite:<5d}\n")
    if StdI.Sz2 is not None:
        fp.write(f"2Sz            {StdI.Sz2:<5d}\n")
    if StdI.ncond is not None:
        fp.write(f"Ncond          {StdI.ncond:<5d}\n")
    fp.write(f"Lanczos_max    {StdI.Lanczos_max:<5d}\n")
    fp.write(f"initial_iv     {StdI.initial_iv:<5d}\n")
    if StdI.nvec is not None:
        fp.write(f"nvec           {StdI.nvec:<5d}\n")
    fp.write(f"exct           {StdI.exct:<5d}\n")
    fp.write(f"LanczosEps     {StdI.LanczosEps:<5d}\n")
    fp.write(f"LanczosTarget  {StdI.LanczosTarget:<5d}\n")
    fp.write(f"LargeValue     {StdI.LargeValue:<25.15e}\n")
    fp.write(f"NumAve         {StdI.NumAve:<5d}\n")
    fp.write(f"ExpecInterval  {StdI.ExpecInterval:<5d}\n")
    fp.write(f"NOmega         {StdI.Nomega:<5d}\n")
    fp.write(
        f"OmegaMax       {StdI.OmegaMax:<25.15e} {StdI.OmegaIm:<25.15e}\n"
    )
    fp.write(
        f"OmegaMin       {StdI.OmegaMin:<25.15e} {StdI.OmegaIm:<25.15e}\n"
    )
    fp.write(
        f"OmegaOrg       {StdI.OmegaOrg:<25.15e} {0.0:<25.15e}\n"
    )
    fp.write(f"PreCG          {0:<5d}\n")
    if StdI.method == MethodType.TIME_EVOLUTION:
        fp.write(f"ExpandCoef     {StdI.ExpandCoef:<5d}\n")


def _write_modpara_mvmc(fp, StdI: StdIntList) -> None:
    """Write the mVMC-specific body of ``modpara.def``.

    Parameters
    ----------
    fp : file object
        Open file handle for ``modpara.def``.
    StdI : StdIntList
        The global parameter structure.
    """
    fp.write("VMC_Cal_Parameters\n")
    fp.write("--------------------\n")
    fp.write(f"CDataFileHead  {StdI.CDataFileHead}\n")
    fp.write(f"CParaFileHead  {StdI.CParaFileHead}\n")
    fp.write("--------------------\n")
    fp.write(f"NVMCCalMode    {StdI.NVMCCalMode}\n")
    fp.write(f"NLanczosMode   {StdI.NLanczosMode}\n")
    fp.write("--------------------\n")
    fp.write(f"NDataIdxStart  {StdI.NDataIdxStart}\n")
    fp.write(f"NDataQtySmp    {StdI.NDataQtySmp}\n")
    fp.write("--------------------\n")
    fp.write(f"Nsite          {StdI.nsite}\n")
    fp.write(f"Ncond          {StdI.ncond:<5d}\n")
    if StdI.Sz2 is not None:
        fp.write(f"2Sz            {StdI.Sz2}\n")
    if StdI.NSPGaussLeg is not None:
        fp.write(f"NSPGaussLeg    {StdI.NSPGaussLeg}\n")
    if StdI.NSPStot is not None:
        fp.write(f"NSPStot        {StdI.NSPStot}\n")
    fp.write(f"NMPTrans       {StdI.NMPTrans}\n")
    fp.write(f"NSROptItrStep  {StdI.NSROptItrStep}\n")
    fp.write(f"NSROptItrSmp   {StdI.NSROptItrSmp}\n")
    fp.write(f"DSROptRedCut   {StdI.DSROptRedCut:.10f}\n")
    fp.write(f"DSROptStaDel   {StdI.DSROptStaDel:.10f}\n")
    fp.write(f"DSROptStepDt   {StdI.DSROptStepDt:.10f}\n")
    fp.write(f"NVMCWarmUp     {StdI.NVMCWarmUp}\n")
    fp.write(f"NVMCInterval   {StdI.NVMCInterval}\n")
    fp.write(f"NVMCSample     {StdI.NVMCSample}\n")
    fp.write(f"NExUpdatePath  {StdI.NExUpdatePath}\n")
    fp.write(f"RndSeed        {StdI.RndSeed}\n")
    fp.write(f"NSplitSize     {StdI.NSplitSize}\n")
    fp.write(f"NStore         {StdI.NStore}\n")
    fp.write(f"NSRCG          {StdI.NSRCG}\n")


def _write_modpara_uhf_hwave(fp, StdI: StdIntList) -> None:
    """Write the UHF / H-wave body of ``modpara.def``.

    UHF and H-wave share identical file content except for the header
    banner line.  The correct banner is selected from
    :data:`_MODPARA_BANNER`.

    Parameters
    ----------
    fp : file object
        Open file handle for ``modpara.def``.
    StdI : StdIntList
        The global parameter structure.
    """
    banner = _MODPARA_BANNER[StdI.solver]
    fp.write(f"{banner}\n")
    fp.write("--------------------\n")
    fp.write(f"CDataFileHead  {StdI.CDataFileHead}\n")
    fp.write("CParaFileHead  zqp\n")
    fp.write("--------------------\n")
    fp.write(f"Nsite          {StdI.nsite}\n")
    if StdI.Sz2 is not None:
        fp.write(f"2Sz            {StdI.Sz2:<5d}\n")
    # UHF/HWAVE emit an unset Ncond as the integer sentinel (legacy format).
    ncond_out = StdI.ncond if StdI.ncond is not None else NaN_i
    fp.write(f"Ncond          {ncond_out:<5d}\n")
    fp.write(f"IterationMax   {StdI.Iteration_max}\n")
    fp.write(f"EPS            {StdI.eps}\n")
    fp.write(f"Mix            {StdI.mix:.10f}\n")
    fp.write(f"RndSeed        {StdI.RndSeed}\n")
    fp.write(f"EpsSlater      {StdI.eps_slater}\n")
    fp.write(f"NMPTrans       {StdI.NMPTrans}\n")


_MODPARA_BANNER: dict[str, str] = {
    SolverType.UHF:   "UHF_Cal_Parameters",
    SolverType.HWAVE: "HWAVE_Cal_Parameters",
}
"""Maps UHF/HWAVE solver type to the ``modpara.def`` banner line."""


def print_mod_para(StdI: StdIntList) -> None:
    """Write ``modpara.def`` containing model / calculation parameters.

    This is the Python translation of the C function ``PrintModPara()``.
    The common header is written here; the solver-specific body is provided
    by :meth:`ExpertModeSolverPlugin.write_modpara_body`.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  A large number of solver-specific
        fields are read; see the individual body-writer functions for details.
    """
    from ..plugin import get_plugin, ExpertModeSolverPlugin
    plugin = get_plugin(StdI.solver)
    assert isinstance(plugin, ExpertModeSolverPlugin)

    with open("modpara.def", "w") as fp:
        fp.write("--------------------\n")
        fp.write("Model_Parameters   0\n")
        fp.write("--------------------\n")
        plugin.write_modpara_body(fp, StdI)

    logger.info("     modpara.def is written.")


class GreenFunctionIndices:
    """Generator for one- and two-body Green-function index lists.

    Encapsulates the site/spin parameters and helper logic needed to
    produce the index tuples written to ``greenone.def`` and
    ``greentwo.def``.

    Parameters
    ----------
    nsite : int
        Total number of sites in the system.
    NsiteUC : int
        Number of sites in the unit cell.
    locspinflag : list of int
        Per-site local-spin flag array (length ``nsite``).
    is_kondo : bool
        Whether the model is Kondo (doubles the UC range).
    is_mvmc : bool
        Whether the solver is mVMC (uses alternative index order for
        two-body off-diagonal spin entries).

    Attributes
    ----------
    nsite : int
    NsiteUC : int
    locspinflag : list of int
    is_kondo : bool
    is_mvmc : bool
    """

    def __init__(
        self,
        nsite: int,
        NsiteUC: int,
        locspinflag: list[int],
        is_kondo: bool,
        is_mvmc: bool = False,
    ) -> None:
        self.nsite = nsite
        self.NsiteUC = NsiteUC
        self.locspinflag = locspinflag
        self.is_kondo = is_kondo
        self.is_mvmc = is_mvmc

        # Pre-compute lookup tables to avoid per-call overhead in tight loops
        self._spin_max = [1 if locspinflag[s] == 0 else locspinflag[s]
                          for s in range(nsite)]
        self._is_local_spin = [locspinflag[s] != 0 for s in range(nsite)]

    # ------------------------------------------------------------------
    #  Low-level helpers
    # ------------------------------------------------------------------

    def spin_max(self, site: int) -> int:
        """Return the maximum spin index for *site*.

        For itinerant sites (``locspinflag == 0``), returns 1.
        For local-spin sites, returns the local-spin value.

        Parameters
        ----------
        site : int
            Site index.

        Returns
        -------
        int
            Maximum spin index (inclusive upper bound).
        """
        return self._spin_max[site]

    def skip_local_spin_pair(self, site_a: int, site_b: int) -> bool:
        """Return True if a pair of distinct local-spin sites should be skipped.

        Parameters
        ----------
        site_a : int
            First site index.
        site_b : int
            Second site index.

        Returns
        -------
        bool
            ``True`` if both sites are local-spin and distinct.
        """
        return (site_a != site_b
                and self._is_local_spin[site_a]
                and self._is_local_spin[site_b])

    def kondo_site(self, isite: int) -> int:
        """Map a Kondo unit-cell index to the physical site index.

        Parameters
        ----------
        isite : int
            Unit-cell site index (``0 .. 2*NsiteUC - 1``).

        Returns
        -------
        int
            Physical site index.
        """
        if isite >= self.NsiteUC:
            return isite - self.NsiteUC + self.nsite // 2
        return isite

    # ------------------------------------------------------------------
    #  One-body index generators
    # ------------------------------------------------------------------

    def green1_corr(self) -> list[tuple[int, int, int, int]]:
        """Generate one-body Green-function indices for correlation mode.

        Only same-spin pairs (``ispin == jspin``) are kept.  For Kondo
        models the unit-cell range is doubled to include localized sites.

        Returns
        -------
        list of tuple[int, int, int, int]
            Index tuples ``(isite, ispin, jsite, jspin)``.
        """
        indices: list[tuple[int, int, int, int]] = []
        xkondo = 2 if self.is_kondo else 1
        is_local = self._is_local_spin

        for isite in range(self.NsiteUC * xkondo):
            isite2 = self.kondo_site(isite)
            for ispin in range(self._spin_max[isite2] + 1):
                for jsite in range(self.nsite):
                    if is_local[isite2] and is_local[jsite] and isite2 != jsite:
                        continue
                    if ispin <= self._spin_max[jsite]:
                        indices.append((isite2, ispin, jsite, ispin))

        return indices

    def green1_raw(self) -> list[tuple[int, int, int, int]]:
        """Generate one-body Green-function indices for raw (full) mode.

        All site-spin combinations are emitted, subject to the constraint
        that pairs of distinct local-spin sites are skipped.

        Returns
        -------
        list of tuple[int, int, int, int]
            Index tuples ``(isite, ispin, jsite, jspin)``.
        """
        site_spins = [(s, sp) for s in range(self.nsite)
                      for sp in range(self._spin_max[s] + 1)]
        is_local = self._is_local_spin
        indices: list[tuple[int, int, int, int]] = [
            (isite, ispin, jsite, jspin)
            for isite, ispin in site_spins
            for jsite, jspin in site_spins
            if not (is_local[isite] and is_local[jsite] and isite != jsite)
        ]
        return indices

    # ------------------------------------------------------------------
    #  Two-body index generators
    # ------------------------------------------------------------------

    def green2_corr(self) -> list[tuple[int, int, int, int, int, int, int, int]]:
        """Generate two-body Green-function indices for correlation mode.

        For each pair ``(site1, site3)`` only spin combinations satisfying
        ``spin1 - spin2 + spin3 - spin4 == 0`` are kept.  For mVMC, when
        ``spin1 != spin2`` or ``spin3 != spin4`` a special index order is
        used; otherwise the standard HPhi order is used.

        Returns
        -------
        list of tuple[int, int, int, int, int, int, int, int]
            Index tuples of 8 elements each.
        """
        indices: list[tuple[int, int, int, int, int, int, int, int]] = []
        xkondo = 2 if self.is_kondo else 1
        is_mvmc = self.is_mvmc
        spin_max = self._spin_max

        for site1 in range(self.NsiteUC * xkondo):
            site1k = self.kondo_site(site1)
            S1Max = spin_max[site1k]

            for spin1 in range(S1Max + 1):
                for spin2 in range(S1Max + 1):

                    for site3 in range(self.nsite):
                        S3Max = spin_max[site3]

                        # spin4 = spin1 - spin2 + spin3
                        for spin3 in range(S3Max + 1):
                            spin4 = spin1 - spin2 + spin3
                            if spin4 < 0 or spin4 > S3Max:
                                continue
                            if is_mvmc and (
                                spin1 != spin2 or spin3 != spin4
                            ):
                                indices.append((
                                    site1k, spin1,
                                    site3, spin4,
                                    site3, spin3,
                                    site1k, spin2,
                                ))
                            else:
                                indices.append((
                                    site1k, spin1,
                                    site1k, spin2,
                                    site3, spin3,
                                    site3, spin4,
                                ))

        return indices

    def green2_raw(self) -> list[tuple[int, int, int, int, int, int, int, int]]:
        """Generate two-body Green-function indices for raw (full) mode.

        All four-site combinations are emitted with local-spin pair
        filtering on ``(site1, site2)`` and ``(site3, site4)``.

        Returns
        -------
        list of tuple[int, int, int, int, int, int, int, int]
            Index tuples of 8 elements each.
        """
        site_spins = [(s, sp) for s in range(self.nsite)
                      for sp in range(self._spin_max[s] + 1)]
        is_local = self._is_local_spin

        # Precompute valid site pairs (no local-spin pair skip)
        indices: list[tuple[int, int, int, int, int, int, int, int]] = []
        for (s1, sp1), (s2, sp2) in product(site_spins, repeat=2):
            if is_local[s1] and is_local[s2] and s1 != s2:
                continue
            for (s3, sp3), (s4, sp4) in product(site_spins, repeat=2):
                if is_local[s3] and is_local[s4] and s3 != s4:
                    continue
                indices.append((s1, sp1, s2, sp2, s3, sp3, s4, sp4))

        return indices


def print_1_green(StdI: StdIntList) -> None:
    """Write ``greenone.def`` listing one-body Green-function indices.

    This is the Python translation of the C function ``Print1Green()``.
    Index generation is delegated to :class:`GreenFunctionIndices`.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``ioutputmode`` : int -- 0 (none), 1 (correlation), or 2 (raw).
        - ``model`` : ModelType -- model type (Kondo triggers doubled UC).
        - ``NsiteUC`` : int -- number of sites in the unit cell.
        - ``nsite`` : int -- total number of sites.
        - ``locspinflag`` : list of int -- per-site local-spin flag.
    """
    if StdI.ioutputmode == 0:
        return

    gf = GreenFunctionIndices(
        StdI.nsite, StdI.NsiteUC, StdI.locspinflag,
        is_kondo=(StdI.model == ModelType.KONDO),
    )

    if StdI.ioutputmode == 1:
        greenindx = gf.green1_corr()
    else:
        greenindx = gf.green1_raw()

    ngreen = len(greenindx)

    with open("greenone.def", "w") as fp:
        lines = [
            "===============================\n",
            f"NCisAjs {ngreen:10d}\n",
            "===============================\n",
            "======== Green functions ======\n",
            "===============================\n",
        ]
        for i0, s0, i1, s1 in greenindx:
            lines.append(f"{i0:5d} {s0:5d} {i1:5d} {s1:5d}\n")
        fp.write("".join(lines))

    logger.info("    greenone.def is written.")


def print_2_green(StdI: StdIntList) -> None:
    """Write ``greentwo.def`` listing two-body Green-function indices.

    This is the Python translation of the C function ``Print2Green()``.
    Index generation is delegated to :class:`GreenFunctionIndices`.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``ioutputmode`` : int
        - ``model`` : ModelType
        - ``solver`` : SolverType
        - ``NsiteUC`` : int
        - ``nsite`` : int
        - ``locspinflag`` : list of int
    """
    if StdI.ioutputmode not in (1, 2):
        return

    gf = GreenFunctionIndices(
        StdI.nsite, StdI.NsiteUC, StdI.locspinflag,
        is_kondo=(StdI.model == ModelType.KONDO),
        is_mvmc=(StdI.solver == SolverType.mVMC),
    )

    if StdI.ioutputmode == 1:
        greenindx = gf.green2_corr()
    else:
        greenindx = gf.green2_raw()

    ngreen = len(greenindx)
    with open("greentwo.def", "w") as fp:
        lines = [
            "=============================================\n",
            f"NCisAjsCktAltDC {ngreen:10d}\n",
            "=============================================\n",
            "======== Green functions for Sq AND Nq ======\n",
            "=============================================\n",
        ]
        for i0, s0, i1, s1, i2, s2, i3, s3 in greenindx:
            lines.append(
                f"{i0:5d} {s0:5d} {i1:5d} {s1:5d} "
                f"{i2:5d} {s2:5d} {i3:5d} {s3:5d}\n"
            )
        fp.write("".join(lines))

    logger.info("    greentwo.def is written.")


def unsupported_system(model: str, lattice: str) -> None:
    """Print an error message and abort for an unsupported model/lattice pair.

    This is the Python translation of the C function
    ``UnsupportedSystem()``.

    Parameters
    ----------
    model : str
        The model name specified by the user.
    lattice : str
        The lattice name specified by the user.

    Raises
    ------
    ValueError
        Always raised after logging the error message.
    """
    msg = (
        f"Unsupported combination in the STANDARD MODE: MODEL = {model}, "
        f"LATTICE = {lattice}. "
        "Please use the EXPERT MODE, or write a NEW FUNCTION and post us."
    )
    logger.error(msg)
    raise ValueError(msg)


def check_output_mode(StdI: StdIntList) -> None:
    """Verify and set the integer output-mode flag from the string keyword.

    This is the Python translation of the C function
    ``CheckOutputMode()``.  The mapping is:

    - ``"non"`` / ``"none"`` / ``"off"`` -> ``ioutputmode = 0``
    - ``"cor"`` / ``"corr"`` / ``"correlation"`` -> ``ioutputmode = 1``
    - ``None`` (unset, default) -> ``ioutputmode = 1``
    - ``"raw"`` / ``"all"`` / ``"full"`` -> ``ioutputmode = 2``
    - anything else -> error and exit

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  ``StdI.outputmode`` is read and
        ``StdI.ioutputmode`` is set **in place**.

    Raises
    ------
    ValueError
        If ``StdI.outputmode`` does not match any recognised keyword.
    """
    if StdI.outputmode is None:
        StdI.ioutputmode = 1
        logger.info(
            "      ioutputmode = %-10d  ######  DEFAULT VALUE IS USED  ######",
            StdI.ioutputmode,
        )
    else:
        mode = OUTPUT_MODE_TO_INT.get(StdI.outputmode)
        if mode is None:
            msg = f"Unsupported OutPutMode : {StdI.outputmode}"
            logger.error(msg)
            raise ValueError(msg)
        StdI.ioutputmode = mode
        logger.info("      ioutputmode = %-10d", StdI.ioutputmode)


def _check_mod_para_hphi(StdI: StdIntList) -> None:
    """Set HPhi-specific default model parameters.

    Handles Lanczos parameters, spectrum frequency grid, and large-value
    cutoff.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure, modified in place.
    """
    StdI.Lanczos_max = print_val_i("Lanczos_max", StdI.Lanczos_max, 2000)
    StdI.initial_iv = print_val_i("initial_iv", StdI.initial_iv, -1)
    # nvec is not given a default here (commented out in C)
    StdI.exct = print_val_i("exct", StdI.exct, 1)
    StdI.LanczosEps = print_val_i("LanczosEps", StdI.LanczosEps, 14)
    StdI.LanczosTarget = print_val_i("LanczosTarget", StdI.LanczosTarget, 2)
    if StdI.LanczosTarget < StdI.exct:
        StdI.LanczosTarget = StdI.exct
    StdI.NumAve = print_val_i("NumAve", StdI.NumAve, 5)
    StdI.ExpecInterval = print_val_i("ExpecInterval", StdI.ExpecInterval, 20)
    StdI.Nomega = print_val_i("NOmega", StdI.Nomega, 200)
    StdI.OmegaMax = print_val_d(
        "OmegaMax", StdI.OmegaMax, StdI.LargeValue * StdI.nsite
    )
    StdI.OmegaMin = print_val_d(
        "OmegaMin", StdI.OmegaMin, -StdI.LargeValue * StdI.nsite
    )
    StdI.OmegaOrg = print_val_d("OmegaOrg", StdI.OmegaOrg, 0.0)
    StdI.OmegaIm = print_val_d(
        "OmegaIm", StdI.OmegaIm, 0.01 * int(StdI.LargeValue)
    )


def _check_mod_para_mvmc(StdI: StdIntList) -> None:
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


def _check_mod_para_uhf(StdI: StdIntList) -> None:
    """Set UHF/H-wave-specific default model parameters.

    Handles random seed, iteration limit, mixing, convergence, and
    symmetry-projection parameters.  UHF and H-wave share identical
    defaults.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure, modified in place.
    """
    StdI.RndSeed = print_val_i("RndSeed", StdI.RndSeed, 123456789)
    StdI.Iteration_max = print_val_i("Iteration_max", StdI.Iteration_max, 1000)
    StdI.mix = print_val_d("Mix", StdI.mix, 0.5)
    StdI.eps = print_val_i("eps", StdI.eps, 8)
    StdI.eps_slater = print_val_i("EpsSlater", StdI.eps_slater, 6)
    StdI.NMPTrans = print_val_i("NMPTrans", StdI.NMPTrans, 0)


# -------------------------------------------------------------------
#  Conserved-quantity validation rules
# -------------------------------------------------------------------
#
class _ConservedQtyRule(NamedTuple):
    """Validation rule for conserved quantities (ncond and 2Sz).

    Attributes
    ----------
    ncond_label : str
        Label used for the ncond check (``"nelec"`` or ``"ncond"``).
    ncond_action : str or None
        Action for ncond: ``"required"``, ``"not_used"``, or ``None``.
    sz2_action : str or None
        Action for 2Sz: ``"required"``, ``"not_used"``, ``"default_0"``,
        or ``None``.
    """

    ncond_label: str
    ncond_action: str | None
    sz2_action: str | None


# The key is ``(ModelType, is_hphi: bool, lGC: int)``.
# For Spin model, is_hphi is ignored (keyed as True and False with same rule).

_CONSERVED_QTY_RULES: dict[tuple, _ConservedQtyRule] = {
    # Hubbard
    (ModelType.HUBBARD, True,  0): _ConservedQtyRule("nelec", "required", None),
    (ModelType.HUBBARD, True,  1): _ConservedQtyRule("nelec", "not_used", "not_used"),
    (ModelType.HUBBARD, False, 0): _ConservedQtyRule("ncond", "required", "default_0"),
    (ModelType.HUBBARD, False, 1): _ConservedQtyRule("ncond", "required", "not_used"),
    # Spin (is_hphi dimension doesn't matter — same rules)
    (ModelType.SPIN,    True,  0): _ConservedQtyRule("ncond", "not_used", "required"),
    (ModelType.SPIN,    True,  1): _ConservedQtyRule("ncond", "not_used", "not_used"),
    (ModelType.SPIN,    False, 0): _ConservedQtyRule("ncond", "not_used", "required"),
    (ModelType.SPIN,    False, 1): _ConservedQtyRule("ncond", "not_used", "not_used"),
    # Kondo
    (ModelType.KONDO,   True,  0): _ConservedQtyRule("ncond", "required", None),
    (ModelType.KONDO,   True,  1): _ConservedQtyRule("nelec", "not_used", "not_used"),
    (ModelType.KONDO,   False, 0): _ConservedQtyRule("ncond", "required", "default_0"),
    (ModelType.KONDO,   False, 1): _ConservedQtyRule("ncond", "required", "not_used"),
}
"""Rules for validating ``ncond`` and ``2Sz`` by (model, is_hphi, lGC)."""

# Dispatch tables for ncond and Sz2 validation actions
_NCOND_ACTION_DISPATCH: dict[str, Callable] = {
    "required": required_val_i,
    "not_used": not_used_i,
}
"""Maps ncond action strings to validation functions (label, value) -> None."""

_SZ2_ACTION_DISPATCH: dict[str, Callable] = {
    "required": lambda StdI: required_val_i("2Sz", StdI.Sz2),
    "not_used": lambda StdI: not_used_i("2Sz", StdI.Sz2),
    "default_0": lambda StdI: setattr(StdI, 'Sz2', print_val_i("2Sz", StdI.Sz2, 0)),
}
"""Maps Sz2 action strings to validation functions (StdI) -> None."""


def _check_conserved_quantities(StdI: StdIntList) -> None:
    """Validate conserved quantities (``ncond`` and ``2Sz``) for the current model.

    The checks depend on ``StdI.model``, ``StdI.solver``, and ``StdI.lGC``
    (grand-canonical flag):

    - **Hubbard** (HPhi, canonical): ``nelec`` required.
    - **Hubbard** (HPhi, GC): ``nelec`` and ``2Sz`` unused.
    - **Hubbard** (other solvers): ``ncond`` required; ``2Sz`` defaults to 0
      in canonical, unused in GC.
    - **Spin**: ``ncond`` unused (set to 0 for mVMC); ``2Sz`` required in
      canonical, unused in GC.
    - **Kondo** (HPhi, canonical): ``ncond`` required.
    - **Kondo** (HPhi, GC): ``nelec`` and ``2Sz`` unused.
    - **Kondo** (other solvers): ``ncond`` required; ``2Sz`` defaults to 0
      in canonical, unused in GC.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure, modified in place.
    """
    is_hphi = (StdI.solver == SolverType.HPhi)
    key = (StdI.model, is_hphi, StdI.lGC)
    rule = _CONSERVED_QTY_RULES.get(key)
    if rule is None:
        return
    ncond_label, ncond_action, sz2_action = rule

    # Apply ncond validation via dispatch
    if ncond_action is not None:
        _NCOND_ACTION_DISPATCH[ncond_action](ncond_label, StdI.ncond)

    # Spin + mVMC: set ncond = 0 (after not_used check, matching C order)
    if StdI.model == ModelType.SPIN and StdI.solver == SolverType.mVMC:
        StdI.ncond = 0

    # Apply Sz2 validation via dispatch
    if sz2_action is not None:
        _SZ2_ACTION_DISPATCH[sz2_action](StdI)


def check_mod_para(StdI: StdIntList) -> None:
    """Validate and set default values for solver-specific model parameters.

    This is the Python translation of the C function ``CheckModPara()``.
    Solver-specific defaults are applied via the plugin's
    :meth:`SolverPlugin.set_defaults`.  Then the conserved quantities
    (``ncond``, ``2Sz``) are checked via :func:`_check_conserved_quantities`.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Many fields are read and
        modified **in place** via the helper functions
        :func:`print_val_i`, :func:`print_val_d`, :func:`not_used_i`,
        and :func:`required_val_i`.
    """
    # ------------------------------------------------------------------
    #  Solver-specific defaults (via plugin)
    # ------------------------------------------------------------------
    from ..plugin import get_plugin
    get_plugin(StdI.solver).set_defaults(StdI)

    # ------------------------------------------------------------------
    #  Conserved quantities: ncond and 2Sz
    # ------------------------------------------------------------------
    _check_conserved_quantities(StdI)
