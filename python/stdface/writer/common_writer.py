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
build_loc_spn / build_trans / build_namelist / build_modpara
    Build the corresponding ``XxxData`` for the common definition files.
build_green_one / build_green_two
    Build the one-/two-body Green-function index data.
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
from dataclasses import dataclass
from itertools import product
from pathlib import Path
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


@dataclass
class LocSpnData:
    """Per-site local-spin flags (``locspn.def``).

    Attributes
    ----------
    flags : list of int
        One flag per site (0 = itinerant electron, nonzero = local spin).
    """

    flags: list

    def write(self, directory: Path = Path(".")) -> None:
        nlocspin = sum(1 for f in self.flags if f != 0)
        lines = ["================================ \n",
                 f"NlocalSpin {nlocspin:5d}  \n",
                 "================================ \n",
                 "========i_1LocSpn_0IteElc ====== \n",
                 "================================ \n"]
        for isite, flag in enumerate(self.flags):
            lines.append(f"{isite:5d}  {flag:5d}\n")
        with open(Path(directory) / "locspn.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("    locspn.def is written.")

    def to_dict(self) -> dict:
        return {"flags": list(self.flags)}

    @classmethod
    def from_dict(cls, data: dict) -> "LocSpnData":
        return cls(flags=list(data["flags"]))


def build_loc_spn(StdI: StdIntList) -> LocSpnData:
    """Build :class:`LocSpnData` from ``StdI.locspinflag`` (first ``nsite``)."""
    return LocSpnData(flags=[int(StdI.locspinflag[i]) for i in range(StdI.nsite)])


@dataclass
class TransData:
    """One-body transfer integrals (``trans.def``).

    Attributes
    ----------
    rows : list of tuple
        ``(i, s_i, j, s_j, re, im)`` per non-zero transfer term, with the
        amplitude already split into real / imaginary parts.
    """

    rows: list

    def write(self, directory: Path = Path(".")) -> None:
        """Write ``trans.def`` to *directory*."""
        lines = ["======================== \n",
                 f"NTransfer {len(self.rows):7d}  \n",
                 "======================== \n",
                 "========i_j_s_tijs====== \n",
                 "======================== \n"]
        for i0, s0, i1, s1, re, im in self.rows:
            lines.append(
                f"{i0:5d} {s0:5d} {i1:5d} {s1:5d} "
                f"{re:25.15f} {im:25.15f}\n"
            )
        with open(Path(directory) / "trans.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("      trans.def is written.")

    def to_dict(self) -> dict:
        return {"rows": [list(r) for r in self.rows]}

    @classmethod
    def from_dict(cls, data: dict) -> "TransData":
        return cls(rows=[tuple(r) for r in data["rows"]])


def build_trans(StdI: StdIntList) -> TransData:
    """Build :class:`TransData` from ``StdI.trans_list``.

    Duplicate index quadruples are merged (amplitudes summed) and entries
    below ``AMPLITUDE_EPS`` are suppressed.  This is the single
    ``complex -> (real, imag)`` conversion boundary for transfers.
    """
    ntrans = len(StdI.trans_list)
    vals = np.array([t[0] for t in StdI.trans_list], dtype=complex)
    indx = np.array([t[1:5] for t in StdI.trans_list], dtype=int).reshape(ntrans, 4)
    _merge_duplicate_terms(indx, vals, ntrans)

    rows = []
    for ktrans in range(ntrans):
        val = vals[ktrans]
        if abs(val) > AMPLITUDE_EPS:
            i0, s0, i1, s1 = indx[ktrans]
            rows.append((int(i0), int(s0), int(i1), int(s1),
                         float(val.real), float(val.imag)))
    return TransData(rows=rows)


def _require_expert_plugin(plugin, caller: str) -> None:
    """Reject non-Expert plugins (e.g. H-wave) reaching Expert-only builders."""
    from ..plugin import ExpertModeSolverPlugin
    if not isinstance(plugin, ExpertModeSolverPlugin):
        msg = (f"{caller} requires an ExpertModeSolverPlugin, "
               f"got {type(plugin).__name__} (solver {plugin.name!r})")
        logger.error(msg)
        raise ValueError(msg)


def _namelist_entries_hphi(StdI: StdIntList) -> list:
    """Return HPhi-specific ``(keyword, filename)`` namelist entries."""
    entries: list = [("CalcMod", "calcmod.def")]
    if StdI.SpectrumBody == 1:
        entries.append(("SingleExcitation", "single.def"))
    else:
        entries.append(("PairExcitation", "pair.def"))
    if StdI.method == MethodType.TIME_EVOLUTION:
        if StdI.PumpBody == 1:
            entries.append(("TEOneBody", "teone.def"))
        elif StdI.PumpBody == 2:
            entries.append(("TETwoBody", "tetwo.def"))
    entries.append(("SpectrumVec", f"{StdI.CDataFileHead}_eigenvec_0"))
    if StdI.lBoost == 1:
        entries.append(("Boost", "boost.def"))
    return entries


def _namelist_entries_mvmc(StdI: StdIntList) -> list:
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


# Interaction output-flag attribute -> (namelist keyword, filename)
_INTERACTION_FLAGS: list[tuple[str, str, str]] = [
    ("LCintra",   "CoulombIntra", "coulombintra.def"),
    ("LCinter",   "CoulombInter", "coulombinter.def"),
    ("LHund",     "Hund",         "hund.def"),
    ("LEx",       "Exchange",     "exchange.def"),
    ("LPairLift", "PairLift",     "pairlift.def"),
    ("LPairHopp", "PairHop",      "pairhopp.def"),
    ("Lintr",     "InterAll",     "interall.def"),
]
"""Maps each interaction output flag to its namelist ``(keyword, filename)``.

When the flag equals 1 the corresponding definition file is listed in
``namelist.def``.
"""


@dataclass
class NamelistData:
    """``namelist.def`` content as ordered ``(keyword, filename)`` entries."""

    entries: list

    def to_text(self) -> str:
        return "".join(f"{kw:>16}  {fn}\n" for kw, fn in self.entries)

    def write(self, directory: Path = Path(".")) -> None:
        with open(Path(directory) / "namelist.def", "w") as fp:
            fp.write(self.to_text())
        logger.info("    namelist.def is written.")

    def to_dict(self) -> dict:
        return {"entries": [list(e) for e in self.entries]}

    @classmethod
    def from_dict(cls, data: dict) -> "NamelistData":
        return cls(entries=[tuple(e) for e in data["entries"]])


def build_namelist(StdI: StdIntList) -> NamelistData:
    """Build :class:`NamelistData` from the current build results on *StdI*.

    The interaction entries are selected from the ``L*`` flags (set by
    :func:`build_interactions`); the Green entries from ``ioutputmode``;
    solver-specific entries from the plugin.  (When the output container
    is introduced, these can be derived from the assembled data objects.)
    """
    from ..plugin import get_plugin
    plugin = get_plugin(StdI.solver)
    _require_expert_plugin(plugin, "build_namelist")

    entries: list = [
        ("ModPara", "modpara.def"),
        ("LocSpin", "locspn.def"),
        ("Trans", "trans.def"),
    ]
    for flag_attr, kw, fn in _INTERACTION_FLAGS:
        if getattr(StdI, flag_attr) == 1:
            entries.append((kw, fn))
    if StdI.ioutputmode != 0:
        entries.append(("OneBodyG", "greenone.def"))
        if plugin.has_two_body_green(StdI):
            entries.append(("TwoBodyG", "greentwo.def"))
    entries += plugin.namelist_entries(StdI)
    return NamelistData(entries=entries)


# ---------------------------------------------------------------------------
#  ModParaData: data / format separation for modpara.def
# ---------------------------------------------------------------------------
#
#  A modpara body is generated as an ordered list of *line descriptors*
#  (tuples).  The descriptor carries the semantic value plus a presentation
#  hint (format string); ``to_text`` renders the ``.def`` byte-for-byte and
#  ``to_dict`` exposes only the semantic ``key -> value`` mapping.
#
#    ("sep",)                         -> "--------------------"
#    ("raw", text)                    -> text  (banner, fixed lines)
#    ("kv",  key, value, fmt)         -> f"{key:<15}{value:{fmt}}"
#    ("kv2", key, v1, f1, v2, f2)     -> f"{key:<15}{v1:{f1}} {v2:{f2}}"

_MODPARA_SEP = "--------------------"


def _render_modpara_line(line: tuple) -> str:
    """Render one modpara line descriptor to its ``.def`` text (with newline)."""
    kind = line[0]
    if kind == "sep":
        return _MODPARA_SEP + "\n"
    if kind == "raw":
        return f"{line[1]}\n"
    if kind == "kv":
        _, key, value, fmt = line
        return f"{key:<15}{value:{fmt}}\n"
    if kind == "kv2":
        _, key, v1, f1, v2, f2 = line
        return f"{key:<15}{v1:{f1}} {v2:{f2}}\n"
    raise ValueError(f"unknown modpara line descriptor: {line!r}")


@dataclass
class ModParaData:
    """``modpara.def`` content as ordered line descriptors.

    ``to_text`` renders the file; ``to_dict`` exposes the semantic
    parameter mapping (separators / banner / format hints dropped).
    """

    lines: list

    def to_text(self) -> str:
        return "".join(_render_modpara_line(ln) for ln in self.lines)

    def write(self, directory: Path = Path(".")) -> None:
        with open(Path(directory) / "modpara.def", "w") as fp:
            fp.write(self.to_text())
        logger.info("     modpara.def is written.")

    def to_dict(self) -> dict:
        params: dict = {}
        for ln in self.lines:
            if ln[0] == "kv":
                params[ln[1]] = ln[2]
            elif ln[0] == "kv2":
                params[ln[1]] = [ln[2], ln[4]]
        return {"params": params}


def build_modpara(StdI: StdIntList) -> ModParaData:
    """Build :class:`ModParaData` (common header + solver-specific body)."""
    from ..plugin import get_plugin
    plugin = get_plugin(StdI.solver)
    _require_expert_plugin(plugin, "build_modpara")
    lines: list = [("sep",), ("raw", "Model_Parameters   0"), ("sep",)]
    lines += plugin.modpara_lines(StdI)
    return ModParaData(lines=lines)


def _modpara_lines_hphi(StdI: StdIntList) -> list:
    """Return the HPhi-specific body line descriptors of ``modpara.def``."""
    lines: list = [
        ("raw", "HPhi_Cal_Parameters"),
        ("sep",),
        ("kv", "CDataFileHead", StdI.CDataFileHead, ""),
        ("kv", "CParaFileHead", "zqp", ""),
        ("sep",),
        ("kv", "Nsite", StdI.nsite, "<5d"),
    ]
    if StdI.Sz2 is not None:
        lines.append(("kv", "2Sz", StdI.Sz2, "<5d"))
    if StdI.ncond is not None:
        lines.append(("kv", "Ncond", StdI.ncond, "<5d"))
    lines += [
        ("kv", "Lanczos_max", StdI.Lanczos_max, "<5d"),
        ("kv", "initial_iv", StdI.initial_iv, "<5d"),
    ]
    if StdI.nvec is not None:
        lines.append(("kv", "nvec", StdI.nvec, "<5d"))
    lines += [
        ("kv", "exct", StdI.exct, "<5d"),
        ("kv", "LanczosEps", StdI.LanczosEps, "<5d"),
        ("kv", "LanczosTarget", StdI.LanczosTarget, "<5d"),
        ("kv", "LargeValue", StdI.LargeValue, "<25.15e"),
        ("kv", "NumAve", StdI.NumAve, "<5d"),
        ("kv", "ExpecInterval", StdI.ExpecInterval, "<5d"),
        ("kv", "NOmega", StdI.Nomega, "<5d"),
        ("kv2", "OmegaMax", StdI.OmegaMax, "<25.15e", StdI.OmegaIm, "<25.15e"),
        ("kv2", "OmegaMin", StdI.OmegaMin, "<25.15e", StdI.OmegaIm, "<25.15e"),
        ("kv2", "OmegaOrg", StdI.OmegaOrg, "<25.15e", 0.0, "<25.15e"),
        ("kv", "PreCG", 0, "<5d"),
    ]
    if StdI.method == MethodType.TIME_EVOLUTION:
        lines.append(("kv", "ExpandCoef", StdI.ExpandCoef, "<5d"))
    return lines


def _modpara_lines_mvmc(StdI: StdIntList) -> list:
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


def _modpara_lines_uhf_hwave(StdI: StdIntList) -> list:
    """Return the UHF / H-wave body line descriptors of ``modpara.def``.

    UHF and H-wave share identical content except the banner line, which
    is selected from :data:`_MODPARA_BANNER`.
    """
    lines: list = [
        ("raw", _MODPARA_BANNER[StdI.solver]),
        ("sep",),
        ("kv", "CDataFileHead", StdI.CDataFileHead, ""),
        ("kv", "CParaFileHead", "zqp", ""),
        ("sep",),
        ("kv", "Nsite", StdI.nsite, ""),
    ]
    if StdI.Sz2 is not None:
        lines.append(("kv", "2Sz", StdI.Sz2, "<5d"))
    # UHF/HWAVE emit an unset Ncond as the integer sentinel (legacy format).
    ncond_out = StdI.ncond if StdI.ncond is not None else NaN_i
    lines += [
        ("kv", "Ncond", ncond_out, "<5d"),
        ("kv", "IterationMax", StdI.Iteration_max, ""),
        ("kv", "EPS", StdI.eps, ""),
        ("kv", "Mix", StdI.mix, ".10f"),
        ("kv", "RndSeed", StdI.RndSeed, ""),
        ("kv", "EpsSlater", StdI.eps_slater, ""),
        ("kv", "NMPTrans", StdI.NMPTrans, ""),
    ]
    return lines


_MODPARA_BANNER: dict[str, str] = {
    SolverType.UHF: "UHF_Cal_Parameters",
}
"""Maps the solver type to the ``modpara.def`` banner line (UHF only:
H-wave never emits ``modpara.def``, so no other entry is reachable)."""


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


@dataclass
class GreenOneData:
    """One-body Green-function indices (``greenone.def``)."""

    rows: list  # (i0, s0, i1, s1)

    def write(self, directory: Path = Path(".")) -> None:
        lines = [
            "===============================\n",
            f"NCisAjs {len(self.rows):10d}\n",
            "===============================\n",
            "======== Green functions ======\n",
            "===============================\n",
        ]
        for i0, s0, i1, s1 in self.rows:
            lines.append(f"{i0:5d} {s0:5d} {i1:5d} {s1:5d}\n")
        with open(Path(directory) / "greenone.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("    greenone.def is written.")

    def to_dict(self) -> dict:
        return {"rows": [list(r) for r in self.rows]}

    @classmethod
    def from_dict(cls, data: dict) -> "GreenOneData":
        return cls(rows=[tuple(r) for r in data["rows"]])


def build_green_one(StdI: StdIntList) -> "GreenOneData | None":
    """Build :class:`GreenOneData`, or ``None`` when output is disabled."""
    if StdI.ioutputmode == 0:
        return None

    gf = GreenFunctionIndices(
        StdI.nsite, StdI.NsiteUC, StdI.locspinflag,
        is_kondo=(StdI.model == ModelType.KONDO),
    )
    greenindx = gf.green1_corr() if StdI.ioutputmode == 1 else gf.green1_raw()
    return GreenOneData(rows=[tuple(int(x) for x in idx) for idx in greenindx])


@dataclass
class GreenTwoData:
    """Two-body Green-function indices (``greentwo.def``)."""

    rows: list  # (i0, s0, i1, s1, i2, s2, i3, s3)

    def write(self, directory: Path = Path(".")) -> None:
        lines = [
            "=============================================\n",
            f"NCisAjsCktAltDC {len(self.rows):10d}\n",
            "=============================================\n",
            "======== Green functions for Sq AND Nq ======\n",
            "=============================================\n",
        ]
        for i0, s0, i1, s1, i2, s2, i3, s3 in self.rows:
            lines.append(
                f"{i0:5d} {s0:5d} {i1:5d} {s1:5d} "
                f"{i2:5d} {s2:5d} {i3:5d} {s3:5d}\n"
            )
        with open(Path(directory) / "greentwo.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("    greentwo.def is written.")

    def to_dict(self) -> dict:
        return {"rows": [list(r) for r in self.rows]}

    @classmethod
    def from_dict(cls, data: dict) -> "GreenTwoData":
        return cls(rows=[tuple(r) for r in data["rows"]])


def build_green_two(StdI: StdIntList) -> "GreenTwoData | None":
    """Build :class:`GreenTwoData`, or ``None`` when output is disabled."""
    if StdI.ioutputmode not in (1, 2):
        return None

    gf = GreenFunctionIndices(
        StdI.nsite, StdI.NsiteUC, StdI.locspinflag,
        is_kondo=(StdI.model == ModelType.KONDO),
        is_mvmc=(StdI.solver == SolverType.mVMC),
    )
    greenindx = gf.green2_corr() if StdI.ioutputmode == 1 else gf.green2_raw()
    return GreenTwoData(rows=[tuple(int(x) for x in idx) for idx in greenindx])


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
