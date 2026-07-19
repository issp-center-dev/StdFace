"""mVMC solver field container (C3).

Holds every field the mVMC solver reads.  The variational / SR-opt /
orbital block is mVMC-unique; the sublattice / symmetry block
(``NMPTrans`` / ``RndSeed`` / ``Lsub`` / ``Wsub`` / ``Hsub`` / ``NCellsub`` /
``boxsub`` / ``rboxsub``) is shared with UHF + H-wave and declared here as a
duplicated copy per the policy in ``dev/solver_config_split.md`` §2-2.

The mVMC-unique fields are read only by mVMC-specific code
(``common_writer._{modpara_lines,namelist_entries,check_mod_para}_mvmc`` and
``solvers/mvmc/``), so they resolve via the attached config for mVMC runs.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class MVMCConfig:
    # --- mVMC-unique --------------------------------------------------------
    CParaFileHead: str | None = None
    NVMCCalMode: int | None = None
    NLanczosMode: int | None = None
    NDataIdxStart: int | None = None
    NDataQtySmp: int | None = None
    NSPGaussLeg: int | None = None
    NSROptItrStep: int | None = None
    NSROptItrSmp: int | None = None
    DSROptRedCut: float | None = None
    DSROptStaDel: float | None = None
    DSROptStepDt: float | None = None
    NVMCWarmUp: int | None = None
    NVMCInterval: int | None = None
    NVMCSample: int | None = None
    NExUpdatePath: int | None = None
    NSplitSize: int | None = None
    NSPStot: int | None = None
    NStore: int | None = None
    NSRCG: int | None = None
    ComplexType: int | None = None
    Orb: None = None
    AntiOrb: None = None
    NOrb: int = 0
    NSym: int = 0

    # --- Sublattice / symmetry — shared with UHF + H-wave (duplicated) ------
    NMPTrans: int | None = None
    RndSeed: int | None = None
    Lsub: int | None = None
    Wsub: int | None = None
    Hsub: int | None = None
    NCellsub: int = 0
    boxsub: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
    rboxsub: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
