"""UHF solver field container (C3).

Holds every field the UHF solver reads.  None of these are *unique* to
UHF: ``mix`` / ``eps`` / ``eps_slater`` / ``Iteration_max`` are shared with
H-wave (see ``common_writer._modpara_lines_uhf_hwave`` / ``_check_mod_para_uhf``)
and the sublattice / symmetry block is shared with mVMC and H-wave.  Per the
duplication policy (``dev/solver_config_split.md`` §2-2) each solver config
declares its own copy; only one config is ever attached at a time, so the
copies never collide at runtime.

Until C3-5 removes these names from ``StdIntList`` the corresponding dataclass
fields shadow this config (attribute access resolves to ``StdIntList`` first),
so this container is dormant — it is exercised live only after the unified
removal.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class UHFConfig:
    # UHF + H-wave shared
    mix: float | None = None
    eps: int | None = None
    eps_slater: int | None = None
    Iteration_max: int | None = None

    # Sublattice / symmetry — shared with mVMC + H-wave (duplicated copy)
    NMPTrans: int | None = None
    RndSeed: int | None = None
    Lsub: int | None = None
    Wsub: int | None = None
    Hsub: int | None = None
    NCellsub: int = 0
    boxsub: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
    rboxsub: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
