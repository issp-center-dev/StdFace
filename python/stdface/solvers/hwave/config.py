"""H-wave solver field container (C3), shared by UHFR and UHFK.

Holds every field the H-wave plugins read.  ``fileprefix`` / ``export_all`` /
``lattice_gp`` are H-wave-unique (UHFK Wannier export);
``mix`` / ``eps`` / ``eps_slater`` / ``Iteration_max`` are shared with UHF and
the sublattice / symmetry block with mVMC + UHF.  Per the duplication policy
(``dev/solver_config_split.md`` §2-2) this declares its own copies.

``calcmode`` is intentionally **not** here: it is solver-selection metadata
read by ``HWavePlugin``'s output-mode dispatch (``_hwave_output_mode``)
alongside ``solver``, so it stays on ``StdIntList``.

This single config is registered under ``HWAVE`` as well as the ``UHFR`` /
``UHFK`` aliases, so all three solver names attach (and share) the same
container.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class HWaveConfig:
    # H-wave-unique (UHFK Wannier export)
    fileprefix: str | None = None
    export_all: int | None = None
    lattice_gp: int | None = None

    # Shared with UHF (duplicated copy)
    mix: float | None = None
    eps: int | None = None
    eps_slater: int | None = None
    Iteration_max: int | None = None

    # Sublattice / symmetry — shared with mVMC + UHF (duplicated copy)
    NMPTrans: int | None = None
    RndSeed: int | None = None
    Lsub: int | None = None
    Wsub: int | None = None
    Hsub: int | None = None
    NCellsub: int = 0
    boxsub: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
    rboxsub: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
