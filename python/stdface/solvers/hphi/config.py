"""HPhi solver field container (C3).

Holds every HPhi-unique field: calculation method / restart / I/O modes,
Lanczos / TPQ / FullDiag parameters, spectrum, real-time (pump) evolution,
and Boost.  HPhi does not use the sublattice / symmetry block, so (unlike
the other solver configs) this carries no duplicated shared fields.

Every read of these fields in shared lattice code is guarded by
``StdI.solver == SolverType.HPhi`` (interaction_builder pump blocks,
chain/ladder Boost via an early ``solver != HPhi`` return), so non-HPhi
runs never reach them and they resolve via the attached config for HPhi.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class HPhiConfig:
    # Method / restart / I/O
    method: str | None = None
    Restart: str | None = None
    InitialVecType: str | None = None
    EigenVecIO: str | None = None
    HamIO: str | None = None
    FlgTemp: int = 0

    # Lanczos / TPQ / FullDiag
    Lanczos_max: int | None = None
    initial_iv: int | None = None
    nvec: int | None = None
    exct: int | None = None
    LanczosEps: int | None = None
    LanczosTarget: int | None = None
    NumAve: int | None = None
    ExpecInterval: int | None = None
    LargeValue: float | None = None
    NGPU: int | None = None
    Scalapack: int | None = None

    # Boost
    list_6spin_pair: None = None
    list_6spin_star: None = None
    num_pivot: int = 0
    ishift_nspin: int = 0

    # Spectrum
    CalcSpec: str | None = None
    SpectrumType: str | None = None
    Nomega: int | None = None
    OmegaMax: float | None = None
    OmegaMin: float | None = None
    OmegaOrg: float | None = None
    OmegaIm: float | None = None
    SpectrumQ: np.ndarray = field(default_factory=lambda: np.zeros(3))
    SpectrumBody: int = 0
    OutputExVec: str | None = None

    # Real-time evolution / pump
    dt: float | None = None
    tshift: float | None = None
    tdump: float | None = None
    freq: float | None = None
    Uquench: float | None = None
    VecPot: np.ndarray = field(default_factory=lambda: np.zeros(3))
    PumpType: str | None = None
    PumpBody: int = 0
    npump: None = None
    pumpindx: None = None
    pump: None = None
    At: None = None
    ExpandCoef: int | None = None
