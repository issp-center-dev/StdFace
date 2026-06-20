"""
Variables used in the Standard mode.

This module defines the StdIntList dataclass which contains all the variables
and parameters used in the Standard mode of HPhi-mVMC-StdFace. In the original
C code these variables are passed as a pointer to the StdIntList structure.

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

from dataclasses import dataclass, field
from enum import Enum

import numpy as np


class ModelType(str, Enum):
    """Canonical model type identifiers.

    Inherits from ``str`` so that ``ModelType.SPIN == "spin"`` is ``True``,
    preserving full backward compatibility with existing string comparisons.

    Attributes
    ----------
    SPIN : str
        Pure spin model (S=1/2 or general S).
    HUBBARD : str
        Fermion Hubbard model.
    KONDO : str
        Kondo lattice model (itinerant + localised spins).
    """

    SPIN = "spin"
    HUBBARD = "hubbard"
    KONDO = "kondo"


class SolverType(str, Enum):
    """Canonical solver type identifiers.

    Inherits from ``str`` so that ``SolverType.HPhi == "HPhi"`` is ``True``,
    preserving full backward compatibility with existing string comparisons.

    Attributes
    ----------
    HPhi : str
        Exact-diagonalisation / Lanczos solver.
    mVMC : str
        Many-variable Variational Monte Carlo solver.
    UHF : str
        Unrestricted Hartree-Fock solver.
    HWAVE : str
        H-wave solver as supplied on the command line; resolved to
        :attr:`UHFR` or :attr:`UHFK` by ``_resolve_solver_name``.
    UHFR : str
        H-wave real-space UHF mode (``calcmode = uhfr``).
    UHFK : str
        H-wave Wannier90 mode (``calcmode = uhfk`` or ``rpa``).
    """

    HPhi = "HPhi"
    mVMC = "mVMC"
    UHF = "UHF"
    HWAVE = "HWAVE"
    UHFR = "UHFR"
    UHFK = "UHFK"


class MethodType(str, Enum):
    """Canonical HPhi calculation method identifiers.

    Inherits from ``str`` so that ``MethodType.LANCZOS == "lanczos"`` is
    ``True``, preserving full backward compatibility with existing string
    comparisons.

    Attributes
    ----------
    LANCZOS : str
        Lanczos diagonalisation.
    LANCZOS_ENERGY : str
        Lanczos (energy-only, skip eigenvector).
    TPQ : str
        Thermal Pure Quantum state method.
    FULLDIAG : str
        Full (direct) diagonalisation.
    CG : str
        Conjugate Gradient method.
    TIME_EVOLUTION : str
        Real-time evolution (pump / quench).
    CTPQ : str
        Canonical Thermal Pure Quantum state method.
    """

    LANCZOS = "lanczos"
    LANCZOS_ENERGY = "lanczosenergy"
    TPQ = "tpq"
    FULLDIAG = "fulldiag"
    CG = "cg"
    TIME_EVOLUTION = "timeevolution"
    CTPQ = "ctpq"


# ---------------------------------------------------------------------------
#  Sentinel constants (matching the C code)
# ---------------------------------------------------------------------------

NaN_i: int = 2147483647
"""Sentinel for an unset integer parameter (same as ``INT_MAX`` in C)."""


def is_unset_or_trivial_d(val: float, trivial: float = 0.0) -> bool:
    """Return True if a float parameter is unset (NaN) or equals *trivial*.

    Helper for :meth:`SolverPlugin.validate` implementations.  An unset
    float parameter is ``None``; a residual ``NaN`` matrix element is
    also treated as unset.
    """
    import math
    return val is None or math.isnan(val) or val == trivial


# ---------------------------------------------------------------------------
#  Numerical tolerances
# ---------------------------------------------------------------------------

AMPLITUDE_EPS: float = 1e-6
"""Threshold for treating an amplitude as non-zero in output files.

Terms with ``abs(coeff) <= AMPLITUDE_EPS`` are suppressed when writing
interaction, transfer, and pump definition files.
"""

ZERO_BODY_EPS: float = 1e-12
"""Threshold for skipping zero-amplitude terms when accumulating
one-body (transfer) and two-body (InterAll) Hamiltonian entries.
"""


# ---------------------------------------------------------------------------
#  Sub-objects of StdIntList (C1: data-structure split)
# ---------------------------------------------------------------------------


def _delegate(sub: str, name: str) -> property:
    """Build a property on ``StdIntList`` delegating to ``self.<sub>.<name>``.

    Used by the C1 façade so existing ``StdI.<name>`` access (read, write,
    and in-place numpy mutation) transparently reaches the sub-object.
    """
    def getter(self):
        return getattr(getattr(self, sub), name)

    def setter(self, value):
        setattr(getattr(self, sub), name, value)

    return property(getter, setter)


@dataclass
class HamiltonianTerms:
    """Hamiltonian term lists and their output flags.

    Holds the transfer / interaction term lists (built during lattice
    setup) and the per-type ``L*`` output flags set when writing the
    ``.def`` files.  ``StdIntList`` delegates to an instance of this class
    via façade properties (see :func:`_delegate`).
    """

    trans_list: list = field(default_factory=list)
    Lintr: int = 0
    intr_list: list = field(default_factory=list)
    LCintra: int = 0
    Cintra_list: list = field(default_factory=list)
    LCinter: int = 0
    Cinter_list: list = field(default_factory=list)
    LHund: int = 0
    Hund_list: list = field(default_factory=list)
    LEx: int = 0
    Ex_list: list = field(default_factory=list)
    LPairLift: int = 0
    PairLift_list: list = field(default_factory=list)
    LPairHopp: int = 0
    PairHopp_list: list = field(default_factory=list)


@dataclass
class LatticeGeometry:
    """Lattice geometry: super-cell box, unit-cell data and boundary phase.

    Holds the lattice name, lattice constants/dimensions, the super-cell
    ``box`` / ``rbox`` / ``direct`` matrices, unit-cell data
    (``NCell`` / ``Cell`` / ``NsiteUC`` / ``tau`` / ``nsite``) and the
    boundary phase vectors.  ``StdIntList`` delegates to an instance of
    this class via façade properties (see :func:`_delegate`).
    """

    lattice: str | None = None
    a: float | None = None
    length: np.ndarray = field(default_factory=lambda: np.zeros(3))
    W: int | None = None
    L: int | None = None
    Height: int | None = None
    direct: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    box: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
    rbox: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
    NCell: int = 0
    Cell: None = None
    NsiteUC: int = 0
    tau: None = None
    nsite: int = 0
    phase: np.ndarray = field(default_factory=lambda: np.zeros(3))
    ExpPhase: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=complex))
    AntiPeriod: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))


@dataclass
class ModelInput:
    """Model Hamiltonian parameters (hoppings, Coulomb, spin couplings, field).

    Holds the user-facing model name and the Hamiltonian coupling
    constants: hoppings ``t*``, Coulomb ``U``/``V*``, spin couplings
    ``J*``/``J*All`` and single-ion ``D``, plus the magnetic-field
    parameters.  ``StdIntList`` delegates to an instance of this class via
    façade properties (see :func:`_delegate`).  Calculation selectors
    (``lGC`` / ``S2`` / ``Sz2`` / ``ncond`` / ``method`` / ``lBoost``)
    remain on ``StdIntList`` for now; they may move here at C2.
    """

    model: str | None = None
    mu: float | None = None

    # Hopping parameters (complex)
    t: complex | None = None
    tp: complex | None = None
    t0: complex | None = None
    t0p: complex | None = None
    t0pp: complex | None = None
    t1: complex | None = None
    t1p: complex | None = None
    t1pp: complex | None = None
    t2: complex | None = None
    t2p: complex | None = None
    t2pp: complex | None = None
    tpp: complex | None = None

    # Coulomb parameters (float)
    U: float | None = None
    V: float | None = None
    Vp: float | None = None
    V0: float | None = None
    V0p: float | None = None
    V0pp: float | None = None
    V1: float | None = None
    V1p: float | None = None
    V1pp: float | None = None
    V2: float | None = None
    V2p: float | None = None
    V2pp: float | None = None
    Vpp: float | None = None

    # Isotropic/anisotropic diagonal spin couplings (float)
    JAll: float | None = None
    JpAll: float | None = None
    J0All: float | None = None
    J0pAll: float | None = None
    J0ppAll: float | None = None
    J1All: float | None = None
    J1pAll: float | None = None
    J1ppAll: float | None = None
    J2All: float | None = None
    J2pAll: float | None = None
    J2ppAll: float | None = None
    JppAll: float | None = None

    # Spin coupling matrices (3x3 float arrays)
    J: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    Jp: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    J0: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    J0p: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    J0pp: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    J1: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    J1p: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    J1pp: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    J2: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    J2p: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    J2pp: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    Jpp: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    D: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))

    # Magnetic field parameters
    h: float | None = None
    Gamma: float | None = None
    Gamma_y: float | None = None
    K: float | None = None


@dataclass
class StdIntList:
    """Main structure containing all parameters and variables for Standard mode.

    Attributes
    ----------
    lattice : str
        Name of lattice. Input parameter.
    a : float
        The lattice constant. Input parameter.
    length : np.ndarray
        Anisotropic lattice constant (shape ``(3,)``), input parameter
        ``wlength``, ``llength``, ``hlength``.
    W : int
        Number of sites along the 1st axis, input parameter.
    L : int
        Number of sites along the 2nd axis, input parameter.
    Height : int
        Number of sites along the 3rd axis, input parameter.
    direct : np.ndarray
        The unit direct lattice vector (shape ``(3, 3)``).
        Set in ``StdFace_InitSite()``.
    box : np.ndarray
        The shape of the super-cell (shape ``(3, 3)``, int).
        Input parameter ``a0W``, ``a0L``, ``a0H``, etc. or defined from
        ``W``, etc. in ``StdFace_InitSite()``.
    rbox : np.ndarray
        The inversion of ``box`` (shape ``(3, 3)``, int).
        Set in ``StdFace_InitSite()``.
    NCell : int
        The number of the unit cell in the super-cell (determinant of ``box``).
        Set in ``StdFace_InitSite()``.
    Cell : np.ndarray or None
        ``[NCell][3]`` The cell position in the fractional coordinate.
        Malloc and set in ``StdFace_InitSite()``.
    NsiteUC : int
        Number of sites in the unit cell. Defined in the beginning of each
        lattice function.
    tau : np.ndarray or None
        Cell-internal site position in the fractional coordinate.
        Defined in the beginning of each lattice function.
    model : str
        Name of model, input parameter.
    mu : float
        Chemical potential, input parameter.
    t : complex
        Nearest-neighbor hopping, input parameter.
    tp : complex
        2nd-nearest hopping, input parameter.
    t0 : complex
        Anisotropic hopping (1st), input parameter.
    t0p : complex
        Anisotropic hopping (2nd), input parameter.
    t0pp : complex
        Anisotropic hopping (3rd), input parameter.
    t1 : complex
        Anisotropic hopping (1st), input parameter.
    t1p : complex
        Anisotropic hopping (2nd), input parameter.
    t1pp : complex
        Anisotropic hopping (3rd), input parameter.
    t2 : complex
        Anisotropic hopping (1st), input parameter.
    t2p : complex
        Anisotropic hopping (2nd), input parameter.
    t2pp : complex
        Anisotropic hopping (3rd), input parameter.
    tpp : complex
        3rd-nearest hopping, input parameter.
    U : float
        On-site Coulomb potential, input parameter.
    V : float
        Off-site Coulomb potential (1st), input parameter.
    Vp : float
        Off-site Coulomb potential (2nd), input parameter.
    V0 : float
        Anisotropic Coulomb potential (1st), input parameter.
    V0p : float
        Anisotropic Coulomb potential (2nd), input parameter.
    V0pp : float
        Anisotropic Coulomb potential (3rd), input parameter.
    V1 : float
        Anisotropic Coulomb potential (1st), input parameter.
    V1p : float
        Anisotropic Coulomb potential (2nd), input parameter.
    V1pp : float
        Anisotropic Coulomb potential (3rd), input parameter.
    V2 : float
        Anisotropic Coulomb potential (1st), input parameter.
    V2p : float
        Anisotropic Coulomb potential (2nd), input parameter.
    V2pp : float
        Anisotropic Coulomb potential (3rd), input parameter.
    Vpp : float
        Off-site Coulomb potential (3rd), input parameter.
    JAll : float
        Isotropic, diagonal spin coupling (1st Near.), input parameter ``J``.
    JpAll : float
        Isotropic, diagonal spin coupling (2nd Near), input parameter ``Jp``.
    J0All : float
        Anisotropic, diagonal spin coupling (1st Near), input parameter ``J0``.
    J0pAll : float
        Anisotropic, diagonal spin coupling (2nd Near), input parameter ``J0'``.
    J0ppAll : float
        Anisotropic, diagonal spin coupling (3rd Near), input parameter ``J0''``.
    J1All : float
        Anisotropic, diagonal spin coupling (1st Near), input parameter ``J1``.
    J1pAll : float
        Anisotropic, diagonal spin coupling (2nd Near), input parameter ``J1'``.
    J1ppAll : float
        Anisotropic, diagonal spin coupling (3rd Near), input parameter ``J1''``.
    J2All : float
        Anisotropic, diagonal spin coupling (1st Near), input parameter ``J2``.
    J2pAll : float
        Anisotropic, diagonal spin coupling (2nd Near), input parameter ``J2'``.
    J2ppAll : float
        Anisotropic, diagonal spin coupling (3rd Near), input parameter ``J2''``.
    JppAll : float
        Isotropic, diagonal spin coupling (3rd Near), input parameter ``J''``.
    J : np.ndarray
        Isotropic, diagonal/off-diagonal spin coupling (1st Near.)
        (shape ``(3, 3)``).
    Jp : np.ndarray
        Isotropic, diagonal/off-diagonal spin coupling (2nd Near.)
        (shape ``(3, 3)``).
    J0 : np.ndarray
        Anisotropic, diagonal/off-diagonal spin coupling (1st Near.)
        (shape ``(3, 3)``).
    J0p : np.ndarray
        Anisotropic, diagonal/off-diagonal spin coupling (2nd Near.)
        (shape ``(3, 3)``).
    J0pp : np.ndarray
        Anisotropic, diagonal/off-diagonal spin coupling (3rd Near.)
        (shape ``(3, 3)``).
    J1 : np.ndarray
        Anisotropic, diagonal/off-diagonal spin coupling (1st Near.)
        (shape ``(3, 3)``).
    J1p : np.ndarray
        Anisotropic, diagonal/off-diagonal spin coupling (2nd Near.)
        (shape ``(3, 3)``).
    J1pp : np.ndarray
        Anisotropic, diagonal/off-diagonal spin coupling (3rd Near.)
        (shape ``(3, 3)``).
    J2 : np.ndarray
        Anisotropic, diagonal/off-diagonal spin coupling (1st Near.)
        (shape ``(3, 3)``).
    J2p : np.ndarray
        Anisotropic, diagonal/off-diagonal spin coupling (2nd Near.)
        (shape ``(3, 3)``).
    J2pp : np.ndarray
        Anisotropic, diagonal/off-diagonal spin coupling (3rd Near.)
        (shape ``(3, 3)``).
    Jpp : np.ndarray
        Isotropic, diagonal/off-diagonal spin coupling (3rd Near.)
        (shape ``(3, 3)``).
    D : np.ndarray
        Coefficient for S_iz S_iz (shape ``(3, 3)``).
        Only ``D[2][2]`` is used.
    h : float
        Longitudinal magnetic field, input parameter.
    Gamma : float
        Transverse magnetic field, input parameter.
    Gamma_y : float
        Transverse magnetic field (y), input parameter.
    K : float
        4-spin term. Not used.
    phase : np.ndarray
        Boundary phase (shape ``(3,)``), input parameter ``phase0``, etc.
    ExpPhase : np.ndarray
        exp(i * pi * phase / 180) (shape ``(3,)``, complex).
    AntiPeriod : np.ndarray
        If corresponding ``phase`` == 180, it becomes 1 (shape ``(3,)``, int).
    nsite : int
        Number of sites, set in each lattice file.
    locspinflag : np.ndarray or None
        ``[nsite]`` LocSpin in Expert mode.
    trans_list : list of tuple
        One-body transfer terms as ``(amp, isite, ispin, jsite, jspin)``
        tuples (``amp`` is complex).  Replaces the former
        ``trans`` / ``transindx`` / ``ntrans`` arrays + counter.
    Lintr : int
        Print ``interall.def`` or not.
    intr_list : list of tuple
        General two-body (InterAll) terms as
        ``(amp, i1, s1, i2, s2, i3, s3, i4, s4)`` tuples.  Replaces the
        former ``intr`` / ``intrindx`` / ``nintr`` arrays + counter.
    LCintra : int
        Print ``coulombintra.def`` or not.
    Cintra_list : list of tuple
        Intra-site Coulomb terms as ``(coeff, isite)`` tuples.
    LCinter : int
        Print ``coulombinter.def`` or not.
    Cinter_list : list of tuple
        Inter-site Coulomb terms as ``(coeff, isite, jsite)`` tuples.
    LHund : int
        Print ``hund.def`` or not.
    Hund_list : list of tuple
        Hund terms as ``(coeff, isite, jsite)`` tuples.
    LEx : int
        Print ``exchange.def`` or not.
    Ex_list : list of tuple
        Exchange terms as ``(coeff, isite, jsite)`` tuples.
    LPairLift : int
        Print ``pairlift.def`` or not.
    PairLift_list : list of tuple
        Pair-lift terms as ``(coeff, isite, jsite)`` tuples.
    LPairHopp : int
        Print ``pairhopp.def`` or not.
    PairHopp_list : list of tuple
        Pair-hopping terms as ``(coeff, isite, jsite)`` tuples.

    The six lists above replace the former ``X`` / ``XIndx`` / ``NX``
    arrays + counters.

    lBoost : int
        Boost flag.
    ncond : int
        Number of electrons, input from file.
    lGC : int
        Switch for computing grand-canonical ensemble (== 1).
    S2 : int
        Total spin |S| of a local spin, input from file.
    outputmode : str
        Select amount of correlation function, input from file.
    CDataFileHead : str
        Header of the output files, input from file.
    Sz2 : int
        Total Sz, input from file.
    ioutputmode : int
        Switch associated to ``outputmode``.
    cutoff_t : float
        Cutoff for the hopping in wannier90, input from file.
    cutoff_u : float
        Cutoff for the Coulomb in wannier90, input from file.
    cutoff_j : float
        Cutoff for the Hund in wannier90, input from file.
    cutoff_length_t : float
        Cutoff for R in wannier90, input from file.
    cutoff_length_U : float
        Cutoff for R in wannier90, input from file.
    cutoff_length_J : float
        Cutoff for R in wannier90, input from file.
    cutoff_tR : np.ndarray
        Cutoff R-vector for hopping (shape ``(3,)``, int).
    cutoff_UR : np.ndarray
        Cutoff R-vector for Coulomb (shape ``(3,)``, int).
    cutoff_JR : np.ndarray
        Cutoff R-vector for Hund (shape ``(3,)``, int).
    cutoff_tVec : np.ndarray
        Cutoff vector for hopping (shape ``(3, 3)``).
    cutoff_UVec : np.ndarray
        Cutoff vector for Coulomb (shape ``(3, 3)``).
    cutoff_JVec : np.ndarray
        Cutoff vector for Hund (shape ``(3, 3)``).
    lambda_ : float
        Tuning parameter of U and J in wannier90, input from file.
        Named ``lambda_`` because ``lambda`` is a Python keyword.
    lambda_U : float
        Tuning parameter of U in wannier90, input from file.
    lambda_J : float
        Tuning parameter of J in wannier90, input from file.
    double_counting_mode : str
        Select mode of double counting, input from file.
    alpha : float
        Tuning parameter of chemical potential correction in wannier90,
        input from file.
    solver : str
        Which solver is active: ``"HPhi"``, ``"mVMC"``, ``"UHF"``,
        or ``"HWAVE"``.
    method : str
        (HPhi) The name of method, input from file.
    Restart : str
        (HPhi) The name of restart mode, input from file.
    InitialVecType : str
        (HPhi) The name of initial-guess type, input from file.
    EigenVecIO : str
        (HPhi) The name of I/O mode for eigenvector, input from file.
    HamIO : str
        (HPhi) The name of I/O mode for Hamiltonian, input from file.
    FlgTemp : int
        (HPhi) Temperature flag.
    Lanczos_max : int
        (HPhi) The maximum number of iterations, input from file.
    initial_iv : int
        (HPhi) The number for generating random number, input from file.
    nvec : int
        (HPhi) Number of vectors.
    exct : int
        (HPhi) The number of eigenvectors to be computed, input from file.
    LanczosEps : int
        (HPhi) Convergence threshold for the Lanczos method.
    LanczosTarget : int
        (HPhi) Which eigenvector is used for the convergence check.
    NumAve : int
        (HPhi) Number of trials for TPQ calculation.
    ExpecInterval : int
        (HPhi) Interval for the iteration when the expectation value is
        computed.
    LargeValue : float
        (HPhi) The shift parameter for the TPQ calculation.
    NGPU : int
        (HPhi) Number of GPU for FullDiag.
    Scalapack : int
        (HPhi) Flag for FullDiag w/ Scalapack.
    list_6spin_pair : None
        (HPhi) Boost 6-spin pair list (dynamic 3D int array).
    list_6spin_star : None
        (HPhi) Boost 6-spin star list (dynamic 2D int array).
    num_pivot : int
        (HPhi) Boost pivot number.
    ishift_nspin : int
        (HPhi) Boost spin shift.
    CalcSpec : str
        (HPhi) The name of mode for spectrum, input from file.
    SpectrumType : str
        (HPhi) The type of mode for spectrum, input from file.
    Nomega : int
        (HPhi) Number of frequencies, input from file.
    OmegaMax : float
        (HPhi) Maximum of frequency for spectrum, input from file.
    OmegaMin : float
        (HPhi) Minimum of frequency for spectrum, input from file.
    OmegaOrg : float
        (HPhi) Origin of frequency for spectrum, input from file.
    OmegaIm : float
        (HPhi) Imaginary part of frequency.
    SpectrumQ : np.ndarray
        (HPhi) Wavenumber (q-vector) in fractional coordinate
        (shape ``(3,)``).
    SpectrumBody : int
        (HPhi) One- or two-body excitation, defined from ``SpectrumType``.
    OutputExVec : str
        (HPhi) The name of output mode for the excited vector, input from file.
    dt : float
        (HPhi) Time step.
    tshift : float
        (HPhi) Shift of time-step of laser.
    tdump : float
        (HPhi) Time scale of dumping.
    freq : float
        (HPhi) Frequency of laser.
    Uquench : float
        (HPhi) Quenched on-site potential.
    VecPot : np.ndarray
        (HPhi) Vector potential (shape ``(3,)``).
    PumpType : str
        (HPhi) The type of pump.
    PumpBody : int
        (HPhi) One- or two-body pumping, defined from ``PumpType``.
    npump : None
        (HPhi) Number of pump transfers (dynamic 1D int array).
    pumpindx : None
        (HPhi) Site/spin indices for pump (dynamic 3D int array).
    pump : None
        (HPhi) Coefficient for pump (dynamic 2D complex array).
    At : None
        (HPhi) Vector potential time series (dynamic 2D float array).
    ExpandCoef : int
        (HPhi) The number of Hamiltonian-vector operations for time evolution.
    CParaFileHead : str
        (mVMC) Header of the optimized wavefunction, input from file.
    NVMCCalMode : int
        (mVMC) Optimization (=0) or compute correlation function (=1),
        input from file.
    NLanczosMode : int
        (mVMC) Power Lanczos (=1), input from file.
    NDataIdxStart : int
        (mVMC) Start index of trials, input from file.
    NDataQtySmp : int
        (mVMC) Number of trials, input from file.
    NSPGaussLeg : int
        (mVMC) Number of Gauss-Legendre points for spin projection,
        input from file.
    NMPTrans : int
        (mVMC/UHF/HWAVE) Number of translation symmetry.
    NSROptItrStep : int
        (mVMC) Number of iterations for stochastic reconfiguration.
    NSROptItrSmp : int
        (mVMC) Number of steps for sampling.
    NSROptFixSmp : int
        (mVMC) Stochastic reconfiguration parameter.
    DSROptRedCut : float
        (mVMC) Stochastic reconfiguration parameter, input from file.
    DSROptStaDel : float
        (mVMC) Stochastic reconfiguration parameter, input from file.
    DSROptStepDt : float
        (mVMC) Stochastic reconfiguration parameter, input from file.
    NVMCWarmUp : int
        (mVMC) VMC warm-up steps.
    NVMCInterval : int
        (mVMC) VMC interval.
    NVMCSample : int
        (mVMC) VMC sample count.
    NExUpdatePath : int
        (mVMC) Exchange update path.
    RndSeed : int
        (mVMC/UHF/HWAVE) Random seed.
    NSplitSize : int
        (mVMC) Split size.
    NSPStot : int
        (mVMC) Total spin S.
    NStore : int
        (mVMC) Store flag.
    NSRCG : int
        (mVMC) SR-CG flag.
    ComplexType : int
        (mVMC) Complex type flag.
    Lsub : int
        (mVMC/UHF/HWAVE) Sublattice L.
    Wsub : int
        (mVMC/UHF/HWAVE) Sublattice W.
    Hsub : int
        (mVMC/UHF/HWAVE) Sublattice H.
    NCellsub : int
        (mVMC/UHF/HWAVE) Number of cells in a sublattice.
    boxsub : np.ndarray
        (mVMC/UHF/HWAVE) Sublattice box (shape ``(3, 3)``, int).
    rboxsub : np.ndarray
        (mVMC/UHF/HWAVE) Sublattice inverse box (shape ``(3, 3)``, int).
    Orb : np.ndarray or None
        (mVMC) ``[nsite][nsite]`` Orbital index (dynamic 2D int array).
    AntiOrb : np.ndarray or None
        (mVMC) ``[nsite][nsite]`` Anti-periodic switch (dynamic 2D int array).
    NOrb : int
        (mVMC) Number of independent orbital index.
    NSym : int
        (mVMC) Number of translation symmetries.
    mix : float
        (UHF/HWAVE) Linear mixing ratio for update.
    eps : int
        (UHF/HWAVE) Convergence threshold for Green's functions.
    eps_slater : int
        (UHF/HWAVE) Convergence threshold for Slater's functions.
    Iteration_max : int
        (UHF/HWAVE) Max number for iterations.
    calcmode : str
        (HWAVE) Calculation mode: UHF, UHFk.
    fileprefix : str
        (HWAVE) Prefix of output filenames.
    export_all : int
        (HWAVE) Output zero elements in UHFk mode.
    lattice_gp : int
        (HWAVE) Create ``lattice.gp``.
    """

    # ------------------------------------------------------------------
    #  Parameters for LATTICE
    # ------------------------------------------------------------------
    # C1: lattice geometry lives in a sub-object; the names below are
    # façade properties delegating to ``self._lattice`` (see _delegate).
    _lattice: LatticeGeometry = field(default_factory=LatticeGeometry)
    lattice = _delegate("_lattice", "lattice")
    a = _delegate("_lattice", "a")
    length = _delegate("_lattice", "length")
    W = _delegate("_lattice", "W")
    L = _delegate("_lattice", "L")
    Height = _delegate("_lattice", "Height")
    direct = _delegate("_lattice", "direct")
    box = _delegate("_lattice", "box")
    rbox = _delegate("_lattice", "rbox")
    NCell = _delegate("_lattice", "NCell")
    Cell = _delegate("_lattice", "Cell")
    NsiteUC = _delegate("_lattice", "NsiteUC")
    tau = _delegate("_lattice", "tau")

    # ------------------------------------------------------------------
    #  Parameters for MODEL (delegated to _model)
    # ------------------------------------------------------------------
    _model: ModelInput = field(default_factory=ModelInput)
    model = _delegate("_model", "model")
    mu = _delegate("_model", "mu")
    t = _delegate("_model", "t")
    tp = _delegate("_model", "tp")
    t0 = _delegate("_model", "t0")
    t0p = _delegate("_model", "t0p")
    t0pp = _delegate("_model", "t0pp")
    t1 = _delegate("_model", "t1")
    t1p = _delegate("_model", "t1p")
    t1pp = _delegate("_model", "t1pp")
    t2 = _delegate("_model", "t2")
    t2p = _delegate("_model", "t2p")
    t2pp = _delegate("_model", "t2pp")
    tpp = _delegate("_model", "tpp")
    U = _delegate("_model", "U")
    V = _delegate("_model", "V")
    Vp = _delegate("_model", "Vp")
    V0 = _delegate("_model", "V0")
    V0p = _delegate("_model", "V0p")
    V0pp = _delegate("_model", "V0pp")
    V1 = _delegate("_model", "V1")
    V1p = _delegate("_model", "V1p")
    V1pp = _delegate("_model", "V1pp")
    V2 = _delegate("_model", "V2")
    V2p = _delegate("_model", "V2p")
    V2pp = _delegate("_model", "V2pp")
    Vpp = _delegate("_model", "Vpp")
    JAll = _delegate("_model", "JAll")
    JpAll = _delegate("_model", "JpAll")
    J0All = _delegate("_model", "J0All")
    J0pAll = _delegate("_model", "J0pAll")
    J0ppAll = _delegate("_model", "J0ppAll")
    J1All = _delegate("_model", "J1All")
    J1pAll = _delegate("_model", "J1pAll")
    J1ppAll = _delegate("_model", "J1ppAll")
    J2All = _delegate("_model", "J2All")
    J2pAll = _delegate("_model", "J2pAll")
    J2ppAll = _delegate("_model", "J2ppAll")
    JppAll = _delegate("_model", "JppAll")
    J = _delegate("_model", "J")
    Jp = _delegate("_model", "Jp")
    J0 = _delegate("_model", "J0")
    J0p = _delegate("_model", "J0p")
    J0pp = _delegate("_model", "J0pp")
    J1 = _delegate("_model", "J1")
    J1p = _delegate("_model", "J1p")
    J1pp = _delegate("_model", "J1pp")
    J2 = _delegate("_model", "J2")
    J2p = _delegate("_model", "J2p")
    J2pp = _delegate("_model", "J2pp")
    Jpp = _delegate("_model", "Jpp")
    D = _delegate("_model", "D")
    h = _delegate("_model", "h")
    Gamma = _delegate("_model", "Gamma")
    Gamma_y = _delegate("_model", "Gamma_y")
    K = _delegate("_model", "K")

    # ------------------------------------------------------------------
    #  Phase for the boundary (delegated to _lattice)
    # ------------------------------------------------------------------
    phase = _delegate("_lattice", "phase")
    ExpPhase = _delegate("_lattice", "ExpPhase")
    AntiPeriod = _delegate("_lattice", "AntiPeriod")

    # ------------------------------------------------------------------
    #  Transfer, Interaction, Locspin
    # ------------------------------------------------------------------
    nsite = _delegate("_lattice", "nsite")
    locspinflag: None = None
    # C1: Hamiltonian terms live in a sub-object; the 15 names below are
    # façade properties delegating to ``self._terms`` (see _delegate).
    _terms: HamiltonianTerms = field(default_factory=HamiltonianTerms)
    trans_list = _delegate("_terms", "trans_list")
    Lintr = _delegate("_terms", "Lintr")
    intr_list = _delegate("_terms", "intr_list")
    LCintra = _delegate("_terms", "LCintra")
    Cintra_list = _delegate("_terms", "Cintra_list")
    LCinter = _delegate("_terms", "LCinter")
    Cinter_list = _delegate("_terms", "Cinter_list")
    LHund = _delegate("_terms", "LHund")
    Hund_list = _delegate("_terms", "Hund_list")
    LEx = _delegate("_terms", "LEx")
    Ex_list = _delegate("_terms", "Ex_list")
    LPairLift = _delegate("_terms", "LPairLift")
    PairLift_list = _delegate("_terms", "PairLift_list")
    LPairHopp = _delegate("_terms", "LPairHopp")
    PairHopp_list = _delegate("_terms", "PairHopp_list")
    lBoost: int = 0

    # ------------------------------------------------------------------
    #  Calculation conditions
    # ------------------------------------------------------------------
    ncond: int | None = None
    lGC: int = 0
    S2: int | None = None
    outputmode: str | None = None
    CDataFileHead: str | None = None
    Sz2: int | None = None
    ioutputmode: int = 0

    # ------------------------------------------------------------------
    #  Wannier90 mode
    # ------------------------------------------------------------------
    cutoff_t: float | None = None
    cutoff_u: float | None = None
    cutoff_j: float | None = None
    cutoff_length_t: float | None = None
    cutoff_length_U: float | None = None
    cutoff_length_J: float | None = None
    cutoff_tR: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))
    cutoff_UR: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))
    cutoff_JR: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))
    cutoff_tVec: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    cutoff_UVec: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    cutoff_JVec: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    lambda_: float | None = None
    lambda_U: float | None = None
    lambda_J: float | None = None
    double_counting_mode: str | None = None
    alpha: float | None = None

    # ------------------------------------------------------------------
    #  Solver selector
    # ------------------------------------------------------------------
    solver: str = ""

    # ------------------------------------------------------------------
    #  Active solver config (C3: dynamic single-active delegation)
    # ------------------------------------------------------------------
    # Set by ``_attach_solver_config`` (top of ``_reset_vals``) from the
    # config registry, keyed on the raw solver name.  While solver-specific
    # fields still live on this dataclass (pre-C3-5) this stays effectively
    # dormant: with no config registered for a solver it remains ``None`` and
    # attribute access falls through to the normal dataclass fields.
    _solver_cfg: object | None = None

    # ------------------------------------------------------------------
    #  HPhi fields
    # ------------------------------------------------------------------
    method: str | None = None
    Restart: str | None = None
    InitialVecType: str | None = None
    EigenVecIO: str | None = None
    HamIO: str | None = None
    FlgTemp: int = 0
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
    list_6spin_pair: None = None
    list_6spin_star: None = None
    num_pivot: int = 0
    ishift_nspin: int = 0
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

    # ------------------------------------------------------------------
    #  mVMC fields
    # ------------------------------------------------------------------
    CParaFileHead: str | None = None
    NVMCCalMode: int | None = None
    NLanczosMode: int | None = None
    NDataIdxStart: int | None = None
    NDataQtySmp: int | None = None
    NSPGaussLeg: int | None = None
    NMPTrans: int | None = None
    NSROptItrStep: int | None = None
    NSROptItrSmp: int | None = None
    NSROptFixSmp: int = 0
    DSROptRedCut: float | None = None
    DSROptStaDel: float | None = None
    DSROptStepDt: float | None = None
    NVMCWarmUp: int | None = None
    NVMCInterval: int | None = None
    NVMCSample: int | None = None
    NExUpdatePath: int | None = None
    RndSeed: int | None = None
    NSplitSize: int | None = None
    NSPStot: int | None = None
    NStore: int | None = None
    NSRCG: int | None = None
    ComplexType: int | None = None
    Lsub: int | None = None
    Wsub: int | None = None
    Hsub: int | None = None
    NCellsub: int = 0
    boxsub: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
    rboxsub: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
    Orb: None = None
    AntiOrb: None = None
    NOrb: int = 0
    NSym: int = 0

    # ------------------------------------------------------------------
    #  UHF / HWAVE fields
    # ------------------------------------------------------------------
    mix: float | None = None
    eps: int | None = None
    eps_slater: int | None = None
    Iteration_max: int | None = None

    # ------------------------------------------------------------------
    #  HWAVE-only fields
    # ------------------------------------------------------------------
    calcmode: str | None = None
    fileprefix: str | None = None
    export_all: int | None = None
    lattice_gp: int | None = None

    # ------------------------------------------------------------------
    #  C3: dynamic delegation to the active solver config
    # ------------------------------------------------------------------
    def __getattr__(self, name: str):
        # Called only when normal lookup fails (dataclass fields and the C1
        # _delegate properties are resolved first).  Route to the active
        # solver config if it owns the name; otherwise behave like a normal
        # missing attribute.
        cfg = self.__dict__.get("_solver_cfg")
        if cfg is not None and hasattr(cfg, name):
            return getattr(cfg, name)
        raise AttributeError(
            f"{type(self).__name__!r} object has no attribute {name!r}")

    def __setattr__(self, name: str, value) -> None:
        cfg = self.__dict__.get("_solver_cfg")
        if (cfg is not None
                and name not in _STDI_OWN_FIELDS
                and hasattr(cfg, name)):
            setattr(cfg, name, value)
        else:
            object.__setattr__(self, name, value)


# Names that belong to StdIntList itself (dataclass fields + C1 _delegate
# façade properties) and must never be routed to the solver config.
# Computed once at import; reflects whatever fields currently live on the
# class, so it shrinks automatically as solver fields move out in C3-2..C3-5.
import dataclasses as _dataclasses  # noqa: E402

_STDI_OWN_FIELDS: frozenset[str] = frozenset(
    {f.name for f in _dataclasses.fields(StdIntList)}
    | {name for name, val in vars(StdIntList).items()
       if isinstance(val, property)}
)
