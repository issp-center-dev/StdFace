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

NaN_d: float = float("nan")
"""Sentinel for an unset float parameter (IEEE NaN)."""

NaN_c: complex = complex(float("nan"), 0.0)
"""Sentinel for an unset complex parameter (real part is NaN)."""

UNSET_STRING: str = "****"
"""Sentinel for an unset string parameter (C convention ``"****"``)."""


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
    ntrans : int
        Number of transfer, counted in each lattice file.
    transindx : np.ndarray or None
        ``[ntrans][4]`` Site/spin indices of one-body term.
    trans : np.ndarray or None
        ``[ntrans]`` Coefficient of one-body term (complex).
    nintr : int
        Number of InterAll, counted in each lattice file.
    Lintr : int
        Print ``interall.def`` or not.
    intrindx : np.ndarray or None
        ``[nintr][8]`` Site/spin indices of two-body term.
    intr : np.ndarray or None
        ``[nintr]`` Coefficient of general two-body term (complex).
    NCintra : int
        Number of intra-site Coulomb interaction.
    LCintra : int
        Print ``coulombintra.def`` or not.
    CintraIndx : np.ndarray or None
        ``[NCintra][1]`` Site indices of intra-site Coulomb term.
    Cintra : np.ndarray or None
        ``[NCintra]`` Coefficient of intra-site Coulomb term.
    NCinter : int
        Number of inter-site Coulomb interaction.
    LCinter : int
        Print ``coulombinter.def`` or not.
    CinterIndx : np.ndarray or None
        ``[NCinter][2]`` Site indices of inter-site Coulomb term.
    Cinter : np.ndarray or None
        ``[NCinter]`` Coefficient of inter-site Coulomb term.
    NHund : int
        Number of Hund term.
    LHund : int
        Print ``hund.def`` or not.
    HundIndx : np.ndarray or None
        ``[NHund][2]`` Site indices of Hund term.
    Hund : np.ndarray or None
        ``[NHund]`` Coefficient of Hund term.
    NEx : int
        Number of exchange term.
    LEx : int
        Print ``exchange.def`` or not.
    ExIndx : np.ndarray or None
        ``[NEx][2]`` Site indices of exchange term.
    Ex : np.ndarray or None
        ``[NEx]`` Coefficient of exchange term.
    NPairLift : int
        Number of pair-lift term.
    LPairLift : int
        Print ``pairlift.def`` or not.
    PLIndx : np.ndarray or None
        ``[NPairLift][2]`` Site indices of pair-lift term.
    PairLift : np.ndarray or None
        ``[NPairLift]`` Coefficient of pair-lift term.
    NPairHopp : int
        Number of pair-hopping term.
    LPairHopp : int
        Print ``pairhopp.def`` or not.
    PHIndx : np.ndarray or None
        ``[NPairHopp][2]`` Site indices of pair-hopping term.
    PairHopp : np.ndarray or None
        ``[NPairHopp]`` Coefficient of pair-hopping term.
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
    lattice: str = UNSET_STRING
    a: float = 0.0
    length: np.ndarray = field(default_factory=lambda: np.zeros(3))
    W: int = 0
    L: int = 0
    Height: int = 0
    direct: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    box: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
    rbox: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=int))
    NCell: int = 0
    Cell: None = None
    NsiteUC: int = 0
    tau: None = None

    # ------------------------------------------------------------------
    #  Parameters for MODEL
    # ------------------------------------------------------------------
    model: str = UNSET_STRING
    mu: float = 0.0

    # Hopping parameters (complex)
    t: complex = 0 + 0j
    tp: complex = 0 + 0j
    t0: complex = 0 + 0j
    t0p: complex = 0 + 0j
    t0pp: complex = 0 + 0j
    t1: complex = 0 + 0j
    t1p: complex = 0 + 0j
    t1pp: complex = 0 + 0j
    t2: complex = 0 + 0j
    t2p: complex = 0 + 0j
    t2pp: complex = 0 + 0j
    tpp: complex = 0 + 0j

    # Coulomb parameters (float)
    U: float = 0.0
    V: float = 0.0
    Vp: float = 0.0
    V0: float = 0.0
    V0p: float = 0.0
    V0pp: float = 0.0
    V1: float = 0.0
    V1p: float = 0.0
    V1pp: float = 0.0
    V2: float = 0.0
    V2p: float = 0.0
    V2pp: float = 0.0
    Vpp: float = 0.0

    # Isotropic/anisotropic diagonal spin couplings (float)
    JAll: float = 0.0
    JpAll: float = 0.0
    J0All: float = 0.0
    J0pAll: float = 0.0
    J0ppAll: float = 0.0
    J1All: float = 0.0
    J1pAll: float = 0.0
    J1ppAll: float = 0.0
    J2All: float = 0.0
    J2pAll: float = 0.0
    J2ppAll: float = 0.0
    JppAll: float = 0.0

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
    h: float = 0.0
    Gamma: float = 0.0
    Gamma_y: float = 0.0
    K: float = 0.0

    # ------------------------------------------------------------------
    #  Phase for the boundary
    # ------------------------------------------------------------------
    phase: np.ndarray = field(default_factory=lambda: np.zeros(3))
    ExpPhase: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=complex))
    AntiPeriod: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))

    # ------------------------------------------------------------------
    #  Transfer, Interaction, Locspin
    # ------------------------------------------------------------------
    nsite: int = 0
    locspinflag: None = None
    ntrans: int = 0
    transindx: None = None
    trans: None = None
    nintr: int = 0
    Lintr: int = 0
    intrindx: None = None
    intr: None = None
    NCintra: int = 0
    LCintra: int = 0
    CintraIndx: None = None
    Cintra: None = None
    NCinter: int = 0
    LCinter: int = 0
    CinterIndx: None = None
    Cinter: None = None
    NHund: int = 0
    LHund: int = 0
    HundIndx: None = None
    Hund: None = None
    NEx: int = 0
    LEx: int = 0
    ExIndx: None = None
    Ex: None = None
    NPairLift: int = 0
    LPairLift: int = 0
    PLIndx: None = None
    PairLift: None = None
    NPairHopp: int = 0
    LPairHopp: int = 0
    PHIndx: None = None
    PairHopp: None = None
    lBoost: int = 0

    # ------------------------------------------------------------------
    #  Calculation conditions
    # ------------------------------------------------------------------
    ncond: int = 0
    lGC: int = 0
    S2: int = 0
    outputmode: str = UNSET_STRING
    CDataFileHead: str = UNSET_STRING
    Sz2: int = 0
    ioutputmode: int = 0

    # ------------------------------------------------------------------
    #  Wannier90 mode
    # ------------------------------------------------------------------
    cutoff_t: float = 0.0
    cutoff_u: float = 0.0
    cutoff_j: float = 0.0
    cutoff_length_t: float = 0.0
    cutoff_length_U: float = 0.0
    cutoff_length_J: float = 0.0
    cutoff_tR: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))
    cutoff_UR: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))
    cutoff_JR: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))
    cutoff_tVec: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    cutoff_UVec: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    cutoff_JVec: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    lambda_: float = 0.0
    lambda_U: float = 0.0
    lambda_J: float = 0.0
    double_counting_mode: str = UNSET_STRING
    alpha: float = 0.0

    # ------------------------------------------------------------------
    #  Solver selector
    # ------------------------------------------------------------------
    solver: str = ""

    # ------------------------------------------------------------------
    #  HPhi fields
    # ------------------------------------------------------------------
    method: str = UNSET_STRING
    Restart: str = UNSET_STRING
    InitialVecType: str = UNSET_STRING
    EigenVecIO: str = UNSET_STRING
    HamIO: str = UNSET_STRING
    FlgTemp: int = 0
    Lanczos_max: int = 0
    initial_iv: int = 0
    nvec: int = 0
    exct: int = 0
    LanczosEps: int = 0
    LanczosTarget: int = 0
    NumAve: int = 0
    ExpecInterval: int = 0
    LargeValue: float = 0.0
    NGPU: int = 0
    Scalapack: int = 0
    list_6spin_pair: None = None
    list_6spin_star: None = None
    num_pivot: int = 0
    ishift_nspin: int = 0
    CalcSpec: str = UNSET_STRING
    SpectrumType: str = UNSET_STRING
    Nomega: int = 0
    OmegaMax: float = 0.0
    OmegaMin: float = 0.0
    OmegaOrg: float = 0.0
    OmegaIm: float = 0.0
    SpectrumQ: np.ndarray = field(default_factory=lambda: np.zeros(3))
    SpectrumBody: int = 0
    OutputExVec: str = UNSET_STRING
    dt: float = 0.0
    tshift: float = 0.0
    tdump: float = 0.0
    freq: float = 0.0
    Uquench: float = 0.0
    VecPot: np.ndarray = field(default_factory=lambda: np.zeros(3))
    PumpType: str = UNSET_STRING
    PumpBody: int = 0
    npump: None = None
    pumpindx: None = None
    pump: None = None
    At: None = None
    ExpandCoef: int = 0

    # ------------------------------------------------------------------
    #  mVMC fields
    # ------------------------------------------------------------------
    CParaFileHead: str = UNSET_STRING
    NVMCCalMode: int = 0
    NLanczosMode: int = 0
    NDataIdxStart: int = 0
    NDataQtySmp: int = 0
    NSPGaussLeg: int = 0
    NMPTrans: int = 0
    NSROptItrStep: int = 0
    NSROptItrSmp: int = 0
    NSROptFixSmp: int = 0
    DSROptRedCut: float = 0.0
    DSROptStaDel: float = 0.0
    DSROptStepDt: float = 0.0
    NVMCWarmUp: int = 0
    NVMCInterval: int = 0
    NVMCSample: int = 0
    NExUpdatePath: int = 0
    RndSeed: int = 0
    NSplitSize: int = 0
    NSPStot: int = 0
    NStore: int = 0
    NSRCG: int = 0
    ComplexType: int = 0
    Lsub: int = 0
    Wsub: int = 0
    Hsub: int = 0
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
    mix: float = 0.0
    eps: int = 0
    eps_slater: int = 0
    Iteration_max: int = 0

    # ------------------------------------------------------------------
    #  HWAVE-only fields
    # ------------------------------------------------------------------
    calcmode: str = UNSET_STRING
    fileprefix: str = UNSET_STRING
    export_all: int = 0
    lattice_gp: int = 0
