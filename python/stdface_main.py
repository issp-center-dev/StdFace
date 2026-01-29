"""
Read input file and write files for Expert mode; initialize variables; check parameters.

This module is the Python translation of ``StdFace_main.c``.  It provides the
top-level entry point for the Standard-mode input generator used by HPhi, mVMC,
UHF, and H-wave.

The following lattices are supported:

- 1D Chain
- 1D Ladder
- 2D Tetragonal
- 2D Triangular
- 2D Honeycomb
- 2D Kagome
- 3D Simple Orthorhombic
- 3D Face Centered Orthorhombic
- 3D Pyrochlore

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

from collections.abc import Callable
from typing import NamedTuple

from stdface_vals import (
    StdIntList, ModelType, SolverType, MethodType,
    NaN_i, NaN_d, NaN_c, UNSET_STRING,
)
from param_check import exit_program
from writer.hphi_writer import (
    large_value as _large_value,
    vector_potential as _vector_potential,
)
from writer.common_writer import (
    unsupported_system as _unsupported_system,
)
from writer.solver_writer import get_solver_writer
from keyword_parser import (
    trim_space_quote as _trim_space_quote,
    parse_common_keyword as _parse_common_keyword,
    parse_solver_keyword as _parse_solver_keyword,
)
from lattice import (
    chain_lattice,
    square_lattice,
    ladder,
    triangular_lattice,
    honeycomb_lattice,
    kagome,
    orthorhombic,
    fc_ortho,
    pyrochlore,
    wannier90 as wannier90_mod,
)

# ---------------------------------------------------------------------------
#  Lattice dispatch tables
# ---------------------------------------------------------------------------
#  Maps every recognised lattice alias to its builder function.
#  Each function has the signature ``(StdI: StdIntList) -> None``.
LATTICE_DISPATCH: dict[str, Callable[[StdIntList], None]] = {
    "chain":                        chain_lattice.chain,
    "chainlattice":                 chain_lattice.chain,
    "face-centeredorthorhombic":    fc_ortho.fc_ortho,
    "fcorthorhombic":               fc_ortho.fc_ortho,
    "fco":                          fc_ortho.fc_ortho,
    "face-centeredcubic":           fc_ortho.fc_ortho,
    "fccubic":                      fc_ortho.fc_ortho,
    "fcc":                          fc_ortho.fc_ortho,
    "honeycomb":                    honeycomb_lattice.honeycomb,
    "honeycomblattice":             honeycomb_lattice.honeycomb,
    "kagome":                       kagome.kagome,
    "kagomelattice":                kagome.kagome,
    "ladder":                       ladder.ladder,
    "ladderlattice":                ladder.ladder,
    "orthorhombic":                 orthorhombic.orthorhombic,
    "simpleorthorhombic":           orthorhombic.orthorhombic,
    "cubic":                        orthorhombic.orthorhombic,
    "simplecubic":                  orthorhombic.orthorhombic,
    "pyrochlore":                   pyrochlore.pyrochlore,
    "tetragonal":                   square_lattice.tetragonal,
    "tetragonallattice":            square_lattice.tetragonal,
    "square":                       square_lattice.tetragonal,
    "squarelattice":                square_lattice.tetragonal,
    "triangular":                   triangular_lattice.triangular,
    "triangularlattice":            triangular_lattice.triangular,
    "wannier90":                    wannier90_mod.wannier90,
}
"""Maps lattice name aliases to their builder functions."""

BOOST_DISPATCH: dict[str, Callable[[StdIntList], None]] = {
    "chain":            chain_lattice.chain_boost,
    "chainlattice":     chain_lattice.chain_boost,
    "honeycomb":        honeycomb_lattice.honeycomb_boost,
    "honeycomblattice": honeycomb_lattice.honeycomb_boost,
    "kagome":           kagome.kagome_boost,
    "kagomelattice":    kagome.kagome_boost,
    "ladder":           ladder.ladder_boost,
    "ladderlattice":    ladder.ladder_boost,
}
"""Maps lattice name aliases to their boost builder functions (HPhi only)."""

# ---------------------------------------------------------------------------
#  Model name normalisation table
# ---------------------------------------------------------------------------


class _ModelConfig(NamedTuple):
    """Configuration resolved from a user-facing model name alias.

    Attributes
    ----------
    model : ModelType
        Canonical model type (HUBBARD, SPIN, or KONDO).
    lGC : int
        Grand-canonical flag (0 = canonical, 1 = grand-canonical).
    lBoost : int
        Boost extension flag (0 = off, 1 = on, HPhi only).
    """

    model: ModelType
    lGC: int
    lBoost: int


#  Maps user-facing model name aliases to a ``_ModelConfig``.
#  Entries that require ``solver == SolverType.HPhi`` are in a separate dict.
MODEL_ALIASES: dict[str, _ModelConfig] = {
    "fermionhubbard":   _ModelConfig(ModelType.HUBBARD, 0, 0),
    "hubbard":          _ModelConfig(ModelType.HUBBARD, 0, 0),
    "fermionhubbardgc": _ModelConfig(ModelType.HUBBARD, 1, 0),
    "hubbardgc":        _ModelConfig(ModelType.HUBBARD, 1, 0),
    "spin":             _ModelConfig(ModelType.SPIN,    0, 0),
    "spingc":           _ModelConfig(ModelType.SPIN,    1, 0),
    "kondolattice":     _ModelConfig(ModelType.KONDO,   0, 0),
    "kondo":            _ModelConfig(ModelType.KONDO,   0, 0),
    "kondolatticegc":   _ModelConfig(ModelType.KONDO,   1, 0),
    "kondogc":          _ModelConfig(ModelType.KONDO,   1, 0),
}
"""Maps model name aliases to a :class:`_ModelConfig`."""

MODEL_ALIASES_HPHI_BOOST: dict[str, _ModelConfig] = {
    "spingcboost":      _ModelConfig(ModelType.SPIN,    1, 1),
    "spingccma":        _ModelConfig(ModelType.SPIN,    1, 1),
}
"""HPhi-only model aliases that enable the Boost extension."""

# ---------------------------------------------------------------------------
#  Method name normalisation table (HPhi only)
# ---------------------------------------------------------------------------
METHOD_ALIASES: dict[str, MethodType] = {
    "direct":          MethodType.FULLDIAG,
    "alldiag":         MethodType.FULLDIAG,
    "te":              MethodType.TIME_EVOLUTION,
    "time-evolution":  MethodType.TIME_EVOLUTION,
}
"""Maps HPhi method name aliases to their canonical forms."""

# Sentinel constants imported from stdface_vals:
# NaN_i, NaN_d, NaN_c, UNSET_STRING


# ===================================================================
#  Solver-specific field reset tables
# ===================================================================

# Each solver maps to a list of (field_name, value) pairs for scalar
# assignments, plus a list of (field_name, value) pairs for array-fill
# assignments (``getattr(StdI, name)[:] = value`` or ``[:, :] = value``).
#
# UHF and HWAVE share a common base; HWAVE adds two extra fields.

# ===================================================================
#  Common (solver-independent) field reset tables
# ===================================================================
#
#  _COMMON_RESET_SCALARS lists (field_name, sentinel_value) pairs for
#  scalar fields that are reset to NaN_d, NaN_i, or NaN_c by every
#  call to _reset_vals().  These replace ~80 individual ``StdI.x = val``
#  assignments.

_COMMON_RESET_SCALARS: list[tuple[str, object]] = [
    # Lattice scalars
    ("a", NaN_d),
    # Magnetic field
    ("Gamma", NaN_d),
    ("Gamma_y", NaN_d),
    ("h", NaN_d),
    # Lattice dimensions
    ("Height", NaN_i),
    ("L", NaN_i),
    ("W", NaN_i),
    # Isotropic scalar spin couplings
    ("JAll", NaN_d),
    ("JpAll", NaN_d),
    ("JppAll", NaN_d),
    ("J0All", NaN_d),
    ("J0pAll", NaN_d),
    ("J0ppAll", NaN_d),
    ("J1All", NaN_d),
    ("J1pAll", NaN_d),
    ("J1ppAll", NaN_d),
    ("J2All", NaN_d),
    ("J2pAll", NaN_d),
    ("J2ppAll", NaN_d),
    # Anisotropy / single-ion
    ("K", NaN_d),
    # Chemical potential / spin
    ("mu", NaN_d),
    ("S2", NaN_i),
    # Hopping parameters (complex)
    ("t", NaN_c),
    ("tp", NaN_c),
    ("tpp", NaN_c),
    ("t0", NaN_c),
    ("t0p", NaN_c),
    ("t0pp", NaN_c),
    ("t1", NaN_c),
    ("t1p", NaN_c),
    ("t1pp", NaN_c),
    ("t2", NaN_c),
    ("t2p", NaN_c),
    ("t2pp", NaN_c),
    # Coulomb parameters
    ("U", NaN_d),
    ("V", NaN_d),
    ("Vp", NaN_d),
    ("Vpp", NaN_d),
    ("V0", NaN_d),
    ("V0p", NaN_d),
    ("V0pp", NaN_d),
    ("V1", NaN_d),
    ("V1p", NaN_d),
    ("V1pp", NaN_d),
    ("V2", NaN_d),
    ("V2p", NaN_d),
    ("V2pp", NaN_d),
    # Calculation conditions
    ("ncond", NaN_i),
    ("Sz2", NaN_i),
    # Wannier90 cutoffs
    ("cutoff_t", NaN_d),
    ("cutoff_u", NaN_d),
    ("cutoff_j", NaN_d),
    ("cutoff_length_t", NaN_d),
    ("cutoff_length_U", NaN_d),
    ("cutoff_length_J", NaN_d),
    ("lambda_", NaN_d),
    ("lambda_U", NaN_d),
    ("lambda_J", NaN_d),
    ("alpha", NaN_d),
]
"""Common scalar field resets — ``setattr(StdI, name, value)``."""


_COMMON_RESET_ARRAYS: list[tuple[str, object]] = [
    # Lattice vectors
    ("length", NaN_d),
    ("box", NaN_i),
    ("direct", NaN_d),
    # 3x3 spin coupling matrices
    ("J", NaN_d),
    ("Jp", NaN_d),
    ("Jpp", NaN_d),
    ("J0", NaN_d),
    ("J0p", NaN_d),
    ("J0pp", NaN_d),
    ("J1", NaN_d),
    ("J1p", NaN_d),
    ("J1pp", NaN_d),
    ("J2", NaN_d),
    ("J2p", NaN_d),
    ("J2pp", NaN_d),
    # Phase / boundary
    ("phase", NaN_d),
    # Wannier90 cutoff vectors
    ("cutoff_tR", NaN_i),
    ("cutoff_UR", NaN_i),
    ("cutoff_JR", NaN_i),
    ("cutoff_tVec", NaN_d),
    ("cutoff_UVec", NaN_d),
    ("cutoff_JVec", NaN_d),
]
"""Common array-fill field resets — ``getattr(StdI, name)[...] = value``."""


_UHF_BASE_SCALARS: list[tuple[str, object]] = [
    ("NMPTrans", NaN_i),
    ("RndSeed", NaN_i),
    ("mix", NaN_d),
    ("eps", NaN_i),
    ("eps_slater", NaN_i),
    ("Iteration_max", NaN_i),
    ("Hsub", NaN_i),
    ("Lsub", NaN_i),
    ("Wsub", NaN_i),
]

_UHF_BASE_ARRAYS: list[tuple[str, object]] = [
    ("boxsub", NaN_i),
]

_SOLVER_RESET_SCALARS: dict[SolverType, list[tuple[str, object]]] = {
    SolverType.HPhi: [
        ("LargeValue", NaN_d),
        ("OmegaMax", NaN_d),
        ("OmegaMin", NaN_d),
        ("OmegaOrg", NaN_d),
        ("OmegaIm", NaN_d),
        ("Nomega", NaN_i),
        ("FlgTemp", 1),
        ("Lanczos_max", NaN_i),
        ("initial_iv", NaN_i),
        ("nvec", NaN_i),
        ("exct", NaN_i),
        ("LanczosEps", NaN_i),
        ("LanczosTarget", NaN_i),
        ("NumAve", NaN_i),
        ("ExpecInterval", NaN_i),
        ("dt", NaN_d),
        ("tdump", NaN_d),
        ("tshift", NaN_d),
        ("freq", NaN_d),
        ("Uquench", NaN_d),
        ("ExpandCoef", NaN_i),
        ("NGPU", NaN_i),
        ("Scalapack", NaN_i),
    ],
    SolverType.mVMC: [
        ("NVMCCalMode", NaN_i),
        ("NLanczosMode", NaN_i),
        ("NDataIdxStart", NaN_i),
        ("NDataQtySmp", NaN_i),
        ("NSPGaussLeg", NaN_i),
        ("NSPStot", NaN_i),
        ("NMPTrans", NaN_i),
        ("NSROptItrStep", NaN_i),
        ("NSROptItrSmp", NaN_i),
        ("DSROptRedCut", NaN_d),
        ("DSROptStaDel", NaN_d),
        ("DSROptStepDt", NaN_d),
        ("NVMCWarmUp", NaN_i),
        ("NVMCInterval", NaN_i),
        ("NVMCSample", NaN_i),
        ("NExUpdatePath", NaN_i),
        ("RndSeed", NaN_i),
        ("NSplitSize", NaN_i),
        ("NStore", NaN_i),
        ("NSRCG", NaN_i),
        ("ComplexType", NaN_i),
        ("Hsub", NaN_i),
        ("Lsub", NaN_i),
        ("Wsub", NaN_i),
    ],
    SolverType.UHF: _UHF_BASE_SCALARS,
    SolverType.HWAVE: _UHF_BASE_SCALARS + [
        ("export_all", NaN_i),
        ("lattice_gp", NaN_i),
    ],
}
"""Scalar field resets for each solver — ``setattr(StdI, name, value)``."""

_SOLVER_RESET_ARRAYS: dict[SolverType, list[tuple[str, object]]] = {
    SolverType.HPhi: [
        ("SpectrumQ", NaN_d),
        ("VecPot", NaN_d),
    ],
    SolverType.mVMC: [
        ("boxsub", NaN_i),
    ],
    SolverType.UHF: _UHF_BASE_ARRAYS,
    SolverType.HWAVE: _UHF_BASE_ARRAYS,
}
"""Array-fill field resets for each solver — ``getattr(StdI, name)[...] = value``."""


def _apply_field_resets(StdI: StdIntList, solver: SolverType) -> None:
    """Apply solver-specific field resets from the data tables.

    Parameters
    ----------
    StdI : StdIntList
        The parameter structure to reset (modified in place).
    solver : SolverType
        The solver whose field-reset tables to apply.
    """
    for name, value in _SOLVER_RESET_SCALARS.get(solver, ()):
        setattr(StdI, name, value)
    for name, value in _SOLVER_RESET_ARRAYS.get(solver, ()):
        arr = getattr(StdI, name)
        arr[...] = value


# ===================================================================
#  _reset_vals
# ===================================================================


def _reset_vals(StdI: StdIntList) -> None:
    """Clear / initialize every field in *StdI* to its sentinel value.

    This is the Python translation of the C function
    ``StdFace_ResetVals()``.  Fields that have not been specified by
    the user are filled with NaN sentinels so that duplicate-input
    detection and default-value assignment work correctly later.

    The common (solver-independent) fields are driven by the
    ``_COMMON_RESET_SCALARS`` and ``_COMMON_RESET_ARRAYS`` tables.
    Solver-specific fields are handled by ``_apply_field_resets()``.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure whose fields are reset **in
        place**.
    """
    # --- Common scalar fields (table-driven) --------------------------------
    for name, value in _COMMON_RESET_SCALARS:
        setattr(StdI, name, value)

    # --- Common array fields (table-driven) ---------------------------------
    for name, value in _COMMON_RESET_ARRAYS:
        getattr(StdI, name)[...] = value

    # --- D matrix: zero everywhere except D[2][2] = NaN_d -------------------
    StdI.D[:, :] = 0.0
    StdI.D[2, 2] = NaN_d

    # --- Solver-specific fields (table-driven) ------------------------------
    _apply_field_resets(StdI, StdI.solver)

    # --- Boost (always zero, not NaN) ---------------------------------------
    StdI.lBoost = 0


# Keyword parsing helpers (_trim_space_quote,
# _store_with_check_dup_s/sl/i/d/c) and keyword parsers
# (_parse_common_keyword, _parse_solver_keyword)
# have been moved to keyword_parser.py


# ===================================================================
#  Model / method resolution
# ===================================================================


def _resolve_model_and_method(StdI: StdIntList, solver: str) -> None:
    """Normalise the model name and HPhi method, updating *StdI* in place.

    Looks up ``StdI.model`` in :data:`MODEL_ALIASES` (and, for HPhi,
    :data:`MODEL_ALIASES_HPHI_BOOST`) to resolve the canonical
    :class:`ModelType`, the grand-canonical flag ``lGC``, and the Boost
    flag ``lBoost``.

    For HPhi, also normalises ``StdI.method`` via :data:`METHOD_ALIASES`
    and, if the method is time-evolution, computes the vector potential.

    Parameters
    ----------
    StdI : StdIntList
        Parameter structure (modified in place).
    solver : str
        Solver name.

    Raises
    ------
    SystemExit
        If the model name is not recognised.
    """
    StdI.lGC = 0
    StdI.lBoost = 0

    model_info = MODEL_ALIASES.get(StdI.model)
    if model_info is None and solver == SolverType.HPhi:
        model_info = MODEL_ALIASES_HPHI_BOOST.get(StdI.model)
    if model_info is not None:
        StdI.model, StdI.lGC, StdI.lBoost = model_info
    else:
        _unsupported_system(StdI.model, StdI.lattice)

    if solver == SolverType.HPhi:
        StdI.method = METHOD_ALIASES.get(StdI.method, StdI.method)

        if StdI.method == MethodType.TIME_EVOLUTION:
            _vector_potential(StdI)


# ===================================================================
#  Lattice construction and Boost
# ===================================================================


def _build_lattice_and_boost(StdI: StdIntList, solver: str) -> None:
    """Dispatch to the lattice builder and, for HPhi, apply LargeValue and Boost.

    Looks up ``StdI.lattice`` in :data:`LATTICE_DISPATCH` to generate
    the Hamiltonian definition files.  For the HPhi solver, also
    computes the large value and optionally runs the Boost builder.

    Parameters
    ----------
    StdI : StdIntList
        Parameter structure (modified in place).
    solver : str
        Solver name.

    Raises
    ------
    SystemExit
        If the lattice is not recognised.
    """
    lattice = StdI.lattice
    lattice_builder = LATTICE_DISPATCH.get(lattice)
    if lattice_builder is not None:
        lattice_builder(StdI)
    else:
        _unsupported_system(StdI.model, StdI.lattice)

    if solver == SolverType.HPhi:
        _large_value(StdI)

        if StdI.lBoost == 1:
            boost_builder = BOOST_DISPATCH.get(lattice)
            if boost_builder is not None:
                boost_builder(StdI)
            else:
                _unsupported_system(StdI.model, StdI.lattice)


# ===================================================================
#  Input file parsing
# ===================================================================


def _parse_input_file(fname: str, StdI: StdIntList, solver: str) -> None:
    """Open and parse a Standard-mode input file into *StdI*.

    Each non-blank, non-comment line must contain ``keyword = value``.
    Common keywords are tried first, then solver-specific keywords.
    The program exits on duplicate keywords, missing ``=``, or
    unrecognised keywords.

    Parameters
    ----------
    fname : str
        Path to the Standard-mode input file.
    StdI : StdIntList
        Parameter structure to populate (modified in place).
    solver : str
        Solver name (``"HPhi"``, ``"mVMC"``, ``"UHF"``, or ``"HWAVE"``).

    Raises
    ------
    SystemExit
        If the file cannot be opened, a line lacks ``=``, or a keyword
        is unrecognised.
    """
    try:
        fp_in = open(fname, "r")
    except OSError:
        print(f"\n  ERROR !  Cannot open input file {fname} !\n")
        exit_program(-1)

    print(f"\n  Open Standard-Mode Inputfile {fname} \n")

    with fp_in:
        for raw_line in fp_in:
            line = _trim_space_quote(raw_line)

            if line.startswith("//") or line == "":
                print("  Skipping a line.")
                continue

            parts = line.split("=", 1)
            if len(parts) < 2:
                print('\n  ERROR !  "=" is NOT found !\n')
                exit_program(-1)

            keyword = parts[0].lower()
            value = parts[1]
            print(f"  KEYWORD : {keyword:<20s} | VALUE : {value} ")

            if not _parse_common_keyword(keyword, value, StdI):
                if not _parse_solver_keyword(keyword, value, StdI, solver):
                    print("ERROR ! Unsupported Keyword in Standard mode!")
                    exit_program(-1)


# ===================================================================
#  stdface_main -- top-level entry point
# ===================================================================


def stdface_main(fname: str, solver: str = "HPhi") -> None:
    """Read a Standard-mode input file and generate Expert-mode definition files.

    This is the Python translation of the C function ``StdFace_main()``
    (``StdFace_main.c``, lines 2455--3067).  It performs the following
    sequence of operations:

    1. Create and initialise a :class:`StdIntList` parameter structure.
    2. Parse the input file, mapping keywords to structure fields.
    3. Validate and normalise the model name, lattice name and method.
    4. Dispatch to the appropriate lattice-generation function.
    5. Write all Expert-mode definition files required by the chosen
       solver.

    Parameters
    ----------
    fname : str
        Path to the Standard-mode input file.
    solver : str, optional
        Name of the solver that will consume the generated files.
        Must be one of ``"HPhi"``, ``"mVMC"``, ``"UHF"`` or ``"HWAVE"``.
        Defaults to ``"HPhi"``.  This replaces the compile-time
        ``#ifdef`` mechanism of the original C code.

    Raises
    ------
    SystemExit
        If the input file cannot be opened, a keyword is duplicated or
        unrecognised, or a required parameter is missing.
    """
    # ------------------------------------------------------------------
    #  Initialise
    # ------------------------------------------------------------------
    StdI = StdIntList()
    StdI.solver = solver

    print("\n######  Input Parameter of Standard Intarface  ######")

    _reset_vals(StdI)
    _parse_input_file(fname, StdI, solver)

    # ------------------------------------------------------------------
    #  Construct Model
    # ------------------------------------------------------------------
    print("")
    print("#######  Construct Model  #######")
    print("")

    # CDataFileHead default
    if StdI.CDataFileHead == UNSET_STRING:
        StdI.CDataFileHead = "zvo"
        print(f"    CDataFileHead = {'zvo':<12s}######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"    CDataFileHead = {StdI.CDataFileHead}")

    _resolve_model_and_method(StdI, solver)

    _build_lattice_and_boost(StdI, solver)

    # ------------------------------------------------------------------
    #  Print Expert input files
    # ------------------------------------------------------------------
    print("")
    print("######  Print Expert input files  ######")
    print("")

    writer = get_solver_writer(solver)
    writer.write(StdI)

    # ------------------------------------------------------------------
    #  Finalise
    # ------------------------------------------------------------------
    print("\n######  Input files are generated.  ######\n")
