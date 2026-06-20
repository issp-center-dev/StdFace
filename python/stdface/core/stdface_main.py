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

import logging

import numpy as np
from collections.abc import Callable
from typing import NamedTuple

from .stdface_vals import (
    StdIntList, ModelType, SolverType, MethodType,
    NaN_i,
)
from ..solvers.hphi.writer import (
    vector_potential as _vector_potential,
)
from ..writer.common_writer import (
    unsupported_system as _unsupported_system,
)
from .keyword_parser import (
    parse_common_keyword as _parse_common_keyword,
    _apply_keyword_table,
)
from .input_source import StanFileSource
from ..lattice import get_lattice as _get_lattice

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
#  Lattice dispatch (via plugin registry)
# ---------------------------------------------------------------------------
#  Backward-compatible dict-like objects that delegate to the lattice plugin
#  registry.  New code should use ``get_lattice(name).setup(StdI)`` and
#  ``get_lattice(name).boost(StdI)`` directly.


class _LatticeDispatchProxy:
    """Dict-like proxy over the lattice plugin registry (setup methods)."""

    def __getitem__(self, key: str) -> Callable[[StdIntList], None]:
        return _get_lattice(key).setup

    def get(self, key: str, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key: str) -> bool:
        try:
            _get_lattice(key)
            return True
        except KeyError:
            return False

    def items(self):
        from ..lattice import get_all_lattices
        result = []
        for plugin in get_all_lattices():
            for alias in plugin.aliases:
                result.append((alias, plugin.setup))
        return result


class _BoostDispatchProxy:
    """Dict-like proxy over the lattice plugin registry (boost methods)."""

    # Only lattices whose boost() is overridden (not the ABC default)
    _BOOST_LATTICES = {"chain", "chainlattice", "honeycomb", "honeycomblattice",
                       "kagome", "kagomelattice", "ladder", "ladderlattice"}

    def __getitem__(self, key: str) -> Callable[[StdIntList], None]:
        if key not in self._BOOST_LATTICES:
            raise KeyError(key)
        return _get_lattice(key).boost

    def get(self, key: str, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key: str) -> bool:
        return key in self._BOOST_LATTICES

    def items(self):
        result = []
        seen = set()
        for alias in self._BOOST_LATTICES:
            try:
                plugin = _get_lattice(alias)
                if id(plugin) not in seen:
                    seen.add(id(plugin))
                    for a in plugin.aliases:
                        if a in self._BOOST_LATTICES:
                            result.append((a, plugin.boost))
            except KeyError:
                pass
        return result


LATTICE_DISPATCH = _LatticeDispatchProxy()
"""Dict-like proxy that maps lattice aliases to setup functions via the plugin registry."""

BOOST_DISPATCH = _BoostDispatchProxy()
"""Dict-like proxy that maps lattice aliases to boost functions via the plugin registry."""

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
# NaN_i (int matrices); float matrices use np.nan


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
#  scalar fields that are reset to None by every
#  call to _reset_vals().  These replace ~80 individual ``StdI.x = val``
#  assignments.

_COMMON_RESET_SCALARS: list[tuple[str, object]] = [
    # Lattice scalars
    ("a", None),
    # Magnetic field
    ("Gamma", None),
    ("Gamma_y", None),
    ("h", None),
    # Lattice dimensions
    ("Height", None),
    ("L", None),
    ("W", None),
    # Isotropic scalar spin couplings
    ("JAll", None),
    ("JpAll", None),
    ("JppAll", None),
    ("J0All", None),
    ("J0pAll", None),
    ("J0ppAll", None),
    ("J1All", None),
    ("J1pAll", None),
    ("J1ppAll", None),
    ("J2All", None),
    ("J2pAll", None),
    ("J2ppAll", None),
    # Anisotropy / single-ion
    ("K", None),
    # Chemical potential / spin
    ("mu", None),
    ("S2", None),
    # Hopping parameters (complex)
    ("t", None),
    ("tp", None),
    ("tpp", None),
    ("t0", None),
    ("t0p", None),
    ("t0pp", None),
    ("t1", None),
    ("t1p", None),
    ("t1pp", None),
    ("t2", None),
    ("t2p", None),
    ("t2pp", None),
    # Coulomb parameters
    ("U", None),
    ("V", None),
    ("Vp", None),
    ("Vpp", None),
    ("V0", None),
    ("V0p", None),
    ("V0pp", None),
    ("V1", None),
    ("V1p", None),
    ("V1pp", None),
    ("V2", None),
    ("V2p", None),
    ("V2pp", None),
    # Calculation conditions
    ("ncond", None),
    ("Sz2", None),
    # Wannier90 cutoffs
    ("cutoff_t", None),
    ("cutoff_u", None),
    ("cutoff_j", None),
    ("cutoff_length_t", None),
    ("cutoff_length_U", None),
    ("cutoff_length_J", None),
    ("lambda_", None),
    ("lambda_U", None),
    ("lambda_J", None),
    ("alpha", None),
]
"""Common scalar field resets — ``setattr(StdI, name, value)``."""


_COMMON_RESET_ARRAYS: list[tuple[str, object]] = [
    # Lattice vectors
    ("length", np.nan),
    ("box", NaN_i),
    ("direct", np.nan),
    # 3x3 spin coupling matrices
    ("J", np.nan),
    ("Jp", np.nan),
    ("Jpp", np.nan),
    ("J0", np.nan),
    ("J0p", np.nan),
    ("J0pp", np.nan),
    ("J1", np.nan),
    ("J1p", np.nan),
    ("J1pp", np.nan),
    ("J2", np.nan),
    ("J2p", np.nan),
    ("J2pp", np.nan),
    # Phase / boundary
    ("phase", np.nan),
    # Wannier90 cutoff vectors
    ("cutoff_tR", NaN_i),
    ("cutoff_UR", NaN_i),
    ("cutoff_JR", NaN_i),
    ("cutoff_tVec", np.nan),
    ("cutoff_UVec", np.nan),
    ("cutoff_JVec", np.nan),
]
"""Common array-fill field resets — ``getattr(StdI, name)[...] = value``."""


def _apply_field_resets(StdI: StdIntList, solver: SolverType) -> None:
    """Apply solver-specific field resets from the plugin registry.

    Parameters
    ----------
    StdI : StdIntList
        The parameter structure to reset (modified in place).
    solver : SolverType
        The solver whose field-reset tables to apply.
    """
    from ..plugin import get_plugin
    plugin = get_plugin(solver)
    for name, value in plugin.reset_scalars:
        setattr(StdI, name, value)
    for name, value in plugin.reset_arrays:
        arr = getattr(StdI, name)
        arr[...] = value


# ===================================================================
#  _reset_vals
# ===================================================================


def _attach_solver_config(StdI: StdIntList) -> None:
    """Attach the active solver's config container to *StdI* (C3).

    Keyed on the *raw* solver name (``StdI.solver``) so it runs before
    ``_resolve_solver_name`` — at which point ``HWAVE`` is not yet split
    into ``UHFR`` / ``UHFK`` (both share the same H-wave config).  A no-op
    when no config is registered for the solver (pre-C3-2 state), leaving
    ``StdI._solver_cfg`` as ``None`` so attribute access falls through to
    the flat dataclass fields.
    """
    # Idempotent: StdIntList.__setattr__ already attaches on solver-set, so
    # this just ensures the config matches the current solver (type-guarded).
    StdI._sync_solver_config(StdI.solver)


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
    # --- Active solver config container (C3; no-op until configs exist) -----
    _attach_solver_config(StdI)

    # --- Common scalar fields (table-driven) --------------------------------
    for name, value in _COMMON_RESET_SCALARS:
        setattr(StdI, name, value)

    # --- Common array fields (table-driven) ---------------------------------
    for name, value in _COMMON_RESET_ARRAYS:
        getattr(StdI, name)[...] = value

    # --- D matrix: zero everywhere except D[2][2] = NaN -------------------
    StdI.D[:, :] = 0.0
    StdI.D[2, 2] = np.nan

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
    ValueError
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


def _build_lattice_and_boost(StdI: StdIntList, solver: str) -> tuple:
    """Run the lattice builder (and HPhi LargeValue/Boost) and return aux data.

    Looks up ``StdI.lattice``, runs ``setup``, then ``post_lattice`` /
    ``validate``.  The lattice-level, solver-independent outputs (gnuplot
    and ``geometry.dat``) are **returned, not written** so the caller
    decides where/whether to emit them.

    Parameters
    ----------
    StdI : StdIntList
        Parameter structure (modified in place).
    solver : str
        Solver name.

    Returns
    -------
    tuple
        ``(gp_data, geo_data, xsf_data)`` — each ``None`` when not
        applicable.  ``xsf_data`` (``lattice.xsf``) is produced only for
        3-D lattices (``ndim == 3``).

    Raises
    ------
    ValueError
        If the lattice is not recognised.
    """
    from ..lattice.geometry_output import build_geometry, build_xsf
    lattice = StdI.lattice
    try:
        lattice_plugin = _get_lattice(lattice)
    except KeyError:
        _unsupported_system(StdI.model, StdI.lattice)  # raises

    gp_data = lattice_plugin.setup(StdI)
    geo_data = build_geometry(StdI)
    xsf_data = build_xsf(StdI) if lattice_plugin.ndim == 3 else None

    from ..plugin import get_plugin
    try:
        plugin = get_plugin(solver)
    except KeyError:
        return gp_data, geo_data, xsf_data
    plugin.post_lattice(StdI)
    plugin.validate(StdI)
    return gp_data, geo_data, xsf_data


# ===================================================================
#  Input file parsing
# ===================================================================


def _parse_solver_keyword_via_plugin(
    keyword: str, value: str, StdI: StdIntList, solver: str
) -> bool:
    """Parse a solver-specific keyword using the plugin registry.

    Parameters
    ----------
    keyword : str
        The lowered keyword string.
    value : str
        The raw value string from the input file.
    StdI : StdIntList
        The global parameter structure, modified in place.
    solver : str
        The solver type.

    Returns
    -------
    bool
        True if the keyword was recognised, False otherwise.
    """
    from ..plugin import get_plugin
    return _apply_keyword_table(get_plugin(solver).keyword_table, keyword, value, StdI)


def _apply_keywords(data: dict, StdI: StdIntList, solver: str) -> None:
    """Apply a ``{keyword: value}`` dict to *StdI*.

    Common keywords are tried first, then solver-specific keywords.
    Alias collisions (two keywords mapping to the same field) are caught
    by the duplicate-checking store helpers.

    Parameters
    ----------
    data : dict
        Keyword/value pairs (keywords already lower-cased).
    StdI : StdIntList
        Parameter structure to populate (modified in place).
    solver : str
        Solver name (``"HPhi"``, ``"mVMC"``, ``"UHF"``, or ``"HWAVE"``).

    Raises
    ------
    ValueError
        If a keyword is unrecognised.
    """
    for keyword, value in data.items():
        keyword = str(keyword).lower()
        if keyword == "solver":
            continue  # solver is selected via the argument, not a keyword
        # The store helpers parse strings (stan.in values are strings); dict /
        # TOML / JSON sources may carry native ints/floats, so normalise here.
        value = str(value)
        logger.info("  KEYWORD : %-20s | VALUE : %s ", keyword, value)
        if not _parse_common_keyword(keyword, value, StdI):
            if not _parse_solver_keyword_via_plugin(keyword, value, StdI, solver):
                msg = f"Unsupported Keyword in Standard mode: {keyword}"
                logger.error(msg)
                raise ValueError(msg)


def _parse_input_file(fname: str, StdI: StdIntList, solver: str) -> None:
    """Open and parse a Standard-mode input file into *StdI*.

    Delegates file reading to :class:`StanFileSource` and applies the
    resulting keyword/value dict to *StdI*.

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
    FileNotFoundError
        If the input file cannot be opened.
    ValueError
        If a line lacks ``=``, a keyword is duplicated, or a keyword is
        unrecognised.
    """
    logger.info("Open Standard-Mode Inputfile %s", fname)
    data = StanFileSource(fname).load()
    _apply_keywords(data, StdI, solver)


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
    FileNotFoundError
        If the input file cannot be opened.
    ValueError
        If a keyword is duplicated or unrecognised, or a required
        parameter is missing.
    """
    # ------------------------------------------------------------------
    #  Initialise
    # ------------------------------------------------------------------
    StdI = StdIntList()
    StdI.solver = solver

    logger.info("######  Input Parameter of Standard Intarface  ######")

    _reset_vals(StdI)
    _parse_input_file(fname, StdI, solver)

    # ------------------------------------------------------------------
    #  Construct Model
    # ------------------------------------------------------------------
    logger.info("#######  Construct Model  #######")

    # CDataFileHead default
    if StdI.CDataFileHead is None:
        StdI.CDataFileHead = "zvo"
        logger.info("    CDataFileHead = %-12s######  DEFAULT VALUE IS USED  ######", "zvo")
    else:
        logger.info("    CDataFileHead = %s", StdI.CDataFileHead)

    _resolve_model_and_method(StdI, solver)

    gp_data, geo_data, xsf_data = _build_lattice_and_boost(StdI, solver)
    # Lattice-level (solver-independent) outputs on their own path.
    if gp_data is not None:
        gp_data.write()
    if geo_data is not None:
        geo_data.write()
    if xsf_data is not None:
        xsf_data.write()

    # ------------------------------------------------------------------
    #  Print Expert input files
    # ------------------------------------------------------------------
    logger.info("######  Print Expert input files  ######")

    from ..plugin import get_plugin
    plugin = get_plugin(solver)
    plugin.write(StdI)

    # ------------------------------------------------------------------
    #  Finalise
    # ------------------------------------------------------------------
    logger.info("######  Input files are generated.  ######")


# ===================================================================
#  generate() — library API (D5)
# ===================================================================


def _make_source(source):
    """Select an :class:`InputSource` from a path / dict.

    ``.toml`` / ``.json`` suffixes pick the matching reader; a ``dict``
    uses :class:`DictSource`; anything else is treated as a ``stan.in``
    file.
    """
    from .input_source import (
        DictSource, TOMLSource, JSONSource, StanFileSource,
    )
    if isinstance(source, dict):
        return DictSource(source)
    from pathlib import Path
    path = Path(source)
    if path.suffix == ".toml":
        return TOMLSource(path)
    if path.suffix == ".json":
        return JSONSource(path)
    return StanFileSource(path)


def _check_solver_conflict(solver_arg: str, data: dict) -> None:
    """Raise if *data* names a solver that differs from *solver_arg*."""
    file_solver = data.get("solver")
    if file_solver is not None and str(file_solver).lower() != str(solver_arg).lower():
        msg = (f"solver conflict: argument={solver_arg!r}, file={file_solver!r}. "
               "Specify solver in only one place.")
        logger.error(msg)
        raise ValueError(msg)


def _build_stdintlist(data: dict, solver: str) -> StdIntList:
    """Build and prepare a :class:`StdIntList` from a keyword dict."""
    StdI = StdIntList()
    StdI.solver = solver
    _reset_vals(StdI)
    _apply_keywords(data, StdI, solver)
    if StdI.CDataFileHead is None:
        StdI.CDataFileHead = "zvo"
    _resolve_model_and_method(StdI, StdI.solver)
    return StdI


def _build_output(StdI: StdIntList):
    """Return the :class:`SolverOutput` for the active solver.

    Each plugin's ``build_output`` returns the right container (the H-wave
    plugin dispatches UHFR ``.def`` vs UHFK Wannier90 by ``calcmode``).
    """
    from ..plugin import get_plugin
    return get_plugin(StdI.solver).build_output(StdI)


def generate(source, solver: str = "HPhi", output_dir=".", output_format=None):
    """Generate solver input from *source* and return a ``SolverOutput``.

    Parameters
    ----------
    source : str | Path | dict
        ``stan.in`` path, ``.toml`` / ``.json`` path, or a keyword dict.
    solver : str
        Solver name (``HPhi`` / ``mVMC`` / ``UHF`` / ``UHFR`` / ``UHFK``).
        Must match any ``solver`` given inside *source*.
    output_dir : str | Path | None
        Where to write the files.  ``None`` writes nothing (data-only).
    output_format : OutputFormat | None
        Output format strategy; defaults to :class:`DefFileFormat`.

    Returns
    -------
    SolverOutput
        The assembled output container (``ExpertModeOutput`` or
        ``WannierModeOutput``).
    """
    import os
    import tempfile
    import contextlib
    from pathlib import Path
    from .output import DefFileFormat

    data = _make_source(source).load()
    _check_solver_conflict(solver, data)

    @contextlib.contextmanager
    def _chdir(path):
        prev = os.getcwd()
        os.chdir(path)
        try:
            yield
        finally:
            os.chdir(prev)

    def _run():
        # The build runs with cwd = the target directory so that the
        # solver-specific files emitted eagerly by ``write_solver_files``
        # (HPhi excitation/calcmod, mVMC variational) land there too -- not
        # just the SolverOutput / gnuplot / geometry written explicitly.
        StdI = _build_stdintlist(data, solver)
        gp_data, geo_data, xsf_data = _build_lattice_and_boost(StdI, StdI.solver)
        out = _build_output(StdI)
        fmt = output_format or DefFileFormat()
        fmt.write_output(out, Path("."))
        if gp_data is not None:
            gp_data.write(Path("."))
        if geo_data is not None:
            geo_data.write(Path("."))
        if xsf_data is not None:
            xsf_data.write(Path("."))
        return out

    if output_dir is None:
        # Data-only: build in a throwaway directory so nothing lands in cwd.
        with tempfile.TemporaryDirectory() as tmp, _chdir(tmp):
            return _run()

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    with _chdir(directory):
        return _run()
