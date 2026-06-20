"""Solver plugin interface and registry.

This module defines the :class:`SolverPlugin` abstract base class that every
solver plugin must implement, and a plugin registry for discovering and
retrieving plugins by name.

New solver plugins should subclass :class:`SolverPlugin` and implement the
required abstract properties and methods.  The :meth:`write` method provides
a template for the common output sequence (locspn → trans → interactions →
modpara → solver-specific → green → namelist); solvers with a different
sequence can override :meth:`write` entirely.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .core.stdface_vals import StdIntList


class SolverPlugin(ABC):
    """Abstract base class for solver plugins.

    Each solver plugin provides:

    - A keyword table for parsing solver-specific input keywords.
    - Field reset tables for initialising solver-specific fields on
      :class:`StdIntList`.
    - A :meth:`write` method to generate Expert-mode definition files.
    - Optional hooks for post-lattice processing and field initialisation.

    Template Method
    ---------------
    The default :meth:`write` implementation calls a sequence of steps
    that is common to most solvers::

        write_locspn → write_trans → write_interactions →
        check_and_write_modpara → write_solver_specific →
        write_green → write_namelist

    Subclasses should override :meth:`write_solver_specific` for
    solver-specific output (e.g. excitation files, variational parameters).
    Solvers that need a completely different sequence (e.g. H-wave in
    Wannier90 export mode) can override :meth:`write` directly.

    Attributes
    ----------
    name : str
        Canonical solver name (e.g. ``"HPhi"``).
    keyword_table : dict[str, tuple]
        Solver-specific keyword dispatch table.
    reset_scalars : list[tuple[str, object]]
        ``(field_name, sentinel_value)`` pairs for scalar field resets.
    reset_arrays : list[tuple[str, object]]
        ``(field_name, sentinel_value)`` pairs for array-fill field resets.
    """

    # ------------------------------------------------------------------
    #  Abstract properties — must be implemented by every plugin
    # ------------------------------------------------------------------

    @property
    @abstractmethod
    def name(self) -> str:
        """Canonical solver name."""

    @property
    @abstractmethod
    def keyword_table(self) -> dict[str, tuple]:
        """Solver-specific keyword dispatch table."""

    @property
    @abstractmethod
    def reset_scalars(self) -> list[tuple[str, object]]:
        """Scalar field reset table."""

    @property
    @abstractmethod
    def reset_arrays(self) -> list[tuple[str, object]]:
        """Array-fill field reset table."""

    # ------------------------------------------------------------------
    #  Output generation
    # ------------------------------------------------------------------

    @abstractmethod
    def write(self, StdI: StdIntList) -> None:
        """Write all Expert-mode definition files for this solver.

        Each plugin implements its own explicit output sequence.  Expert-mode
        solvers (HPhi/mVMC/UHF) build the common files via
        :meth:`ExpertModeSolverPlugin._write_common_files`; other solvers
        (H-wave) write their own file set.

        Parameters
        ----------
        StdI : StdIntList
            The fully-populated parameter structure.
        """

    # ------------------------------------------------------------------
    #  Optional lifecycle hooks
    # ------------------------------------------------------------------

    def post_lattice(self, StdI: StdIntList) -> None:
        """Optional hook called after lattice construction.

        Override this to perform solver-specific post-lattice processing
        (e.g. computing LargeValue for HPhi, running Boost).
        The default implementation does nothing.

        Parameters
        ----------
        StdI : StdIntList
            The parameter structure (modified in place).
        """

    def init_fields(self, StdI: StdIntList) -> None:
        """Optional hook to initialise solver-specific attributes on StdI.

        Override this to set up dynamic attributes that don't exist on the
        base :class:`StdIntList` dataclass. The default implementation
        does nothing.

        Parameters
        ----------
        StdI : StdIntList
            The parameter structure (modified in place).
        """

    def set_defaults(self, StdI: StdIntList) -> None:
        """Apply solver-specific default model parameters (default: no-op).

        Called from :func:`check_mod_para`.  Lives on the base class (not
        :class:`ExpertModeSolverPlugin`) because non-Expert solvers such as
        H-wave (and UHFR after B4) also run ``check_mod_para``.

        Parameters
        ----------
        StdI : StdIntList
            The parameter structure (modified in place).
        """

    def validate(self, StdI: StdIntList) -> None:
        """Solver-specific validation after lattice construction (default: no-op).

        Called between ``post_lattice`` and ``write``.  Override to reject
        unsupported parameter combinations; raise :class:`ValueError` on
        invalid input.

        Parameters
        ----------
        StdI : StdIntList
            The fully-constructed parameter structure (read-only here).
        """


class ExpertModeSolverPlugin(SolverPlugin):
    """Base class for solvers that emit ``modpara.def`` / ``namelist.def``.

    HPhi, mVMC and UHF inherit from this class.  The solver-specific body
    of ``modpara.def`` and the extra entries of ``namelist.def`` are
    provided by :meth:`modpara_lines` and :meth:`write_namelist_body`;
    :meth:`has_two_body_green` controls whether ``greentwo.def`` is listed.

    These methods replace the former ``_MODPARA_BODY_DISPATCH`` /
    ``_NAMELIST_BODY_DISPATCH`` dispatch tables.
    """

    @abstractmethod
    def modpara_lines(self, StdI: StdIntList) -> list:
        """Return the solver-specific body line descriptors of ``modpara.def``.

        See :func:`stdface.writer.common_writer.build_modpara` for the
        descriptor format.
        """

    def namelist_entries(self, StdI: StdIntList) -> list:
        """Return solver-specific ``(keyword, filename)`` namelist entries.

        Default: no extra entries (used by UHF).
        """
        return []

    def has_two_body_green(self, StdI: StdIntList) -> bool:
        """Whether ``greentwo.def`` is listed in ``namelist.def`` (default True)."""
        return True

    def write_solver_files(self, StdI: StdIntList) -> None:
        """Hook: emit solver-specific files (default: none, used by UHF).

        Runs after the common data is built and before ``namelist`` is
        built, because HPhi's excitation / calcmod set ``SpectrumBody`` /
        ``PumpBody`` which ``build_namelist`` reads.
        """

    def build_output(self, StdI: StdIntList) -> "ExpertModeOutput":
        """Assemble the Expert-mode output container for *StdI*.

        Orchestrates the prep / build sequence so that namelist (built
        last) sees the flags set by the interaction build and the
        solver-specific hook.
        """
        from .writer.common_writer import (
            check_mod_para, build_loc_spn, build_trans, build_modpara,
            check_output_mode, build_green_one, build_green_two, build_namelist,
        )
        from .writer.interaction_writer import build_interactions
        from .core.output import ExpertModeOutput

        check_mod_para(StdI)
        locspn = build_loc_spn(StdI)
        trans = build_trans(StdI)
        interactions = build_interactions(StdI)   # sets L* flags
        modpara = build_modpara(StdI)
        check_output_mode(StdI)
        green_one = build_green_one(StdI)
        green_two = (build_green_two(StdI)
                     if self.has_two_body_green(StdI) else None)
        self.write_solver_files(StdI)             # sets SpectrumBody/PumpBody
        namelist = build_namelist(StdI)
        return ExpertModeOutput(
            locspn=locspn, trans=trans, interactions=interactions,
            modpara=modpara, namelist=namelist,
            green_one=green_one, green_two=green_two,
        )

    def write(self, StdI: StdIntList) -> None:
        """Write all Expert-mode files via the output container."""
        self.build_output(StdI).write()


class WannierModeSolverPlugin(SolverPlugin):
    """Base class for solvers that emit Wannier90-format files (UHFK / RPA).

    Writes ``geometry``/``transfer``/``coulombintra`` files via the
    Wannier90 writer instead of the Expert-mode ``.def`` files.
    """

    def write(self, StdI: StdIntList) -> None:
        from .core.output import build_wannier_output
        build_wannier_output(StdI).write()
        self.write_wannier_extras(StdI)

    def write_wannier_extras(self, StdI: StdIntList) -> None:
        """Write additional Wannier-mode files (default: none)."""


# ---------------------------------------------------------------------------
#  Plugin registry
# ---------------------------------------------------------------------------

_plugins: dict[str, SolverPlugin] = {}


def register(plugin: SolverPlugin, *aliases: str) -> None:
    """Register a solver plugin under its name and any extra *aliases*.

    Parameters
    ----------
    plugin : SolverPlugin
        The plugin instance to register.
    *aliases : str
        Additional names the same instance answers to (e.g. ``HWavePlugin``
        is registered under ``HWAVE`` plus ``UHFR`` / ``UHFK``).

    Raises
    ------
    ValueError
        If any of the names is already registered.
    """
    for name in (plugin.name, *aliases):
        if name in _plugins:
            raise ValueError(
                f"Solver plugin {name!r} is already registered"
            )
        _plugins[name] = plugin


def get_plugin(name: str) -> SolverPlugin:
    """Retrieve a registered solver plugin by name.

    Parameters
    ----------
    name : str
        The solver name (e.g. ``"HPhi"``).

    Returns
    -------
    SolverPlugin
        The registered plugin instance.

    Raises
    ------
    KeyError
        If no plugin with that name is registered.
    """
    if name not in _plugins:
        _discover_plugins()
    if name not in _plugins:
        raise KeyError(
            f"No solver plugin registered for {name!r}. "
            f"Available: {', '.join(_plugins) or '(none)'}"
        )
    return _plugins[name]


def _discover_plugins() -> None:
    """Load built-in plugins by importing the solvers package.

    This function is called lazily on first access to ensure plugins
    are registered even if the solvers package hasn't been imported yet.
    """
    try:
        import stdface.solvers  # noqa: F401 — triggers auto-registration
    except ImportError:
        # Solvers package not available (e.g., minimal install, missing
        # dependencies, or corrupted package). Plugin lookup will fail
        # with a helpful KeyError listing available plugins.
        pass


# ---------------------------------------------------------------------------
#  Solver config registry (C3: per-solver field containers)
# ---------------------------------------------------------------------------
#
# Maps a *raw* solver name (incl. the pre-resolution ``HWAVE`` alias, which
# shares ``UHFRPlugin``/``UHFKPlugin``'s config) to a zero-arg factory.  This
# is separate from the plugin registry because the config must be attached to
# StdIntList *before* solver-name resolution (``_reset_vals`` runs before
# ``_resolve_solver_name``), when ``HWAVE`` has no registered plugin yet.
# Solver modules populate it on import via :func:`register_config`.

_config_factories: dict[str, object] = {}


def register_config(name: str, factory) -> None:
    """Register a zero-arg config factory for a (raw) solver name."""
    _config_factories[name] = factory


def get_config_factory(name: str):
    """Return the config factory for *name*, or ``None`` if none registered."""
    if name not in _config_factories:
        _discover_plugins()
    return _config_factories.get(name)
