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
    provided by :meth:`write_modpara_body` and :meth:`write_namelist_body`;
    :meth:`has_two_body_green` controls whether ``greentwo.def`` is listed.

    These methods replace the former ``_MODPARA_BODY_DISPATCH`` /
    ``_NAMELIST_BODY_DISPATCH`` dispatch tables.
    """

    @abstractmethod
    def write_modpara_body(self, fp, StdI: StdIntList) -> None:
        """Write the solver-specific body of ``modpara.def``."""

    def write_namelist_body(self, fp, StdI: StdIntList) -> None:
        """Write solver-specific extra entries in ``namelist.def``.

        Default: no extra entries (used by UHF).
        """

    def has_two_body_green(self, StdI: StdIntList) -> bool:
        """Whether ``greentwo.def`` is listed in ``namelist.def`` (default True)."""
        return True

    def _write_common_files(self, StdI: StdIntList) -> None:
        """Write the files common to all Expert-mode solvers.

        Replaces the former ``write()`` template + per-step wrapper methods.
        ``namelist.def`` is intentionally **not** written here: each solver
        writes its solver-specific files first (HPhi's excitation/calcmod set
        ``SpectrumBody`` / ``PumpBody`` that ``print_namelist`` reads) and
        calls ``print_namelist`` last.
        """
        from .writer.common_writer import (
            print_loc_spin, print_trans, check_mod_para, print_mod_para,
            check_output_mode, print_1_green, print_2_green,
        )
        from .writer.interaction_writer import print_interactions

        print_loc_spin(StdI)
        print_trans(StdI)
        print_interactions(StdI)
        check_mod_para(StdI)
        print_mod_para(StdI)
        check_output_mode(StdI)
        print_1_green(StdI)
        if self.has_two_body_green(StdI):
            print_2_green(StdI)


class WannierModeSolverPlugin(SolverPlugin):
    """Base class for solvers that emit Wannier90-format files (UHFK / RPA).

    Writes ``geometry``/``transfer``/``coulombintra`` files via the
    Wannier90 writer instead of the Expert-mode ``.def`` files.
    """

    def write(self, StdI: StdIntList) -> None:
        from .writer.wannier90_writer import export_geometry, export_interaction
        export_geometry(StdI)
        export_interaction(StdI)
        self.write_wannier_extras(StdI)

    def write_wannier_extras(self, StdI: StdIntList) -> None:
        """Write additional Wannier-mode files (default: none)."""


# ---------------------------------------------------------------------------
#  Plugin registry
# ---------------------------------------------------------------------------

_plugins: dict[str, SolverPlugin] = {}


def register(plugin: SolverPlugin) -> None:
    """Register a solver plugin.

    Parameters
    ----------
    plugin : SolverPlugin
        The plugin instance to register.

    Raises
    ------
    ValueError
        If a plugin with the same name is already registered.
    """
    if plugin.name in _plugins:
        raise ValueError(
            f"Solver plugin {plugin.name!r} is already registered"
        )
    _plugins[plugin.name] = plugin


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
