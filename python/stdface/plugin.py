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

try:
    from importlib.metadata import entry_points
except ImportError:
    # Python < 3.10
    from importlib_metadata import entry_points  # type: ignore[import-not-found]


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
    #  Template method for writing output files
    # ------------------------------------------------------------------

    def write(self, StdI: StdIntList) -> None:
        """Write all Expert-mode definition files for this solver.

        This is a template method that calls the common output steps in
        order.  Override individual steps or this entire method as needed.

        Parameters
        ----------
        StdI : StdIntList
            The fully-populated parameter structure.
        """
        self.write_locspn(StdI)
        self.write_trans(StdI)
        self.write_interactions(StdI)
        self.check_and_write_modpara(StdI)
        self.write_solver_specific(StdI)
        self.check_output_mode(StdI)
        self.write_green(StdI)
        self.write_namelist(StdI)

    # ------------------------------------------------------------------
    #  Default step implementations — override as needed
    # ------------------------------------------------------------------

    def write_locspn(self, StdI: StdIntList) -> None:
        """Write locspn.def."""
        from .writer.common_writer import print_loc_spin
        print_loc_spin(StdI)

    def write_trans(self, StdI: StdIntList) -> None:
        """Write trans.def."""
        from .writer.common_writer import print_trans
        print_trans(StdI)

    def write_interactions(self, StdI: StdIntList) -> None:
        """Write interaction definition files."""
        from .writer.interaction_writer import print_interactions
        print_interactions(StdI)

    def check_and_write_modpara(self, StdI: StdIntList) -> None:
        """Check and write modpara.def."""
        from .writer.common_writer import check_mod_para, print_mod_para
        check_mod_para(StdI)
        print_mod_para(StdI)

    def write_solver_specific(self, StdI: StdIntList) -> None:
        """Write solver-specific output files.

        Override this method to add solver-specific output (e.g. excitation
        files for HPhi, variational parameter files for mVMC).
        The default implementation does nothing.

        Parameters
        ----------
        StdI : StdIntList
            The fully-populated parameter structure.
        """

    def check_output_mode(self, StdI: StdIntList) -> None:
        """Check and validate the output mode setting."""
        from .writer.common_writer import check_output_mode
        check_output_mode(StdI)

    def write_green(self, StdI: StdIntList) -> None:
        """Write Green's function definition files.

        Default writes both greenone.def and greentwo.def.
        Override to write only a subset.

        Parameters
        ----------
        StdI : StdIntList
            The fully-populated parameter structure.
        """
        from .writer.common_writer import print_1_green, print_2_green
        print_1_green(StdI)
        print_2_green(StdI)

    def write_namelist(self, StdI: StdIntList) -> None:
        """Write namelist.def."""
        from .writer.common_writer import print_namelist
        print_namelist(StdI)

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
    """Load built-in and external plugins.

    This function:
    1. Imports the built-in stdface.solvers package to auto-register built-in plugins
    2. Discovers and loads external plugins via the "stdface.solvers" entry point

    Called lazily on first access to ensure plugins are registered even if the
    solvers package hasn't been imported yet.
    """
    # Load built-in plugins
    try:
        import stdface.solvers  # noqa: F401 — triggers auto-registration
    except ImportError:
        pass

    # Load external plugins from entry points
    eps = entry_points()
    # Handle both old (dict) and new (SelectableGroups) API
    if hasattr(eps, 'select'):
        # Python 3.10+ / importlib.metadata >= 3.6
        solver_eps = eps.select(group='stdface.solvers')
    else:
        # Python 3.9 / older importlib_metadata
        solver_eps = eps.get('stdface.solvers', [])

    for ep in solver_eps:
        try:
            plugin_class = ep.load()
            # Instantiate and register if not already registered
            # (built-in plugins are already registered via import above)
            if callable(plugin_class):
                plugin_instance = plugin_class()
                if plugin_instance.name not in _plugins:
                    register(plugin_instance)
        except Exception:
            # Silently skip plugins that fail to load
            # This prevents one broken plugin from breaking the entire system
            pass
