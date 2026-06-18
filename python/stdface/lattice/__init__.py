"""Lattice construction subpackage.

This package contains modules for building lattice geometries and their
associated interaction terms.  Each lattice type is implemented in its own
module; shared utilities (site initialisation, geometry output, interaction
building, input parameter helpers) are also included.

The :class:`LatticePlugin` ABC and plugin registry provide a uniform
interface for lattice discovery and dispatch.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..core.stdface_vals import StdIntList
    from .site_util import GnuplotData


class LatticePlugin(ABC):
    """Abstract base class for lattice plugins.

    Each lattice plugin declares its geometry metadata and provides
    a :meth:`setup` method that builds the lattice Hamiltonian terms.
    Lattices that support HPhi Boost mode override :meth:`boost`.

    Attributes
    ----------
    name : str
        Canonical lattice name (e.g. ``"chain"``).
    aliases : list[str]
        All recognised name aliases (including the canonical name).
    ndim : int
        Number of spatial dimensions (1, 2, or 3).
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Canonical lattice name."""

    @property
    @abstractmethod
    def aliases(self) -> list[str]:
        """All recognised name aliases for this lattice."""

    @property
    @abstractmethod
    def ndim(self) -> int:
        """Number of spatial dimensions (1, 2, or 3)."""

    @abstractmethod
    def setup(self, StdI: StdIntList) -> "GnuplotData | None":
        """Build the lattice geometry and Hamiltonian terms.

        Parameters
        ----------
        StdI : StdIntList
            The parameter structure (modified in place).

        Returns
        -------
        GnuplotData or None
            The ``lattice.gp`` data for 2-D lattices (``None`` for 3-D
            lattices or when gnuplot output is suppressed).
        """

    def boost(self, StdI: StdIntList) -> None:
        """Build HPhi Boost-mode data for this lattice.

        The default implementation does nothing.  Override in lattices
        that support Boost (chain, honeycomb, kagome, ladder).

        Parameters
        ----------
        StdI : StdIntList
            The parameter structure (modified in place).
        """


# ---------------------------------------------------------------------------
#  Lattice plugin registry
# ---------------------------------------------------------------------------

_lattice_plugins: dict[str, LatticePlugin] = {}
"""Maps every alias to its :class:`LatticePlugin` instance."""


def register_lattice(plugin: LatticePlugin) -> None:
    """Register a lattice plugin under all its aliases.

    Parameters
    ----------
    plugin : LatticePlugin
        The plugin instance to register.

    Raises
    ------
    ValueError
        If any alias is already registered to a *different* plugin.
    """
    for alias in plugin.aliases:
        existing = _lattice_plugins.get(alias)
        if existing is not None and existing is not plugin:
            raise ValueError(
                f"Lattice alias {alias!r} is already registered to "
                f"{existing.name!r}, cannot register {plugin.name!r}"
            )
        _lattice_plugins[alias] = plugin


def get_lattice(name: str) -> LatticePlugin:
    """Retrieve a registered lattice plugin by alias.

    Parameters
    ----------
    name : str
        A lattice name alias (e.g. ``"chain"``, ``"squarelattice"``).

    Returns
    -------
    LatticePlugin
        The registered plugin instance.

    Raises
    ------
    KeyError
        If no plugin is registered under *name*.
    """
    if name not in _lattice_plugins:
        _discover_lattices()
    if name not in _lattice_plugins:
        raise KeyError(
            f"No lattice plugin registered for {name!r}. "
            f"Available: {', '.join(sorted(set(p.name for p in _lattice_plugins.values()))) or '(none)'}"
        )
    return _lattice_plugins[name]


def get_all_lattices() -> list[LatticePlugin]:
    """Return all unique registered lattice plugins.

    Returns
    -------
    list[LatticePlugin]
        Unique plugin instances (deduplicated by identity).
    """
    if not _lattice_plugins:
        _discover_lattices()
    seen: set[int] = set()
    result: list[LatticePlugin] = []
    for plugin in _lattice_plugins.values():
        pid = id(plugin)
        if pid not in seen:
            seen.add(pid)
            result.append(plugin)
    return result


def _discover_lattices() -> None:
    """Import all lattice modules to trigger auto-registration."""
    from . import chain_lattice, square_lattice, ladder, triangular_lattice  # noqa: F401
    from . import honeycomb_lattice, kagome, orthorhombic, fc_ortho, pyrochlore  # noqa: F401
    from . import wannier90  # noqa: F401
