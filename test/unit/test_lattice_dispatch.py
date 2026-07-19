"""Unit tests for the lattice plugin registry."""
from __future__ import annotations

import pytest

from stdface.lattice import get_lattice, get_all_lattices, LatticePlugin


class TestLatticeRegistry:
    """Tests for the lattice plugin registry."""

    def test_get_chain(self):
        plugin = get_lattice("chain")
        assert plugin.name == "chain"
        assert plugin.ndim == 1

    def test_get_chain_alias(self):
        assert get_lattice("chain") is get_lattice("chainlattice")

    def test_get_square(self):
        plugin = get_lattice("tetragonal")
        assert plugin.name == "tetragonal"
        assert plugin.ndim == 2

    def test_square_aliases(self):
        p = get_lattice("tetragonal")
        assert get_lattice("tetragonallattice") is p
        assert get_lattice("square") is p
        assert get_lattice("squarelattice") is p

    def test_get_triangular(self):
        plugin = get_lattice("triangular")
        assert plugin.name == "triangular"
        assert plugin.ndim == 2

    def test_get_honeycomb(self):
        plugin = get_lattice("honeycomb")
        assert plugin.name == "honeycomb"
        assert plugin.ndim == 2

    def test_get_kagome(self):
        plugin = get_lattice("kagome")
        assert plugin.name == "kagome"
        assert plugin.ndim == 2

    def test_get_ladder(self):
        plugin = get_lattice("ladder")
        assert plugin.name == "ladder"
        assert plugin.ndim == 1

    def test_get_orthorhombic(self):
        plugin = get_lattice("orthorhombic")
        assert plugin.name == "orthorhombic"
        assert plugin.ndim == 3

    def test_orthorhombic_aliases(self):
        p = get_lattice("orthorhombic")
        assert get_lattice("simpleorthorhombic") is p
        assert get_lattice("cubic") is p
        assert get_lattice("simplecubic") is p

    def test_get_fco(self):
        plugin = get_lattice("fco")
        assert plugin.name == "fco"
        assert plugin.ndim == 3

    def test_fco_aliases(self):
        p = get_lattice("fco")
        assert get_lattice("face-centeredorthorhombic") is p
        assert get_lattice("fcorthorhombic") is p
        assert get_lattice("face-centeredcubic") is p
        assert get_lattice("fccubic") is p
        assert get_lattice("fcc") is p

    def test_get_pyrochlore(self):
        plugin = get_lattice("pyrochlore")
        assert plugin.name == "pyrochlore"
        assert plugin.ndim == 3

    def test_get_wannier90(self):
        plugin = get_lattice("wannier90")
        assert plugin.name == "wannier90"
        assert plugin.ndim == 3

    def test_unknown_raises(self):
        with pytest.raises(KeyError):
            get_lattice("nosuchlattice")

    def test_all_plugins_are_lattice_plugins(self):
        for plugin in get_all_lattices():
            assert isinstance(plugin, LatticePlugin)

    def test_all_have_setup(self):
        for plugin in get_all_lattices():
            assert callable(plugin.setup)

    def test_all_have_boost(self):
        for plugin in get_all_lattices():
            assert callable(plugin.boost)

    def test_get_all_lattices_count(self):
        """There should be 10 unique lattice plugins."""
        assert len(get_all_lattices()) == 10

    def test_boost_overridden_only_by_supported_lattices(self):
        """chain / honeycomb / kagome / ladder override boost(); the rest
        inherit the no-op default (formerly encoded in BOOST_DISPATCH)."""
        overriding = {p.name for p in get_all_lattices()
                      if type(p).boost is not LatticePlugin.boost}
        assert overriding == {"chain", "honeycomb", "kagome", "ladder"}


class TestTriangularPluginSetup:
    """``TriangularPlugin.setup`` delegates to :func:`triangular`."""

    def test_setup_calls_triangular_function(self):
        from unittest.mock import patch

        from stdface.lattice import get_lattice

        StdI = object()
        with patch("stdface.lattice.triangular_lattice.triangular") as fn:
            get_lattice("triangular").setup(StdI)
        fn.assert_called_once_with(StdI)
