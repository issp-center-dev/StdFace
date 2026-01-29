"""Unit tests for lattice dispatch tables.

Tests for the ``LATTICE_DISPATCH`` and ``BOOST_DISPATCH`` dict tables
in ``stdface_main``.
"""
from __future__ import annotations

import pytest

from stdface_main import LATTICE_DISPATCH, BOOST_DISPATCH

from lattice import chain_lattice
from lattice import square_lattice
from lattice import ladder
from lattice import triangular_lattice
from lattice import honeycomb_lattice
from lattice import kagome
from lattice import orthorhombic
from lattice import fc_ortho
from lattice import pyrochlore
from lattice import wannier90 as wannier90_mod


class TestLatticeDispatch:
    """Tests for the LATTICE_DISPATCH table."""

    def test_chain_aliases(self):
        """Test that chain aliases all resolve to chain_lattice.chain."""
        assert LATTICE_DISPATCH["chain"] is chain_lattice.chain
        assert LATTICE_DISPATCH["chainlattice"] is chain_lattice.chain

    def test_square_aliases(self):
        """Test that square/tetragonal aliases resolve correctly."""
        assert LATTICE_DISPATCH["tetragonal"] is square_lattice.tetragonal
        assert LATTICE_DISPATCH["tetragonallattice"] is square_lattice.tetragonal
        assert LATTICE_DISPATCH["square"] is square_lattice.tetragonal
        assert LATTICE_DISPATCH["squarelattice"] is square_lattice.tetragonal

    def test_ladder_aliases(self):
        """Test that ladder aliases resolve correctly."""
        assert LATTICE_DISPATCH["ladder"] is ladder.ladder
        assert LATTICE_DISPATCH["ladderlattice"] is ladder.ladder

    def test_triangular_aliases(self):
        """Test that triangular aliases resolve correctly."""
        assert LATTICE_DISPATCH["triangular"] is triangular_lattice.triangular
        assert LATTICE_DISPATCH["triangularlattice"] is triangular_lattice.triangular

    def test_honeycomb_aliases(self):
        """Test that honeycomb aliases resolve correctly."""
        assert LATTICE_DISPATCH["honeycomb"] is honeycomb_lattice.honeycomb
        assert LATTICE_DISPATCH["honeycomblattice"] is honeycomb_lattice.honeycomb

    def test_kagome_aliases(self):
        """Test that kagome aliases resolve correctly."""
        assert LATTICE_DISPATCH["kagome"] is kagome.kagome
        assert LATTICE_DISPATCH["kagomelattice"] is kagome.kagome

    def test_orthorhombic_aliases(self):
        """Test that orthorhombic/cubic aliases resolve correctly."""
        assert LATTICE_DISPATCH["orthorhombic"] is orthorhombic.orthorhombic
        assert LATTICE_DISPATCH["simpleorthorhombic"] is orthorhombic.orthorhombic
        assert LATTICE_DISPATCH["cubic"] is orthorhombic.orthorhombic
        assert LATTICE_DISPATCH["simplecubic"] is orthorhombic.orthorhombic

    def test_fco_aliases(self):
        """Test that face-centered orthorhombic aliases resolve correctly."""
        assert LATTICE_DISPATCH["face-centeredorthorhombic"] is fc_ortho.fc_ortho
        assert LATTICE_DISPATCH["fcorthorhombic"] is fc_ortho.fc_ortho
        assert LATTICE_DISPATCH["fco"] is fc_ortho.fc_ortho
        assert LATTICE_DISPATCH["face-centeredcubic"] is fc_ortho.fc_ortho
        assert LATTICE_DISPATCH["fccubic"] is fc_ortho.fc_ortho
        assert LATTICE_DISPATCH["fcc"] is fc_ortho.fc_ortho

    def test_pyrochlore(self):
        """Test that pyrochlore resolves correctly."""
        assert LATTICE_DISPATCH["pyrochlore"] is pyrochlore.pyrochlore

    def test_wannier90(self):
        """Test that wannier90 resolves correctly."""
        assert LATTICE_DISPATCH["wannier90"] is wannier90_mod.wannier90

    def test_unknown_returns_none(self):
        """Test that unknown lattice returns None via .get()."""
        assert LATTICE_DISPATCH.get("nosuchlattice") is None

    def test_all_entries_callable(self):
        """Test that every entry in the dispatch table is callable."""
        for name, func in LATTICE_DISPATCH.items():
            assert callable(func), f"LATTICE_DISPATCH[{name!r}] is not callable"


class TestBoostDispatch:
    """Tests for the BOOST_DISPATCH table."""

    def test_chain_boost(self):
        """Test chain boost alias."""
        assert BOOST_DISPATCH["chain"] is chain_lattice.chain_boost
        assert BOOST_DISPATCH["chainlattice"] is chain_lattice.chain_boost

    def test_honeycomb_boost(self):
        """Test honeycomb boost alias."""
        assert BOOST_DISPATCH["honeycomb"] is honeycomb_lattice.honeycomb_boost
        assert BOOST_DISPATCH["honeycomblattice"] is honeycomb_lattice.honeycomb_boost

    def test_kagome_boost(self):
        """Test kagome boost alias."""
        assert BOOST_DISPATCH["kagome"] is kagome.kagome_boost
        assert BOOST_DISPATCH["kagomelattice"] is kagome.kagome_boost

    def test_ladder_boost(self):
        """Test ladder boost alias."""
        assert BOOST_DISPATCH["ladder"] is ladder.ladder_boost
        assert BOOST_DISPATCH["ladderlattice"] is ladder.ladder_boost

    def test_unsupported_lattice_not_in_boost(self):
        """Test that lattices without boost are not in BOOST_DISPATCH."""
        assert "square" not in BOOST_DISPATCH
        assert "triangular" not in BOOST_DISPATCH
        assert "pyrochlore" not in BOOST_DISPATCH
        assert "wannier90" not in BOOST_DISPATCH

    def test_all_entries_callable(self):
        """Test that every entry in the boost table is callable."""
        for name, func in BOOST_DISPATCH.items():
            assert callable(func), f"BOOST_DISPATCH[{name!r}] is not callable"
