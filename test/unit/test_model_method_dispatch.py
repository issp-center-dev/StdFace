"""Unit tests for model and method normalisation dispatch tables.

Tests for the ``MODEL_ALIASES``, ``MODEL_ALIASES_HPHI_BOOST``, and
``METHOD_ALIASES`` dict tables in ``stdface_main``.
"""
from __future__ import annotations

import pytest

from stdface.core.stdface_main import MODEL_ALIASES, MODEL_ALIASES_HPHI_BOOST, METHOD_ALIASES


class TestModelAliases:
    """Tests for the MODEL_ALIASES table."""

    def test_hubbard_aliases(self):
        """Test that hubbard aliases normalise correctly."""
        assert MODEL_ALIASES["fermionhubbard"] == ("hubbard", 0, 0)
        assert MODEL_ALIASES["hubbard"] == ("hubbard", 0, 0)

    def test_hubbard_gc_aliases(self):
        """Test that hubbard GC aliases normalise correctly."""
        assert MODEL_ALIASES["fermionhubbardgc"] == ("hubbard", 1, 0)
        assert MODEL_ALIASES["hubbardgc"] == ("hubbard", 1, 0)

    def test_spin_aliases(self):
        """Test that spin aliases normalise correctly."""
        assert MODEL_ALIASES["spin"] == ("spin", 0, 0)

    def test_spin_gc_alias(self):
        """Test that spingc normalises correctly."""
        assert MODEL_ALIASES["spingc"] == ("spin", 1, 0)

    def test_kondo_aliases(self):
        """Test that kondo aliases normalise correctly."""
        assert MODEL_ALIASES["kondolattice"] == ("kondo", 0, 0)
        assert MODEL_ALIASES["kondo"] == ("kondo", 0, 0)

    def test_kondo_gc_aliases(self):
        """Test that kondo GC aliases normalise correctly."""
        assert MODEL_ALIASES["kondolatticegc"] == ("kondo", 1, 0)
        assert MODEL_ALIASES["kondogc"] == ("kondo", 1, 0)

    def test_unknown_returns_none(self):
        """Test that unknown model returns None via .get()."""
        assert MODEL_ALIASES.get("nosuchmodel") is None

    def test_all_entries_are_tuples(self):
        """Test that every entry is a 3-tuple of (str, int, int)."""
        for alias, info in MODEL_ALIASES.items():
            assert isinstance(info, tuple), f"MODEL_ALIASES[{alias!r}] is not a tuple"
            assert len(info) == 3, f"MODEL_ALIASES[{alias!r}] has length {len(info)}"
            name, lgc, lboost = info
            assert isinstance(name, str)
            assert isinstance(lgc, int)
            assert isinstance(lboost, int)

    def test_gc_flags_are_zero_or_one(self):
        """Test that lGC and lBoost are 0 or 1."""
        for alias, (_, lgc, lboost) in MODEL_ALIASES.items():
            assert lgc in (0, 1), f"MODEL_ALIASES[{alias!r}] lGC={lgc}"
            assert lboost in (0, 1), f"MODEL_ALIASES[{alias!r}] lBoost={lboost}"

    def test_boost_not_in_common_aliases(self):
        """Test that boost models are NOT in the common aliases table."""
        assert "spingcboost" not in MODEL_ALIASES
        assert "spingccma" not in MODEL_ALIASES


class TestModelAliasesHPhiBoost:
    """Tests for the MODEL_ALIASES_HPHI_BOOST table."""

    def test_spingcboost(self):
        """Test spingcboost normalises to spin with GC+Boost."""
        assert MODEL_ALIASES_HPHI_BOOST["spingcboost"] == ("spin", 1, 1)

    def test_spingccma(self):
        """Test spingccma normalises to spin with GC+Boost."""
        assert MODEL_ALIASES_HPHI_BOOST["spingccma"] == ("spin", 1, 1)

    def test_all_have_boost_flag(self):
        """Test that every entry has lBoost=1."""
        for alias, (_, _, lboost) in MODEL_ALIASES_HPHI_BOOST.items():
            assert lboost == 1, f"MODEL_ALIASES_HPHI_BOOST[{alias!r}] lBoost={lboost}"

    def test_all_have_gc_flag(self):
        """Test that every entry has lGC=1."""
        for alias, (_, lgc, _) in MODEL_ALIASES_HPHI_BOOST.items():
            assert lgc == 1, f"MODEL_ALIASES_HPHI_BOOST[{alias!r}] lGC={lgc}"

    def test_no_overlap_with_common(self):
        """Test that HPhi-boost aliases don't overlap with common aliases."""
        overlap = set(MODEL_ALIASES_HPHI_BOOST) & set(MODEL_ALIASES)
        assert overlap == set(), f"Overlap: {overlap}"


class TestMethodAliases:
    """Tests for the METHOD_ALIASES table."""

    def test_direct_to_fulldiag(self):
        """Test that 'direct' maps to 'fulldiag'."""
        assert METHOD_ALIASES["direct"] == "fulldiag"

    def test_alldiag_to_fulldiag(self):
        """Test that 'alldiag' maps to 'fulldiag'."""
        assert METHOD_ALIASES["alldiag"] == "fulldiag"

    def test_te_to_timeevolution(self):
        """Test that 'te' maps to 'timeevolution'."""
        assert METHOD_ALIASES["te"] == "timeevolution"

    def test_time_evolution_to_timeevolution(self):
        """Test that 'time-evolution' maps to 'timeevolution'."""
        assert METHOD_ALIASES["time-evolution"] == "timeevolution"

    def test_canonical_names_not_aliased(self):
        """Test that canonical method names are not in the alias table."""
        assert "fulldiag" not in METHOD_ALIASES
        assert "lanczos" not in METHOD_ALIASES
        assert "timeevolution" not in METHOD_ALIASES
        assert "cg" not in METHOD_ALIASES

    def test_unknown_returns_none(self):
        """Test that unknown method returns None via .get()."""
        assert METHOD_ALIASES.get("nosuchmethod") is None

    def test_all_values_are_strings(self):
        """Test that every value in the table is a string."""
        for alias, canonical in METHOD_ALIASES.items():
            assert isinstance(canonical, str), f"METHOD_ALIASES[{alias!r}] is not a str"
