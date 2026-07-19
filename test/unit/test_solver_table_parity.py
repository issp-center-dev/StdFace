"""Cross-solver parity tests for intentionally duplicated definitions.

Policy: each ``solvers/<name>/`` package is self-contained so that it can
eventually be split out as an independent package.  Overlapping keyword
tables, reset tables and config fields are therefore *deliberate copies*
(the C3 duplication policy, ``dev/solver_config_split.md`` §2-2), not
shared modules.  While the solvers co-reside in this repository these
tests keep the copies in sync — a missing entry on one side (like the
``fileprefix`` reset gap fixed in PB-4) is caught here instead of
surfacing as a runtime leak.

When a solver is extracted into its own package, deleting the parity
assertions that involve it is the deliberate act of making it
independent.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from stdface.core.keyword_parser import _grid3x3_keywords
from stdface.plugin import get_plugin
from stdface.solvers.hwave.config import HWaveConfig
from stdface.solvers.mvmc.config import MVMCConfig
from stdface.solvers.uhf.config import UHFConfig


def _dummy_store(keyword, value, current):
    raise AssertionError("never called")


#: Keys of the sublattice / symmetry keyword block shared by mVMC / UHF /
#: H-wave: the scalar keys plus the 3x3 boxsub grid (a0wsub, ...).
SUBLATTICE_KEYS = frozenset(
    {"nmptrans", "rndseed", "hsub", "lsub", "wsub"}
    | set(_grid3x3_keywords("{a}{c}sub", "boxsub", _dummy_store, int))
)

#: Reset-scalar names of the sublattice block (all reset to None).
SUBLATTICE_RESET_NAMES = frozenset(
    {"NMPTrans", "RndSeed", "Hsub", "Lsub", "Wsub"})

#: Config fields of the sublattice block shared by all three configs.
SUBLATTICE_CONFIG_FIELDS = frozenset(
    {"NMPTrans", "RndSeed", "Lsub", "Wsub", "Hsub",
     "NCellsub", "boxsub", "rboxsub"})

#: Config fields of the UHF iteration block shared by UHF and H-wave.
UHF_BASE_CONFIG_FIELDS = frozenset(
    {"mix", "eps", "eps_slater", "Iteration_max"})


def _plugins():
    return {name: get_plugin(name) for name in ("mVMC", "UHF", "HWAVE")}


class TestKeywordTableParity:
    def test_uhf_table_is_contained_in_hwave(self):
        """H-wave accepts the whole UHF keyword surface, entry-identical.

        A keyword added to UHF but not H-wave (or diverging in store
        function / target field) fails here.
        """
        uhf = get_plugin("UHF").keyword_table
        hwave = get_plugin("HWAVE").keyword_table
        for key, entry in uhf.items():
            assert key in hwave, f"keyword {key!r} missing from H-wave"
            assert hwave[key] == entry, f"keyword {key!r} diverged"

    def test_sublattice_block_identical_across_three(self):
        tables = {n: p.keyword_table for n, p in _plugins().items()}
        for key in sorted(SUBLATTICE_KEYS):
            entries = {n: t.get(key) for n, t in tables.items()}
            assert None not in entries.values(), (
                f"sublattice keyword {key!r} missing: {entries}")
            assert len(set(entries.values())) == 1, (
                f"sublattice keyword {key!r} diverged: {entries}")


class TestResetTableParity:
    def test_uhf_reset_scalars_contained_in_hwave(self):
        uhf = dict(get_plugin("UHF").reset_scalars)
        hwave = dict(get_plugin("HWAVE").reset_scalars)
        for name, value in uhf.items():
            assert name in hwave, f"reset for {name!r} missing from H-wave"
            assert hwave[name] == value, f"reset for {name!r} diverged"

    def test_sublattice_reset_scalars_across_three(self):
        for pname, plugin in _plugins().items():
            resets = dict(plugin.reset_scalars)
            for name in sorted(SUBLATTICE_RESET_NAMES):
                assert name in resets, f"{pname}: no reset for {name!r}"
                assert resets[name] is None, f"{pname}: {name!r} reset diverged"

    def test_boxsub_reset_array_across_three(self):
        from stdface.core.stdface_vals import NaN_i
        for pname, plugin in _plugins().items():
            arrays = dict(plugin.reset_arrays)
            assert arrays.get("boxsub") == NaN_i, (
                f"{pname}: boxsub reset missing or diverged")

    def test_reset_scalars_cover_all_config_scalars(self):
        """Every scalar config field has a reset entry (the PB-4 gap class).

        Array fields (covered by reset_arrays) and non-parameter state are
        exempted explicitly.
        """
        exempt = {
            "boxsub", "rboxsub",       # arrays (reset_arrays)
            "NCellsub",                # derived, not an input parameter
        }
        for solver, cfg_cls in [("mVMC", MVMCConfig), ("UHF", UHFConfig),
                                ("HWAVE", HWaveConfig)]:
            plugin = get_plugin(solver)
            reset_names = {n for n, _ in plugin.reset_scalars}
            reset_names |= {n for n, _ in plugin.reset_arrays}
            cfg_names = {f.name for f in dataclasses.fields(cfg_cls)}
            missing = cfg_names - reset_names - exempt
            # mVMC / HPhi carry computed state in their configs; only flag
            # plain parameter-like scalars (annotation "<type> | None").
            missing = {
                n for n in missing
                if "| None" in str(next(f.type for f in
                                        dataclasses.fields(cfg_cls)
                                        if f.name == n))
            }
            assert not missing, (
                f"{solver}: config fields without reset entry: {sorted(missing)}")


class TestConfigFieldParity:
    @staticmethod
    def _field_map(cfg_cls):
        return {f.name: f for f in dataclasses.fields(cfg_cls)}

    @staticmethod
    def _default_of(f):
        if f.default_factory is not dataclasses.MISSING:
            return f.default_factory()
        return f.default

    def _assert_fields_match(self, names, classes):
        maps = {cls.__name__: self._field_map(cls) for cls in classes}
        for name in sorted(names):
            fields = {}
            for cls_name, fmap in maps.items():
                assert name in fmap, f"{cls_name} lacks shared field {name!r}"
                fields[cls_name] = fmap[name]
            types = {str(f.type) for f in fields.values()}
            assert len(types) == 1, f"field {name!r} type diverged: {types}"
            defaults = [self._default_of(f) for f in fields.values()]
            first = defaults[0]
            for other in defaults[1:]:
                if isinstance(first, np.ndarray):
                    assert (other.shape == first.shape
                            and other.dtype == first.dtype
                            and np.array_equal(other, first)), (
                        f"field {name!r} array default diverged")
                else:
                    assert other == first, f"field {name!r} default diverged"

    def test_sublattice_fields_across_three_configs(self):
        self._assert_fields_match(
            SUBLATTICE_CONFIG_FIELDS, [MVMCConfig, UHFConfig, HWaveConfig])

    def test_uhf_base_fields_between_uhf_and_hwave(self):
        self._assert_fields_match(
            UHF_BASE_CONFIG_FIELDS, [UHFConfig, HWaveConfig])
