"""Edge-case tests for solver/lattice plugin registries.

Covers duplicate registration, lazy discovery failure, ``get_all_lattices``
with an empty cache, HPhi ``post_lattice`` KeyError, and keyword-table
properties that may be unused in normal parsing paths.
"""
from __future__ import annotations

import builtins
from unittest.mock import patch

import pytest

import stdface.plugin as plugin_mod
from stdface.core.stdface_vals import StdIntList, SolverType
from stdface.plugin import get_plugin, register
from stdface.solvers.hphi._plugin import HPhiPlugin


class TestSolverPluginRegister:
    def test_register_same_name_raises_valueerror(self):
        with pytest.raises(ValueError, match="already registered"):
            register(HPhiPlugin())


class TestSolverPluginDiscoverImportError:
    def test_discover_plugins_swallows_import_error(self):
        saved = dict(plugin_mod._plugins)
        plugin_mod._plugins.clear()
        real_import = builtins.__import__

        def _import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "stdface.solvers":
                raise ImportError("simulated missing solvers package")
            return real_import(name, globals, locals, fromlist, level)

        try:
            with patch.object(builtins, "__import__", _import):
                plugin_mod._discover_plugins()
        finally:
            plugin_mod._plugins.clear()
            plugin_mod._plugins.update(saved)

        assert SolverType.HPhi in plugin_mod._plugins


class TestHPhiPostLatticeBoostKeyError:
    def test_unknown_lattice_with_lboost_exits(self):
        pl = get_plugin(SolverType.HPhi)
        StdI = StdIntList()
        StdI.solver = SolverType.HPhi
        StdI.model = "hubbard"
        StdI.lattice = "__no_such_lattice__"
        StdI.lBoost = 1
        with patch("stdface.solvers.hphi._plugin.large_value"):
            with pytest.raises(SystemExit):
                pl.post_lattice(StdI)


class TestSolverKeywordTableProperties:
    def test_hwave_mvmc_uhf_expose_keyword_tables(self):
        for name in (SolverType.HWAVE, SolverType.mVMC, SolverType.UHF):
            kt = get_plugin(name).keyword_table
            assert isinstance(kt, dict)
            assert kt
