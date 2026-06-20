"""C3-1: dynamic single-active-config delegation on StdIntList.

These tests exercise the delegation machinery (``__getattr__`` /
``__setattr__`` + config registry) in isolation, using a throwaway dummy
config.  No solver fields have moved yet, so production behaviour is
unchanged; this only verifies the routing contract that C3-2..C3-5 rely on.
"""
from dataclasses import dataclass

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList, _STDI_OWN_FIELDS
from stdface import plugin


@dataclass
class _DummyConfig:
    foo: int | None = None
    bar: str | None = None


def test_no_config_missing_attr_raises():
    """With no config attached, a genuinely missing attr raises AttributeError."""
    s = StdIntList()
    assert s._solver_cfg is None
    with pytest.raises(AttributeError):
        _ = s.this_attr_does_not_exist


def test_getset_routes_to_attached_config():
    s = StdIntList()
    s._solver_cfg = _DummyConfig()
    # name owned by the config (not an StdIntList field) -> routed
    assert "foo" not in _STDI_OWN_FIELDS
    s.foo = 42
    assert s.foo == 42
    assert s._solver_cfg.foo == 42
    # and the value really lives on the config, not as an StdI instance attr
    assert "foo" not in s.__dict__


def test_own_fields_never_route_to_config():
    """A C1 façade property and a flat field stay on StdIntList even with a config."""
    s = StdIntList()
    s._solver_cfg = _DummyConfig()
    # flat own field (avoid 'solver', which would auto-attach a real config)
    s.lGC = 1
    assert s.lGC == 1
    assert s._solver_cfg.bar is None  # untouched
    # C1 façade property (delegates to _lattice); in-place numpy still works
    s.box[0, 0] = 7
    assert s.box[0, 0] == 7
    assert np.array_equal(s._lattice.box, s.box)


def test_cell_map_not_routed_to_config():
    """Dynamic cache attrs the config does not declare stay on StdIntList."""
    s = StdIntList()
    s._solver_cfg = _DummyConfig()
    s._cell_map = {"x": 1}
    assert s._cell_map == {"x": 1}
    assert not hasattr(s._solver_cfg, "_cell_map")


def test_config_registry_register_and_get():
    name = "_c3_test_solver"
    plugin.register_config(name, _DummyConfig)
    try:
        factory = plugin.get_config_factory(name)
        assert factory is _DummyConfig
        assert isinstance(factory(), _DummyConfig)
        assert plugin.get_config_factory("_c3_unregistered_solver") is None
    finally:
        plugin._config_factories.pop(name, None)


# --- C3-2: UHFConfig attaches in the live flow (dormant / shadowed) ---------

def test_uhf_config_attached_by_reset_vals():
    """A UHF run attaches a UHFConfig; its fields now resolve via the config."""
    from stdface.core.stdface_main import _reset_vals
    from stdface.solvers.uhf.config import UHFConfig

    s = StdIntList()
    s.solver = "UHF"
    _reset_vals(s)
    assert isinstance(s._solver_cfg, UHFConfig)
    # mix moved off StdIntList (C3-5): writes and reads route to the config.
    s.mix = 0.25
    assert s.mix == 0.25
    assert s._solver_cfg.mix == 0.25
    assert "mix" not in s.__dict__


def test_unregistered_solver_attaches_no_config():
    """A solver name with no registered config keeps _solver_cfg None."""
    from stdface.core.stdface_main import _attach_solver_config

    s = StdIntList()
    s.solver = "_no_such_solver"
    _attach_solver_config(s)
    assert s._solver_cfg is None


# --- Duplication-drift guard (dev/solver_config_split.md §2-2) ---------------
#
# The sublattice / symmetry block is duplicated across MVMCConfig / UHFConfig /
# HWaveConfig on purpose (per-solver separability).  No shared constant binds
# them, so this invariant test catches any future drift in field set or
# defaults instead.

_SUBLATTICE_FIELDS = (
    "NMPTrans", "RndSeed", "Lsub", "Wsub", "Hsub", "NCellsub", "boxsub", "rboxsub",
)


def test_sublattice_block_consistent_across_configs():
    from stdface.solvers.mvmc.config import MVMCConfig
    from stdface.solvers.uhf.config import UHFConfig
    from stdface.solvers.hwave.config import HWaveConfig

    configs = [MVMCConfig(), UHFConfig(), HWaveConfig()]
    for name in _SUBLATTICE_FIELDS:
        values = [getattr(c, name) for c in configs]
        for v in values:
            assert np.array_equal(v, values[0]), f"{name} default drift: {values}"
    # boxsub / rboxsub must be (3,3) int arrays on every config
    for c in configs:
        for arr_name in ("boxsub", "rboxsub"):
            arr = getattr(c, arr_name)
            assert arr.shape == (3, 3)
            assert np.issubdtype(arr.dtype, np.integer)


def test_hphi_config_has_no_sublattice_block():
    """HPhi does not use the sublattice block, so HPhiConfig must not carry it."""
    from stdface.solvers.hphi.config import HPhiConfig

    cfg = HPhiConfig()
    for name in _SUBLATTICE_FIELDS:
        assert not hasattr(cfg, name)
