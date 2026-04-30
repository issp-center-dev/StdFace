"""Unit tests for all lattice modules.

Shared tests that verify each translated lattice module can be imported
and that the expected lattice functions exist and are callable.
"""
from __future__ import annotations

import importlib

import pytest


# ---------------------------------------------------------------------------
#  Module / function pairs to test
# ---------------------------------------------------------------------------

# Each entry is (module_name, list_of_expected_functions).
LATTICE_MODULES = [
    ("stdface.lattice.square_lattice", ["tetragonal"]),
    ("stdface.lattice.ladder", ["ladder", "ladder_boost"]),
    ("stdface.lattice.triangular_lattice", ["triangular"]),
    ("stdface.lattice.honeycomb_lattice", ["honeycomb", "honeycomb_boost"]),
    ("stdface.lattice.kagome", ["kagome", "kagome_boost"]),
    ("stdface.lattice.orthorhombic", ["orthorhombic"]),
    ("stdface.lattice.fc_ortho", ["fc_ortho"]),
    ("stdface.lattice.pyrochlore", ["pyrochlore"]),
    ("stdface.lattice.chain_lattice", ["chain", "chain_boost"]),
]

# Flat list of (module_name, function_name) for parametrized tests.
_MODULE_FUNC_PAIRS = [
    (mod, fn)
    for mod, fns in LATTICE_MODULES
    for fn in fns
]


# ---------------------------------------------------------------------------
#  Tests: modules import without error
# ---------------------------------------------------------------------------


class TestLatticeImports:
    """Verify that every lattice module can be imported."""

    @pytest.mark.parametrize(
        "module_name",
        [mod for mod, _ in LATTICE_MODULES],
        ids=[mod for mod, _ in LATTICE_MODULES],
    )
    def test_import(self, module_name: str):
        """Each lattice module should import without raising."""
        mod = importlib.import_module(module_name)
        assert mod is not None


# ---------------------------------------------------------------------------
#  Tests: expected functions exist and are callable
# ---------------------------------------------------------------------------


class TestLatticeFunctionsExist:
    """Verify that every expected lattice function exists and is callable."""

    @pytest.mark.parametrize(
        "module_name,func_name",
        _MODULE_FUNC_PAIRS,
        ids=[f"{m}.{f}" for m, f in _MODULE_FUNC_PAIRS],
    )
    def test_function_exists(self, module_name: str, func_name: str):
        """The function should be an attribute of the module."""
        mod = importlib.import_module(module_name)
        assert hasattr(mod, func_name), (
            f"{module_name} is missing function '{func_name}'"
        )

    @pytest.mark.parametrize(
        "module_name,func_name",
        _MODULE_FUNC_PAIRS,
        ids=[f"{m}.{f}" for m, f in _MODULE_FUNC_PAIRS],
    )
    def test_function_callable(self, module_name: str, func_name: str):
        """The function should be callable."""
        mod = importlib.import_module(module_name)
        fn = getattr(mod, func_name)
        assert callable(fn), (
            f"{module_name}.{func_name} is not callable"
        )


# ---------------------------------------------------------------------------
#  Tests: modules expose StdIntList via their imports
# ---------------------------------------------------------------------------


class TestLatticeModuleContents:
    """Additional sanity checks on module contents."""

    @pytest.mark.parametrize(
        "module_name",
        [mod for mod, _ in LATTICE_MODULES],
        ids=[mod for mod, _ in LATTICE_MODULES],
    )
    def test_module_has_docstring(self, module_name: str):
        """Each lattice module should have a module-level docstring."""
        mod = importlib.import_module(module_name)
        assert mod.__doc__ is not None, (
            f"{module_name} is missing a module docstring"
        )

    @pytest.mark.parametrize(
        "module_name,func_name",
        _MODULE_FUNC_PAIRS,
        ids=[f"{m}.{f}" for m, f in _MODULE_FUNC_PAIRS],
    )
    def test_function_has_docstring(self, module_name: str, func_name: str):
        """Each lattice function should have a docstring."""
        mod = importlib.import_module(module_name)
        fn = getattr(mod, func_name)
        assert fn.__doc__ is not None, (
            f"{module_name}.{func_name} is missing a docstring"
        )

    @pytest.mark.parametrize(
        "module_name,func_name",
        _MODULE_FUNC_PAIRS,
        ids=[f"{m}.{f}" for m, f in _MODULE_FUNC_PAIRS],
    )
    def test_function_accepts_StdI_argument(self, module_name: str, func_name: str):
        """Each lattice function should accept at least one argument (StdI)."""
        import inspect
        mod = importlib.import_module(module_name)
        fn = getattr(mod, func_name)
        sig = inspect.signature(fn)
        assert len(sig.parameters) >= 1, (
            f"{module_name}.{func_name} should accept at least one parameter (StdI)"
        )
