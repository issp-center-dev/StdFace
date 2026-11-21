"""HPhi Boost is available on four lattices only; the others must be
rejected as the C generator does, not answered with a plain SpinGC input."""
from __future__ import annotations

import pytest

from stdface import generate
from stdface.lattice import LatticePlugin, get_lattice

_SUPPORTED = {
    "chain": {"L": 8, "J0": 1.0},
    "honeycomb": {"W": 3, "L": 2, "J0": 1.0, "J1": 0.5, "J2": 0.3},
    "kagome": {"W": 3, "L": 2, "J0": 1.0, "J1": 0.5, "J2": 0.3},
    "ladder": {"L": 8, "W": 2, "J0": 1.0, "J1": 1.0},
}
_UNSUPPORTED = {
    "square": {"W": 2, "L": 2, "J": 1.0},
    "triangular": {"W": 2, "L": 2, "J": 1.0},
    "orthorhombic": {"W": 2, "L": 2, "Height": 2, "J": 1.0},
    "pyrochlore": {"W": 1, "L": 1, "Height": 1, "J": 1.0},
    "fco": {"W": 2, "L": 2, "Height": 2, "J": 1.0},
}
_BOOST_MODELS = ("SpinGCBoost", "SpinGCCMA")


def _run(tmp_path, lattice, params, model):
    data = {"lattice": lattice, "model": model, "method": "Lanczos", **params}
    out = tmp_path / "out"
    generate(data, solver="HPhi", output_dir=out)
    return {p.name for p in out.iterdir()}


@pytest.mark.parametrize("model", _BOOST_MODELS)
@pytest.mark.parametrize("lattice", sorted(_UNSUPPORTED))
def test_boost_is_rejected_on_unsupported_lattice(tmp_path, lattice, model):
    with pytest.raises(ValueError, match=f"MODEL = spin, LATTICE = {lattice}"):
        _run(tmp_path, lattice, _UNSUPPORTED[lattice], model)
    out = tmp_path / "out"
    assert not out.exists() or not any(out.glob("*.def"))


@pytest.mark.parametrize("lattice", sorted(_UNSUPPORTED))
def test_plain_spingc_is_still_accepted(tmp_path, lattice):
    assert "modpara.def" in _run(tmp_path, lattice, _UNSUPPORTED[lattice], "SpinGC")


@pytest.mark.parametrize("model", _BOOST_MODELS)
@pytest.mark.parametrize("lattice", sorted(_SUPPORTED))
def test_boost_is_generated_on_supported_lattice(tmp_path, lattice, model):
    assert "boost.def" in _run(tmp_path, lattice, _SUPPORTED[lattice], model)


@pytest.mark.parametrize("lattice", sorted({**_SUPPORTED, **_UNSUPPORTED}) + ["wannier90"])
def test_support_flag_matches_the_boost_override(lattice):
    """A lattice declares Boost support exactly when it implements boost()."""
    plugin = get_lattice(lattice)
    overrides = type(plugin).boost is not LatticePlugin.boost
    assert plugin.supports_boost == overrides == (lattice in _SUPPORTED)
