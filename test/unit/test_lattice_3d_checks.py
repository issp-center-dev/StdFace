"""Behaviour checks for the 3D lattices that mirror upstream C fixes."""
from __future__ import annotations

import pytest

from stdface import generate


def _trans(path):
    """Return trans.def as {(site_i, spin_i, site_j, spin_j): complex}."""
    terms = {}
    for ln in (path / "trans.def").read_text().splitlines()[5:]:
        i, si, j, sj, re, im = ln.split()
        key = (int(i), int(si), int(j), int(sj))
        terms[key] = terms.get(key, 0) + complex(float(re), float(im))
    return terms


class TestFccKondoLocalField:
    """The magnetic field acts on the localized spins as well (C: 39bd64f)."""

    _INPUT = {"model": "KondoGC", "lattice": "fcc", "L": 1, "W": 1, "Height": 1,
              "t": 0.0, "J": 1.0, "method": "CG"}

    @pytest.mark.parametrize("cells", [1, 2])
    def test_every_site_gets_the_same_field_terms(self, tmp_path, cells):
        data = {**self._INPUT, "L": cells, "h": 0.2, "Gamma": 0.4, "Gamma_y": 0.6}
        generate(data, solver="HPhi", output_dir=tmp_path)
        terms = _trans(tmp_path)
        # localized spins are sites 0..n-1, itinerant sites n..2n-1
        for site in range(2 * cells):
            assert terms[(site, 0, site, 0)] == pytest.approx(0.1)
            assert terms[(site, 1, site, 1)] == pytest.approx(-0.1)
            assert terms[(site, 1, site, 0)] == pytest.approx(0.2 + 0.3j)
            assert terms[(site, 0, site, 1)] == pytest.approx(0.2 - 0.3j)
        assert len(terms) == 4 * 2 * cells

    def test_no_field_no_local_terms(self, tmp_path):
        generate(self._INPUT, solver="HPhi", output_dir=tmp_path)
        assert not (tmp_path / "trans.def").exists() or _trans(tmp_path) == {}


_BASE = {"L": 2, "W": 1, "Height": 1, "method": "CG"}
_MODEL_INPUT = {
    "Spin": {"model": "SpinGC", "J": 1.0},
    "Hubbard": {"model": "HubbardGC", "t": 1.0},
    "Kondo": {"model": "KondoGC", "t": 1.0, "J": 1.0},
}
_DOUBLE_PRIME_DIRECTIONAL = ["t0''", "t1''", "t2''", "V0''", "V1''", "V2''",
                             "J0''", "J1''", "J2''", "J0''xy"]


def _rejected_keys(lattice, model):
    """Keys the lattice accepts syntactically but never wires (C: fcec403, 8634d27)."""
    keys = list(_DOUBLE_PRIME_DIRECTIONAL)
    if lattice == "pyrochlore":
        keys += ["t'", "t''", "V'", "V''", "J'", "J''"]
    elif lattice == "fcc":
        keys += ["t''", "V''", "J''"]
        if model != "Spin":
            keys.append("J'")
    elif model != "Spin":
        keys += ["J'", "J''"]
    if model == "Spin":
        keys += ["V2", "V0'", "V1'", "V2'"]
        if lattice == "orthorhombic":
            keys.append("V''")
    return keys


def _active_keys(lattice, model):
    """Keys that must keep changing the Hamiltonian."""
    if model == "Spin":
        keys = ["J0'", "J1'", "J2'"]
        if lattice == "orthorhombic":
            keys += ["J'", "J''"]
        elif lattice == "fcc":
            keys.append("J'")
    else:
        keys = ["t0'", "t1'", "t2'", "V0'", "V1'", "V2'"]
        if lattice == "orthorhombic":
            keys += ["t'", "t''", "V'", "V''"]
        elif lattice == "fcc":
            keys += ["t'", "V'"]
    return keys


def _input(lattice, model, key=None, value=None):
    data = {**_BASE, "lattice": lattice, **_MODEL_INPUT[model]}
    if lattice == "pyrochlore" and key is not None:
        # the scalar fallback and a directional override are mutually exclusive
        if model == "Spin" and key in ("J0'", "J1'", "J2'"):
            del data["J"]
            data["h"] = 0.125
        elif model != "Spin" and key in ("t0'", "t1'", "t2'"):
            del data["t"]
            data["mu"] = 0.125
    if key is not None:
        data[key] = value
    return data


def _hamiltonian(tmp_path, data):
    generate(data, solver="HPhi", output_dir=tmp_path)
    names = ("trans", "interall", "coulombintra", "coulombinter", "hund",
             "exchange", "pairhop", "pairlift")
    return {n: (tmp_path / f"{n}.def").read_text()
            for n in names if (tmp_path / f"{n}.def").exists()}


_LATTICES = ["orthorhombic", "pyrochlore", "fcc"]
_MODELS = ["Spin", "Hubbard", "Kondo"]


class TestUnused3DParameters:
    @pytest.mark.parametrize("lattice", _LATTICES)
    @pytest.mark.parametrize("model", _MODELS)
    def test_unwired_parameters_are_rejected(self, lattice, model):
        for key in _rejected_keys(lattice, model):
            for value in (0.375, 0):
                with pytest.raises(ValueError, match="will NOT be USED"):
                    generate(_input(lattice, model, key, value),
                             solver="HPhi", output_dir=None)

    @pytest.mark.parametrize("lattice", _LATTICES)
    @pytest.mark.parametrize("model", _MODELS)
    def test_wired_parameters_still_change_the_hamiltonian(self, tmp_path, lattice, model):
        for n, key in enumerate(_active_keys(lattice, model)):
            zero = _hamiltonian(tmp_path / f"z{n}", _input(lattice, model, key, 0))
            nonzero = _hamiltonian(tmp_path / f"n{n}", _input(lattice, model, key, 0.375))
            assert zero != nonzero, key
