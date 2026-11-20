"""Wannier90 double counting: density at an unlisted R is zero (C: a4eb46e)."""
from __future__ import annotations

import pytest

from stdface import generate


def _write(path, coordinates, values):
    lines = ["Wannier regression fixture", "1", str(len(coordinates)),
             " ".join("1" for _ in coordinates)]
    lines += [f"{r} 1 1 {v}" for r, v in zip(coordinates, values)]
    path.write_text("\n".join(lines) + "\n")


_R3 = ["-1 0 0", "0 0 0", "1 0 0"]


def _run(tmp_path, monkeypatch, name, density_r, density_values, mode):
    work = tmp_path / name
    work.mkdir()
    (work / "zvo_geom.dat").write_text("1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n")
    _write(work / "zvo_hr.dat", _R3, ["-1 0", "0 0", "-1 0"])
    _write(work / "zvo_ur.dat", _R3, ["1 0", "4 0", "1 0"])
    _write(work / "zvo_jr.dat", _R3, ["0.25 0", "0 0", "0.25 0"])
    _write(work / "zvo_dr.dat", density_r, density_values)
    monkeypatch.chdir(work)
    data = {"lattice": "wannier90", "model": "HubbardGC", "W": 3, "L": 1,
            "Height": 1, "method": "CG", "doublecounting": mode,
            "cutoff_length_U": -1, "cutoff_length_J": -1}
    out = work / "out"
    generate(data, solver="HPhi", output_dir=out)
    return {p.name: p.read_text() for p in out.iterdir() if p.suffix == ".def"}


@pytest.mark.parametrize("mode", ["hartree", "hartree_u", "full"])
def test_unlisted_density_equals_explicit_zero(tmp_path, monkeypatch, mode):
    """Listing only R=0 gives the same files as listing zeros at R=+-1."""
    sparse = _run(tmp_path, monkeypatch, "sparse", ["0 0 0"], ["0.5 0"], mode)
    dense = _run(tmp_path, monkeypatch, "dense", _R3, ["0 0", "0.5 0", "0 0"], mode)
    assert sparse == dense
    assert "trans.def" in sparse


def test_listed_density_still_enters_the_fock_term(tmp_path, monkeypatch):
    zero = _run(tmp_path, monkeypatch, "zero", _R3, ["0 0", "0.5 0", "0 0"], "full")
    listed = _run(tmp_path, monkeypatch, "listed", _R3, ["0.1 0", "0.5 0", "0.1 0"], "full")
    assert zero["trans.def"] != listed["trans.def"]
