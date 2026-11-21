"""Wannier90 spin model: the superexchange needs a usable on-site U (C: 73221b6)."""
from __future__ import annotations

import pytest

from stdface import generate

_HOP = "-0.5 0;0 0;0 0;-0.5 0"
_LOCAL = "-0.2 0;0 0;0 0;-0.2 0"
_NO_HOP = "0 0;0 0;0 0;0 0"


def _write(path, n, coordinates, blocks):
    """Write a Wannier90-format matrix file: one n*n block per R vector."""
    lines = ["Wannier regression fixture", str(n), str(len(coordinates)),
             " ".join("1" for _ in coordinates)]
    for r, block in zip(coordinates, blocks):
        values = block.split(";")
        assert len(values) == n * n
        k = 0
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                lines.append(f"{r} {i} {j} {values[k]}")
                k += 1
    path.write_text("\n".join(lines) + "\n")


def _case(tmp_path, monkeypatch, *, u0="4 0", u1="8 0", hop=_HOP, ur="0 0 0",
          write_ur=True, write_jr=True, **extra):
    (tmp_path / "zvo_geom.dat").write_text("1 0 0\n0 1 0\n0 0 1\n2\n0 0 0\n0 0 0\n")
    _write(tmp_path / "zvo_hr.dat", 2, ["-1 0 0", "0 0 0", "1 0 0"], [hop, _LOCAL, hop])
    if write_ur:
        _write(tmp_path / "zvo_ur.dat", 2, [ur], [f"{u0};0 0;0 0;{u1}"])
    if write_jr:
        _write(tmp_path / "zvo_jr.dat", 2, ["0 0 0"], [_NO_HOP])
    monkeypatch.chdir(tmp_path)
    data = {"lattice": "wannier90", "model": "SpinGC", "W": 3, "L": 1, "Height": 1,
            "method": "CG", "h": 0.125, **extra}
    generate(data, solver="HPhi", output_dir=tmp_path)
    path = tmp_path / "exchange.def"
    if not path.exists():
        return {}
    return {(int(i), int(j)): float(v)
            for i, j, v in (ln.split() for ln in path.read_text().splitlines()[5:])}


_NEEDS_U = "requires finite, nonzero on-site U"


class TestSuperexchange:
    def test_values_follow_each_orbitals_u(self, tmp_path, monkeypatch):
        ex = _case(tmp_path, monkeypatch)
        assert len(ex) == 6
        for (i, j), value in ex.items():
            assert i % 2 == j % 2
            assert value == pytest.approx(0.125 if i % 2 == 0 else 0.0625)

    def test_negative_u_is_allowed(self, tmp_path, monkeypatch):
        ex = _case(tmp_path, monkeypatch, u0="-4 0", u1="-8 0")
        assert sorted(set(ex.values())) == pytest.approx([-0.125, -0.0625])

    def test_orbital_without_hopping_needs_no_u(self, tmp_path, monkeypatch):
        ex = _case(tmp_path, monkeypatch, u1="0 0", hop="-0.5 0;0 0;0 0;0 0")
        assert len(ex) == 3 and all(i % 2 == 0 for i, _ in ex)

    def test_no_hopping_needs_no_u_at_all(self, tmp_path, monkeypatch):
        assert _case(tmp_path, monkeypatch, u0="0 0", hop=_NO_HOP) == {}

    def test_missing_jr_file_is_fine(self, tmp_path, monkeypatch):
        assert len(_case(tmp_path, monkeypatch, write_jr=False)) == 6


class TestRejectedU:
    @pytest.mark.parametrize("kwargs", [
        {"u0": "0 0"},                                   # zero on the first orbital
        {"u1": "0 0"},                                   # zero on the second
        {"u1": "0 1"},                                   # only an imaginary part
        {"u1": "inf 0"},
        {"u1": "nan 0"},
        {"u1": "0.001 0", "cutoff_u": 0.01},             # removed by the cutoff
        {"ur": "1 0 0", "cutoff_length_U": -1},          # only an off-site U is listed
        {"lambda_U": 0},                                 # every U scaled to zero
        {"write_ur": False},                             # no _ur.dat at all
    ])
    def test_unusable_on_site_u_is_an_error(self, tmp_path, monkeypatch, kwargs):
        with pytest.raises(ValueError, match=_NEEDS_U):
            _case(tmp_path, monkeypatch, **kwargs)
