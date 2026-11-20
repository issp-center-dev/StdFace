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
