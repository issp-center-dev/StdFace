"""Unit tests for the spinless fermion model (HPhi, chain lattice only)."""
from __future__ import annotations

import pytest

from stdface import generate

_SPINLESS = {"model": "SpinlessFermion", "lattice": "chain", "L": 4,
             "t": 1.0, "nelec": 2, "method": "Lanczos"}


def _generate(tmp_path, data):
    generate(data, solver="HPhi", output_dir=tmp_path)
    return tmp_path


def _trans_rows(path):
    """Return trans.def rows as (i, si, j, sj, re, im)."""
    rows = []
    for ln in (path / "trans.def").read_text().splitlines()[5:]:
        i, si, j, sj, re, im = ln.split()
        rows.append((int(i), int(si), int(j), int(sj), float(re), float(im)))
    return rows


class TestSpinlessFermionOutput:
    def test_hopping_uses_a_single_spin_channel(self, tmp_path):
        rows = _trans_rows(_generate(tmp_path, _SPINLESS))
        assert len(rows) == 2 * 4
        assert all(si == 0 and sj == 0 for _, si, _, sj, _, _ in rows)
        assert {(i, j) for i, _, j, _, _, _ in rows} == {
            (i, (i + d) % 4) for i in range(4) for d in (1, -1)}
        assert all(re == 1.0 and im == 0.0 for *_, re, im in rows)

    def test_further_neighbour_hopping(self, tmp_path):
        data = {**_SPINLESS, "L": 8, "nelec": 4, "t'": 0.5, "t''": 0.25}
        rows = _trans_rows(_generate(tmp_path, data))
        assert len(rows) == 3 * 2 * 8
        by_distance = {min((i - j) % 8, (j - i) % 8): re
                       for i, _, j, _, re, _ in rows}
        assert by_distance == {1: 1.0, 2: 0.5, 3: 0.25}

    def test_no_two_body_files_and_no_local_spins(self, tmp_path):
        out = _generate(tmp_path, _SPINLESS)
        names = {p.name for p in out.iterdir()}
        assert not names & {"coulombintra.def", "coulombinter.def", "hund.def",
                            "exchange.def", "pairlift.def", "interall.def"}
        locspin = (out / "locspn.def").read_text().splitlines()[5:]
        assert [ln.split() for ln in locspin] == [[str(i), "0"] for i in range(4)]

    def test_calcmodel_and_conserved_quantities(self, tmp_path):
        out = _generate(tmp_path, _SPINLESS)
        assert "CalcModel   7" in (out / "calcmod.def").read_text().splitlines()
        modpara = (out / "modpara.def").read_text()
        assert "Ncond          2" in modpara
        assert "2Sz" not in modpara

    def test_spinless_alias_matches_canonical_name(self, tmp_path):
        a = _generate(tmp_path / "a", _SPINLESS)
        b = _generate(tmp_path / "b", {**_SPINLESS, "model": "Spinless"})
        for name in ("trans.def", "calcmod.def", "modpara.def", "namelist.def"):
            assert (a / name).read_text() == (b / name).read_text()

    def test_momentum_sector(self, tmp_path):
        out = _generate(tmp_path, {**_SPINLESS, "MomentumIndex": 1})
        assert (out / "qptransidx.def").read_text().startswith("# MomentumIndex 1\n")


class TestSpinlessFermionRejection:
    @pytest.mark.parametrize("solver", ["mVMC", "UHF", "HWAVE"])
    def test_only_hphi_accepts_the_model(self, solver):
        data = {k: v for k, v in _SPINLESS.items() if k != "method"}
        with pytest.raises(ValueError, match="Unsupported combination"):
            generate(data, solver=solver, output_dir=None)

    @pytest.mark.parametrize("lattice", ["square", "ladder", "nosuchlattice"])
    def test_only_chain_lattice(self, lattice):
        data = {**_SPINLESS, "lattice": lattice, "W": 2, "L": 2}
        with pytest.raises(ValueError, match="supports only chain lattice"):
            generate(data, solver="HPhi", output_dir=None)

    def test_gc_variant_does_not_exist(self):
        with pytest.raises(ValueError, match="Unsupported combination"):
            generate({**_SPINLESS, "model": "SpinlessFermionGC"},
                     solver="HPhi", output_dir=None)

    def test_nelec_is_required(self):
        data = {k: v for k, v in _SPINLESS.items() if k != "nelec"}
        with pytest.raises(ValueError, match="ncond is NOT specified"):
            generate(data, solver="HPhi", output_dir=None)

    @pytest.mark.parametrize("key, value", [
        ("2Sz", 0), ("2S", 1), ("U", 4.0), ("V", 1.0), ("V'", 1.0),
        ("mu", 0.5), ("h", 0.5), ("Gamma", 0.5), ("J", 1.0), ("D", 1.0),
    ])
    def test_unused_parameters_are_rejected(self, key, value):
        with pytest.raises(ValueError, match="will NOT be USED"):
            generate({**_SPINLESS, key: value}, solver="HPhi", output_dir=None)

    def test_momentum_method_restriction(self):
        data = {**_SPINLESS, "method": "TPQ", "MomentumIndex": 1}
        with pytest.raises(ValueError, match="only Lanczos and CG"):
            generate(data, solver="HPhi", output_dir=None)


class TestOtherModelsUnaffected:
    """The h/Gamma defaults moved into the per-model branches of chain()."""

    def test_hubbard_chain_field_terms_still_generated(self, tmp_path):
        data = {"model": "Hubbard", "lattice": "chain", "L": 4, "t": 1.0,
                "U": 4.0, "h": 0.3, "nelec": 4, "2Sz": 0, "method": "Lanczos"}
        rows = _trans_rows(_generate(tmp_path, data))
        onsite = [(i, si, re) for i, si, j, sj, re, _ in rows
                  if i == j and si == sj]
        assert sorted(onsite) == sorted(
            [(i, 0, 0.15) for i in range(4)] + [(i, 1, -0.15) for i in range(4)])
