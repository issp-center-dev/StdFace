"""Unit tests for the HPhi MomentumIndex support (qptransidx.def / TransSym)."""
from __future__ import annotations

import pytest

from stdface import generate

_SPIN_XY = {"model": "Spin", "lattice": "chain", "L": 4,
            "Jx": 1.0, "Jy": 1.0, "2Sz": 0, "method": "CG"}
_HUBBARD = {"model": "Hubbard", "lattice": "chain", "L": 4,
            "t": 1.0, "U": 4.0, "nelec": 4, "2Sz": 0, "method": "Lanczos"}


def _generate(tmp_path, data):
    generate(data, solver="HPhi", output_dir=tmp_path)
    return tmp_path


class TestTransSymOutput:
    def test_not_generated_without_momentum_index(self, tmp_path):
        out = _generate(tmp_path, _SPIN_XY)
        assert not (out / "qptransidx.def").exists()
        assert "TransSym" not in (out / "namelist.def").read_text()

    def test_qptransidx_content(self, tmp_path):
        out = _generate(tmp_path, {**_SPIN_XY, "MomentumIndex": 1})
        lines = (out / "qptransidx.def").read_text().splitlines()
        sep = "=" * 45
        assert lines[:6] == [
            "# MomentumIndex 1",
            sep,
            "NQPTrans          4",
            sep,
            "======== TrIdx_TrWeight_and_TrIdx_i_xi ======",
            sep,
        ]
        assert lines[6:10] == [
            "0         1.000000000000000         0.000000000000000",
            "1         0.000000000000000        -1.000000000000000",
            "2        -1.000000000000000         0.000000000000000",
            "3         0.000000000000000         1.000000000000000",
        ]
        perm = [tuple(int(x) for x in ln.split()) for ln in lines[10:]]
        assert perm == [(op, site, (site + op) % 4, 1)
                        for op in range(4) for site in range(4)]
        assert lines[10] == "    0      0      0      1"

    @pytest.mark.parametrize("L, k", [(4, 1), (4, 2), (4, 3), (6, 3), (8, 2)])
    def test_weights_have_no_negative_zero(self, tmp_path, L, k):
        out = _generate(tmp_path, {**_SPIN_XY, "L": L, "MomentumIndex": k})
        weights = (out / "qptransidx.def").read_text().splitlines()[6:6 + L]
        assert not any("-0.000000000000000" in ln for ln in weights)

    def test_namelist_lists_transsym_before_spectrumvec(self, tmp_path):
        out = _generate(tmp_path, {**_HUBBARD, "MomentumIndex": 2})
        keys = [ln.split()[0]
                for ln in (out / "namelist.def").read_text().splitlines()]
        assert "        TransSym  qptransidx.def" in (
            out / "namelist.def").read_text().splitlines()
        assert keys.index("TransSym") + 1 == keys.index("SpectrumVec")

    def test_chainlattice_alias_and_zero_sector(self, tmp_path):
        data = {**_HUBBARD, "lattice": "chainlattice", "MomentumIndex": 0}
        out = _generate(tmp_path, data)
        weights = (out / "qptransidx.def").read_text().splitlines()[6:10]
        assert all(ln.split()[1:] == ["1.000000000000000", "0.000000000000000"]
                   for ln in weights)


_REJECTED = [
    ({**_SPIN_XY, "lattice": "square", "W": 2, "L": 2}, "only chain lattice"),
    ({**_SPIN_XY, "MomentumIndex": 4}, "0 <= MomentumIndex < L"),
    ({**_SPIN_XY, "MomentumIndex": -1}, "0 <= MomentumIndex < L"),
    ({**_SPIN_XY, "phase0": 30.0}, "boundary phase"),
    ({**_SPIN_XY, "method": "TPQ"}, "only Lanczos and CG"),
    ({**_SPIN_XY, "model": "SpinGC"}, "only canonical Spin model"),
    ({**_SPIN_XY, "2S": 2}, "only Spin-1/2"),
    ({**_SPIN_XY, "Jz": 1.0}, "exchange-only Spin-1/2 chain"),
    ({**_SPIN_XY, "h": 0.5}, "exchange-only Spin-1/2 chain"),
    ({**_SPIN_XY, "Jy": 0.5}, "exchange-only Spin-1/2 chain"),
    ({**_HUBBARD, "model": "HubbardGC"}, "only canonical Hubbard model"),
    ({k: v for k, v in _HUBBARD.items() if k != "2Sz"}, "requires nelec and 2Sz"),
    ({**_HUBBARD, "mu": 0.3}, "mu/h/Gamma/Gamma_y"),
    ({**_HUBBARD, "t'": 0.2}, "next-neighbor hopping"),
    ({**_HUBBARD, "V": 1.0}, "V/CoulombInter"),
    ({**_HUBBARD, "method": "FullDiag"}, "only Lanczos and CG"),
    ({"model": "Kondo", "lattice": "chain", "L": 4, "t": 1.0, "J": 1.0,
      "nelec": 4, "2Sz": 0, "method": "CG"}, "only Spin, SpinlessFermion, and Hubbard"),
]


class TestMomentumValidation:
    @pytest.mark.parametrize("data, message", _REJECTED)
    def test_unsupported_hamiltonian_is_rejected(self, data, message):
        data = {"MomentumIndex": 1, **data}
        with pytest.raises(ValueError, match=message):
            generate(data, solver="HPhi", output_dir=None)

    @pytest.mark.parametrize("data", [_SPIN_XY, _HUBBARD])
    def test_supported_hamiltonian_is_accepted(self, data):
        out = generate({**data, "MomentumIndex": 1}, solver="HPhi",
                       output_dir=None)
        assert any(f.fname == "qptransidx.def" for f in out.solver_files)

    def test_keyword_is_hphi_only(self):
        with pytest.raises(ValueError, match="Unsupported Keyword"):
            generate({"model": "Hubbard", "lattice": "chain", "L": 4, "t": 1.0,
                      "U": 4.0, "nelec": 4, "MomentumIndex": 1},
                     solver="mVMC", output_dir=None)
