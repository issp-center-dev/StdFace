"""Unit tests for the public generate() library API (D5)."""
from __future__ import annotations

import json
import os

import pytest

from stdface import generate, ExpertModeOutput, WannierModeOutput, JSONFormat


def _write_minimal_w90_inputs(path, prefix="zvo", NsiteUC=1):
    """Write a minimal *_geom.dat / *_hr.dat pair for the wannier90 lattice."""
    with open(os.path.join(path, f"{prefix}_geom.dat"), "w") as f:
        f.write("  1.0  0.0  0.0\n  0.0  1.0  0.0\n  0.0  0.0  1.0\n")
        f.write(f"  {NsiteUC}\n")
        for i in range(NsiteUC):
            f.write("  0.0  0.0  0.0\n")
    with open(os.path.join(path, f"{prefix}_hr.dat"), "w") as f:
        f.write("  written by test\n")
        f.write(f"  {NsiteUC}\n  1\n  1\n")
        for i in range(NsiteUC):
            for j in range(NsiteUC):
                re_val = -1.0 if i == j else -0.5
                f.write(f"  0  0  0  {i + 1}  {j + 1}  {re_val:12.6f}  0.000000\n")
from stdface.core.stdface_main import _make_source, _check_solver_conflict
from stdface.core.input_source import (
    DictSource, StanFileSource, TOMLSource, JSONSource,
)

_HUB = {"model": "Hubbard", "lattice": "Tetragonal",
        "W": 2, "L": 2, "t": 1.0, "U": 4.0, "nelec": 4}
_HPHI = {**_HUB, "method": "Lanczos"}


class TestMakeSource:
    def test_dict(self):
        assert isinstance(_make_source({"a": 1}), DictSource)

    def test_suffixes(self, tmp_path):
        assert isinstance(_make_source(tmp_path / "x.toml"), TOMLSource)
        assert isinstance(_make_source(tmp_path / "x.json"), JSONSource)
        assert isinstance(_make_source(tmp_path / "x.in"), StanFileSource)


class TestSolverConflict:
    def test_no_conflict_when_absent(self):
        _check_solver_conflict("HPhi", {"model": "Hubbard"})  # no raise

    def test_match_is_ok_case_insensitive(self):
        _check_solver_conflict("HPhi", {"solver": "hphi"})  # no raise

    def test_conflict_raises(self):
        with pytest.raises(ValueError):
            _check_solver_conflict("HPhi", {"solver": "mVMC"})


class TestGenerate:
    def test_returns_expert_output(self):
        out = generate(_HPHI, solver="HPhi", output_dir=None)
        assert isinstance(out, ExpertModeOutput)
        assert out.to_dict()["modpara"]["params"]["Nsite"] == 4

    def test_output_dir_none_is_file_free(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        generate(_HPHI, solver="HPhi", output_dir=None)
        assert os.listdir(".") == []

    def test_output_dir_collects_all_files(self, tmp_path):
        generate(_HPHI, solver="HPhi", output_dir=tmp_path / "o")
        fs = set(os.listdir(tmp_path / "o"))
        # data-backed + lattice-level + HPhi solver-specific all land together
        assert {"modpara.def", "namelist.def", "trans.def", "greenone.def",
                "geometry.dat", "lattice.gp", "calcmod.def"} <= fs

    def test_uhfr_partial_container(self, tmp_path):
        out = generate(_HUB, solver="UHFR", output_dir=tmp_path / "r")
        assert isinstance(out, ExpertModeOutput)
        assert out.locspn is None and out.modpara is None and out.namelist is None
        fs = set(os.listdir(tmp_path / "r"))
        assert "trans.def" in fs and "modpara.def" not in fs

    def test_uhfk_returns_wannier_output(self, tmp_path):
        # UHFK is reached via solver=HWAVE + calcmode (pre-resolution path).
        s = {"model": "Hubbard", "lattice": "Tetragonal", "W": 2, "L": 2,
             "calcmode": "uhfk"}
        out = generate(s, solver="HWAVE", output_dir=tmp_path / "k")
        assert isinstance(out, WannierModeOutput)
        assert (tmp_path / "k" / "geom.dat").exists()
        assert not (tmp_path / "k" / "geometry.dat").exists()  # suppressed

    def test_wannier90_inputs_resolved_against_caller_cwd(self, tmp_path,
                                                          monkeypatch):
        """Wannier90 data files come from the caller's cwd, not output_dir (#P8)."""
        monkeypatch.chdir(tmp_path)
        _write_minimal_w90_inputs(tmp_path, prefix="zvo")

        s = {"model": "Hubbard", "lattice": "wannier90",
             "W": 2, "L": 2, "Height": 1, "calcmode": "uhfr",
             "nelec": 2, "2Sz": 0}
        out = generate(s, solver="HWAVE", output_dir=tmp_path / "out")
        assert isinstance(out, ExpertModeOutput)
        assert (tmp_path / "out" / "trans.def").exists()

    def test_json_format(self, tmp_path):
        generate(_HPHI, solver="HPhi", output_dir=tmp_path / "j",
                 output_format=JSONFormat())
        doc = json.loads((tmp_path / "j" / "stdface_output.json").read_text())
        assert doc["modpara"]["params"]["Nsite"] == 4

    def test_conflict_propagates(self):
        with pytest.raises(ValueError):
            generate({**_HPHI, "solver": "mVMC"}, solver="HPhi", output_dir=None)


class TestUHFkPhaseGuard:
    """UHFk / RPA assume translational symmetry, so a boundary twist
    (``phase``) must be rejected; the real-space UHFR mode is unaffected."""

    _UHFK = {**_HUB, "calcmode": "uhfk"}

    def test_uhfk_without_phase_ok(self):
        out = generate(self._UHFK, solver="HWAVE", output_dir=None)
        assert isinstance(out, WannierModeOutput)

    def test_uhfk_with_phase_raises(self):
        with pytest.raises(ValueError, match="phase"):
            generate({**self._UHFK, "phase0": 180}, solver="HWAVE",
                     output_dir=None)

    def test_rpa_with_phase_raises(self):
        with pytest.raises(ValueError, match="phase"):
            generate({**self._UHFK, "calcmode": "rpa", "phase0": 180},
                     solver="HWAVE", output_dir=None)

    def test_uhfr_with_phase_ok(self):
        # UHFR is real-space (.def); a boundary twist is allowed.
        out = generate({**self._UHFK, "calcmode": "uhfr", "phase0": 180},
                       solver="HWAVE", output_dir=None)
        assert isinstance(out, ExpertModeOutput)
