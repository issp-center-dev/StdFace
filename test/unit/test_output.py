"""Unit tests for the D2 solver output containers (core/output.py)."""
from __future__ import annotations

import inspect
import os

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.core.output import (
    SolverOutput, WannierModeOutput, build_wannier_output,
)
import stdface.writer.wannier90_writer as ew


def _make_uhfk_stdi():
    s = StdIntList()
    s.solver = "HWAVE"
    s.NsiteUC = 1
    s.NCell = 2
    s.fileprefix = ""
    s.export_all = None
    s.box = np.array([[2, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
    s.rbox = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 2]], dtype=int)
    s.Cell = np.zeros((2, 3), dtype=int)
    s.Cell[1, 0] = 1
    s.direct = np.eye(3)
    s.tau = np.zeros((1, 3))
    s.trans_list = [(-1.0 + 0j, 0, 0, 1, 0)]
    s.Cintra_list = [(4.0, 0), (4.0, 1)]
    for a in ("Cinter_list", "Hund_list", "Ex_list",
              "PairLift_list", "PairHopp_list"):
        setattr(s, a, [])
    return s


class TestSolverOutputABC:
    def test_is_abstract(self):
        assert inspect.isabstract(SolverOutput)
        with pytest.raises(TypeError):
            SolverOutput()  # cannot instantiate

    def test_wannier_is_subclass(self):
        assert issubclass(WannierModeOutput, SolverOutput)


class TestWannierModeOutput:
    def test_build_shape(self):
        out = build_wannier_output(_make_uhfk_stdi())
        assert isinstance(out, WannierModeOutput)
        assert out.geom_fname == "geom.dat"
        assert out.interactions  # at least transfer + coulombintra

    def test_prefix_applied(self):
        s = _make_uhfk_stdi()
        s.fileprefix = "pre"
        out = build_wannier_output(s)
        assert out.geom_fname == "pre_geom.dat"

    def test_to_dict_keys(self):
        d = build_wannier_output(_make_uhfk_stdi()).to_dict()
        assert set(d) == {"geometry", "geom_fname", "interactions"}

    def test_write_parity_with_direct_export(self, tmp_path):
        via_direct = tmp_path / "direct"
        via_cont = tmp_path / "cont"
        via_direct.mkdir()
        orig = os.getcwd()
        os.chdir(via_direct)
        try:
            ew.export_geometry(_make_uhfk_stdi())
            ew.export_interaction(_make_uhfk_stdi())
        finally:
            os.chdir(orig)
        build_wannier_output(_make_uhfk_stdi()).write(via_cont)

        files = sorted(p.name for p in via_direct.iterdir())
        assert files == sorted(p.name for p in via_cont.iterdir())
        assert files  # non-empty
        for f in files:
            assert (via_direct / f).read_text() == (via_cont / f).read_text(), f

    def test_write_creates_directory(self, tmp_path):
        target = tmp_path / "nested" / "out"
        build_wannier_output(_make_uhfk_stdi()).write(target)
        assert (target / "geom.dat").exists()
