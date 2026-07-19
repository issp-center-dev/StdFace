"""Unit tests for the D2 solver output containers (core/output.py)."""
from __future__ import annotations

import inspect

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.core.output import (
    SolverOutput, WannierModeOutput, build_wannier_output, ExpertModeOutput,
    OutputFormat, DefFileFormat, JSONFormat,
)
from stdface.writer.common_writer import (
    LocSpnData, TransData, NamelistData, ModParaData, GreenOneData,
)
from stdface.writer.interaction_writer import InteractionData


def _make_expert_output(green_two=None):
    return ExpertModeOutput(
        locspn=LocSpnData([0, 1]),
        trans=TransData([(0, 0, 1, 1, 1.0, 0.0)]),
        interactions=[InteractionData(
            "coulombintra.def", "NCoulombIntra", "banner", [(0, 4.0)])],
        modpara=ModParaData([("raw", "HPhi_Cal_Parameters"),
                             ("kv", "Nsite", 4, "<5d")]),
        namelist=NamelistData([("ModPara", "modpara.def")]),
        green_one=GreenOneData([(0, 0, 0, 0)]),
        green_two=green_two,
    )


class TestExpertModeOutput:
    """D2b-2: ExpertModeOutput container (write / to_dict)."""

    def test_is_solver_output(self):
        assert issubclass(ExpertModeOutput, SolverOutput)

    def test_write_emits_files(self, tmp_path):
        _make_expert_output().write(tmp_path)
        for fn in ("locspn.def", "trans.def", "coulombintra.def",
                   "modpara.def", "greenone.def", "namelist.def"):
            assert (tmp_path / fn).exists(), fn
        # green_two is None -> not written
        assert not (tmp_path / "greentwo.def").exists()

    def test_to_dict_keys(self):
        d = _make_expert_output().to_dict()
        assert set(d) == {"locspn", "trans", "interactions", "modpara",
                          "namelist", "green_one", "green_two",
                          "solver_files"}
        assert d["green_two"] is None
        assert d["interactions"][0]["count_label"] == "NCoulombIntra"

    def test_creates_nested_directory(self, tmp_path):
        target = tmp_path / "a" / "b"
        _make_expert_output().write(target)
        assert (target / "modpara.def").exists()


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

    def test_write_creates_directory(self, tmp_path):
        target = tmp_path / "nested" / "out"
        build_wannier_output(_make_uhfk_stdi()).write(target)
        assert (target / "geom.dat").exists()


class TestOutputFormats:
    """D3: OutputFormat strategies (DefFileFormat / JSONFormat)."""

    def test_outputformat_is_abstract(self):
        with pytest.raises(TypeError):
            OutputFormat()

    def test_deffileformat_writes_def_files(self, tmp_path):
        DefFileFormat().write_output(_make_expert_output(), tmp_path)
        assert (tmp_path / "modpara.def").exists()
        assert (tmp_path / "namelist.def").exists()
        assert (tmp_path / "trans.def").exists()

    def test_jsonformat_writes_single_json(self, tmp_path):
        out = _make_expert_output()
        JSONFormat().write_output(out, tmp_path)
        path = tmp_path / "stdface_output.json"
        assert path.exists()
        import json
        assert json.loads(path.read_text()) == out.to_dict()
        # no .def files emitted by the JSON format
        assert not (tmp_path / "modpara.def").exists()

    def test_jsonformat_custom_filename_and_nested_dir(self, tmp_path):
        target = tmp_path / "sub"
        JSONFormat("out.json").write_output(_make_expert_output(), target)
        assert (target / "out.json").exists()

    def test_jsonformat_on_wannier_output(self, tmp_path):
        out = build_wannier_output(_make_uhfk_stdi())
        JSONFormat().write_output(out, tmp_path)
        import json
        assert json.loads((tmp_path / "stdface_output.json").read_text()) == out.to_dict()
