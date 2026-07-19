"""Unit tests for geometry_output module.

Tests for ``build_xsf`` / ``build_geometry`` and their data classes.
"""
from __future__ import annotations

import math
import os
import tempfile

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.lattice.geometry_output import (
    _cell_diff, build_geometry, GeometryData, build_xsf, XsfData,
)


# The production print_xsf / print_geometry entry points were removed (they
# were thin wrappers over build_*().write()).  These local shims preserve the
# behaviour the content tests below exercise.
def print_xsf(StdI):
    build_xsf(StdI).write()


def print_geometry(StdI):
    data = build_geometry(StdI)
    if data is not None:
        data.write()


def _make_stdi(
    lattice: str = "chain",
    model: str = "hubbard",
    solver: str = "HPhi",
    L: int = 4,
) -> StdIntList:
    """Create a minimal StdIntList for geometry tests.

    Sets up a 1D chain with identity direct vectors (a=1) by default.
    """
    StdI = StdIntList()
    StdI.solver = solver
    StdI.model = model
    StdI.lattice = lattice

    # Box = [[L,0,0],[0,1,0],[0,0,1]]
    StdI.box[:, :] = 0
    StdI.box[0, 0] = L
    StdI.box[1, 1] = 1
    StdI.box[2, 2] = 1

    StdI.NCell = L
    StdI.NsiteUC = 1
    StdI.nsite = L

    # Direct lattice vectors (identity)
    StdI.direct[:, :] = 0.0
    StdI.direct[0, 0] = 1.0
    StdI.direct[1, 1] = 1.0
    StdI.direct[2, 2] = 1.0

    # Cell array
    StdI.Cell = np.zeros((L, 3), dtype=float)
    for i in range(L):
        StdI.Cell[i, 0] = i

    # tau
    StdI.tau = np.zeros((1, 3))

    # Phase
    StdI.phase[:] = 0.0

    # length (for CONVVEC)
    StdI.length[:] = 1.0

    return StdI


class TestPrintGeometry:
    """Tests for print_geometry."""

    def test_creates_geometry_dat(self):
        """Test that geometry.dat is created."""
        StdI = _make_stdi()
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_geometry(StdI)
                assert os.path.exists("geometry.dat")
            finally:
                os.chdir(orig)

    def test_contains_direct_vectors(self):
        """Test that geometry.dat contains the direct lattice vectors."""
        StdI = _make_stdi()
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_geometry(StdI)
                with open("geometry.dat") as f:
                    content = f.read()
                # First 3 lines are the direct vectors
                lines = content.strip().split("\n")
                assert len(lines) >= 3
                # First line should contain 1.0 0.0 0.0 (identity row 0)
                vals = [float(x) for x in lines[0].split()]
                assert vals[0] == pytest.approx(1.0)
                assert vals[1] == pytest.approx(0.0)
                assert vals[2] == pytest.approx(0.0)
            finally:
                os.chdir(orig)

    def test_contains_phase(self):
        """Test that geometry.dat contains the phase line."""
        StdI = _make_stdi()
        StdI.phase[0] = 90.0
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_geometry(StdI)
                with open("geometry.dat") as f:
                    lines = f.readlines()
                # Line 4 (index 3) is phase
                phase_vals = [float(x) for x in lines[3].split()]
                assert phase_vals[0] == pytest.approx(90.0)
            finally:
                os.chdir(orig)

    def test_cell_count(self):
        """Test that the correct number of cell lines are written."""
        StdI = _make_stdi(L=4)
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_geometry(StdI)
                with open("geometry.dat") as f:
                    lines = f.readlines()
                # 3 direct + 1 phase + 3 box + 4 cells = 11
                assert len(lines) == 11
            finally:
                os.chdir(orig)

    def test_kondo_doubles_cells(self):
        """Test that kondo model doubles the cell lines."""
        StdI = _make_stdi(L=4, model="kondo")
        StdI.nsite = 8  # kondo doubles
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_geometry(StdI)
                with open("geometry.dat") as f:
                    lines = f.readlines()
                # 3 direct + 1 phase + 3 box + 4 cells + 4 kondo = 15
                assert len(lines) == 15
            finally:
                os.chdir(orig)

    def test_hwave_uhfk_skips(self):
        """Test that HWAVE with uhfk mode skips geometry output."""
        StdI = _make_stdi(solver="HWAVE")
        StdI.calcmode = "uhfk"
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_geometry(StdI)
                assert not os.path.exists("geometry.dat")
            finally:
                os.chdir(orig)


class TestPrintXsf:
    """Tests for print_xsf."""

    def test_creates_lattice_xsf(self):
        """Test that lattice.xsf is created."""
        StdI = _make_stdi()
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_xsf(StdI)
                assert os.path.exists("lattice.xsf")
            finally:
                os.chdir(orig)

    def test_contains_crystal_header(self):
        """Test that lattice.xsf contains CRYSTAL and PRIMVEC headers."""
        StdI = _make_stdi()
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_xsf(StdI)
                with open("lattice.xsf") as f:
                    content = f.read()
                assert "CRYSTAL" in content
                assert "PRIMVEC" in content
                assert "PRIMCOORD" in content
            finally:
                os.chdir(orig)

    def test_no_convvec_for_chain(self):
        """Test that chain lattice does not write CONVVEC."""
        StdI = _make_stdi(lattice="chain")
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_xsf(StdI)
                with open("lattice.xsf") as f:
                    content = f.read()
                assert "CONVVEC" not in content
            finally:
                os.chdir(orig)

    def test_convvec_for_orthorhombic(self):
        """Test that orthorhombic lattice writes CONVVEC."""
        StdI = _make_stdi(lattice="orthorhombic")
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_xsf(StdI)
                with open("lattice.xsf") as f:
                    content = f.read()
                assert "CONVVEC" in content
            finally:
                os.chdir(orig)

    def test_primcoord_count(self):
        """Test that PRIMCOORD section has correct atom count."""
        StdI = _make_stdi(L=4)
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_xsf(StdI)
                with open("lattice.xsf") as f:
                    lines = f.readlines()
                # Find PRIMCOORD line, next line has "N 1"
                for i, line in enumerate(lines):
                    if "PRIMCOORD" in line:
                        count_line = lines[i + 1].strip()
                        parts = count_line.split()
                        assert int(parts[0]) == 4  # NCell * NsiteUC
                        break
                else:
                    pytest.fail("PRIMCOORD not found")
            finally:
                os.chdir(orig)


class TestCellDiff:
    """Tests for _cell_diff helper."""

    def test_basic_difference(self):
        """Test basic cell coordinate difference."""
        Cell = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ])
        diff = _cell_diff(Cell, 2, 0)
        assert diff == [2, 1, 0]

    def test_returns_integers(self):
        """Test that result is a list of integers."""
        Cell = np.array([
            [0.0, 0.0, 0.0],
            [3.0, 2.0, 1.0],
        ])
        diff = _cell_diff(Cell, 1, 0)
        assert all(isinstance(x, int) for x in diff)
        assert diff == [3, 2, 1]

    def test_negative_difference(self):
        """Test negative cell difference (reversed indices)."""
        Cell = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 3.0],
        ])
        diff = _cell_diff(Cell, 0, 1)
        assert diff == [-1, -2, -3]

    def test_same_cell(self):
        """Test that same cell gives zero difference."""
        Cell = np.array([
            [5.0, 5.0, 5.0],
            [1.0, 2.0, 3.0],
        ])
        diff = _cell_diff(Cell, 0, 0)
        assert diff == [0, 0, 0]

    def test_float_to_int_conversion(self):
        """Test that float coordinates are converted to int."""
        Cell = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 2.7, 3.9],  # floats that should truncate
        ])
        diff = _cell_diff(Cell, 1, 0)
        # int() truncates towards zero
        assert diff == [1, 2, 3]


class TestGeometryDataBuildWrite:
    """D5-1: build_geometry / GeometryData (independent lattice output)."""

    def test_build_returns_data_with_native_types(self):
        d = build_geometry(_make_stdi(L=4))
        assert isinstance(d, GeometryData)
        # JSON-native: plain int/float only
        assert all(isinstance(x, float) for r in d.direct for x in r)
        assert all(isinstance(x, int) for r in d.box for x in r)
        assert all(isinstance(x, int) for s in d.sites for x in s)
        assert len(d.sites) == 4  # NCell * NsiteUC

    def test_suppressed_for_explicit_wannier_modes(self):
        """Suppression is the H-wave plugin's call (wants_geometry_file)."""
        s = _make_stdi(L=4, solver="HWAVE")
        s.calcmode = "uhfk"
        assert build_geometry(s) is None
        s.calcmode = "rpa"
        assert build_geometry(s) is None

    def test_suppressed_for_uhfk_alias_without_calcmode(self):
        """The UHFK alias suppresses geometry.dat even with calcmode unset
        (formerly inconsistent with the HWAVE + calcmode=uhfk spelling)."""
        s = _make_stdi(L=4, solver="UHFK")
        assert build_geometry(s) is None

    def test_emitted_for_hwave_default_calcmode(self):
        """HWAVE with calcmode unset keeps emitting geometry.dat
        (compatibility; see the uhfr_first_step integration case)."""
        s = _make_stdi(L=4, solver="HWAVE")
        assert build_geometry(s) is not None

    def test_emitted_for_expert_solvers(self):
        assert build_geometry(_make_stdi(L=4, solver="HPhi")) is not None

    def test_to_from_dict_roundtrip(self):
        d = build_geometry(_make_stdi(L=4))
        assert GeometryData.from_dict(d.to_dict()).to_dict() == d.to_dict()

    def test_write_parity_with_print(self, tmp_path, monkeypatch):
        s = _make_stdi(L=4)
        monkeypatch.chdir(tmp_path)
        print_geometry(s)
        via_print = (tmp_path / "geometry.dat").read_text()
        build_geometry(s).write(tmp_path)
        via_build = (tmp_path / "geometry.dat").read_text()
        assert via_print == via_build

    def test_kondo_doubles_sites(self):
        d = build_geometry(_make_stdi(model="kondo", L=4))
        assert len(d.sites) == 8  # doubled


class TestXsfDataBuildWrite:
    """xsf-1: build_xsf / XsfData (independent lattice output)."""

    def test_build_returns_native_types(self):
        d = build_xsf(_make_stdi(L=4))
        assert isinstance(d, XsfData)
        assert all(isinstance(x, float) for r in d.primvec for x in r)
        assert all(isinstance(x, float) for c in d.coords for x in c)
        assert len(d.coords) == 4  # NCell * NsiteUC

    def test_convvec_present_for_orthorhombic(self):
        assert build_xsf(_make_stdi(lattice="orthorhombic")).convvec is not None

    def test_no_convvec_for_chain(self):
        assert build_xsf(_make_stdi(lattice="chain")).convvec is None

    def test_to_from_dict_roundtrip(self):
        d = build_xsf(_make_stdi(lattice="orthorhombic", L=4))
        assert XsfData.from_dict(d.to_dict()).to_dict() == d.to_dict()

    def test_write_parity_with_print(self, tmp_path, monkeypatch):
        for lat in ("chain", "orthorhombic"):
            s = _make_stdi(lattice=lat, L=4)
            monkeypatch.chdir(tmp_path)
            print_xsf(s)
            via_print = (tmp_path / "lattice.xsf").read_text()
            build_xsf(s).write(tmp_path)
            via_build = (tmp_path / "lattice.xsf").read_text()
            assert via_print == via_build, lat
