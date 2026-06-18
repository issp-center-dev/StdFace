"""Unit tests for GnuplotBuffer / GnuplotData (A4)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.lattice.site_util import (
    GnuplotBuffer,
    GnuplotData,
    _LATTICE_GP_FOOTER,
)


def _make_stdi_2d() -> StdIntList:
    """Minimal StdIntList with the fields _write_gnuplot_header reads."""
    StdI = StdIntList()
    StdI.box[:, :] = 0
    StdI.box[0, 0] = 2
    StdI.box[1, 1] = 2
    StdI.box[2, 2] = 1
    StdI.direct[:, :] = 0.0
    StdI.direct[0, 0] = 1.0
    StdI.direct[1, 1] = 1.0
    StdI.direct[2, 2] = 1.0
    return StdI


class TestGnuplotBuffer:
    def test_empty_is_falsy(self):
        assert not GnuplotBuffer()

    def test_nonempty_is_truthy(self):
        buf = GnuplotBuffer()
        buf.add(0, 1, 0.0, 0.0, 1.0, 0.0, 1)
        assert buf

    def test_build_returns_gnuplotdata(self):
        buf = GnuplotBuffer()
        buf.add(0, 1, 0.0, 0.0, 1.0, 0.0, 1)
        data = buf.build(_make_stdi_2d())
        assert isinstance(data, GnuplotData)

    def test_build_includes_header_bond_footer(self):
        buf = GnuplotBuffer()
        buf.add(0, 1, 0.0, 0.0, 1.0, 0.0, 1)
        content = buf.build(_make_stdi_2d()).content
        assert "set xrange" in content          # header
        assert 'set label "0"' in content       # bond endpoint label
        assert "set arrow from" in content       # bond arrow (connect < 3)
        assert content.endswith(_LATTICE_GP_FOOTER)

    def test_no_arrow_when_connect_ge_3(self):
        buf = GnuplotBuffer()
        buf.add(0, 1, 0.0, 0.0, 1.0, 0.0, 3)
        content = buf.build(_make_stdi_2d()).content
        # labels present, but no bond arrow for connect >= 3
        assert 'set label' in content
        assert "set arrow from 0.000000, 0.000000 to 1.000000" not in content

    def test_build_is_pure_no_file(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        buf = GnuplotBuffer()
        buf.add(0, 1, 0.0, 0.0, 1.0, 0.0, 1)
        buf.build(_make_stdi_2d())
        assert not (tmp_path / "lattice.gp").exists()


class TestGnuplotData:
    def test_write_creates_file(self, tmp_path):
        data = GnuplotData(content="hello gp\n")
        data.write(tmp_path)
        assert (tmp_path / "lattice.gp").read_text() == "hello gp\n"

    def test_write_defaults_to_cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        GnuplotData(content="x\n").write()
        assert (tmp_path / "lattice.gp").read_text() == "x\n"

    def test_to_dict_from_dict_roundtrip(self):
        data = GnuplotData(content="some\nscript\n")
        restored = GnuplotData.from_dict(data.to_dict())
        assert restored == data
        assert restored.content == "some\nscript\n"
