"""Unit tests for boost_output module.

Tests for ``write_boost_mag_field``, ``write_boost_j_full``,
``write_boost_j_symmetric``, ``write_boost_6spin_star``, and
``write_boost_6spin_pair``.
"""
from __future__ import annotations

import io

import numpy as np
import pytest

from stdface_vals import StdIntList
from lattice.boost_output import (
    write_boost_mag_field,
    write_boost_j_full,
    write_boost_j_symmetric,
    write_boost_6spin_star,
    write_boost_6spin_pair,
)


# ---------------------------------------------------------------------------
# write_boost_mag_field
# ---------------------------------------------------------------------------
class TestWriteBoostMagField:
    """Tests for write_boost_mag_field."""

    def test_output_format(self):
        """Output contains header and three scaled field components."""
        StdI = StdIntList()
        StdI.Gamma = 2.0
        StdI.Gamma_y = 4.0
        StdI.h = 6.0
        buf = io.StringIO()
        write_boost_mag_field(buf, StdI)
        lines = buf.getvalue().splitlines()
        assert lines[0] == "# Magnetic field"
        vals = lines[1].split()
        assert len(vals) == 3
        assert float(vals[0]) == pytest.approx(-1.0)   # -0.5 * 2.0
        assert float(vals[1]) == pytest.approx(-2.0)   # -0.5 * 4.0
        assert float(vals[2]) == pytest.approx(-3.0)   # -0.5 * 6.0

    def test_zero_fields(self):
        """All zeros when fields are zero."""
        StdI = StdIntList()
        StdI.Gamma = 0.0
        StdI.Gamma_y = 0.0
        StdI.h = 0.0
        buf = io.StringIO()
        write_boost_mag_field(buf, StdI)
        lines = buf.getvalue().splitlines()
        vals = lines[1].split()
        for v in vals:
            assert float(v) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# write_boost_j_full
# ---------------------------------------------------------------------------
class TestWriteBoostJFull:
    """Tests for write_boost_j_full."""

    def test_identity_matrix(self):
        """Identity matrix scaled by 0.25 produces diagonal 0.25 values."""
        J = np.eye(3)
        buf = io.StringIO()
        write_boost_j_full(buf, J)
        lines = buf.getvalue().splitlines()
        assert len(lines) == 3
        for i, line in enumerate(lines):
            vals = [float(v) for v in line.split()]
            for j, v in enumerate(vals):
                expected = 0.25 if i == j else 0.0
                assert v == pytest.approx(expected)

    def test_custom_scale(self):
        """Custom scale factor is applied correctly."""
        J = np.ones((3, 3))
        buf = io.StringIO()
        write_boost_j_full(buf, J, scale=0.5)
        lines = buf.getvalue().splitlines()
        for line in lines:
            vals = [float(v) for v in line.split()]
            assert all(v == pytest.approx(0.5) for v in vals)

    def test_default_scale(self):
        """Default scale is 0.25."""
        J = np.ones((3, 3))
        buf = io.StringIO()
        write_boost_j_full(buf, J)
        lines = buf.getvalue().splitlines()
        for line in lines:
            vals = [float(v) for v in line.split()]
            assert all(v == pytest.approx(0.25) for v in vals)

    def test_three_rows(self):
        """Output always has exactly three rows."""
        J = np.zeros((3, 3))
        buf = io.StringIO()
        write_boost_j_full(buf, J)
        lines = buf.getvalue().splitlines()
        assert len(lines) == 3


# ---------------------------------------------------------------------------
# write_boost_j_symmetric
# ---------------------------------------------------------------------------
class TestWriteBoostJSymmetric:
    """Tests for write_boost_j_symmetric."""

    def test_identity_matrix(self):
        """Identity matrix produces diagonal 0.25 values."""
        J = np.eye(3)
        buf = io.StringIO()
        write_boost_j_symmetric(buf, J)
        lines = buf.getvalue().splitlines()
        assert len(lines) == 3
        for i, line in enumerate(lines):
            vals = [float(v) for v in line.split()]
            for j, v in enumerate(vals):
                expected = 0.25 if i == j else 0.0
                assert v == pytest.approx(expected)

    def test_uses_upper_triangle(self):
        """Only upper-triangle elements of J are used."""
        J = np.zeros((3, 3))
        J[0, 1] = 2.0
        J[1, 0] = 999.0  # should be ignored
        buf = io.StringIO()
        write_boost_j_symmetric(buf, J)
        lines = buf.getvalue().splitlines()
        # Row 0: J[0,0]=0, J[0,1]=2, J[0,2]=0
        row0 = [float(v) for v in lines[0].split()]
        assert row0[1] == pytest.approx(0.5)  # 0.25 * 2.0
        # Row 1: J[0,1]=2, J[1,1]=0, J[1,2]=0
        row1 = [float(v) for v in lines[1].split()]
        assert row1[0] == pytest.approx(0.5)  # 0.25 * 2.0 (from J[0,1])

    def test_three_rows(self):
        """Output always has exactly three rows."""
        J = np.zeros((3, 3))
        buf = io.StringIO()
        write_boost_j_symmetric(buf, J)
        lines = buf.getvalue().splitlines()
        assert len(lines) == 3

    def test_default_scale(self):
        """Default scale is 0.25."""
        J = np.ones((3, 3))
        buf = io.StringIO()
        write_boost_j_symmetric(buf, J)
        lines = buf.getvalue().splitlines()
        for line in lines:
            vals = [float(v) for v in line.split()]
            assert all(v == pytest.approx(0.25) for v in vals)


# ---------------------------------------------------------------------------
# write_boost_6spin_star
# ---------------------------------------------------------------------------
class TestWriteBoost6spinStar:
    """Tests for write_boost_6spin_star."""

    def _make_stdi(self, num_pivot: int = 2) -> StdIntList:
        StdI = StdIntList()
        StdI.num_pivot = num_pivot
        StdI.list_6spin_star = np.zeros((num_pivot, 7), dtype=int)
        for ip in range(num_pivot):
            StdI.list_6spin_star[ip, :] = [8, 1, 1, 1, 1, 1, 1]
        return StdI

    def test_writes_header(self):
        """Output starts with the header comment."""
        StdI = self._make_stdi(1)
        buf = io.StringIO()
        write_boost_6spin_star(buf, StdI)
        lines = buf.getvalue().splitlines()
        assert lines[0] == "# StdI->list_6spin_star"

    def test_writes_pivot_labels(self):
        """Each pivot has a '# pivot N' label."""
        StdI = self._make_stdi(3)
        buf = io.StringIO()
        write_boost_6spin_star(buf, StdI)
        text = buf.getvalue()
        for i in range(3):
            assert f"# pivot {i}" in text

    def test_writes_correct_values(self):
        """Values from list_6spin_star are written correctly."""
        StdI = self._make_stdi(1)
        StdI.list_6spin_star[0, :] = [7, 2, 3, 4, 5, 6, 7]
        buf = io.StringIO()
        write_boost_6spin_star(buf, StdI)
        lines = buf.getvalue().splitlines()
        # line 0: header, line 1: "# pivot 0", line 2: values
        vals = [int(v) for v in lines[2].split()]
        assert vals == [7, 2, 3, 4, 5, 6, 7]

    def test_num_pivot_lines(self):
        """Number of value lines equals num_pivot."""
        StdI = self._make_stdi(4)
        buf = io.StringIO()
        write_boost_6spin_star(buf, StdI)
        lines = buf.getvalue().splitlines()
        # 1 header + 4*(1 label + 1 value) = 9 lines
        assert len(lines) == 1 + 4 * 2


# ---------------------------------------------------------------------------
# write_boost_6spin_pair
# ---------------------------------------------------------------------------
class TestWriteBoost6spinPair:
    """Tests for write_boost_6spin_pair."""

    def _make_stdi(self, num_pivot: int = 1, num_intr: int = 3) -> StdIntList:
        StdI = StdIntList()
        StdI.num_pivot = num_pivot
        StdI.list_6spin_star = np.zeros((num_pivot, 7), dtype=int)
        for ip in range(num_pivot):
            StdI.list_6spin_star[ip, 0] = num_intr
        StdI.list_6spin_pair = np.zeros((num_pivot, 7, num_intr), dtype=int)
        return StdI

    def test_writes_header(self):
        """Output starts with the header comment."""
        StdI = self._make_stdi()
        buf = io.StringIO()
        write_boost_6spin_pair(buf, StdI)
        lines = buf.getvalue().splitlines()
        assert lines[0] == "# StdI->list_6spin_pair"

    def test_writes_correct_interaction_count(self):
        """Number of interaction rows per pivot equals list_6spin_star[ip, 0]."""
        StdI = self._make_stdi(num_pivot=1, num_intr=5)
        buf = io.StringIO()
        write_boost_6spin_pair(buf, StdI)
        lines = buf.getvalue().splitlines()
        # 1 header + 1 pivot label + 5 interaction rows = 7
        assert len(lines) == 7

    def test_writes_correct_values(self):
        """Values from list_6spin_pair are written correctly."""
        StdI = self._make_stdi(num_pivot=1, num_intr=2)
        StdI.list_6spin_pair[0, :, 0] = [0, 1, 2, 3, 4, 5, 1]
        StdI.list_6spin_pair[0, :, 1] = [1, 2, 0, 3, 4, 5, 2]
        buf = io.StringIO()
        write_boost_6spin_pair(buf, StdI)
        lines = buf.getvalue().splitlines()
        row0 = [int(v) for v in lines[2].split()]
        row1 = [int(v) for v in lines[3].split()]
        assert row0 == [0, 1, 2, 3, 4, 5, 1]
        assert row1 == [1, 2, 0, 3, 4, 5, 2]

    def test_multiple_pivots(self):
        """Multiple pivots each get their own section."""
        StdI = self._make_stdi(num_pivot=2, num_intr=2)
        buf = io.StringIO()
        write_boost_6spin_pair(buf, StdI)
        text = buf.getvalue()
        assert "# pivot 0" in text
        assert "# pivot 1" in text
