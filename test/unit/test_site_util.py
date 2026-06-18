"""Unit tests for site_util module.

Tests for ``_cell_vector``, ``_fold_to_cell``, ``_fold_site``,
``_validate_box_params``, ``init_site``, ``find_site``, ``set_label``,
and ``lattice_gp``.
"""
from __future__ import annotations

import io
import math

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.lattice.site_util import (
    _cell_vector,
    _fold_to_cell, _fold_site, _find_cell_index, _write_gnuplot_header,
    _write_gnuplot_bond,
    _validate_box_params,
    _det_and_cofactor, _compute_reciprocal_box, _enumerate_cells,
    _LATTICE_GP_FOOTER,
    lattice_gp, close_lattice_xsf,
    init_site, find_site, set_label,
    set_local_spin_flags,
)
from stdface.core.stdface_vals import ModelType, SolverType, NaN_i


def _make_stdi_chain(L: int = 4) -> StdIntList:
    """Create a minimal StdIntList with a 1D chain super-cell already set up.

    Bypasses ``init_site`` by manually setting the fields that ``init_site``
    would compute (box, rbox, NCell, Cell, tau, ExpPhase, etc.).
    """
    StdI = StdIntList()
    StdI.solver = "HPhi"
    StdI.model = "hubbard"
    StdI.pi = math.acos(-1.0)
    StdI.pi180 = StdI.pi / 180.0

    # 1D chain: box = [[L,0,0],[0,1,0],[0,0,1]]
    StdI.box[:, :] = 0
    StdI.box[0, 0] = L
    StdI.box[1, 1] = 1
    StdI.box[2, 2] = 1

    StdI.NCell = L
    StdI.NsiteUC = 1

    # Direct lattice vectors (identity)
    StdI.direct[:, :] = 0.0
    StdI.direct[0, 0] = 1.0
    StdI.direct[1, 1] = 1.0
    StdI.direct[2, 2] = 1.0

    # Reciprocal box (cofactor matrix of box)
    StdI.rbox[:, :] = 0
    StdI.rbox[0, 0] = 1
    StdI.rbox[1, 1] = L
    StdI.rbox[2, 2] = L

    # Cell array
    StdI.Cell = np.zeros((L, 3), dtype=int)
    for i in range(L):
        StdI.Cell[i, 0] = i

    # tau
    StdI.tau = np.zeros((1, 3))

    # Phase = 0 → ExpPhase = 1
    StdI.phase[:] = 0.0
    for ii in range(3):
        StdI.ExpPhase[ii] = 1.0 + 0j
        StdI.AntiPeriod[ii] = 0

    return StdI


# ===================================================================
#  _cell_vector
# ===================================================================


class TestCellVector:
    """Tests for _cell_vector helper."""

    def test_basic_extraction(self):
        """Extract a row from a float Cell array as int list."""
        Cell = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        assert _cell_vector(Cell, 0) == [1, 2, 3]
        assert _cell_vector(Cell, 1) == [4, 5, 6]

    def test_returns_int(self):
        """Values are int, not float."""
        Cell = np.array([[1.0, 0.0, 0.0]])
        result = _cell_vector(Cell, 0)
        assert all(isinstance(v, int) for v in result)

    def test_negative_values(self):
        """Negative coordinates handled correctly."""
        Cell = np.array([[-1.0, -2.0, -3.0]])
        assert _cell_vector(Cell, 0) == [-1, -2, -3]

    def test_zero_vector(self):
        """Zero coordinate row."""
        Cell = np.zeros((3, 3))
        assert _cell_vector(Cell, 0) == [0, 0, 0]

    def test_returns_list(self):
        """Return type is a plain list, not numpy array."""
        Cell = np.array([[1.0, 2.0, 3.0]])
        result = _cell_vector(Cell, 0)
        assert isinstance(result, list)
        assert len(result) == 3

    def test_equality_comparison(self):
        """Lists from _cell_vector support == comparison."""
        Cell = np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]])
        assert _cell_vector(Cell, 0) == _cell_vector(Cell, 1)

    def test_used_in_find_cell_index(self):
        """_cell_vector result is compatible with _find_cell_index."""
        StdI = _make_stdi_chain(4)
        cv = _cell_vector(StdI.Cell, 2)
        assert _find_cell_index(StdI, cv) == 2


# ===================================================================
#  _fold_to_cell
# ===================================================================


class TestFoldToCell:
    """Tests for the _fold_to_cell pure function."""

    def test_identity_cell_inside(self):
        """Coordinate inside a 4×1×1 cell stays unchanged."""
        box = np.array([[4, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
        rbox = np.array([[1, 0, 0], [0, 4, 0], [0, 0, 4]], dtype=int)
        ncell = 4
        nBox, fold = _fold_to_cell(rbox, ncell, box, [2, 0, 0])
        assert nBox == [0, 0, 0]
        assert fold == [2, 0, 0]

    def test_wraps_at_boundary(self):
        """Coordinate at cell boundary wraps to origin."""
        box = np.array([[4, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
        rbox = np.array([[1, 0, 0], [0, 4, 0], [0, 0, 4]], dtype=int)
        ncell = 4
        nBox, fold = _fold_to_cell(rbox, ncell, box, [4, 0, 0])
        assert nBox[0] == 1
        assert fold == [0, 0, 0]

    def test_negative_wraps(self):
        """Negative coordinate wraps to the end of the cell."""
        box = np.array([[4, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
        rbox = np.array([[1, 0, 0], [0, 4, 0], [0, 0, 4]], dtype=int)
        ncell = 4
        nBox, fold = _fold_to_cell(rbox, ncell, box, [-1, 0, 0])
        assert nBox[0] == -1
        assert fold == [3, 0, 0]

    def test_2d_square_cell(self):
        """Fold in a 3×3 square cell."""
        box = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 1]], dtype=int)
        rbox = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 9]], dtype=int)
        ncell = 9
        nBox, fold = _fold_to_cell(rbox, ncell, box, [3, 1, 0])
        assert nBox[0] == 1
        assert fold == [0, 1, 0]

    def test_matches_fold_site(self):
        """_fold_to_cell with main-cell params matches _fold_site."""
        StdI = _make_stdi_chain(6)
        coord = [7, 0, 0]
        expected = _fold_site(StdI, coord)
        result = _fold_to_cell(StdI.rbox, StdI.NCell, StdI.box, coord)
        assert result == expected

    def test_subcell_params(self):
        """_fold_to_cell works with sub-lattice parameters."""
        # 2×1×1 sub-cell inside a 4×1×1 cell
        boxsub = np.array([[2, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
        rboxsub = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 2]], dtype=int)
        ncellsub = 2
        nBox, fold = _fold_to_cell(rboxsub, ncellsub, boxsub, [3, 0, 0])
        assert nBox[0] == 1
        assert fold == [1, 0, 0]

    def test_origin(self):
        """Origin coordinate produces zero nBox and zero fold."""
        box = np.array([[4, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
        rbox = np.array([[1, 0, 0], [0, 4, 0], [0, 0, 4]], dtype=int)
        ncell = 4
        nBox, fold = _fold_to_cell(rbox, ncell, box, [0, 0, 0])
        assert nBox == [0, 0, 0]
        assert fold == [0, 0, 0]


# ===================================================================
#  _fold_site
# ===================================================================


class TestFoldSite:
    """Tests for _fold_site."""

    def test_inside_cell_stays(self):
        """Test that a site already inside the cell is unchanged."""
        StdI = _make_stdi_chain(4)
        nBox, fold = _fold_site(StdI, [0, 0, 0])
        assert nBox == [0, 0, 0]
        assert fold == [0, 0, 0]

    def test_site_at_boundary(self):
        """Test that a site at the cell boundary wraps around."""
        StdI = _make_stdi_chain(4)
        nBox, fold = _fold_site(StdI, [4, 0, 0])
        # [4,0,0] is outside the L=4 cell → should fold to [0,0,0] with nBox=[1,0,0]
        assert nBox[0] == 1
        assert fold == [0, 0, 0]

    def test_negative_site_wraps(self):
        """Test that a negative coordinate wraps correctly."""
        StdI = _make_stdi_chain(4)
        nBox, fold = _fold_site(StdI, [-1, 0, 0])
        # [-1,0,0] wraps to [3,0,0] with nBox=[-1,0,0]
        assert nBox[0] == -1
        assert fold == [3, 0, 0]

    def test_interior_site(self):
        """Test that an interior site (not at boundary) stays."""
        StdI = _make_stdi_chain(4)
        nBox, fold = _fold_site(StdI, [2, 0, 0])
        assert nBox == [0, 0, 0]
        assert fold == [2, 0, 0]


# ===================================================================
#  init_site
# ===================================================================


class TestInitSite:
    """Tests for init_site."""

    @staticmethod
    def _make_init_stdi(L: int = 4, W: int | None = None,
                        Height: int | None = None,
                        NsiteUC: int = 1) -> StdIntList:
        """Helper to create a StdIntList ready for init_site.

        Sets box to NaN_i (as _reset_vals does in the real code) so that
        init_site sees the LWH path rather than the box path.
        """
        NaN_i = 2147483647
        StdI = StdIntList()
        StdI.pi = math.acos(-1.0)
        StdI.pi180 = StdI.pi / 180.0
        StdI.L = L
        StdI.W = W if W is not None else NaN_i
        StdI.Height = Height if Height is not None else NaN_i
        StdI.NsiteUC = NsiteUC
        StdI.direct[:, :] = 0.0
        StdI.direct[0, 0] = 1.0
        StdI.direct[1, 1] = 1.0
        StdI.direct[2, 2] = 1.0
        StdI.phase[:] = 0.0
        # box must be NaN_i so init_site uses LWH path (not box path)
        StdI.box[:, :] = NaN_i
        return StdI

    def test_lwh_sets_box(self):
        """Test that L/W/Height sets the box matrix."""
        StdI = self._make_init_stdi(L=4, W=2)
        init_site(StdI, None, 2)
        assert StdI.box[0, 0] == 2  # W
        assert StdI.box[1, 1] == 4  # L
        assert StdI.NCell == 8

    def test_computes_ncell(self):
        """Test that NCell is computed as det(box)."""
        StdI = self._make_init_stdi(L=3, W=3)
        init_site(StdI, None, 2)
        assert StdI.NCell == 9

    def test_allocates_cell_array(self):
        """Test that Cell array has correct shape."""
        StdI = self._make_init_stdi(L=4)
        init_site(StdI, None, 2)
        assert StdI.Cell.shape == (4, 3)

    def test_allocates_tau(self):
        """Test that tau array has correct shape."""
        StdI = self._make_init_stdi(L=4, NsiteUC=2)
        init_site(StdI, None, 2)
        assert StdI.tau.shape == (2, 3)

    def test_anti_periodic_phase(self):
        """Test that phase=180 sets AntiPeriod flag."""
        StdI = self._make_init_stdi(L=4)
        StdI.phase[0] = 180.0
        init_site(StdI, None, 2)
        assert StdI.AntiPeriod[0] == 1
        assert StdI.AntiPeriod[1] == 0

    def test_writes_gnuplot_header(self):
        """Test that 2D init_site writes gnuplot lattice.gp header."""
        StdI = self._make_init_stdi(L=2, W=2)
        fp = io.StringIO()
        init_site(StdI, fp, 2)
        content = fp.getvalue()
        assert "set xrange" in content
        assert "set yrange" in content
        assert "set arrow" in content


# ===================================================================
#  _validate_box_params
# ===================================================================


class TestValidateBoxParams:
    """Tests for _validate_box_params helper."""

    def _nan_box(self) -> np.ndarray:
        """Return a 3x3 box filled with NaN_i."""
        return np.full((3, 3), NaN_i, dtype=int)

    def test_lwh_specified_fills_diagonal(self):
        """When L/W/Height are set, box becomes diag(W, L, Height)."""
        box = self._nan_box()
        L, W, H = _validate_box_params(4, 3, 2, box)
        assert (L, W, H) == (4, 3, 2)
        expected = np.array([[3, 0, 0], [0, 4, 0], [0, 0, 2]], dtype=int)
        np.testing.assert_array_equal(box, expected)

    def test_lwh_partial_defaults_to_one(self):
        """Unset L/W/Height entries default to 1."""
        box = self._nan_box()
        L, W, H = _validate_box_params(6, NaN_i, NaN_i, box)
        assert (L, W, H) == (6, 1, 1)
        expected = np.array([[1, 0, 0], [0, 6, 0], [0, 0, 1]], dtype=int)
        np.testing.assert_array_equal(box, expected)

    def test_box_specified_uses_defaults_identity(self):
        """When box entries set and no defaults, identity is used as default."""
        box = self._nan_box()
        box[0, 0] = 5
        # Others remain NaN_i → should be replaced by identity defaults
        L, W, H = _validate_box_params(NaN_i, NaN_i, NaN_i, box)
        assert (L, W, H) == (NaN_i, NaN_i, NaN_i)
        # box[0,0]=5 kept; off-diag default=0, diag default=1
        assert box[0, 0] == 5
        assert box[1, 1] == 1
        assert box[2, 2] == 1
        assert box[0, 1] == 0

    def test_box_specified_uses_custom_defaults(self):
        """When defaults are given, they are used for unset box entries."""
        box = self._nan_box()
        defaults = np.array([[10, 20, 30], [40, 50, 60], [70, 80, 90]])
        L, W, H = _validate_box_params(NaN_i, NaN_i, NaN_i, box, defaults=defaults)
        np.testing.assert_array_equal(box, defaults)

    def test_neither_specified_fills_from_defaults(self):
        """When nothing is specified, box is filled from defaults."""
        box = self._nan_box()
        L, W, H = _validate_box_params(NaN_i, NaN_i, NaN_i, box)
        # Default is identity
        expected = np.eye(3, dtype=int)
        np.testing.assert_array_equal(box, expected)

    def test_conflict_raises_system_exit(self):
        """Both L/W/Height and box entries specified → exit."""
        box = self._nan_box()
        box[0, 0] = 5  # box entry set
        with pytest.raises(SystemExit):
            _validate_box_params(4, NaN_i, NaN_i, box)  # L set too

    def test_conflict_with_suffix(self):
        """Conflict error uses suffix in label."""
        box = self._nan_box()
        box[1, 1] = 3
        with pytest.raises(SystemExit):
            _validate_box_params(NaN_i, 2, NaN_i, box, suffix="sub")

    def test_suffix_affects_labels(self, caplog):
        """Suffix parameter changes the logged parameter names."""
        import logging

        box = self._nan_box()
        with caplog.at_level(logging.INFO):
            L, W, H = _validate_box_params(4, 3, 2, box, suffix="sub")
        assert "Lsub" in caplog.text
        assert "Wsub" in caplog.text
        assert "Hsub" in caplog.text

    def test_no_suffix_height_label(self, caplog):
        """Without suffix, height label is 'Height'."""
        import logging

        box = self._nan_box()
        with caplog.at_level(logging.INFO):
            L, W, H = _validate_box_params(NaN_i, NaN_i, 5, box)
        assert "Height" in caplog.text

    def test_box_modified_in_place(self):
        """Box array is modified in-place, not replaced."""
        box = self._nan_box()
        original_id = id(box)
        _validate_box_params(4, 3, 2, box)
        assert id(box) == original_id


# ===================================================================
#  _det_and_cofactor
# ===================================================================


class TestDetAndCofactor:
    """Tests for _det_and_cofactor pure function."""

    def test_identity_matrix(self):
        """Identity matrix → det=1, cofactor=identity."""
        box = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        det, cof = _det_and_cofactor(box)
        assert det == 1
        np.testing.assert_array_equal(cof, np.eye(3))

    def test_diagonal_matrix(self):
        """Diagonal matrix → det=product of diagonals."""
        box = np.array([[4, 0, 0], [0, 3, 0], [0, 0, 1]], dtype=float)
        det, cof = _det_and_cofactor(box)
        assert det == 12
        assert cof[0, 0] == 3  # cofactor of (0,0) = 3*1 = 3
        assert cof[1, 1] == 4  # cofactor of (1,1) = 4*1 = 4
        assert cof[2, 2] == 12  # cofactor of (2,2) = 4*3 = 12

    def test_negative_det_flipped(self):
        """Negative determinant is flipped to positive."""
        # det of [[0,1,0],[2,0,0],[0,0,1]] is -2, should be flipped to 2
        box = np.array([[0, 1, 0], [2, 0, 0], [0, 0, 1]], dtype=float)
        det, cof = _det_and_cofactor(box)
        assert det == 2
        # cofactor should also be sign-flipped

    def test_zero_det(self):
        """Degenerate matrix returns det=0 (no exit — caller decides)."""
        box = np.array([[1, 0, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
        det, cof = _det_and_cofactor(box)
        assert det == 0

    def test_1d_chain(self):
        """1D chain L=4 → det=4, cofactor[0,0]=1, cofactor[1,1]=4."""
        box = np.array([[4, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        det, cof = _det_and_cofactor(box)
        assert det == 4
        assert cof[0, 0] == 1
        assert cof[1, 1] == 4

    def test_returns_float_array(self):
        """Cofactor matrix should be a numpy float array."""
        box = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]], dtype=float)
        det, cof = _det_and_cofactor(box)
        assert isinstance(cof, np.ndarray)
        assert cof.shape == (3, 3)

    def test_off_diagonal_box(self):
        """Non-diagonal box with known determinant."""
        # [[2,1,0],[0,3,0],[0,0,1]] → det = 2*3*1 - 0 = 6
        box = np.array([[2, 1, 0], [0, 3, 0], [0, 0, 1]], dtype=float)
        det, cof = _det_and_cofactor(box)
        assert det == 6


# ===================================================================
#  _compute_reciprocal_box
# ===================================================================


class TestComputeReciprocalBox:
    """Tests for _compute_reciprocal_box."""

    @staticmethod
    def _make_stdi_with_box(box: list[list[int]]) -> StdIntList:
        """Create StdIntList with a given box matrix."""
        StdI = StdIntList()
        for i in range(3):
            for j in range(3):
                StdI.box[i, j] = box[i][j]
        return StdI

    def test_identity_box(self):
        """Identity box → NCell=1, rbox=identity."""
        StdI = self._make_stdi_with_box([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        _compute_reciprocal_box(StdI)
        assert StdI.NCell == 1
        np.testing.assert_array_equal(StdI.rbox, np.eye(3, dtype=int))

    def test_1d_chain(self):
        """1D chain with L=4 → NCell=4."""
        StdI = self._make_stdi_with_box([[4, 0, 0], [0, 1, 0], [0, 0, 1]])
        _compute_reciprocal_box(StdI)
        assert StdI.NCell == 4
        assert StdI.rbox[0, 0] == 1
        assert StdI.rbox[1, 1] == 4

    def test_2d_square(self):
        """2D square W=3,L=3 → NCell=9."""
        StdI = self._make_stdi_with_box([[3, 0, 0], [0, 3, 0], [0, 0, 1]])
        _compute_reciprocal_box(StdI)
        assert StdI.NCell == 9

    def test_negative_det_flips_sign(self):
        """Negative determinant is corrected to positive NCell."""
        # box with det=-2
        StdI = self._make_stdi_with_box([[0, 1, 0], [2, 0, 0], [0, 0, 1]])
        _compute_reciprocal_box(StdI)
        assert StdI.NCell == 2

    def test_zero_det_exits(self):
        """Degenerate box (det=0) causes SystemExit."""
        StdI = self._make_stdi_with_box([[1, 0, 0], [1, 0, 0], [0, 0, 1]])
        with pytest.raises(SystemExit):
            _compute_reciprocal_box(StdI)


# ===================================================================
#  _enumerate_cells
# ===================================================================


class TestEnumerateCells:
    """Tests for _enumerate_cells."""

    @staticmethod
    def _make_stdi_for_enum(box: list[list[int]], NCell: int) -> StdIntList:
        """Create StdIntList with box and NCell ready for enumeration."""
        StdI = StdIntList()
        StdI.pi = math.acos(-1.0)
        StdI.pi180 = StdI.pi / 180.0
        for i in range(3):
            for j in range(3):
                StdI.box[i, j] = box[i][j]
        StdI.NCell = NCell
        # Compute rbox (needed by _fold_site)
        for ii in range(3):
            for jj in range(3):
                StdI.rbox[ii, jj] = (int(StdI.box[(ii + 1) % 3, (jj + 1) % 3])
                                     * int(StdI.box[(ii + 2) % 3, (jj + 2) % 3])
                                     - int(StdI.box[(ii + 1) % 3, (jj + 2) % 3])
                                     * int(StdI.box[(ii + 2) % 3, (jj + 1) % 3]))
        return StdI

    def test_1d_chain_cells(self):
        """1D chain L=4 produces 4 cells: [0,0,0]..[3,0,0]."""
        StdI = self._make_stdi_for_enum([[4, 0, 0], [0, 1, 0], [0, 0, 1]], 4)
        _enumerate_cells(StdI)
        assert StdI.Cell.shape == (4, 3)
        # All cells should have unique x-coordinates 0..3
        xs = sorted(StdI.Cell[:, 0].tolist())
        assert xs == [0, 1, 2, 3]

    def test_2d_square_cells(self):
        """2D square 2x2 produces 4 cells."""
        StdI = self._make_stdi_for_enum([[2, 0, 0], [0, 2, 0], [0, 0, 1]], 4)
        _enumerate_cells(StdI)
        assert StdI.Cell.shape == (4, 3)

    def test_single_cell(self):
        """Identity box produces exactly 1 cell at origin."""
        StdI = self._make_stdi_for_enum([[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1)
        _enumerate_cells(StdI)
        assert StdI.Cell.shape == (1, 3)
        np.testing.assert_array_equal(StdI.Cell[0], [0, 0, 0])

    def test_tilted_box(self):
        """Tilted box [[2,1,0],[0,2,0],[0,0,1]] → NCell=4."""
        StdI = self._make_stdi_for_enum([[2, 1, 0], [0, 2, 0], [0, 0, 1]], 4)
        _enumerate_cells(StdI)
        assert StdI.Cell.shape == (4, 3)
        # All 4 cells should be distinct
        rows = set(tuple(r) for r in StdI.Cell.tolist())
        assert len(rows) == 4


# ===================================================================
#  find_site
# ===================================================================


class TestFindSite:
    """Tests for find_site."""

    def test_same_site(self):
        """Test that zero displacement gives same site."""
        StdI = _make_stdi_chain(4)
        isite, jsite, Cphase, dR = find_site(StdI, 0, 0, 0, 0, 0, 0, 0, 0)
        assert isite == 0
        assert jsite == 0
        assert Cphase == pytest.approx(1.0 + 0j)

    def test_nearest_neighbour(self):
        """Test nearest-neighbour displacement."""
        StdI = _make_stdi_chain(4)
        isite, jsite, Cphase, dR = find_site(StdI, 0, 0, 0, 1, 0, 0, 0, 0)
        assert isite == 0
        assert jsite == 1
        assert Cphase == pytest.approx(1.0 + 0j)

    def test_boundary_wraps(self):
        """Test that displacement beyond boundary wraps with phase=1."""
        StdI = _make_stdi_chain(4)
        # Site 3, displacement +1 → wraps to site 0
        isite, jsite, Cphase, dR = find_site(StdI, 3, 0, 0, 1, 0, 0, 0, 0)
        assert isite == 3
        assert jsite == 0
        assert Cphase == pytest.approx(1.0 + 0j)

    def test_kondo_offset(self):
        """Test that kondo model adds NCell*NsiteUC offset."""
        StdI = _make_stdi_chain(4)
        StdI.model = "kondo"
        isite, jsite, Cphase, dR = find_site(StdI, 0, 0, 0, 1, 0, 0, 0, 0)
        # kondo adds NCell*NsiteUC = 4*1 = 4
        assert isite == 4
        assert jsite == 5

    def test_dr_vector(self):
        """Test that dR vector is computed correctly."""
        StdI = _make_stdi_chain(4)
        _, _, _, dR = find_site(StdI, 0, 0, 0, 1, 0, 0, 0, 0)
        assert dR[0] == pytest.approx(-1.0)
        assert dR[1] == pytest.approx(0.0)
        assert dR[2] == pytest.approx(0.0)

    def test_anti_periodic_phase(self):
        """Test that anti-periodic boundary gives Cphase=-1."""
        StdI = _make_stdi_chain(4)
        StdI.ExpPhase[0] = -1.0 + 0j  # anti-periodic in W direction
        StdI.AntiPeriod[0] = 1
        # Site 3, displacement +1 → wraps across boundary
        _, _, Cphase, _ = find_site(StdI, 3, 0, 0, 1, 0, 0, 0, 0)
        assert Cphase == pytest.approx(-1.0 + 0j)


# ===================================================================
#  _write_gnuplot_bond
# ===================================================================


class TestWriteGnuplotBond:
    """Tests for _write_gnuplot_bond helper."""

    def test_single_digit_sites(self):
        """Labels use 1-digit format for sites < 10."""
        fp = io.StringIO()
        _write_gnuplot_bond(fp, 3, 7, 1.0, 2.0, 3.0, 4.0, 1)
        content = fp.getvalue()
        assert 'set label "3"' in content
        assert 'set label "7"' in content

    def test_double_digit_sites(self):
        """Labels use 2-digit format for sites >= 10."""
        fp = io.StringIO()
        _write_gnuplot_bond(fp, 12, 25, 0.0, 0.0, 1.0, 1.0, 1)
        content = fp.getvalue()
        assert 'set label "12"' in content
        assert 'set label "25"' in content

    def test_arrow_written_for_connect_lt_3(self):
        """Arrow command is written when connect < 3."""
        for connect in (1, 2):
            fp = io.StringIO()
            _write_gnuplot_bond(fp, 0, 1, 0.0, 0.0, 1.0, 1.0, connect)
            content = fp.getvalue()
            assert "set arrow" in content
            assert f"nohead ls {connect}" in content

    def test_no_arrow_for_connect_ge_3(self):
        """No arrow command when connect >= 3."""
        fp = io.StringIO()
        _write_gnuplot_bond(fp, 0, 1, 0.0, 0.0, 1.0, 1.0, 3)
        content = fp.getvalue()
        assert "set arrow" not in content

    def test_positions_in_output(self):
        """Site positions appear in the label commands."""
        fp = io.StringIO()
        _write_gnuplot_bond(fp, 0, 1, 1.5, 2.5, 3.5, 4.5, 1)
        content = fp.getvalue()
        assert "1.500000" in content
        assert "2.500000" in content
        assert "3.500000" in content
        assert "4.500000" in content

    def test_exactly_three_lines_with_arrow(self):
        """Output has exactly 3 non-empty lines: 2 labels + 1 arrow."""
        fp = io.StringIO()
        _write_gnuplot_bond(fp, 0, 1, 0.0, 0.0, 1.0, 1.0, 1)
        lines = [l for l in fp.getvalue().splitlines() if l.strip()]
        assert len(lines) == 3

    def test_exactly_two_lines_without_arrow(self):
        """Output has exactly 2 non-empty lines when connect >= 3."""
        fp = io.StringIO()
        _write_gnuplot_bond(fp, 0, 1, 0.0, 0.0, 1.0, 1.0, 3)
        lines = [l for l in fp.getvalue().splitlines() if l.strip()]
        assert len(lines) == 2

    def test_mixed_digit_sites(self):
        """One site < 10, one >= 10."""
        fp = io.StringIO()
        _write_gnuplot_bond(fp, 5, 15, 0.0, 0.0, 1.0, 1.0, 1)
        content = fp.getvalue()
        assert 'set label "5"' in content
        assert 'set label "15"' in content


# ===================================================================
#  set_label
# ===================================================================


class TestSetLabel:
    """Tests for set_label."""

    def test_returns_site_indices(self):
        """Test that set_label returns correct site indices."""
        StdI = _make_stdi_chain(4)
        isite, jsite, Cphase, dR = set_label(
            StdI, None, 0, 0, 1, 0, 0, 0, 1)
        assert isite == 0
        assert jsite == 1

    def test_writes_gnuplot_labels(self):
        """Test that gnuplot labels are written when fp is provided."""
        StdI = _make_stdi_chain(4)
        fp = io.StringIO()
        set_label(StdI, fp, 0, 0, 1, 0, 0, 0, 1)
        content = fp.getvalue()
        assert "set label" in content
        assert "set arrow" in content

    def test_no_arrow_for_connect_ge_3(self):
        """Test that connect >= 3 suppresses arrow output."""
        StdI = _make_stdi_chain(4)
        fp = io.StringIO()
        set_label(StdI, fp, 0, 0, 1, 0, 0, 0, 3)
        content = fp.getvalue()
        assert "set label" in content
        assert "set arrow" not in content

    def test_fp_none_no_error(self):
        """Test that fp=None doesn't cause an error."""
        StdI = _make_stdi_chain(4)
        isite, jsite, Cphase, dR = set_label(
            StdI, None, 0, 0, 1, 0, 0, 0, 1)
        # Should complete without error
        assert isinstance(isite, (int, np.integer))


# ===================================================================
#  Backward compatibility
# ===================================================================


class TestBackwardCompatibility:
    """Test that functions are importable from stdface.core.stdface_model_util."""

    def test_import_from_stdface_model_util(self):
        """Test that all 4 functions are re-exported."""
        from stdface.core.stdface_model_util import (
            _fold_site as fs,
            init_site as iis,
            find_site as fis,
            set_label as sl,
        )
        from stdface.lattice.site_util import _fold_site, init_site, find_site, set_label
        assert fs is _fold_site
        assert iis is init_site
        assert fis is find_site
        assert sl is set_label


# ===================================================================
#  _find_cell_index
# ===================================================================


class TestFindCellIndex:
    """Tests for _find_cell_index helper."""

    def test_finds_first_cell(self):
        """Test finding cell index 0 in a 1D chain."""
        StdI = _make_stdi_chain(4)
        assert _find_cell_index(StdI, [0, 0, 0]) == 0

    def test_finds_middle_cell(self):
        """Test finding cell index 2 in a 1D chain."""
        StdI = _make_stdi_chain(4)
        assert _find_cell_index(StdI, [2, 0, 0]) == 2

    def test_finds_last_cell(self):
        """Test finding cell index L-1 in a 1D chain."""
        StdI = _make_stdi_chain(6)
        assert _find_cell_index(StdI, [5, 0, 0]) == 5

    def test_all_cells_found(self):
        """Test that every cell in the lattice can be found."""
        L = 8
        StdI = _make_stdi_chain(L)
        for i in range(L):
            assert _find_cell_index(StdI, [i, 0, 0]) == i

    def test_no_match_returns_zero(self):
        """Test that a non-existent coordinate returns 0."""
        StdI = _make_stdi_chain(4)
        assert _find_cell_index(StdI, [99, 0, 0]) == 0

    def test_2d_square(self):
        """Test finding cells in a 2D square lattice."""
        StdI = StdIntList()
        StdI.NCell = 4
        StdI.Cell = np.array([
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
        ], dtype=int)
        assert _find_cell_index(StdI, [0, 1, 0]) == 2
        assert _find_cell_index(StdI, [1, 1, 0]) == 3

    def test_3d_cell(self):
        """Test finding a cell with nonzero H component."""
        StdI = StdIntList()
        StdI.NCell = 2
        StdI.Cell = np.array([
            [0, 0, 0],
            [0, 0, 1],
        ], dtype=int)
        assert _find_cell_index(StdI, [0, 0, 1]) == 1

    def test_early_return_on_match(self):
        """Test that first match is returned (not last)."""
        StdI = StdIntList()
        StdI.NCell = 3
        StdI.Cell = np.array([
            [1, 0, 0],
            [0, 0, 0],
            [1, 0, 0],  # duplicate of cell 0
        ], dtype=int)
        # Should return first match (index 0), not last (index 2)
        assert _find_cell_index(StdI, [1, 0, 0]) == 0


# ===================================================================
#  set_local_spin_flags
# ===================================================================


class TestSetLocalSpinFlags:
    """Tests for set_local_spin_flags."""

    @staticmethod
    def _make_stdi(model: ModelType, S2: int = 1) -> StdIntList:
        """Create a minimal StdIntList for set_local_spin_flags tests."""
        StdI = StdIntList()
        StdI.model = model
        StdI.S2 = S2
        return StdI

    def test_spin_model(self):
        """Test that SPIN model sets all flags to S2."""
        StdI = self._make_stdi(ModelType.SPIN, S2=3)
        set_local_spin_flags(StdI, 6)
        assert StdI.nsite == 6
        assert len(StdI.locspinflag) == 6
        assert all(StdI.locspinflag == 3)

    def test_hubbard_model(self):
        """Test that HUBBARD model sets all flags to 0."""
        StdI = self._make_stdi(ModelType.HUBBARD)
        set_local_spin_flags(StdI, 8)
        assert StdI.nsite == 8
        assert len(StdI.locspinflag) == 8
        assert all(StdI.locspinflag == 0)

    def test_kondo_model(self):
        """Test that KONDO model doubles nsite and splits flags."""
        StdI = self._make_stdi(ModelType.KONDO, S2=1)
        set_local_spin_flags(StdI, 4)
        assert StdI.nsite == 8
        assert len(StdI.locspinflag) == 8
        # First half: S2, second half: 0
        np.testing.assert_array_equal(StdI.locspinflag[:4], 1)
        np.testing.assert_array_equal(StdI.locspinflag[4:], 0)

    def test_kondo_with_s2_3(self):
        """Test KONDO with S2=3 (spin-3/2)."""
        StdI = self._make_stdi(ModelType.KONDO, S2=3)
        set_local_spin_flags(StdI, 6)
        assert StdI.nsite == 12
        np.testing.assert_array_equal(StdI.locspinflag[:6], 3)
        np.testing.assert_array_equal(StdI.locspinflag[6:], 0)

    def test_single_site(self):
        """Test with a single base site."""
        StdI = self._make_stdi(ModelType.SPIN, S2=1)
        set_local_spin_flags(StdI, 1)
        assert StdI.nsite == 1
        assert StdI.locspinflag[0] == 1

    def test_kondo_single_site(self):
        """Test KONDO with a single base site doubles to 2."""
        StdI = self._make_stdi(ModelType.KONDO, S2=1)
        set_local_spin_flags(StdI, 1)
        assert StdI.nsite == 2
        assert StdI.locspinflag[0] == 1
        assert StdI.locspinflag[1] == 0


# ===================================================================
#  _write_gnuplot_header
# ===================================================================


class TestWriteGnuplotHeader:
    """Tests for _write_gnuplot_header extracted from init_site."""

    @staticmethod
    def _make_stdi_2d(W: int = 2, L: int = 2) -> StdIntList:
        """Create a minimal StdIntList with 2D direct/box for gnuplot tests."""
        StdI = StdIntList()
        StdI.direct[:, :] = 0.0
        StdI.direct[0, 0] = 1.0
        StdI.direct[1, 1] = 1.0
        StdI.direct[2, 2] = 1.0
        StdI.box[:, :] = 0
        StdI.box[0, 0] = W
        StdI.box[1, 1] = L
        StdI.box[2, 2] = 1
        return StdI

    def test_writes_xrange(self):
        """Test that xrange line is present."""
        StdI = self._make_stdi_2d()
        fp = io.StringIO()
        _write_gnuplot_header(fp, StdI)
        assert "set xrange" in fp.getvalue()

    def test_writes_yrange(self):
        """Test that yrange line is present."""
        StdI = self._make_stdi_2d()
        fp = io.StringIO()
        _write_gnuplot_header(fp, StdI)
        assert "set yrange" in fp.getvalue()

    def test_writes_four_arrows(self):
        """Test that exactly four boundary arrows are written."""
        StdI = self._make_stdi_2d()
        fp = io.StringIO()
        _write_gnuplot_header(fp, StdI)
        content = fp.getvalue()
        arrow_lines = [l for l in content.splitlines() if l.startswith("set arrow")]
        assert len(arrow_lines) == 4

    def test_writes_style_lines(self):
        """Test that gnuplot style lines are present."""
        StdI = self._make_stdi_2d()
        fp = io.StringIO()
        _write_gnuplot_header(fp, StdI)
        content = fp.getvalue()
        assert "set style line 1" in content
        assert "set style line 2" in content
        assert "set style line 3" in content

    def test_corner_positions_identity(self):
        """Test corner positions for W=2, L=3 with identity direct."""
        StdI = self._make_stdi_2d(W=2, L=3)
        fp = io.StringIO()
        _write_gnuplot_header(fp, StdI)
        content = fp.getvalue()
        # pos[0] = (0,0), pos[1] = (W,0)=(2,0), pos[2] = (0,L)=(0,3), pos[3]=(2,3)
        # Arrow 0→1: from 0,0 to 2,0
        assert "from 0.000000, 0.000000 to 2.000000, 0.000000" in content
        # Arrow 1→3: from 2,0 to 2,3
        assert "from 2.000000, 0.000000 to 2.000000, 3.000000" in content

    def test_size_square_and_unset(self):
        """Test that layout commands are present."""
        StdI = self._make_stdi_2d()
        fp = io.StringIO()
        _write_gnuplot_header(fp, StdI)
        content = fp.getvalue()
        assert "set size square" in content
        assert "unset key" in content
        assert "unset tics" in content
        assert "unset border" in content


# ===================================================================
#  lattice_gp context manager
# ===================================================================


class TestLatticeGp:
    """Tests for the lattice_gp context manager."""

    def test_yields_file_for_hphi(self, tmp_path, monkeypatch):
        """Test that a file handle is yielded for HPhi solver."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.nsite = 4
        StdI.locspinflag = np.zeros(4, dtype=int)
        with lattice_gp(StdI) as fp:
            assert fp is not None

    def test_yields_file_for_mvmc(self, tmp_path, monkeypatch):
        """Test that a file handle is yielded for mVMC solver."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.solver = SolverType.mVMC
        StdI.nsite = 4
        StdI.locspinflag = np.zeros(4, dtype=int)
        with lattice_gp(StdI) as fp:
            assert fp is not None

    def test_yields_none_for_hwave(self, tmp_path, monkeypatch):
        """Test that None is yielded for H-wave when lattice_gp=0."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.solver = SolverType.HWAVE
        StdI.lattice_gp = 0
        StdI.nsite = 4
        StdI.locspinflag = np.zeros(4, dtype=int)
        with lattice_gp(StdI) as fp:
            assert fp is None

    def test_yields_file_for_hwave_with_lattice_gp(self, tmp_path, monkeypatch):
        """Test that a file handle is yielded for H-wave when lattice_gp=1."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.solver = SolverType.HWAVE
        StdI.lattice_gp = 1
        StdI.nsite = 4
        StdI.locspinflag = np.zeros(4, dtype=int)
        with lattice_gp(StdI) as fp:
            assert fp is not None

    def test_creates_lattice_gp_file(self, tmp_path, monkeypatch):
        """Test that lattice.gp is created on disk."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.nsite = 4
        StdI.locspinflag = np.zeros(4, dtype=int)
        with lattice_gp(StdI) as fp:
            pass
        assert (tmp_path / "lattice.gp").exists()

    def test_writes_footer_on_exit(self, tmp_path, monkeypatch):
        """Test that footer is written when the context manager exits."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.nsite = 4
        StdI.locspinflag = np.zeros(4, dtype=int)
        with lattice_gp(StdI) as fp:
            fp.write("# test header\n")
        content = (tmp_path / "lattice.gp").read_text()
        assert "# test header" in content
        assert "plot '-' w d lc 7" in content
        assert "pause -1" in content

    def test_none_fp_no_error(self, tmp_path, monkeypatch):
        """Test that None fp context manager exits without error."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.solver = SolverType.HWAVE
        StdI.lattice_gp = 0
        StdI.nsite = 4
        StdI.locspinflag = np.zeros(4, dtype=int)
        # Should not raise
        with lattice_gp(StdI) as fp:
            assert fp is None

    def test_calls_print_geometry(self, tmp_path, monkeypatch):
        """Test that geometry.dat is created (proves print_geometry was called)."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.solver = SolverType.HWAVE
        StdI.lattice_gp = 0
        StdI.nsite = 4
        StdI.locspinflag = np.zeros(4, dtype=int)
        with lattice_gp(StdI) as fp:
            pass
        assert (tmp_path / "geometry.dat").exists()

    def test_footer_content_matches_constant(self):
        """Test that _LATTICE_GP_FOOTER has expected content."""
        assert "plot '-' w d lc 7" in _LATTICE_GP_FOOTER
        assert "0.0 0.0" in _LATTICE_GP_FOOTER
        assert "end" in _LATTICE_GP_FOOTER
        assert "pause -1" in _LATTICE_GP_FOOTER


# ===================================================================
#  close_lattice_xsf
# ===================================================================


class TestCloseLatticeXsf:
    """Tests for close_lattice_xsf (3D counterpart of close_lattice_gp)."""

    def test_creates_lattice_xsf(self, tmp_path, monkeypatch):
        """Test that lattice.xsf is created on disk."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.nsite = 4
        StdI.NsiteUC = 1
        StdI.locspinflag = np.zeros(4, dtype=int)
        close_lattice_xsf(StdI)
        assert (tmp_path / "lattice.xsf").exists()

    def test_creates_geometry_dat(self, tmp_path, monkeypatch):
        """Test that geometry.dat is created on disk."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.nsite = 4
        StdI.NsiteUC = 1
        StdI.locspinflag = np.zeros(4, dtype=int)
        close_lattice_xsf(StdI)
        assert (tmp_path / "geometry.dat").exists()

    def test_lattice_xsf_contains_crystal(self, tmp_path, monkeypatch):
        """Test that lattice.xsf has XCrySDen CRYSTAL header."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.nsite = 4
        StdI.NsiteUC = 1
        StdI.locspinflag = np.zeros(4, dtype=int)
        close_lattice_xsf(StdI)
        content = (tmp_path / "lattice.xsf").read_text()
        assert "CRYSTAL" in content

    def test_geometry_dat_has_content(self, tmp_path, monkeypatch):
        """Test that geometry.dat is non-empty."""
        monkeypatch.chdir(tmp_path)
        StdI = _make_stdi_chain(4)
        StdI.nsite = 4
        StdI.NsiteUC = 1
        StdI.locspinflag = np.zeros(4, dtype=int)
        close_lattice_xsf(StdI)
        content = (tmp_path / "geometry.dat").read_text()
        assert len(content) > 0


