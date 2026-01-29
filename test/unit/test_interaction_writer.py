"""Unit tests for the interaction_writer module.

Tests for the :func:`print_interactions` function which was extracted from
``common_writer.py`` into its own module.
"""
from __future__ import annotations

import os
import tempfile

from stdface_vals import StdIntList
from writer.interaction_writer import (
    print_interactions,
    _merge_1idx,
    _merge_2idx,
    _count_nonzero,
    _process_interaction,
    _INTERACTION_TYPES,
    _merge_interall_equivalent,
    _reorder_interall_hermitian,
    _remove_interall_diagonal,
    _write_interall,
    _write_interaction_file,
)


# Sentinel values matching what _reset_vals sets at runtime
NaN_i = 2147483647


def _make_stdi_base(**overrides) -> StdIntList:
    """Create a minimal StdIntList with common fields."""
    StdI = StdIntList()
    StdI.nsite = overrides.get("nsite", 4)
    StdI.solver = overrides.get("solver", "HPhi")
    StdI.model = overrides.get("model", "hubbard")
    StdI.lGC = overrides.get("lGC", 0)
    StdI.lBoost = overrides.get("lBoost", 0)
    return StdI


def _make_empty_interactions(StdI: StdIntList) -> None:
    """Initialise all interaction arrays to empty."""
    StdI.NCintra = 0
    StdI.CintraIndx = []
    StdI.Cintra = []
    StdI.NCinter = 0
    StdI.CinterIndx = []
    StdI.Cinter = []
    StdI.NHund = 0
    StdI.HundIndx = []
    StdI.Hund = []
    StdI.NEx = 0
    StdI.ExIndx = []
    StdI.Ex = []
    StdI.NPairLift = 0
    StdI.PLIndx = []
    StdI.PairLift = []
    StdI.NPairHopp = 0
    StdI.PHIndx = []
    StdI.PairHopp = []
    StdI.nintr = 0
    StdI.intrindx = []
    StdI.intr = []


class TestImportFromInteractionWriter:
    """Verify the module can be imported directly."""

    def test_import_print_interactions(self):
        assert callable(print_interactions)


class TestMergeDuplicateCoulombIntra:
    """Test that duplicate CoulombIntra terms are merged."""

    def test_duplicate_site_merged(self):
        StdI = _make_stdi_base()
        _make_empty_interactions(StdI)
        StdI.NCintra = 3
        StdI.CintraIndx = [[0], [1], [0]]  # site 0 appears twice
        StdI.Cintra = [2.0, 3.0, 5.0]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_interactions(StdI)
                assert StdI.Cintra[0] == 7.0  # merged: 2.0 + 5.0
                assert StdI.Cintra[2] == 0.0  # zeroed
                assert StdI.LCintra == 1
            finally:
                os.chdir(orig)


class TestMergeDuplicateCoulombInter:
    """Test that duplicate CoulombInter terms are merged (symmetric)."""

    def test_symmetric_pair_merged(self):
        StdI = _make_stdi_base()
        _make_empty_interactions(StdI)
        StdI.NCinter = 2
        StdI.CinterIndx = [[0, 1], [1, 0]]  # same pair, reversed
        StdI.Cinter = [1.5, 2.5]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_interactions(StdI)
                assert StdI.Cinter[0] == 4.0  # merged
                assert StdI.Cinter[1] == 0.0
                assert StdI.LCinter == 1
            finally:
                os.chdir(orig)


class TestBoostSuppressesOutput:
    """Test that lBoost=1 sets all flags to 0."""

    def test_boost_suppresses_coulomb_intra(self):
        StdI = _make_stdi_base(lBoost=1)
        _make_empty_interactions(StdI)
        StdI.NCintra = 1
        StdI.CintraIndx = [[0]]
        StdI.Cintra = [4.0]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_interactions(StdI)
                assert StdI.LCintra == 0
                assert not os.path.exists("coulombintra.def")
            finally:
                os.chdir(orig)


class TestExchangeWritten:
    """Test that exchange.def is written for non-zero Exchange terms."""

    def test_exchange_file_created(self):
        StdI = _make_stdi_base()
        _make_empty_interactions(StdI)
        StdI.NEx = 1
        StdI.ExIndx = [[0, 1]]
        StdI.Ex = [0.5]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_interactions(StdI)
                assert StdI.LEx == 1
                assert os.path.exists("exchange.def")
                content = open("exchange.def").read()
                assert "NExchange          1" in content
            finally:
                os.chdir(orig)


class TestInterAllWritten:
    """Test that interall.def is written for non-zero InterAll terms."""

    def test_single_interall_term(self):
        StdI = _make_stdi_base()
        _make_empty_interactions(StdI)
        StdI.nintr = 1
        StdI.intrindx = [[0, 0, 0, 0, 0, 0, 0, 0]]
        StdI.intr = [1.0 + 0j]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_interactions(StdI)
                assert StdI.Lintr == 1
                assert os.path.exists("interall.def")
                content = open("interall.def").read()
                assert "NInterAll       1" in content
            finally:
                os.chdir(orig)


class TestReExportFromCommonWriter:
    """Verify backward-compatible re-export from common_writer."""

    def test_import_from_common_writer_still_works(self):
        from writer.common_writer import print_interactions as pi
        assert pi is print_interactions


class TestMerge1Idx:
    """Tests for the _merge_1idx helper."""

    def test_no_duplicates(self):
        coeff = [1.0, 2.0, 3.0]
        indx = [[0], [1], [2]]
        _merge_1idx(3, indx, coeff)
        assert coeff == [1.0, 2.0, 3.0]

    def test_duplicate_merged(self):
        coeff = [1.0, 2.0, 3.0]
        indx = [[0], [1], [0]]
        _merge_1idx(3, indx, coeff)
        assert coeff[0] == 4.0
        assert coeff[2] == 0.0

    def test_triple_duplicate(self):
        coeff = [1.0, 2.0, 3.0]
        indx = [[5], [5], [5]]
        _merge_1idx(3, indx, coeff)
        assert coeff[0] == 6.0
        assert coeff[1] == 0.0
        assert coeff[2] == 0.0

    def test_empty(self):
        coeff = []
        indx = []
        _merge_1idx(0, indx, coeff)
        assert coeff == []


class TestMerge2Idx:
    """Tests for the _merge_2idx helper."""

    def test_no_duplicates(self):
        coeff = [1.0, 2.0]
        indx = [[0, 1], [2, 3]]
        _merge_2idx(2, indx, coeff)
        assert coeff == [1.0, 2.0]

    def test_same_order_merged(self):
        coeff = [1.0, 2.0]
        indx = [[0, 1], [0, 1]]
        _merge_2idx(2, indx, coeff)
        assert coeff[0] == 3.0
        assert coeff[1] == 0.0

    def test_reversed_order_merged(self):
        coeff = [1.0, 2.0]
        indx = [[0, 1], [1, 0]]
        _merge_2idx(2, indx, coeff)
        assert coeff[0] == 3.0
        assert coeff[1] == 0.0

    def test_empty(self):
        coeff = []
        indx = []
        _merge_2idx(0, indx, coeff)
        assert coeff == []


class TestCountNonzero:
    """Tests for the _count_nonzero helper."""

    def test_all_nonzero(self):
        assert _count_nonzero(3, [1.0, 2.0, 3.0]) == 3

    def test_some_zero(self):
        assert _count_nonzero(4, [1.0, 0.0, 3.0, 0.0]) == 2

    def test_all_below_threshold(self):
        assert _count_nonzero(2, [0.0000001, 0.0000005]) == 0

    def test_at_threshold(self):
        assert _count_nonzero(1, [0.0000011]) == 1

    def test_empty(self):
        assert _count_nonzero(0, []) == 0


class TestInteractionTypesTable:
    """Tests for the _INTERACTION_TYPES metadata table."""

    def test_six_entries(self):
        assert len(_INTERACTION_TYPES) == 6

    def test_coulomb_intra_is_1idx(self):
        assert _INTERACTION_TYPES[0].n_indices == 1
        assert _INTERACTION_TYPES[0].filename == "coulombintra.def"

    def test_remaining_five_are_2idx(self):
        for spec in _INTERACTION_TYPES[1:]:
            assert spec.n_indices == 2

    def test_all_entries_have_required_fields(self):
        required = {"nterms_attr", "indx_attr", "coeff_attr", "flag_attr",
                     "filename", "count_label", "banner", "n_indices"}
        for spec in _INTERACTION_TYPES:
            assert set(spec._fields) == required


class TestMergeInterallEquivalent:
    """Tests for _merge_interall_equivalent (Pass 1)."""

    def test_empty(self):
        intr = []
        intrindx = []
        _merge_interall_equivalent(0, intrindx, intr)
        assert intr == []

    def test_no_match(self):
        intr = [1.0 + 0j, 2.0 + 0j]
        intrindx = [
            [0, 0, 1, 0, 2, 0, 3, 0],
            [4, 0, 5, 0, 6, 0, 7, 0],
        ]
        _merge_interall_equivalent(2, intrindx, intr)
        assert intr[0] == 1.0 + 0j
        assert intr[1] == 2.0 + 0j

    def test_case_a_exact_match(self):
        """All 8 indices identical → coefficients added."""
        intr = [1.0 + 0j, 2.0 + 0j]
        idx = [0, 0, 1, 0, 2, 0, 3, 0]
        intrindx = [list(idx), list(idx)]
        _merge_interall_equivalent(2, intrindx, intr)
        assert intr[0] == 3.0 + 0j
        assert intr[1] == 0.0

    def test_case_b_hermitian_conjugate_add(self):
        """Indices 0-3 <-> 4-7 swap (with exclusion conditions met) → add."""
        # j = [0,0, 1,0, 2,0, 3,0]
        # k = [2,0, 3,0, 0,0, 1,0]  (swapped halves)
        # Exclusion 1: not (j[0]==j[6] and j[1]==j[7]) → 0!=3 or 0!=0 → True
        # Exclusion 2: not (j[2]==j[4] and j[3]==j[5]) → 1!=2 or 0!=0 → True
        intr = [1.0 + 0j, 4.0 + 0j]
        intrindx = [
            [0, 0, 1, 0, 2, 0, 3, 0],
            [2, 0, 3, 0, 0, 0, 1, 0],
        ]
        _merge_interall_equivalent(2, intrindx, intr)
        assert intr[0] == 5.0 + 0j
        assert intr[1] == 0.0

    def test_case_c_partial_swap_subtract(self):
        """Partial swap pattern → subtract."""
        # j = [0,0, 1,1, 2,0, 3,0]
        # k = [2,0, 3,0, 1,1, 0,0]  (j[4-5]↔k[0-1], j[0-1]↔k[4-5])
        # Pattern: j[0]==k[4] (0==0), j[1]==k[5] (0==0),
        #          j[2]==k[2] (1==1), j[3]==k[3] (1==1),
        #          j[4]==k[0] (2==2), j[5]==k[1] (0==0),
        #          j[6]==k[6] (3==3), j[7]==k[7] (0==0)
        # Exclusion: not(j[2]==j[0] and j[3]==j[1]) → 1!=0 → True
        #            not(j[2]==j[4] and j[3]==j[5]) → 1!=2 → True
        intr = [5.0 + 0j, 2.0 + 0j]
        intrindx = [
            [0, 0, 1, 1, 2, 0, 3, 0],
            [2, 0, 1, 1, 0, 0, 3, 0],
        ]
        _merge_interall_equivalent(2, intrindx, intr)
        assert intr[0] == 3.0 + 0j  # 5 - 2
        assert intr[1] == 0.0

    def test_three_terms_two_duplicates(self):
        """Three terms, first and third are exact duplicates."""
        idx = [0, 0, 1, 0, 2, 0, 3, 0]
        intr = [1.0 + 0j, 10.0 + 0j, 2.0 + 0j]
        intrindx = [
            list(idx),
            [4, 0, 5, 0, 6, 0, 7, 0],
            list(idx),
        ]
        _merge_interall_equivalent(3, intrindx, intr)
        assert intr[0] == 3.0 + 0j
        assert intr[1] == 10.0 + 0j
        assert intr[2] == 0.0


class TestReorderInterallHermitian:
    """Tests for _reorder_interall_hermitian (Pass 2)."""

    def test_empty(self):
        _reorder_interall_hermitian(0, [], [])

    def test_no_reorder_needed(self):
        intr = [1.0 + 0j, 2.0 + 0j]
        intrindx = [
            [0, 0, 1, 0, 2, 0, 3, 0],
            [4, 0, 5, 0, 6, 0, 7, 0],
        ]
        orig_k = list(intrindx[1])
        _reorder_interall_hermitian(2, intrindx, intr)
        assert intrindx[1] == orig_k
        assert intr[1] == 2.0 + 0j

    def test_pattern1_direct_hermitian(self):
        """Pattern 1: direct Hermitian conjugate → reorder indices."""
        # j = [0,0, 1,1, 2,0, 3,1]
        # For k to match pattern 1:
        # j[6]==k[4], j[7]==k[5], j[4]==k[6], j[5]==k[7],
        # j[2]==k[0], j[3]==k[1], j[0]==k[2], j[1]==k[3]
        # So k = [1,1, 0,0, 3,1, 2,0]
        # Exclusion: not(k[0]==k[6] and k[1]==k[7]) → 1!=2 → True
        #            not(k[2]==k[4] and k[3]==k[5]) → 0!=3 → True
        intr = [1.0 + 0j, 2.0 + 0j]
        intrindx = [
            [0, 0, 1, 1, 2, 0, 3, 1],
            [1, 1, 0, 0, 3, 1, 2, 0],
        ]
        _reorder_interall_hermitian(2, intrindx, intr)
        # After reorder: k becomes [j[6],j[7], j[4],j[5], j[2],j[3], j[0],j[1]]
        # = [3,1, 2,0, 1,1, 0,0]
        assert intrindx[1] == [3, 1, 2, 0, 1, 1, 0, 0]
        assert intr[1] == 2.0 + 0j  # no sign change

    def test_pattern2_sign_flip(self):
        """Pattern 2: reorder with sign flip."""
        # j = [0,0, 1,1, 2,0, 3,1]
        # For pattern 2 sub-case 1:
        # j[6]==k[4], j[7]==k[5] → k[4]=3, k[5]=1
        # j[4]==k[2], j[5]==k[3] → k[2]=2, k[3]=0
        # j[2]==k[0], j[3]==k[1] → k[0]=1, k[1]=1
        # j[0]==k[6], j[1]==k[7] → k[6]=0, k[7]=0
        # So k = [1,1, ?, 2,0, 3,1, 0,0]  wait, need [k0..k7]
        # k = [1,1, 2,0, 3,1, 0,0]
        # Exclusion: not(k[2]==k[0] and k[3]==k[1]) → 2!=1 → True
        #            not(k[2]==k[4] and k[3]==k[5]) → 2!=3 → True
        intr = [1.0 + 0j, 2.0 + 0j]
        intrindx = [
            [0, 0, 1, 1, 2, 0, 3, 1],
            [1, 1, 2, 0, 3, 1, 0, 0],
        ]
        _reorder_interall_hermitian(2, intrindx, intr)
        assert intrindx[1] == [3, 1, 2, 0, 1, 1, 0, 0]
        assert intr[1] == -2.0 + 0j  # sign flipped


class TestRemoveInterallDiagonal:
    """Tests for _remove_interall_diagonal (Pass 3)."""

    def test_empty(self):
        _remove_interall_diagonal(0, [], [])

    def test_no_diagonal_pair_kept(self):
        """Term without a diagonal pair is kept."""
        intr = [1.0 + 0j]
        intrindx = [[0, 0, 1, 0, 2, 0, 3, 0]]
        _remove_interall_diagonal(1, intrindx, intr)
        assert intr[0] == 1.0 + 0j

    def test_diagonal_pair_04_removed(self):
        """(site0,spin0)==(site4,spin4) with no matching diagonal → zeroed."""
        # idx[0]==idx[4] and idx[1]==idx[5] → diagonal pair
        # None of the four matching-diagonal conditions hold
        intr = [1.0 + 0j]
        intrindx = [[0, 0, 1, 0, 0, 0, 2, 0]]
        # diagonal: (0,0)==(0,0) → True
        # matching diag: (0,0)==(1,0)? No. (0,0)==(2,0)? No.
        #                (0,0)==(1,0)? No. (0,0)==(2,0)? No.
        _remove_interall_diagonal(1, intrindx, intr)
        assert intr[0] == 0.0

    def test_diagonal_pair_26_removed(self):
        """(site2,spin2)==(site6,spin6) with no matching diagonal → zeroed."""
        intr = [1.0 + 0j]
        intrindx = [[0, 0, 1, 0, 2, 0, 1, 0]]
        # diagonal pair: idx[2]==idx[6] (1==1) and idx[3]==idx[7] (0==0)
        # matching diag: (0,0)==(1,0)? No. (0,0)==(1,0)? No.
        #                (2,0)==(1,0)? No. (2,0)==(1,0)? No.
        _remove_interall_diagonal(1, intrindx, intr)
        assert intr[0] == 0.0

    def test_diagonal_pair_with_matching_diagonal_kept(self):
        """Diagonal pair present, but matching diagonal → kept."""
        intr = [1.0 + 0j]
        # Make idx[0]==idx[4] (diagonal pair) AND idx[0]==idx[2] (matching diag)
        intrindx = [[0, 0, 0, 0, 0, 0, 1, 0]]
        _remove_interall_diagonal(1, intrindx, intr)
        assert intr[0] == 1.0 + 0j  # kept

    def test_mixed_terms(self):
        """Multiple terms: one removed, one kept."""
        intr = [1.0 + 0j, 2.0 + 0j]
        intrindx = [
            [0, 0, 1, 0, 0, 0, 2, 0],  # diagonal pair (0,0)==(0,0), no match → zero
            [0, 0, 1, 0, 2, 0, 3, 0],  # no diagonal pair → kept
        ]
        _remove_interall_diagonal(2, intrindx, intr)
        assert intr[0] == 0.0
        assert intr[1] == 2.0 + 0j


class TestWriteInteractionFile:
    """Tests for _write_interaction_file (unified 1-idx / 2-idx writer)."""

    def test_1idx_basic(self, tmp_path, monkeypatch):
        """Write a 1-index file with mixed zero/non-zero terms."""
        monkeypatch.chdir(tmp_path)
        indx = [[0], [1], [2]]
        coeff = [2.5, 0.0, 3.0]
        _write_interaction_file(
            "coulombintra.def", "NCoulombIntra",
            "================== CoulombIntra ================",
            3, indx, coeff, n_indices=1,
        )
        content = (tmp_path / "coulombintra.def").read_text()
        assert "NCoulombIntra          2" in content
        assert "2.500000000000000" in content
        assert "3.000000000000000" in content
        # The zero term should NOT appear in output
        lines = [l for l in content.splitlines() if "." in l and "===" not in l]
        assert len(lines) == 2

    def test_2idx_basic(self, tmp_path, monkeypatch):
        """Write a 2-index file with one non-zero term."""
        monkeypatch.chdir(tmp_path)
        indx = [[0, 1], [2, 3]]
        coeff = [1.5, 0.0]
        _write_interaction_file(
            "exchange.def", "NExchange",
            "====== ExchangeCoupling coupling ============",
            2, indx, coeff, n_indices=2,
        )
        content = (tmp_path / "exchange.def").read_text()
        assert "NExchange          1" in content
        assert "1.500000000000000" in content
        lines = [l for l in content.splitlines() if "." in l and "===" not in l]
        assert len(lines) == 1

    def test_all_zero_still_writes_file(self, tmp_path, monkeypatch):
        """File is always written (flag check is in _process_interaction)."""
        monkeypatch.chdir(tmp_path)
        indx = [[0]]
        coeff = [0.0]
        _write_interaction_file(
            "coulombintra.def", "NCoulombIntra",
            "================== CoulombIntra ================",
            1, indx, coeff, n_indices=1,
        )
        content = (tmp_path / "coulombintra.def").read_text()
        assert "NCoulombIntra          0" in content

    def test_header_format(self, tmp_path, monkeypatch):
        """Verify the 5-line header structure."""
        monkeypatch.chdir(tmp_path)
        _write_interaction_file(
            "hund.def", "NHund",
            "=============== Hund coupling ===============",
            0, [], [], n_indices=2,
        )
        lines = (tmp_path / "hund.def").read_text().splitlines()
        assert lines[0] == "============================================="
        assert lines[1].startswith("NHund")
        assert lines[2] == "============================================="
        assert lines[3] == "=============== Hund coupling ==============="
        assert lines[4] == "============================================="


class TestWriteInterall:
    """Tests for _write_interall."""

    def test_zero_terms_no_file(self):
        StdI = _make_stdi_base()
        _make_empty_interactions(StdI)
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_interall(StdI)
                assert StdI.Lintr == 0
                assert not os.path.exists("interall.def")
            finally:
                os.chdir(orig)

    def test_boost_suppresses_interall(self):
        StdI = _make_stdi_base(lBoost=1)
        _make_empty_interactions(StdI)
        StdI.nintr = 1
        StdI.intrindx = [[0, 0, 1, 0, 2, 0, 3, 0]]
        StdI.intr = [1.0 + 0j]
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_interall(StdI)
                assert StdI.Lintr == 0
                assert not os.path.exists("interall.def")
            finally:
                os.chdir(orig)

    def test_nonzero_writes_file(self):
        StdI = _make_stdi_base()
        _make_empty_interactions(StdI)
        StdI.nintr = 2
        StdI.intrindx = [
            [0, 0, 1, 0, 2, 0, 3, 0],
            [4, 0, 5, 0, 6, 0, 7, 0],
        ]
        StdI.intr = [1.5 + 0.5j, 0.0 + 0j]  # one non-zero, one zero
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_interall(StdI)
                assert StdI.Lintr == 1
                assert os.path.exists("interall.def")
                content = open("interall.def").read()
                assert "NInterAll       1" in content
                assert "1.500000000000000" in content
                assert "0.500000000000000" in content
            finally:
                os.chdir(orig)

    def test_all_zero_no_file(self):
        StdI = _make_stdi_base()
        _make_empty_interactions(StdI)
        StdI.nintr = 2
        StdI.intrindx = [
            [0, 0, 1, 0, 2, 0, 3, 0],
            [4, 0, 5, 0, 6, 0, 7, 0],
        ]
        StdI.intr = [0.0 + 0j, 0.0 + 0j]
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_interall(StdI)
                assert StdI.Lintr == 0
                assert not os.path.exists("interall.def")
            finally:
                os.chdir(orig)
