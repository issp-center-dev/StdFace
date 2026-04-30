"""Unit tests for mvmc_writer module.

Tests for the mVMC solver-specific output functions extracted from
``stdface_main.py``.
"""
from __future__ import annotations

import os
import tempfile

import pytest

from stdface.core.stdface_vals import StdIntList
import numpy as np

from stdface.solvers.mvmc.writer import (
    print_orb,
    print_orb_para,
    print_gutzwiller,
    _has_anti_period,
    _compute_parallel_orbitals,
    _write_orbitalidxpara,
    _write_orbitalidxgen,
    _gutzwiller_momentum_projected,
    _gutzwiller_global_optimization,
    _write_gutzwiller_file,
)

# Sentinel values matching what _reset_vals sets at runtime
NaN_i = 2147483647

# The .def file header is 5 lines:
#   line 0: =============================================
#   line 1: NXxxIdx  <value>
#   line 2: ComplexType  <value>
#   line 3: =============================================
#   line 4: =============================================
HEADER_LINES = 5


def _make_stdi_for_orb(nsite: int = 4, **overrides) -> StdIntList:
    """Create a StdIntList with fields needed by print_orb / print_orb_para.

    Sets up a simple diagonal orbital matrix where Orb[i][j] = 0 for all
    pairs (single orbital) and AntiOrb[i][j] = 1 (no anti-periodic signs).
    """
    StdI = StdIntList()
    StdI.nsite = nsite
    StdI.NOrb = overrides.get("NOrb", 1)
    StdI.ComplexType = overrides.get("ComplexType", 0)

    # Default: all sites share orbital 0, sign = +1
    StdI.Orb = overrides.get(
        "Orb", [[0] * nsite for _ in range(nsite)]
    )
    StdI.AntiOrb = overrides.get(
        "AntiOrb", [[1] * nsite for _ in range(nsite)]
    )
    StdI.AntiPeriod = overrides.get("AntiPeriod", [0, 0, 0])
    return StdI


def _make_stdi_for_gutzwiller(nsite: int = 4, **overrides) -> StdIntList:
    """Create a StdIntList with fields needed by print_gutzwiller."""
    StdI = StdIntList()
    StdI.nsite = nsite
    StdI.NMPTrans = overrides.get("NMPTrans", NaN_i)
    StdI.model = overrides.get("model", "hubbard")

    # Default: diagonal orbital = 0 for all sites
    StdI.Orb = overrides.get(
        "Orb", [[0] * nsite for _ in range(nsite)]
    )
    StdI.locspinflag = overrides.get("locspinflag", [0] * nsite)
    StdI.NsiteUC = overrides.get("NsiteUC", 1)
    StdI.NCell = overrides.get("NCell", nsite)
    return StdI


class TestPrintOrb:
    """Tests for the print_orb function."""

    def test_writes_orbitalidx_def(self):
        """Test that orbitalidx.def is created with correct header."""
        StdI = _make_stdi_for_orb(nsite=2, NOrb=1, ComplexType=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_orb(StdI)
                assert os.path.exists("orbitalidx.def")
                content = open("orbitalidx.def").read()
                assert "NOrbitalIdx          1" in content
                assert "ComplexType          0" in content
            finally:
                os.chdir(orig)

    def test_no_anti_periodic_3_columns(self):
        """Test output format without anti-periodic boundaries (3 columns)."""
        StdI = _make_stdi_for_orb(nsite=2, NOrb=1, AntiPeriod=[0, 0, 0])

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_orb(StdI)
                lines = open("orbitalidx.def").readlines()
                # Data starts after 5-line header
                first_data = lines[HEADER_LINES].split()
                assert len(first_data) == 3  # isite, jsite, Orb (no anti-periodic)
            finally:
                os.chdir(orig)

    def test_anti_periodic_4_columns(self):
        """Test output format with anti-periodic boundaries (4 columns)."""
        StdI = _make_stdi_for_orb(nsite=2, NOrb=1, AntiPeriod=[1, 0, 0])

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_orb(StdI)
                lines = open("orbitalidx.def").readlines()
                first_data = lines[HEADER_LINES].split()
                assert len(first_data) == 4  # isite, jsite, Orb, AntiOrb
            finally:
                os.chdir(orig)

    def test_optimization_line_count(self):
        """Test that NOrb optimization lines are appended."""
        NOrb = 3
        StdI = _make_stdi_for_orb(nsite=2, NOrb=NOrb)
        # Set distinct orbital indices
        StdI.Orb = [[0, 1], [2, 0]]

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_orb(StdI)
                lines = open("orbitalidx.def").readlines()
                # Last NOrb lines should be optimization flags
                opt_lines = lines[-NOrb:]
                for i, line in enumerate(opt_lines):
                    parts = line.split()
                    assert parts[0] == str(i)
                    assert parts[1] == "1"
            finally:
                os.chdir(orig)


class TestPrintOrbPara:
    """Tests for the print_orb_para function."""

    def test_writes_both_files(self):
        """Test that both orbitalidxpara.def and orbitalidxgen.def are created."""
        StdI = _make_stdi_for_orb(nsite=2, NOrb=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_orb_para(StdI)
                assert os.path.exists("orbitalidxpara.def")
                assert os.path.exists("orbitalidxgen.def")
            finally:
                os.chdir(orig)

    def test_para_header(self):
        """Test orbitalidxpara.def header contains NOrbitalIdx and ComplexType."""
        StdI = _make_stdi_for_orb(nsite=2, NOrb=1, ComplexType=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_orb_para(StdI)
                content = open("orbitalidxpara.def").read()
                assert "NOrbitalIdx" in content
                assert "ComplexType          1" in content
            finally:
                os.chdir(orig)

    def test_gen_header_norb_count(self):
        """Test orbitalidxgen.def NOrbitalIdx = NOrb + 2*NOrbGC."""
        StdI = _make_stdi_for_orb(nsite=2, NOrb=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_orb_para(StdI)
                content = open("orbitalidxgen.def").read()
                assert "NOrbitalIdx" in content
                # With single orbital (all 0), NOrbGC = 1 (one pair in lower triangle)
                # Total = NOrb + 2*NOrbGC = 1 + 2*1 = 3
                assert "NOrbitalIdx          3" in content
            finally:
                os.chdir(orig)


class TestPrintGutzwiller:
    """Tests for the print_gutzwiller function."""

    def test_writes_gutzwilleridx_def(self):
        """Test that gutzwilleridx.def is created."""
        StdI = _make_stdi_for_gutzwiller(nsite=4, model="hubbard")

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_gutzwiller(StdI)
                assert os.path.exists("gutzwilleridx.def")
                content = open("gutzwilleridx.def").read()
                assert "NGutzwillerIdx" in content
                assert "ComplexType          0" in content
            finally:
                os.chdir(orig)

    def test_momentum_projected_hubbard(self):
        """Test momentum-projected mode for Hubbard model.

        With NMPTrans=NaN_i (unset) and all sites sharing orbital 0,
        NGutzwiller should be 1 (one unique orbital on diagonal).
        """
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="hubbard", NMPTrans=NaN_i,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_gutzwiller(StdI)
                content = open("gutzwilleridx.def").read()
                assert "NGutzwillerIdx          1" in content
            finally:
                os.chdir(orig)

    def test_momentum_projected_spin(self):
        """Test momentum-projected mode for spin model.

        For spin model, NGutzwiller starts at -1 and all local-spin sites
        are excluded (set to -1). With all locspinflag=1 (spin sites),
        no Gutzwiller indices are counted, so NGutzwiller = 1 (from -(-1)).
        """
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="spin", NMPTrans=1,
            locspinflag=[1, 1, 1, 1],  # all local spin
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_gutzwiller(StdI)
                content = open("gutzwilleridx.def").read()
                # NGutzwiller = -(-1) = 1
                assert "NGutzwillerIdx          1" in content
                # All Gutz[i] mapped to: -1 - (-1) = 0
                lines = open("gutzwilleridx.def").readlines()
                data_lines = lines[HEADER_LINES:HEADER_LINES + 4]
                for line in data_lines:
                    parts = line.split()
                    assert parts[1] == "0"  # all map to Gutzwiller index 0
            finally:
                os.chdir(orig)

    def test_global_optimisation_hubbard(self):
        """Test global-optimisation mode for Hubbard model.

        With NMPTrans=2, NGutzwiller = NsiteUC. Each site maps to
        isite % NsiteUC.
        """
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="hubbard", NMPTrans=2,
            NsiteUC=2, NCell=2,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_gutzwiller(StdI)
                content = open("gutzwilleridx.def").read()
                assert "NGutzwillerIdx          2" in content
                lines = open("gutzwilleridx.def").readlines()
                # Sites: 0->0, 1->1, 2->0, 3->1
                data_lines = lines[HEADER_LINES:HEADER_LINES + 4]
                expected_gutz = [0, 1, 0, 1]
                for i, line in enumerate(data_lines):
                    parts = line.split()
                    assert int(parts[1]) == expected_gutz[i]
            finally:
                os.chdir(orig)

    def test_global_optimisation_spin(self):
        """Test global-optimisation mode for spin model.

        With NMPTrans=2, NGutzwiller=1, all sites map to index 0.
        """
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="spin", NMPTrans=2,
            NsiteUC=2, NCell=2,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_gutzwiller(StdI)
                content = open("gutzwilleridx.def").read()
                assert "NGutzwillerIdx          1" in content
            finally:
                os.chdir(orig)

    def test_hubbard_opt_flags_all_one(self):
        """Test that Hubbard model optimization flags are all 1."""
        StdI = _make_stdi_for_gutzwiller(
            nsite=2, model="hubbard", NMPTrans=NaN_i,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_gutzwiller(StdI)
                lines = open("gutzwilleridx.def").readlines()
                # After header (5 lines) and site lines (nsite=2),
                # remaining lines are optimization flags
                opt_lines = lines[HEADER_LINES + 2:]
                assert len(opt_lines) > 0
                for line in opt_lines:
                    parts = line.split()
                    assert parts[1] == "1"
            finally:
                os.chdir(orig)

    def test_non_hubbard_opt_flag_zero_for_first(self):
        """Test that non-Hubbard model has optimization flag 0 for index 0."""
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="spin", NMPTrans=1,
            locspinflag=[1, 1, 1, 1],
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                print_gutzwiller(StdI)
                lines = open("gutzwilleridx.def").readlines()
                # After header (5) and site lines (4), optimization lines
                opt_lines = lines[HEADER_LINES + 4:]
                assert len(opt_lines) >= 1
                # First opt line: index 0 -> flag 0 (non-hubbard)
                parts = opt_lines[0].split()
                assert parts[0] == "0"
                assert parts[1] == "0"
            finally:
                os.chdir(orig)


# ===================================================================
#  _has_anti_period
# ===================================================================


class TestHasAntiPeriod:
    """Tests for _has_anti_period helper."""

    def test_no_anti_period(self):
        """Test that all-zero AntiPeriod returns False."""
        StdI = StdIntList()
        StdI.AntiPeriod[:] = 0
        assert not _has_anti_period(StdI)

    def test_first_direction_anti(self):
        """Test that AntiPeriod[0]==1 returns True."""
        StdI = StdIntList()
        StdI.AntiPeriod[:] = 0
        StdI.AntiPeriod[0] = 1
        assert _has_anti_period(StdI)

    def test_second_direction_anti(self):
        """Test that AntiPeriod[1]==1 returns True."""
        StdI = StdIntList()
        StdI.AntiPeriod[:] = 0
        StdI.AntiPeriod[1] = 1
        assert _has_anti_period(StdI)

    def test_third_direction_anti(self):
        """Test that AntiPeriod[2]==1 returns True."""
        StdI = StdIntList()
        StdI.AntiPeriod[:] = 0
        StdI.AntiPeriod[2] = 1
        assert _has_anti_period(StdI)

    def test_all_directions_anti(self):
        """Test that all AntiPeriod==1 returns True."""
        StdI = StdIntList()
        StdI.AntiPeriod[:] = 1
        assert _has_anti_period(StdI)

    def test_non_one_value_returns_false(self):
        """Test that AntiPeriod values other than 1 return False."""
        StdI = StdIntList()
        StdI.AntiPeriod[:] = 0
        StdI.AntiPeriod[0] = 2  # not 1
        assert not _has_anti_period(StdI)


# ===================================================================
#  _compute_parallel_orbitals
# ===================================================================


class TestComputeParallelOrbitals:
    """Tests for _compute_parallel_orbitals helper."""

    @staticmethod
    def _identity_orb(n: int):
        """Create identity-like orbital: Orb[i,j] = 0 for all pairs."""
        Orb = np.zeros((n, n), dtype=int)
        AntiOrb = np.ones((n, n), dtype=int)
        return Orb, AntiOrb

    @staticmethod
    def _distance_orb(n: int):
        """Create distance-based orbital: Orb[i,j] = |i-j|."""
        Orb = np.zeros((n, n), dtype=int)
        AntiOrb = np.ones((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                Orb[i, j] = abs(i - j)
        NOrb = n  # max distance is n-1, so indices 0..n-1
        return Orb, AntiOrb, NOrb

    def test_returns_three_values(self):
        """Test that function returns (OrbGC, reverse, NOrbGC)."""
        Orb, AntiOrb = self._identity_orb(3)
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(3, 1, Orb, AntiOrb)
        assert isinstance(OrbGC, list)
        assert isinstance(reverse, list)
        assert isinstance(NOrbGC, int)

    def test_orbgc_dimensions(self):
        """Test that OrbGC has correct dimensions."""
        n = 4
        Orb, AntiOrb = self._identity_orb(n)
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(n, 1, Orb, AntiOrb)
        assert len(OrbGC) == n
        assert all(len(row) == n for row in OrbGC)
        assert len(reverse) == n
        assert all(len(row) == n for row in reverse)

    def test_norbgc_non_negative(self):
        """Test that NOrbGC is non-negative."""
        Orb, AntiOrb, NOrb = self._distance_orb(4)
        _, _, NOrbGC = _compute_parallel_orbitals(4, NOrb, Orb, AntiOrb)
        assert NOrbGC >= 0

    def test_off_diagonal_non_negative(self):
        """Test that off-diagonal renumbered OrbGC values are non-negative."""
        Orb, AntiOrb, NOrb = self._distance_orb(4)
        OrbGC, _, NOrbGC = _compute_parallel_orbitals(4, NOrb, Orb, AntiOrb)
        for i in range(4):
            for j in range(4):
                if i != j:
                    assert OrbGC[i][j] >= 0

    def test_off_diagonal_bounded_by_norbgc(self):
        """Test that off-diagonal OrbGC values are < NOrbGC."""
        Orb, AntiOrb, NOrb = self._distance_orb(4)
        OrbGC, _, NOrbGC = _compute_parallel_orbitals(4, NOrb, Orb, AntiOrb)
        for i in range(4):
            for j in range(4):
                if i != j:
                    assert OrbGC[i][j] < NOrbGC

    def test_symmetrised(self):
        """Test that OrbGC is symmetric after symmetrisation step."""
        Orb, AntiOrb, NOrb = self._distance_orb(4)
        OrbGC, _, _ = _compute_parallel_orbitals(4, NOrb, Orb, AntiOrb)
        for i in range(4):
            for j in range(i + 1, 4):
                assert OrbGC[i][j] == OrbGC[j][i]

    def test_reverse_antisymmetric(self):
        """Test that reverse matrix is antisymmetric after symmetrisation."""
        Orb, AntiOrb, NOrb = self._distance_orb(4)
        _, reverse, _ = _compute_parallel_orbitals(4, NOrb, Orb, AntiOrb)
        for i in range(4):
            for j in range(i + 1, 4):
                assert reverse[i][j] == -reverse[j][i]

    def test_uniform_input_single_orbital(self):
        """Test that all-zero Orb produces a single parallel orbital."""
        n = 3
        Orb, AntiOrb = self._identity_orb(n)
        OrbGC, _, NOrbGC = _compute_parallel_orbitals(n, 1, Orb, AntiOrb)
        # All entries were orbital 0 → symmetrised → renumbered to single index
        assert NOrbGC == 1

    def test_does_not_modify_input(self):
        """Test that input Orb and AntiOrb arrays are not modified."""
        Orb, AntiOrb, NOrb = self._distance_orb(4)
        Orb_copy = Orb.copy()
        AntiOrb_copy = AntiOrb.copy()
        _compute_parallel_orbitals(4, NOrb, Orb, AntiOrb)
        assert np.array_equal(Orb, Orb_copy)
        assert np.array_equal(AntiOrb, AntiOrb_copy)


# ===================================================================
#  _gutzwiller_momentum_projected
# ===================================================================


class TestGutzwillerMomentumProjected:
    """Tests for the _gutzwiller_momentum_projected helper."""

    def test_hubbard_uniform_orbital(self):
        """Hubbard with all diagonal orbitals = 0 → NGutzwiller = 1."""
        StdI = _make_stdi_for_gutzwiller(nsite=4, model="hubbard")
        Gutz = [0] * 4
        NGutz = _gutzwiller_momentum_projected(StdI, Gutz)
        assert NGutz == 1

    def test_hubbard_distinct_orbitals(self):
        """Hubbard with 2 distinct diagonal orbitals → NGutzwiller = 2."""
        Orb = [[0]*4 for _ in range(4)]
        Orb[0][0] = 0; Orb[1][1] = 0; Orb[2][2] = 1; Orb[3][3] = 1
        StdI = _make_stdi_for_gutzwiller(nsite=4, model="hubbard", Orb=Orb)
        Gutz = [0] * 4
        NGutz = _gutzwiller_momentum_projected(StdI, Gutz)
        assert NGutz == 2
        # Sites 0,1 share one index, sites 2,3 share another
        assert Gutz[0] == Gutz[1]
        assert Gutz[2] == Gutz[3]
        assert Gutz[0] != Gutz[2]

    def test_spin_all_local(self):
        """Spin model with all local-spin sites → NGutzwiller = 1."""
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="spin", locspinflag=[1, 1, 1, 1],
        )
        Gutz = [0] * 4
        NGutz = _gutzwiller_momentum_projected(StdI, Gutz)
        assert NGutz == 1  # -(-1) = 1
        # All sites mapped to 0
        assert all(g == 0 for g in Gutz)

    def test_local_spin_sites_excluded(self):
        """Mixed model: local-spin sites excluded from Gutzwiller counting."""
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="kondo", locspinflag=[0, 0, 1, 1],
        )
        Gutz = [0] * 4
        NGutz = _gutzwiller_momentum_projected(StdI, Gutz)
        # Sites 2,3 are local spin → excluded (Gutz = 0 after inversion)
        # Sites 0,1 share orbital 0 → 1 unique index
        assert NGutz >= 1


# ===================================================================
#  _gutzwiller_global_optimization
# ===================================================================


class TestGutzwillerGlobalOptimization:
    """Tests for the _gutzwiller_global_optimization helper."""

    def test_hubbard_nsite_uc(self):
        """Hubbard: NGutzwiller = NsiteUC, sites map to isite mod NsiteUC."""
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="hubbard", NMPTrans=2,
            NsiteUC=2, NCell=2,
        )
        Gutz = [0] * 4
        NGutz = _gutzwiller_global_optimization(StdI, Gutz)
        assert NGutz == 2
        assert Gutz[0] == 0  # cell 0, site 0
        assert Gutz[1] == 1  # cell 0, site 1
        assert Gutz[2] == 0  # cell 1, site 0
        assert Gutz[3] == 1  # cell 1, site 1

    def test_spin_all_zero(self):
        """Spin: NGutzwiller = 1, all sites map to 0."""
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="spin", NMPTrans=2,
            NsiteUC=2, NCell=2,
        )
        Gutz = [0] * 4
        NGutz = _gutzwiller_global_optimization(StdI, Gutz)
        assert NGutz == 1
        assert all(g == 0 for g in Gutz)

    def test_kondo_conduction_and_localized(self):
        """Kondo: NGutzwiller = NsiteUC + 1, conduction→0, localized→isite+1."""
        StdI = _make_stdi_for_gutzwiller(
            nsite=4, model="kondo", NMPTrans=2,
            NsiteUC=1, NCell=2,
        )
        Gutz = [0] * 4
        NGutz = _gutzwiller_global_optimization(StdI, Gutz)
        assert NGutz == 2  # NsiteUC + 1 = 1 + 1
        # Conduction sites (first half)
        assert Gutz[0] == 0
        assert Gutz[1] == 0
        # Localized sites (second half)
        assert Gutz[2] == 1
        assert Gutz[3] == 1


# ===================================================================
#  _write_gutzwiller_file
# ===================================================================


class TestWriteGutzwillerFile:
    """Tests for the _write_gutzwiller_file output helper."""

    def test_writes_gutzwilleridx(self):
        """Test that gutzwilleridx.def is created."""
        StdI = _make_stdi_for_gutzwiller(nsite=2, model="hubbard")
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_gutzwiller_file(StdI, NGutzwiller=1, Gutz=[0, 0])
                assert os.path.exists("gutzwilleridx.def")
                content = open("gutzwilleridx.def").read()
                assert "NGutzwillerIdx          1" in content
            finally:
                os.chdir(orig)

    def test_hubbard_all_opt_flags_one(self):
        """Hubbard: all optimization flags are 1."""
        StdI = _make_stdi_for_gutzwiller(nsite=2, model="hubbard")
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_gutzwiller_file(StdI, NGutzwiller=2, Gutz=[0, 1])
                content = open("gutzwilleridx.def").read()
                lines = content.strip().split("\n")
                # Last 2 lines are opt flags: both should be "   x      1"
                opt_lines = lines[-2:]
                for line in opt_lines:
                    assert line.split()[1] == "1"
            finally:
                os.chdir(orig)

    def test_non_hubbard_first_flag_zero(self):
        """Non-Hubbard: first optimization flag is 0, rest are 1."""
        StdI = _make_stdi_for_gutzwiller(nsite=2, model="spin")
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_gutzwiller_file(StdI, NGutzwiller=2, Gutz=[0, 1])
                content = open("gutzwilleridx.def").read()
                lines = content.strip().split("\n")
                opt_lines = lines[-2:]
                assert opt_lines[0].split()[1] == "0"  # first → 0
                assert opt_lines[1].split()[1] == "1"  # rest → 1
            finally:
                os.chdir(orig)


# ===================================================================
#  _write_orbitalidxpara
# ===================================================================


class TestWriteOrbitalidxpara:
    """Tests for the _write_orbitalidxpara file-writing helper."""

    def test_creates_file(self):
        """Test that orbitalidxpara.def is created."""
        n = 3
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(
            n, 1, [[0]*n for _ in range(n)], [[1]*n for _ in range(n)])
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_orbitalidxpara(n, 0, OrbGC, reverse, NOrbGC)
                assert os.path.exists("orbitalidxpara.def")
            finally:
                os.chdir(orig)

    def test_header_contains_norbitalidx(self):
        """Test that header has NOrbitalIdx count."""
        n = 2
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(
            n, 1, [[0]*n for _ in range(n)], [[1]*n for _ in range(n)])
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_orbitalidxpara(n, 0, OrbGC, reverse, NOrbGC)
                content = open("orbitalidxpara.def").read()
                assert "NOrbitalIdx" in content
                assert "ComplexType" in content
            finally:
                os.chdir(orig)

    def test_upper_triangle_data_lines(self):
        """Test that only upper-triangle pairs (i<j) appear in data."""
        n = 3
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(
            n, 1, [[0]*n for _ in range(n)], [[1]*n for _ in range(n)])
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_orbitalidxpara(n, 0, OrbGC, reverse, NOrbGC)
                lines = open("orbitalidxpara.def").readlines()
                data_lines = lines[HEADER_LINES:]
                # Upper triangle of 3×3 has 3 pairs: (0,1), (0,2), (1,2)
                # Plus NOrbGC optimization lines
                pair_lines = [l for l in data_lines if len(l.split()) == 4]
                assert len(pair_lines) == 3
            finally:
                os.chdir(orig)

    def test_complex_type_passed_through(self):
        """Test that ComplexType value is correctly written."""
        n = 2
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(
            n, 1, [[0]*n for _ in range(n)], [[1]*n for _ in range(n)])
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_orbitalidxpara(n, 1, OrbGC, reverse, NOrbGC)
                content = open("orbitalidxpara.def").read()
                assert "ComplexType          1" in content
            finally:
                os.chdir(orig)


# ===================================================================
#  _write_orbitalidxgen
# ===================================================================


class TestWriteOrbitalidxgen:
    """Tests for the _write_orbitalidxgen file-writing helper."""

    def test_creates_file(self):
        """Test that orbitalidxgen.def is created."""
        StdI = _make_stdi_for_orb(nsite=2)
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(
            StdI.nsite, StdI.NOrb, StdI.Orb, StdI.AntiOrb)
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_orbitalidxgen(StdI, OrbGC, reverse, NOrbGC)
                assert os.path.exists("orbitalidxgen.def")
            finally:
                os.chdir(orig)

    def test_header_norb_count(self):
        """Test that header has NOrb + 2*NOrbGC as NOrbitalIdx."""
        StdI = _make_stdi_for_orb(nsite=2, NOrb=3)
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(
            StdI.nsite, StdI.NOrb, StdI.Orb, StdI.AntiOrb)
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_orbitalidxgen(StdI, OrbGC, reverse, NOrbGC)
                content = open("orbitalidxgen.def").read()
                expected = StdI.NOrb + 2 * NOrbGC
                assert f"NOrbitalIdx {expected:10d}" in content
            finally:
                os.chdir(orig)

    def test_no_anti_period_uses_sign_one(self):
        """Without anti-periodic boundaries, anti-parallel sign column is 1."""
        StdI = _make_stdi_for_orb(nsite=2, AntiPeriod=[0, 0, 0])
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(
            StdI.nsite, StdI.NOrb, StdI.Orb, StdI.AntiOrb)
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_orbitalidxgen(StdI, OrbGC, reverse, NOrbGC)
                lines = open("orbitalidxgen.def").readlines()
                # Anti-parallel section: 6-column lines with spin=0,1
                anti_lines = [l for l in lines[HEADER_LINES:]
                              if len(l.split()) == 6 and l.split()[1] == "0"
                              and l.split()[3] == "1"]
                for line in anti_lines:
                    parts = line.split()
                    assert parts[5] == "1"  # sign column = 1
            finally:
                os.chdir(orig)

    def test_anti_period_uses_anti_orb(self):
        """With anti-periodic boundary, anti-parallel sign uses AntiOrb."""
        n = 2
        AntiOrb = [[-1, 1], [1, -1]]
        StdI = _make_stdi_for_orb(nsite=n, AntiPeriod=[1, 0, 0], AntiOrb=AntiOrb)
        OrbGC, reverse, NOrbGC = _compute_parallel_orbitals(
            StdI.nsite, StdI.NOrb, StdI.Orb, StdI.AntiOrb)
        with tempfile.TemporaryDirectory() as tmpdir:
            orig = os.getcwd()
            os.chdir(tmpdir)
            try:
                _write_orbitalidxgen(StdI, OrbGC, reverse, NOrbGC)
                lines = open("orbitalidxgen.def").readlines()
                # Anti-parallel section: look for negative sign
                anti_lines = [l for l in lines[HEADER_LINES:]
                              if len(l.split()) == 6 and l.split()[1] == "0"
                              and l.split()[3] == "1"]
                signs = [int(l.split()[5]) for l in anti_lines]
                assert -1 in signs  # AntiOrb has -1 entries
            finally:
                os.chdir(orig)
