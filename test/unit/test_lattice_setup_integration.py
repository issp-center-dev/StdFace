"""Integration-style unit tests for lattice modules under ``python/stdface/lattice/``.

Runs each lattice ``setup`` and checks geometry, site counts, and interaction arrays.
Import/docstring smoke tests live in ``test_lattices.py``; the plugin registry is
covered in ``test_lattice_dispatch.py``.
"""
from __future__ import annotations

import math
from typing import Callable

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.lattice import (
    LatticePlugin,
    get_lattice,
    register_lattice,
)
from stdface.lattice.fc_ortho import fc_ortho
from stdface.lattice.honeycomb_lattice import honeycomb, honeycomb_boost
from stdface.lattice.kagome import kagome, kagome_boost
from stdface.lattice.ladder import ladder, ladder_boost
from stdface.lattice.orthorhombic import orthorhombic
from stdface.lattice.pyrochlore import pyrochlore
from stdface.lattice.square_lattice import tetragonal
from stdface.lattice.triangular_lattice import triangular

# ---------------------------------------------------------------------------
#  Constants and shared helpers
# ---------------------------------------------------------------------------

NaN_d = float("nan")
NaN_i = 2147483647
NaN_c = complex(float("nan"), 0.0)


def _fill_unused_fermion_fields(s: StdIntList) -> None:
    """Mark fermion / density parameters unused for spin runs."""
    s.mu = NaN_d
    s.U = NaN_d
    s.t = NaN_c
    s.t0 = NaN_c
    s.tp = NaN_c
    s.t1 = NaN_c
    s.t2 = NaN_c
    s.t0p = NaN_c
    s.t1p = NaN_c
    s.t2p = NaN_c
    s.tpp = NaN_c
    s.t0pp = NaN_c
    s.t1pp = NaN_c
    s.t2pp = NaN_c
    s.V = NaN_d
    s.V0 = NaN_d
    s.V1 = NaN_d
    s.V2 = NaN_d
    s.Vp = NaN_d
    s.V0p = NaN_d
    s.V1p = NaN_d
    s.V2p = NaN_d
    s.Vpp = NaN_d
    s.V0pp = NaN_d
    s.V1pp = NaN_d
    s.V2pp = NaN_d
    s.K = NaN_d


def _fill_j_mats_nan(s: StdIntList) -> None:
    """Set all 3×3 exchange matrices to NaN."""
    for name in (
        "J", "Jp", "Jpp",
        "J0", "J0p", "J0pp",
        "J1", "J1p", "J1pp",
        "J2", "J2p", "J2pp",
    ):
        getattr(s, name)[:, :] = NaN_d


def _fill_j_all_nan(s: StdIntList) -> None:
    """Set scalar J*All fields to NaN (for unused-parameter conflict checks)."""
    for attr in (
        "JAll", "JpAll", "JppAll",
        "J0All", "J0pAll", "J0ppAll",
        "J1All", "J1pAll", "J1ppAll",
        "J2All", "J2pAll", "J2ppAll",
    ):
        setattr(s, attr, NaN_d)


def make_spin_stdint_2d(
    *,
    lattice: str,
    L: int,
    W: int,
    use_shared_jall: bool = True,
) -> StdIntList:
    """Minimal spin ``StdIntList`` for 2D lattices (``init_site`` uses the L/W path)."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "spin"
    s.solver = "HPhi"
    s.lattice = lattice

    s.a = NaN_d
    s.length[:] = NaN_d
    s.direct[:, :] = NaN_d
    s.phase[:] = NaN_d

    s.L = L
    s.W = W
    s.Height = None
    s.box[:, :] = NaN_i

    s.S2 = None
    s.h = NaN_d
    s.Gamma = NaN_d
    s.Gamma_y = NaN_d
    s.D[:, :] = NaN_d

    _fill_j_mats_nan(s)
    _fill_j_all_nan(s)
    _fill_unused_fermion_fields(s)

    if use_shared_jall:
        s.JAll = 1.0
    return s


def make_spin_stdint_3d(*, lattice: str, L: int, W: int, H: int) -> StdIntList:
    """Spin ``StdIntList`` for 3D lattices (orthorhombic, fcc, pyrochlore)."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "spin"
    s.solver = "HPhi"
    s.lattice = lattice

    s.a = NaN_d
    s.length[:] = NaN_d
    s.direct[:, :] = NaN_d
    s.phase[:] = NaN_d

    s.L = L
    s.W = W
    s.Height = H
    s.box[:, :] = NaN_i

    s.S2 = None
    s.h = NaN_d
    s.Gamma = NaN_d
    s.Gamma_y = NaN_d
    s.D[:, :] = NaN_d

    _fill_j_mats_nan(s)
    _fill_j_all_nan(s)
    s.JAll = 1.0
    _fill_unused_fermion_fields(s)
    return s


def make_spin_ladder(L: int, W_legs: int) -> StdIntList:
    """Ladder only: ``input_spin`` paths use per-matrix JkAll values."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "spin"
    s.solver = "HPhi"
    s.lattice = "ladder"

    s.a = NaN_d
    s.length[:] = NaN_d
    s.direct[:, :] = NaN_d
    s.phase[:] = NaN_d

    s.L = L
    s.W = W_legs
    s.Height = None
    s.box[:, :] = NaN_i

    s.S2 = None
    s.h = NaN_d
    s.Gamma = NaN_d
    s.Gamma_y = NaN_d
    s.D[:, :] = NaN_d

    _fill_j_mats_nan(s)
    _fill_j_all_nan(s)
    _fill_unused_fermion_fields(s)

    for attr in ("J0All", "J1All", "J2All", "J1pAll", "J2pAll"):
        setattr(s, attr, 1.0)
    return s


def _assert_post_spin_setup(s: StdIntList, *, nsite_uc: int) -> None:
    """Shared assertions after a spin lattice setup."""
    assert s.NsiteUC == nsite_uc
    assert s.nsite == s.NCell * nsite_uc
    assert s.tau.shape == (nsite_uc, 3)
    assert isinstance(s.trans_list, list) and isinstance(s.intr_list, list)
    assert s.locspinflag is not None and len(s.locspinflag) == s.nsite


# ---------------------------------------------------------------------------
#  2D lattices
# ---------------------------------------------------------------------------


class TestTetragonalSpin:
    """``tetragonal``: square-lattice spin model."""

    def test_small_cell_geometry(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_stdint_2d(lattice="square", L=2, W=2)
        tetragonal(s)
        assert s.NCell == 4
        _assert_post_spin_setup(s, nsite_uc=1)
        np.testing.assert_array_equal(s.box, np.diag([2, 2, 1]))


class TestTriangularSpin:
    """``triangular``: triangular lattice."""

    def test_runs_and_cell_count(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_stdint_2d(lattice="triangular", L=2, W=2)
        triangular(s)
        assert s.NCell == 4
        _assert_post_spin_setup(s, nsite_uc=1)
        # Default a=1 → Ly on the second basis vector along L is Llength * sqrt(3)/2
        assert s.direct[1, 1] == pytest.approx(0.5 * math.sqrt(3.0))


class TestLadderSpin:
    """``ladder``: ladder (NsiteUC equals number of legs)."""

    def test_two_leg_ladder(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_ladder(L=4, W_legs=2)
        ladder(s)
        assert s.NsiteUC == 2
        assert s.W == 1
        assert s.NCell == 4
        _assert_post_spin_setup(s, nsite_uc=2)
        assert s.tau.shape == (2, 3)


class TestHoneycombSpin:
    """``honeycomb``: honeycomb (two sites per unit cell)."""

    def test_two_site_unit_cell(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_stdint_2d(lattice="honeycomb", L=2, W=2)
        honeycomb(s)
        assert s.NsiteUC == 2
        assert s.NCell == 4
        assert s.nsite == 8
        _assert_post_spin_setup(s, nsite_uc=2)


class TestKagomeSpin:
    """``kagome``: kagome (three sites per unit cell)."""

    def test_three_site_unit_cell(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_stdint_2d(lattice="kagome", L=2, W=2)
        kagome(s)
        assert s.NsiteUC == 3
        assert s.NCell == 4
        assert s.nsite == 12
        _assert_post_spin_setup(s, nsite_uc=3)


# ---------------------------------------------------------------------------
#  3D lattices
# ---------------------------------------------------------------------------


class TestOrthorhombicSpin:
    """``orthorhombic``: simple orthorhombic lattice."""

    def test_2x2x2(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_stdint_3d(lattice="orthorhombic", L=2, W=2, H=2)
        orthorhombic(s)
        assert s.NCell == 8
        assert s.NsiteUC == 1
        _assert_post_spin_setup(s, nsite_uc=1)
        assert (tmp_path / "lattice.xsf").exists()


class TestFcOrthoSpin:
    """``fc_ortho``: face-centered orthorhombic lattice."""

    def test_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_stdint_3d(lattice="fco", L=2, W=2, H=2)
        fc_ortho(s)
        assert s.NsiteUC == 1
        assert s.nsite == s.NCell
        _assert_post_spin_setup(s, nsite_uc=1)


class TestPyrochloreSpin:
    """``pyrochlore``: pyrochlore (four sites per unit cell)."""

    def test_four_site_unit_cell(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_stdint_3d(lattice="pyrochlore", L=2, W=2, H=2)
        pyrochlore(s)
        assert s.NsiteUC == 4
        assert s.tau.shape == (4, 3)
        assert s.nsite == s.NCell * 4
        _assert_post_spin_setup(s, nsite_uc=4)


# ---------------------------------------------------------------------------
#  Hubbard (example: square lattice only)
# ---------------------------------------------------------------------------


class TestTetragonalHubbard:
    """Square-lattice Hubbard with spin exchange parameters left unused."""

    def _make_hubbard_square(self, L: int, W: int) -> StdIntList:
        s = StdIntList()
        s.pi = math.acos(-1.0)
        s.pi180 = s.pi / 180.0
        s.model = "hubbard"
        s.solver = "HPhi"
        s.lattice = "square"

        s.a = NaN_d
        s.length[:] = NaN_d
        s.direct[:, :] = NaN_d
        s.phase[:] = NaN_d

        s.L = L
        s.W = W
        s.Height = None
        s.box[:, :] = NaN_i

        s.S2 = None
        s.mu = 0.0
        s.U = 0.0
        s.h = NaN_d
        s.Gamma = NaN_d
        s.Gamma_y = NaN_d
        s.D[:, :] = NaN_d

        _fill_j_mats_nan(s)
        _fill_j_all_nan(s)

        s.t = NaN_c
        s.t0 = complex(1.0, 0.0)
        s.t1 = NaN_c
        s.t2 = NaN_c
        s.tp = NaN_c
        s.t0p = NaN_c
        s.t1p = NaN_c
        s.t2p = NaN_c
        s.t0pp = NaN_c
        s.t1pp = NaN_c
        s.tpp = NaN_c
        s.V = NaN_d
        s.V0 = NaN_d
        s.V1 = NaN_d
        s.V2 = NaN_d
        s.Vp = NaN_d
        s.V0p = NaN_d
        s.V1p = NaN_d
        s.V2p = NaN_d
        s.V0pp = NaN_d
        s.V1pp = NaN_d
        s.K = NaN_d
        return s

    def test_nn_hopping_only(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = self._make_hubbard_square(2, 2)
        tetragonal(s)
        assert s.model == "hubbard"
        assert s.NCell == 4
        assert s.nsite == 4
        assert len(s.trans_list) > 0


# ---------------------------------------------------------------------------
#  Boost (HPhi representative cases)
# ---------------------------------------------------------------------------


class TestLadderBoost:
    """``ladder_boost``: writes ``boost.def`` when constraints are satisfied."""

    def test_creates_boost_def(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_ladder(L=8, W_legs=2)
        ladder(s)
        s.L = 8
        s.S2 = 1
        ladder_boost(s)
        assert (tmp_path / "boost.def").exists()


class TestHoneycombBoost:
    """``honeycomb_boost``: runnable with W=3, L>=2, S2=1."""

    def test_creates_boost_def(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_stdint_2d(lattice="honeycomb", L=2, W=3)
        honeycomb(s)
        s.S2 = 1
        honeycomb_boost(s)
        assert (tmp_path / "boost.def").exists()


class TestKagomeBoost:
    """``kagome_boost``: when W=3, W divisible by 3, and other constraints hold."""

    def test_creates_boost_def(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_spin_stdint_2d(lattice="kagome", L=2, W=3)
        kagome(s)
        s.S2 = 1
        kagome_boost(s)
        assert (tmp_path / "boost.def").exists()


# ---------------------------------------------------------------------------
#  stdface.lattice package API
# ---------------------------------------------------------------------------


class TestLatticePackageAPI:
    """Behavior of ``LatticePlugin`` and ``register_lattice``."""

    def test_cannot_instantiate_abstract_plugin(self):
        class Incomplete(LatticePlugin):
            @property
            def name(self) -> str:
                return "x"

            @property
            def aliases(self) -> list[str]:
                return ["x"]

            @property
            def ndim(self) -> int:
                return 1

        with pytest.raises(TypeError):
            Incomplete()  # type: ignore[misc]

    def test_register_duplicate_alias_raises(self):
        existing = get_lattice("chain")

        class Clash(LatticePlugin):
            @property
            def name(self) -> str:
                return "fakechain"

            @property
            def aliases(self) -> list[str]:
                return ["chain"]

            @property
            def ndim(self) -> int:
                return 1

            def setup(self, StdI: StdIntList) -> None:
                pass

        with pytest.raises(ValueError, match="already registered"):
            register_lattice(Clash())

        assert get_lattice("chain") is existing

    def test_plugin_setup_is_callable(self):
        for name in (
            "tetragonal",
            "triangular",
            "honeycomb",
            "kagome",
            "orthorhombic",
            "fco",
            "pyrochlore",
        ):
            p = get_lattice(name)
            assert isinstance(p, LatticePlugin)
            assert callable(p.setup)


class TestLatticeRoundTrip:
    """Smoke test: ``setup`` completes without error via the registry."""

    @staticmethod
    def _run(plugin_name: str, factory: Callable[[], StdIntList], tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        p = get_lattice(plugin_name)
        s = factory()
        p.setup(s)

    def test_tetragonal_via_registry(self, tmp_path, monkeypatch):
        self._run(
            "square",
            lambda: make_spin_stdint_2d(lattice="square", L=2, W=2),
            tmp_path,
            monkeypatch,
        )

    def test_ladder_via_registry(self, tmp_path, monkeypatch):
        self._run(
            "ladder",
            lambda: make_spin_ladder(L=4, W_legs=2),
            tmp_path,
            monkeypatch,
        )

    def test_orthorhombic_via_registry(self, tmp_path, monkeypatch):
        self._run(
            "cubic",
            lambda: make_spin_stdint_3d(lattice="cubic", L=2, W=2, H=2),
            tmp_path,
            monkeypatch,
        )
