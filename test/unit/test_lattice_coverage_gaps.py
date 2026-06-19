"""Targeted tests for lattice code paths that were under-covered in ``cov.log``.

Focus: Hubbard/Kondo branches, Kondo site-offset loops, and related helpers.
Complements ``test_lattice_setup_integration.py``.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from stdface.core.stdface_vals import StdIntList
from stdface.lattice.chain_lattice import chain
from stdface.lattice.fc_ortho import fc_ortho
from stdface.lattice.honeycomb_lattice import honeycomb
from stdface.lattice.kagome import kagome
from stdface.lattice.ladder import ladder
from stdface.lattice.orthorhombic import orthorhombic
from stdface.lattice.pyrochlore import pyrochlore
from stdface.lattice.square_lattice import tetragonal
from stdface.lattice.triangular_lattice import triangular

NaN_d = float("nan")
NaN_i = 2147483647
NaN_c = complex(float("nan"), 0.0)


def _fill_j_mats_nan(s: StdIntList) -> None:
    for name in (
        "J", "Jp", "Jpp",
        "J0", "J0p", "J0pp",
        "J1", "J1p", "J1pp",
        "J2", "J2p", "J2pp",
    ):
        getattr(s, name)[:, :] = NaN_d


def _fill_j_all_nan(s: StdIntList) -> None:
    for attr in (
        "JAll", "JpAll", "JppAll",
        "J0All", "J0pAll", "J0ppAll",
        "J1All", "J1pAll", "J1ppAll",
        "J2All", "J2pAll", "J2ppAll",
    ):
        setattr(s, attr, NaN_d)


def _base_chain_stdi(L: int) -> StdIntList:
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.solver = "HPhi"
    s.lattice = "chain"
    s.a = NaN_d
    s.length[:] = NaN_d
    s.direct[:, :] = NaN_d
    s.phase[0] = NaN_d
    s.phase[1] = NaN_d
    s.L = L
    s.W = None
    s.Height = None
    s.box[:, :] = NaN_i
    s.h = NaN_d
    s.Gamma = NaN_d
    s.Gamma_y = NaN_d
    s.D[:, :] = NaN_d
    _fill_j_mats_nan(s)
    _fill_j_all_nan(s)
    # chain ``not_used`` before spin/fermion branch (lines 101–113)
    for attr in ("t1", "t2", "t1p", "t2p"):
        setattr(s, attr, NaN_c)
    for attr in ("V1", "V2", "V1p", "V2p"):
        setattr(s, attr, NaN_d)
    s.K = NaN_d
    return s


def make_hubbard_chain(L: int = 4) -> StdIntList:
    """Hubbard chain: fermion branch with only nn hopping ``t0``."""
    s = _base_chain_stdi(L)
    s.model = "hubbard"
    s.S2 = None
    s.mu = 0.0
    s.U = 0.0
    s.t = NaN_c
    s.t0 = complex(1.0, 0.0)
    s.tp = NaN_c
    s.t0p = NaN_c
    s.t0pp = NaN_c
    s.V = NaN_d
    s.V0 = NaN_d
    s.Vp = NaN_d
    s.V0p = NaN_d
    s.V0pp = NaN_d
    return s


def make_kondo_chain(L: int = 4) -> StdIntList:
    """Kondo chain: ``input_spin`` on ``J`` and enlarged interaction counts."""
    s = _base_chain_stdi(L)
    s.model = "kondo"
    s.S2 = None
    s.mu = 0.0
    s.U = 0.0
    s.t = NaN_c
    s.t0 = complex(0.5, 0.0)
    s.tp = NaN_c
    s.t0p = NaN_c
    s.t0pp = NaN_c
    s.V = NaN_d
    s.V0 = NaN_d
    s.Vp = NaN_d
    s.V0p = NaN_d
    s.V0pp = NaN_d
    s.JAll = 1.0
    return s


def make_kondo_square(L: int = 2, W: int = 2) -> StdIntList:
    """Square lattice Kondo: exercises ``input_spin`` and Kondo ``isite`` offset."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "kondo"
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
    s.h = NaN_d
    s.Gamma = NaN_d
    s.Gamma_y = NaN_d
    s.D[:, :] = NaN_d
    _fill_j_mats_nan(s)
    _fill_j_all_nan(s)
    s.mu = 0.0
    s.U = 0.0
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
    s.JAll = 1.0
    return s


def make_hubbard_triangular(L: int = 2, W: int = 2) -> StdIntList:
    """Triangular Hubbard: full fermion branch (hoppings + Coulombs)."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "hubbard"
    s.solver = "HPhi"
    s.lattice = "triangular"
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
    s.mu = 0.0
    s.U = 0.0
    s.t = NaN_c
    s.t0 = complex(1.0, 0.0)
    for attr in (
        "t1", "t2", "tp", "t0p", "t1p", "t2p",
        "t0pp", "t1pp", "t2pp", "tpp",
    ):
        setattr(s, attr, NaN_c)
    for attr in (
        "V", "V0", "V1", "V2", "Vp", "V0p", "V1p", "V2p",
        "Vpp", "V0pp", "V1pp", "V2pp",
    ):
        setattr(s, attr, NaN_d)
    s.K = NaN_d
    return s


def make_hubbard_orthorhombic(L: int = 2, W: int = 2, H: int = 2) -> StdIntList:
    """Simple orthorhombic Hubbard (3D fermion branch including ``tpp`` / ``Vpp``)."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "hubbard"
    s.solver = "HPhi"
    s.lattice = "orthorhombic"
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
    s.mu = 0.0
    s.U = 0.0
    s.t = NaN_c
    s.t0 = complex(1.0, 0.0)
    for attr in ("t1", "t2", "tp", "t0p", "t1p", "t2p", "tpp"):
        setattr(s, attr, NaN_c)
    for attr in ("V", "V0", "V1", "V2", "Vp", "V0p", "V1p", "V2p", "Vpp"):
        setattr(s, attr, NaN_d)
    s.K = NaN_d
    return s


def make_hubbard_fc_ortho(L: int = 2, W: int = 2, H: int = 2) -> StdIntList:
    """Face-centered orthorhombic Hubbard (``J''`` / ``Jpp`` only, no separate tpp line)."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "hubbard"
    s.solver = "HPhi"
    s.lattice = "fco"
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
    s.mu = 0.0
    s.U = 0.0
    s.t = NaN_c
    s.t0 = complex(1.0, 0.0)
    for attr in ("t1", "t2", "tp", "t0p", "t1p", "t2p", "tpp"):
        setattr(s, attr, NaN_c)
    for attr in ("V", "V0", "V1", "V2", "Vp", "V0p", "V1p", "V2p"):
        setattr(s, attr, NaN_d)
    s.K = NaN_d
    return s


def make_hubbard_pyrochlore(L: int = 2, W: int = 2, H: int = 2) -> StdIntList:
    """Pyrochlore Hubbard (no ``t0''`` / ``Vpp`` in the fermion branch)."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "hubbard"
    s.solver = "HPhi"
    s.lattice = "pyrochlore"
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
    s.mu = 0.0
    s.U = 0.0
    s.t = NaN_c
    s.t0 = complex(1.0, 0.0)
    for attr in ("t1", "t2", "tp", "t0p", "t1p", "t2p", "tpp"):
        setattr(s, attr, NaN_c)
    for attr in ("V", "V0", "V1", "V2", "Vp", "V0p", "V1p", "V2p"):
        setattr(s, attr, NaN_d)
    s.K = NaN_d
    return s


def make_hubbard_honeycomb(L: int = 2, W: int = 2) -> StdIntList:
    """Honeycomb Hubbard (same fermion pattern as triangular)."""
    s = make_hubbard_triangular(L, W)
    s.lattice = "honeycomb"
    return s


def make_hubbard_kagome(L: int = 2, W: int = 2) -> StdIntList:
    """Kagome Hubbard: ``t0``–``t2p`` and ``V0``–``V2p`` only (plus ``not_used_j`` on ``J'``)."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "hubbard"
    s.solver = "HPhi"
    s.lattice = "kagome"
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
    s.mu = 0.0
    s.U = 0.0
    s.t = NaN_c
    s.t0 = complex(1.0, 0.0)
    for attr in ("t1", "t2", "tp", "t0p", "t1p", "t2p", "tpp"):
        setattr(s, attr, NaN_c)
    for attr in ("V", "V0", "V1", "V2", "Vp", "V0p", "V1p", "V2p"):
        setattr(s, attr, NaN_d)
    s.K = NaN_d
    return s


def make_hubbard_ladder(L: int = 4, W_legs: int = 2) -> StdIntList:
    """Ladder Hubbard: ``input_hopp`` uses ``t`` (not ``tp`` / ``t0p``)."""
    s = StdIntList()
    s.pi = math.acos(-1.0)
    s.pi180 = s.pi / 180.0
    s.model = "hubbard"
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
    s.t = NaN_c
    s.tp = NaN_c
    s.V = NaN_d
    s.Vp = NaN_d
    s.K = NaN_d
    s.mu = 0.0
    s.U = 0.0
    s.t0 = complex(1.0, 0.0)
    for attr in ("t1", "t2", "t1p", "t2p"):
        setattr(s, attr, NaN_c)
    for attr in ("V0", "V1", "V2", "V1p", "V2p"):
        setattr(s, attr, NaN_d)
    return s


def make_kondo_triangular(L: int = 2, W: int = 2) -> StdIntList:
    """Triangular Kondo: ``input_spin`` on ``J`` and Kondo ``isite`` offset."""
    s = make_hubbard_triangular(L, W)
    s.model = "kondo"
    s.JAll = 1.0
    return s


def make_kondo_honeycomb(L: int = 2, W: int = 2) -> StdIntList:
    """Honeycomb Kondo: same hoppings/Coulombs as Hubbard, ``J`` from ``input_spin``."""
    s = make_hubbard_honeycomb(L, W)
    s.model = "kondo"
    s.JAll = 1.0
    return s


class TestChainHubbardAndKondo:
    """``chain_lattice.chain``: Hubbard and Kondo fermion branches."""

    def test_hubbard_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_hubbard_chain(4)
        chain(s)
        assert s.model == "hubbard"
        assert s.nsite == 4
        assert len(s.trans_list) > 0

    def test_kondo_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_kondo_chain(4)
        chain(s)
        assert s.model == "kondo"
        assert s.nsite == 8
        assert len(s.trans_list) > 0


class TestSquareKondo:
    """``square_lattice.tetragonal``: Kondo branch (``input_spin``, ``isite += NCell``)."""

    def test_kondo_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_kondo_square(2, 2)
        tetragonal(s)
        assert s.model == "kondo"
        assert s.nsite == 8
        assert len(s.trans_list) > 0


class TestTriangularHubbard:
    """``triangular_lattice.triangular``: Hubbard fermion branch."""

    def test_hubbard_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_hubbard_triangular(2, 2)
        triangular(s)
        assert s.model == "hubbard"
        assert s.NCell == 4
        assert s.nsite == 4
        assert len(s.trans_list) > 0


class TestOrthorhombicPyrochloreFCOHoneycombHubbard:
    """3D orthorhombic / FCO / pyrochlore and 2D honeycomb Hubbard paths."""

    def test_orthorhombic_hubbard_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_hubbard_orthorhombic(2, 2, 2)
        orthorhombic(s)
        assert s.model == "hubbard"
        assert s.nsite == 8
        assert len(s.trans_list) > 0

    def test_fc_ortho_hubbard_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_hubbard_fc_ortho(2, 2, 2)
        fc_ortho(s)
        assert s.model == "hubbard"
        assert s.nsite == 8
        assert len(s.trans_list) > 0

    def test_pyrochlore_hubbard_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_hubbard_pyrochlore(2, 2, 2)
        pyrochlore(s)
        assert s.model == "hubbard"
        # ``NsiteUC == 4`` tetrahedral sites per cell × ``L*W*H`` cells
        assert s.nsite == 32
        assert len(s.trans_list) > 0

    def test_honeycomb_hubbard_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_hubbard_honeycomb(2, 2)
        honeycomb(s)
        assert s.model == "hubbard"
        assert s.NCell == 4
        # Honeycomb has ``NsiteUC == 2`` (two sublattice sites) × ``NCell`` cells
        assert s.nsite == 8
        assert len(s.trans_list) > 0

    def test_honeycomb_kondo_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_kondo_honeycomb(2, 2)
        honeycomb(s)
        assert s.model == "kondo"
        assert s.NCell == 4
        # Kondo doubles physical sites: ``2 * NsiteUC * NCell``
        assert s.nsite == 16
        assert len(s.trans_list) > 0


class TestKagomeLadderTriangularKondoHubbard:
    """Kagome / ladder Hubbard and triangular Kondo."""

    def test_kagome_hubbard_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_hubbard_kagome(2, 2)
        kagome(s)
        assert s.model == "hubbard"
        assert s.nsite == 12
        assert len(s.trans_list) > 0

    def test_ladder_hubbard_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_hubbard_ladder(4, 2)
        ladder(s)
        assert s.model == "hubbard"
        assert s.nsite == 8
        assert len(s.trans_list) > 0

    def test_triangular_kondo_runs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        s = make_kondo_triangular(2, 2)
        triangular(s)
        assert s.model == "kondo"
        assert s.nsite == 8
        assert len(s.trans_list) > 0
