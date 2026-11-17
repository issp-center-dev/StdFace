"""HPhi momentum-sector support (``MomentumIndex`` keyword).

When ``MomentumIndex`` is given, HPhi diagonalises within one translation
(momentum) sector of the chain.  This module validates that the generated
Hamiltonian is one HPhi supports there and builds ``qptransidx.def``, the
TransSym input listing the translation operators and their weights.
"""
from __future__ import annotations

import logging
import math

from ...core.stdface_vals import StdIntList, ModelType, MethodType
from ...core.output import SolverFileData

logger = logging.getLogger(__name__)

_EPS = 1.0e-12
_CHAIN_NAMES = ("chain", "chainlattice")
_SUPPORTED_METHODS = (None, MethodType.LANCZOS, MethodType.LANCZOS_ENERGY,
                      MethodType.CG)


def uses_momentum_symmetry(StdI: StdIntList) -> bool:
    return StdI.MomentumIndex is not None


def _fail(msg: str) -> None:
    msg = f"ERROR ! {msg}"
    logger.error(msg)
    raise ValueError(msg)


def _any_nonzero(terms: list) -> bool:
    return any(abs(term[0]) > _EPS for term in terms)


def _nonzero(value) -> bool:
    return value is not None and abs(value) > _EPS


def _check_method(StdI: StdIntList) -> None:
    if StdI.method not in _SUPPORTED_METHODS:
        _fail("MomentumIndex currently supports only Lanczos and CG.")


def _check_spin(StdI: StdIntList) -> None:
    if StdI.lGC != 0:
        _fail("MomentumIndex currently supports only canonical Spin model.")
    if StdI.S2 != 1:
        _fail("MomentumIndex currently supports only Spin-1/2.")
    _check_method(StdI)
    if (_any_nonzero(StdI.trans_list) or _any_nonzero(StdI.intr_list)
            or _any_nonzero(StdI.Cinter_list) or _any_nonzero(StdI.Hund_list)
            or _any_nonzero(StdI.PairLift_list)
            or _any_nonzero(StdI.PairHopp_list)):
        _fail("MomentumIndex currently supports only exchange-only Spin-1/2 "
              "chain with Jz = 0 and no field/general/pair terms.")


def _check_hubbard(StdI: StdIntList) -> None:
    if StdI.lGC != 0:
        _fail("MomentumIndex currently supports only canonical Hubbard model.")
    if StdI.ncond is None or StdI.Sz2 is None:
        _fail("MomentumIndex for Hubbard requires nelec and 2Sz.")
    _check_method(StdI)
    if (_nonzero(StdI.mu) or _nonzero(StdI.h) or _nonzero(StdI.Gamma)
            or _nonzero(StdI.Gamma_y)):
        _fail("MomentumIndex for Hubbard does not support "
              "mu/h/Gamma/Gamma_y terms.")
    if _nonzero(StdI.t0p) or _nonzero(StdI.t0pp):
        _fail("MomentumIndex for Hubbard does not support "
              "next-neighbor hopping terms.")
    if (_nonzero(StdI.V0) or _nonzero(StdI.V0p) or _nonzero(StdI.V0pp)
            or _any_nonzero(StdI.Cinter_list)):
        _fail("MomentumIndex for Hubbard does not support "
              "V/CoulombInter terms.")
    if (_any_nonzero(StdI.intr_list) or _any_nonzero(StdI.Hund_list)
            or _any_nonzero(StdI.PairLift_list)
            or _any_nonzero(StdI.PairHopp_list)):
        _fail("MomentumIndex for Hubbard supports only nearest-neighbor "
              "hopping and onsite U.")


_MODEL_CHECKS = {
    ModelType.SPIN: _check_spin,
    ModelType.HUBBARD: _check_hubbard,
}


def check_momentum_symmetry(StdI: StdIntList) -> None:
    """Reject Hamiltonians HPhi cannot handle in a momentum sector.

    Must run after lattice setup (the term lists are inspected) and is a
    no-op unless ``MomentumIndex`` was given.
    """
    if not uses_momentum_symmetry(StdI):
        return
    if StdI.lattice not in _CHAIN_NAMES:
        _fail("MomentumIndex currently supports only chain lattice.")
    if StdI.L <= 0:
        _fail("MomentumIndex requires positive L.")
    if not 0 <= StdI.MomentumIndex < StdI.L:
        _fail("MomentumIndex must satisfy 0 <= MomentumIndex < L.")
    for idim in range(3):
        if abs(StdI.phase[idim]) > _EPS or StdI.AntiPeriod[idim] != 0:
            _fail("MomentumIndex does not support boundary phase in v1.")
    check = _MODEL_CHECKS.get(StdI.model)
    if check is None:
        _fail("MomentumIndex currently supports only Spin, SpinlessFermion, "
              "and Hubbard models.")
    check(StdI)


def build_momentum_trans_sym(StdI: StdIntList) -> SolverFileData:
    """Build ``qptransidx.def``: the L translations and their phase weights."""
    L = StdI.L
    k = StdI.MomentumIndex
    sep = "=============================================\n"
    lines = [
        # HPhi's reader skips comment lines; this records the sector.
        f"# MomentumIndex {k}\n",
        sep,
        f"NQPTrans {L:10d}\n",
        sep,
        "======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n",
        sep,
    ]
    for op in range(L):
        angle = -2.0 * math.pi * float(k) * float(op) / float(L)
        re = math.cos(angle)
        im = math.sin(angle)
        # Snap rounding noise to +0.0 so "-0.000..." never appears.
        if abs(re) < 1.0e-14:
            re = 0.0
        if abs(im) < 1.0e-14:
            im = 0.0
        lines.append(f"{op} {re:25.15f} {im:25.15f}\n")
    for op in range(L):
        for site in range(L):
            lines.append(f"{op:5d}  {site:5d}  {(site + op) % L:5d}  {1:5d}\n")
    return SolverFileData("qptransidx.def", "".join(lines))
