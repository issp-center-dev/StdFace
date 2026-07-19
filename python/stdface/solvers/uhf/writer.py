"""UHF Expert-mode output bodies (moved from ``writer/common_writer.py``).

The defaults in :func:`set_modpara_defaults` are shared with the H-wave
family, which carries its own copy per the solver self-containment policy
(see CLAUDE.md); behavioural parity is enforced by
``test/unit/test_solver_table_parity.py``.
"""
from __future__ import annotations

import logging

from ...core.stdface_vals import StdIntList, NaN_i
from ...core.param_check import print_val_d, print_val_i

logger = logging.getLogger(__name__)


def modpara_lines(StdI: StdIntList) -> list:
    """Return the UHF body line descriptors of ``modpara.def``."""
    lines: list = [
        ("raw", "UHF_Cal_Parameters"),
        ("sep",),
        ("kv", "CDataFileHead", StdI.CDataFileHead, ""),
        ("kv", "CParaFileHead", "zqp", ""),
        ("sep",),
        ("kv", "Nsite", StdI.nsite, ""),
    ]
    if StdI.Sz2 is not None:
        lines.append(("kv", "2Sz", StdI.Sz2, "<5d"))
    # UHF emits an unset Ncond as the integer sentinel (legacy format).
    ncond_out = StdI.ncond if StdI.ncond is not None else NaN_i
    lines += [
        ("kv", "Ncond", ncond_out, "<5d"),
        ("kv", "IterationMax", StdI.Iteration_max, ""),
        ("kv", "EPS", StdI.eps, ""),
        ("kv", "Mix", StdI.mix, ".10f"),
        ("kv", "RndSeed", StdI.RndSeed, ""),
        ("kv", "EpsSlater", StdI.eps_slater, ""),
        ("kv", "NMPTrans", StdI.NMPTrans, ""),
    ]
    return lines


def set_modpara_defaults(StdI: StdIntList) -> None:
    """Set UHF-specific default model parameters.

    Handles random seed, iteration limit, mixing, convergence, and
    symmetry-projection parameters.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure, modified in place.
    """
    StdI.RndSeed = print_val_i("RndSeed", StdI.RndSeed, 123456789)
    StdI.Iteration_max = print_val_i("Iteration_max", StdI.Iteration_max, 1000)
    StdI.mix = print_val_d("Mix", StdI.mix, 0.5)
    StdI.eps = print_val_i("eps", StdI.eps, 8)
    StdI.eps_slater = print_val_i("EpsSlater", StdI.eps_slater, 6)
    StdI.NMPTrans = print_val_i("NMPTrans", StdI.NMPTrans, 0)
