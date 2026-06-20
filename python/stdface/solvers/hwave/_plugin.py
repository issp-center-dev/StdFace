"""H-wave solver plugins (UHFR / UHFK).

The user passes ``--solver HWAVE``; ``_resolve_solver_name`` then maps it to
:class:`UHFRPlugin` (real-space ``.def`` output) or :class:`UHFKPlugin`
(Wannier90 export) based on ``calcmode``.  UHFR is **not** an
``ExpertModeSolverPlugin`` (it writes no modpara/namelist/locspn).
"""
from __future__ import annotations

from ...plugin import SolverPlugin, WannierModeSolverPlugin, register, register_config
from ...core.stdface_vals import StdIntList, SolverType, NaN_i
from .config import HWaveConfig
from ...core.keyword_parser import (
    store_with_check_dup_i, store_with_check_dup_d, store_with_check_dup_sl,
    _grid3x3_keywords,
)


# ---------------------------------------------------------------------------
#  Output strategies (C4-1: extracted so HWavePlugin can dispatch by calcmode)
# ---------------------------------------------------------------------------

def build_uhfr_output(StdI: StdIntList):
    """Assemble the UHFR real-space ``.def`` output (trans / interactions / green1)."""
    from ...writer.common_writer import (
        build_trans, check_output_mode, check_mod_para, build_green_one,
    )
    from ...writer.interaction_writer import build_interactions
    from ...core.output import ExpertModeOutput

    trans = build_trans(StdI)
    interactions = build_interactions(StdI)
    check_mod_para(StdI)
    check_output_mode(StdI)
    green_one = build_green_one(StdI)
    return ExpertModeOutput(
        trans=trans, interactions=interactions, green_one=green_one)


class UHFRPlugin(SolverPlugin):
    """H-wave real-space UHF mode (``uhfr``).

    Writes ``.def`` files but, unlike Expert mode, no ``modpara.def`` /
    ``namelist.def`` / ``locspn.def`` — so it inherits ``SolverPlugin``
    directly rather than ``ExpertModeSolverPlugin``.
    """

    @property
    def name(self) -> str:
        return SolverType.UHFR

    @property
    def keyword_table(self) -> dict[str, tuple]:
        return _UHFR_KEYWORDS

    @property
    def reset_scalars(self) -> list[tuple[str, object]]:
        return _RESET_SCALARS

    @property
    def reset_arrays(self) -> list[tuple[str, object]]:
        return _RESET_ARRAYS

    def set_defaults(self, StdI: StdIntList) -> None:
        from ...writer.common_writer import _check_mod_para_uhf
        _check_mod_para_uhf(StdI)

    def build_output(self, StdI: StdIntList):
        """Assemble the UHFR partial output (trans / interactions / green1)."""
        return build_uhfr_output(StdI)

    def write(self, StdI: StdIntList) -> None:
        self.build_output(StdI).write()


class UHFKPlugin(WannierModeSolverPlugin):
    """H-wave Wannier90 mode (``uhfk`` / ``rpa``); writes geom/transfer files."""

    @property
    def name(self) -> str:
        return SolverType.UHFK

    @property
    def keyword_table(self) -> dict[str, tuple]:
        return _UHFK_KEYWORDS

    @property
    def reset_scalars(self) -> list[tuple[str, object]]:
        return _RESET_SCALARS

    @property
    def reset_arrays(self) -> list[tuple[str, object]]:
        return _RESET_ARRAYS


def _hwave_output_mode(StdI: StdIntList) -> str:
    """Resolve the H-wave output mode from the solver name and ``calcmode``.

    Equivalent to the legacy ``_resolve_solver_name`` mapping:
    ``UHFR`` iff the solver is UHFR, or HWAVE with ``calcmode == "uhfr"``;
    everything else (uhfk / rpa / unset) selects the UHFK Wannier export.
    """
    if (StdI.solver == SolverType.UHFR
            or (StdI.solver == SolverType.HWAVE and StdI.calcmode == "uhfr")):
        return SolverType.UHFR
    return SolverType.UHFK


class HWavePlugin(SolverPlugin):
    """H-wave family plugin (C4).

    Owns the H-wave *input* (keyword / reset / config = the UHFR/UHFK union)
    and dispatches *output* by ``calcmode`` to the UHFR (.def) or UHFK
    (Wannier90) strategy.  Registered under ``HWAVE`` (and, from C4-2b, also
    ``UHFR`` / ``UHFK``) so parsing and reset run before solver-name
    resolution without a core fallback.
    """

    @property
    def name(self) -> str:
        return SolverType.HWAVE

    @property
    def keyword_table(self) -> dict[str, tuple]:
        return _HWAVE_KEYWORDS

    @property
    def reset_scalars(self) -> list[tuple[str, object]]:
        return _RESET_SCALARS

    @property
    def reset_arrays(self) -> list[tuple[str, object]]:
        return _RESET_ARRAYS

    def set_defaults(self, StdI: StdIntList) -> None:
        from ...writer.common_writer import _check_mod_para_uhf
        _check_mod_para_uhf(StdI)

    def build_output(self, StdI: StdIntList):
        if _hwave_output_mode(StdI) == SolverType.UHFR:
            return build_uhfr_output(StdI)
        from ...core.output import build_wannier_output
        return build_wannier_output(StdI)

    def write(self, StdI: StdIntList) -> None:
        self.build_output(StdI).write()


# -----------------------------------------------------------------------
#  Keyword table
# -----------------------------------------------------------------------

_BOXSUB_KEYWORDS: dict[str, tuple] = _grid3x3_keywords(
    "{a}{c}sub", "boxsub", store_with_check_dup_i, int
)

_UHF_BASE_KEYWORDS: dict[str, tuple] = {
    "iteration_max": (store_with_check_dup_i, "Iteration_max"),
    "rndseed":       (store_with_check_dup_i, "RndSeed"),
    "nmptrans":      (store_with_check_dup_i, "NMPTrans"),
    **_BOXSUB_KEYWORDS,
    "hsub":          (store_with_check_dup_i, "Hsub"),
    "lsub":          (store_with_check_dup_i, "Lsub"),
    "wsub":          (store_with_check_dup_i, "Wsub"),
    "eps":           (store_with_check_dup_i, "eps"),
    "epsslater":     (store_with_check_dup_i, "eps_slater"),
    "mix":           (store_with_check_dup_d, "mix"),
}

# UHFR accepts the same keywords as UHF; UHFK only the Wannier-export keys.
# (``calcmode`` is parsed via the pre-resolution HWAVE table in keyword_parser.)
_UHFR_KEYWORDS: dict[str, tuple] = _UHF_BASE_KEYWORDS

_UHFK_KEYWORDS: dict[str, tuple] = {
    "fileprefix": (store_with_check_dup_sl, "fileprefix"),
    "exportall":  (store_with_check_dup_i,  "export_all"),
    "lattice_gp": (store_with_check_dup_i,  "lattice_gp"),
}

# H-wave family union: UHF base + calcmode + the UHFK Wannier-export keys.
# Used by HWavePlugin (parse runs while solver == "HWAVE", before resolution).
_HWAVE_KEYWORDS: dict[str, tuple] = {
    **_UHF_BASE_KEYWORDS,
    "calcmode":   (store_with_check_dup_sl, "calcmode"),
    **_UHFK_KEYWORDS,
}

# -----------------------------------------------------------------------
#  Reset tables
# -----------------------------------------------------------------------
#
# Field resets actually run pre-resolution (solver == HWAVE) via the legacy
# fallback in stdface_main; these tables satisfy the plugin interface and
# cover the union of UHFR/UHFK fields.

_UHF_BASE_SCALARS: list[tuple[str, object]] = [
    ("NMPTrans", None),
    ("RndSeed", None),
    ("mix", None),
    ("eps", None),
    ("eps_slater", None),
    ("Iteration_max", None),
    ("Hsub", None),
    ("Lsub", None),
    ("Wsub", None),
]

_RESET_SCALARS: list[tuple[str, object]] = _UHF_BASE_SCALARS + [
    ("export_all", None),
    ("lattice_gp", None),
]

_RESET_ARRAYS: list[tuple[str, object]] = [
    ("boxsub", NaN_i),
]


# Auto-register on import
register(UHFRPlugin())
register(UHFKPlugin())
# C4-2a: HWavePlugin owns H-wave input (keyword/reset).  Registering it under
# HWAVE means parse/reset (which run while solver == "HWAVE", before
# _resolve_solver_name) go through the plugin instead of the core fallback.
# Output still flows through UHFR/UHFK after resolution until C4-2b.
register(HWavePlugin())
# One config shared by UHFR/UHFK; also registered under the raw HWAVE alias so
# it attaches before _resolve_solver_name splits HWAVE into UHFR/UHFK.
register_config(SolverType.HWAVE, HWaveConfig)
register_config(SolverType.UHFR, HWaveConfig)
register_config(SolverType.UHFK, HWaveConfig)
