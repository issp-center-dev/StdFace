"""H-wave solver plugins (UHFR / UHFK).

The user passes ``--solver HWAVE``; ``_resolve_solver_name`` then maps it to
:class:`UHFRPlugin` (real-space ``.def`` output) or :class:`UHFKPlugin`
(Wannier90 export) based on ``calcmode``.  UHFR is **not** an
``ExpertModeSolverPlugin`` (it writes no modpara/namelist/locspn).
"""
from __future__ import annotations

from ...plugin import SolverPlugin, WannierModeSolverPlugin, register
from ...core.stdface_vals import StdIntList, SolverType, NaN_i
from ...core.keyword_parser import (
    store_with_check_dup_i, store_with_check_dup_d, store_with_check_dup_sl,
    _grid3x3_keywords,
)


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

    def write(self, StdI: StdIntList) -> None:
        from ...writer.common_writer import (
            print_trans, print_1_green, check_output_mode, check_mod_para,
        )
        from ...writer.interaction_writer import print_interactions

        print_trans(StdI)
        print_interactions(StdI)
        check_mod_para(StdI)
        check_output_mode(StdI)
        print_1_green(StdI)


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
