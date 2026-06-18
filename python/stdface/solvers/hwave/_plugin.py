"""H-wave solver plugin.

Encapsulates all H-wave-specific keyword parsing, field reset tables,
and Expert-mode file writing.
"""
from __future__ import annotations

from ...plugin import SolverPlugin, register
from ...core.stdface_vals import StdIntList, SolverType, NaN_i, NaN_d
from ...core.keyword_parser import (
    store_with_check_dup_i, store_with_check_dup_d, store_with_check_dup_sl,
    _grid3x3_keywords,
)
from ...writer.wannier90_writer import export_geometry, export_interaction


class HWavePlugin(SolverPlugin):
    """Plugin for the H-wave solver."""

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

    def write(self, StdI: StdIntList) -> None:
        """Write H-wave output files.

        Overrides the template method entirely because H-wave has two
        completely different output modes (uhfr vs wannier90 export).
        """
        from ...writer.common_writer import (
            print_trans, print_1_green, check_output_mode, check_mod_para,
        )
        from ...writer.interaction_writer import print_interactions

        if StdI.calcmode == "uhfr":
            print_trans(StdI)
            print_interactions(StdI)
            check_mod_para(StdI)
            check_output_mode(StdI)
            print_1_green(StdI)
        else:
            export_geometry(StdI)
            export_interaction(StdI)


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

_HWAVE_KEYWORDS: dict[str, tuple] = {
    **_UHF_BASE_KEYWORDS,
    "calcmode":   (store_with_check_dup_sl, "calcmode"),
    "fileprefix": (store_with_check_dup_sl, "fileprefix"),
    "exportall":  (store_with_check_dup_i,  "export_all"),
    "lattice_gp": (store_with_check_dup_i,  "lattice_gp"),
}

# -----------------------------------------------------------------------
#  Reset tables
# -----------------------------------------------------------------------

_UHF_BASE_SCALARS: list[tuple[str, object]] = [
    ("NMPTrans", NaN_i),
    ("RndSeed", NaN_i),
    ("mix", NaN_d),
    ("eps", NaN_i),
    ("eps_slater", NaN_i),
    ("Iteration_max", NaN_i),
    ("Hsub", NaN_i),
    ("Lsub", NaN_i),
    ("Wsub", NaN_i),
]

_RESET_SCALARS: list[tuple[str, object]] = _UHF_BASE_SCALARS + [
    ("export_all", NaN_i),
    ("lattice_gp", NaN_i),
]

_RESET_ARRAYS: list[tuple[str, object]] = [
    ("boxsub", NaN_i),
]


# Auto-register on import
register(HWavePlugin())
