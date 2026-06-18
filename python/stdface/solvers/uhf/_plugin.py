"""UHF solver plugin.

Encapsulates all UHF-specific keyword parsing, field reset tables,
and Expert-mode file writing.
"""
from __future__ import annotations

from ...plugin import ExpertModeSolverPlugin, register
from ...core.stdface_vals import StdIntList, SolverType, NaN_i, NaN_d
from ...core.keyword_parser import (
    store_with_check_dup_i, store_with_check_dup_d,
    _grid3x3_keywords,
)


class UHFPlugin(ExpertModeSolverPlugin):
    """Plugin for the UHF (unrestricted Hartree-Fock) solver."""

    @property
    def name(self) -> str:
        return SolverType.UHF

    @property
    def keyword_table(self) -> dict[str, tuple]:
        return _UHF_KEYWORDS

    @property
    def reset_scalars(self) -> list[tuple[str, object]]:
        return _RESET_SCALARS

    @property
    def reset_arrays(self) -> list[tuple[str, object]]:
        return _RESET_ARRAYS

    def set_defaults(self, StdI: StdIntList) -> None:
        from ...writer.common_writer import _check_mod_para_uhf
        _check_mod_para_uhf(StdI)

    def write_modpara_body(self, fp, StdI: StdIntList) -> None:
        from ...writer.common_writer import _write_modpara_uhf_hwave
        _write_modpara_uhf_hwave(fp, StdI)

    def has_two_body_green(self, StdI: StdIntList) -> bool:
        return False

    def write_green(self, StdI: StdIntList) -> None:
        """Write only greenone.def (UHF does not use greentwo)."""
        from ...writer.common_writer import print_1_green
        print_1_green(StdI)


# -----------------------------------------------------------------------
#  Keyword table
# -----------------------------------------------------------------------

_BOXSUB_KEYWORDS: dict[str, tuple] = _grid3x3_keywords(
    "{a}{c}sub", "boxsub", store_with_check_dup_i, int
)

_UHF_KEYWORDS: dict[str, tuple] = {
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

# -----------------------------------------------------------------------
#  Reset tables
# -----------------------------------------------------------------------

_RESET_SCALARS: list[tuple[str, object]] = [
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

_RESET_ARRAYS: list[tuple[str, object]] = [
    ("boxsub", NaN_i),
]


# Auto-register on import
register(UHFPlugin())
