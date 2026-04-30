"""mVMC solver plugin.

Encapsulates all mVMC-specific keyword parsing, field reset tables,
and Expert-mode file writing.
"""
from __future__ import annotations

from ...plugin import SolverPlugin, register
from ...core.stdface_vals import StdIntList, SolverType, NaN_i, NaN_d
from ...core.keyword_parser import (
    store_with_check_dup_s, store_with_check_dup_i, store_with_check_dup_d,
    _grid3x3_keywords,
)
from ...core.param_check import print_val_i
from .variational import generate_orb, proj, print_jastrow
from .writer import print_orb, print_orb_para, print_gutzwiller


class MVMCPlugin(SolverPlugin):
    """Plugin for the mVMC variational Monte Carlo solver."""

    @property
    def name(self) -> str:
        return SolverType.mVMC

    @property
    def keyword_table(self) -> dict[str, tuple]:
        return _MVMC_KEYWORDS

    @property
    def reset_scalars(self) -> list[tuple[str, object]]:
        return _RESET_SCALARS

    @property
    def reset_arrays(self) -> list[tuple[str, object]]:
        return _RESET_ARRAYS

    def write_solver_specific(self, StdI: StdIntList) -> None:
        """Write mVMC-specific variational parameter files."""
        if StdI.lGC == 0 and (StdI.Sz2 == 0 or StdI.Sz2 == NaN_i):
            StdI.ComplexType = print_val_i("ComplexType", StdI.ComplexType, 0)
        else:
            StdI.ComplexType = print_val_i("ComplexType", StdI.ComplexType, 1)

        generate_orb(StdI)
        proj(StdI)
        print_jastrow(StdI)
        if StdI.lGC == 1 or (StdI.Sz2 != 0 and StdI.Sz2 != NaN_i):
            print_orb_para(StdI)
        print_gutzwiller(StdI)
        print_orb(StdI)


# -----------------------------------------------------------------------
#  Keyword table
# -----------------------------------------------------------------------

_BOXSUB_KEYWORDS: dict[str, tuple] = _grid3x3_keywords(
    "{a}{c}sub", "boxsub", store_with_check_dup_i, int
)

_MVMC_KEYWORDS: dict[str, tuple] = {
    **_BOXSUB_KEYWORDS,
    "complextype":    (store_with_check_dup_i, "ComplexType"),
    "cparafilehead":  (store_with_check_dup_s, "CParaFileHead"),
    "dsroptredcut":   (store_with_check_dup_d, "DSROptRedCut"),
    "dsroptstadel":   (store_with_check_dup_d, "DSROptStaDel"),
    "dsroptstepdt":   (store_with_check_dup_d, "DSROptStepDt"),
    "hsub":           (store_with_check_dup_i, "Hsub"),
    "lsub":           (store_with_check_dup_i, "Lsub"),
    "nvmccalmode":    (store_with_check_dup_i, "NVMCCalMode"),
    "ndataidxstart":  (store_with_check_dup_i, "NDataIdxStart"),
    "ndataqtysmp":    (store_with_check_dup_i, "NDataQtySmp"),
    "nlanczosmode":   (store_with_check_dup_i, "NLanczosMode"),
    "nmptrans":       (store_with_check_dup_i, "NMPTrans"),
    "nspgaussleg":    (store_with_check_dup_i, "NSPGaussLeg"),
    "nsplitsize":     (store_with_check_dup_i, "NSplitSize"),
    "nspstot":        (store_with_check_dup_i, "NSPStot"),
    "nsroptitrsmp":   (store_with_check_dup_i, "NSROptItrSmp"),
    "nsroptitrstep":  (store_with_check_dup_i, "NSROptItrStep"),
    "nstore":         (store_with_check_dup_i, "NStore"),
    "nsrcg":          (store_with_check_dup_i, "NSRCG"),
    "nvmcinterval":   (store_with_check_dup_i, "NVMCInterval"),
    "nvmcsample":     (store_with_check_dup_i, "NVMCSample"),
    "nvmcwarmup":     (store_with_check_dup_i, "NVMCWarmUp"),
    "rndseed":        (store_with_check_dup_i, "RndSeed"),
    "wsub":           (store_with_check_dup_i, "Wsub"),
}

# -----------------------------------------------------------------------
#  Reset tables
# -----------------------------------------------------------------------

_RESET_SCALARS: list[tuple[str, object]] = [
    ("NVMCCalMode", NaN_i),
    ("NLanczosMode", NaN_i),
    ("NDataIdxStart", NaN_i),
    ("NDataQtySmp", NaN_i),
    ("NSPGaussLeg", NaN_i),
    ("NSPStot", NaN_i),
    ("NMPTrans", NaN_i),
    ("NSROptItrStep", NaN_i),
    ("NSROptItrSmp", NaN_i),
    ("DSROptRedCut", NaN_d),
    ("DSROptStaDel", NaN_d),
    ("DSROptStepDt", NaN_d),
    ("NVMCWarmUp", NaN_i),
    ("NVMCInterval", NaN_i),
    ("NVMCSample", NaN_i),
    ("NExUpdatePath", NaN_i),
    ("RndSeed", NaN_i),
    ("NSplitSize", NaN_i),
    ("NStore", NaN_i),
    ("NSRCG", NaN_i),
    ("ComplexType", NaN_i),
    ("Hsub", NaN_i),
    ("Lsub", NaN_i),
    ("Wsub", NaN_i),
]

_RESET_ARRAYS: list[tuple[str, object]] = [
    ("boxsub", NaN_i),
]


# Auto-register on import
register(MVMCPlugin())
