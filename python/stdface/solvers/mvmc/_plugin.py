"""mVMC solver plugin.

Encapsulates all mVMC-specific keyword parsing, field reset tables,
and Expert-mode file writing.
"""
from __future__ import annotations

from ...plugin import ExpertModeSolverPlugin, register, register_config
from ...core.stdface_vals import StdIntList, SolverType, NaN_i
from .config import MVMCConfig
from ...core.keyword_parser import (
    store_with_check_dup_s, store_with_check_dup_i, store_with_check_dup_d,
    _grid3x3_keywords,
)
from ...core.param_check import print_val_i
from .variational import generate_orb, proj, build_jastrow
from .writer import build_orb, build_orb_para, build_gutzwiller


class MVMCPlugin(ExpertModeSolverPlugin):
    """Plugin for the mVMC variational Monte Carlo solver."""

    @property
    def name(self) -> str:
        return SolverType.mVMC

    @property
    def keyword_table(self) -> dict[str, tuple]:
        return _MVMC_KEYWORDS

    def set_defaults(self, StdI: StdIntList) -> None:
        from .writer import set_modpara_defaults
        set_modpara_defaults(StdI)

    def modpara_lines(self, StdI: StdIntList) -> list:
        from .writer import modpara_lines
        return modpara_lines(StdI)

    def namelist_entries(self, StdI: StdIntList) -> list:
        from .writer import namelist_entries
        return namelist_entries(StdI)

    @property
    def reset_scalars(self) -> list[tuple[str, object]]:
        return _RESET_SCALARS

    @property
    def reset_arrays(self) -> list[tuple[str, object]]:
        return _RESET_ARRAYS

    def build_solver_files(self, StdI: StdIntList) -> list:
        """Build the mVMC-specific variational files."""
        if StdI.lGC == 0 and (StdI.Sz2 == 0 or StdI.Sz2 is None):
            StdI.ComplexType = print_val_i("ComplexType", StdI.ComplexType, 0)
        else:
            StdI.ComplexType = print_val_i("ComplexType", StdI.ComplexType, 1)

        generate_orb(StdI)
        files = [proj(StdI), build_jastrow(StdI)]
        if StdI.lGC == 1 or (StdI.Sz2 != 0 and StdI.Sz2 is not None):
            files.extend(build_orb_para(StdI))
        files.append(build_gutzwiller(StdI))
        files.append(build_orb(StdI))
        return files


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
    ("CParaFileHead", None),
    ("NVMCCalMode", None),
    ("NLanczosMode", None),
    ("NDataIdxStart", None),
    ("NDataQtySmp", None),
    ("NSPGaussLeg", None),
    ("NSPStot", None),
    ("NMPTrans", None),
    ("NSROptItrStep", None),
    ("NSROptItrSmp", None),
    ("DSROptRedCut", None),
    ("DSROptStaDel", None),
    ("DSROptStepDt", None),
    ("NVMCWarmUp", None),
    ("NVMCInterval", None),
    ("NVMCSample", None),
    ("NExUpdatePath", None),
    ("RndSeed", None),
    ("NSplitSize", None),
    ("NStore", None),
    ("NSRCG", None),
    ("ComplexType", None),
    ("Hsub", None),
    ("Lsub", None),
    ("Wsub", None),
]

_RESET_ARRAYS: list[tuple[str, object]] = [
    ("boxsub", NaN_i),
]


# Auto-register on import
register(MVMCPlugin())
register_config(SolverType.mVMC, MVMCConfig)
