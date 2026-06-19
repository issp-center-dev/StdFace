"""HPhi solver plugin.

Encapsulates all HPhi-specific keyword parsing, field reset tables,
post-lattice hooks (LargeValue, Boost), and Expert-mode file writing.
"""
from __future__ import annotations

import numpy as np

from ...plugin import ExpertModeSolverPlugin, register
from ...core.stdface_vals import StdIntList, SolverType, MethodType
from ...core.keyword_parser import (
    store_with_check_dup_s, store_with_check_dup_sl,
    store_with_check_dup_i, store_with_check_dup_d,
)
from .writer import large_value, print_calc_mod, print_excitation, print_pump


class HPhiPlugin(ExpertModeSolverPlugin):
    """Plugin for the HPhi exact-diagonalisation / Lanczos solver."""

    @property
    def name(self) -> str:
        return SolverType.HPhi

    @property
    def keyword_table(self) -> dict[str, tuple]:
        return _HPHI_KEYWORDS

    def set_defaults(self, StdI: StdIntList) -> None:
        from ...writer.common_writer import _check_mod_para_hphi
        _check_mod_para_hphi(StdI)

    def modpara_lines(self, StdI: StdIntList) -> list:
        from ...writer.common_writer import _modpara_lines_hphi
        return _modpara_lines_hphi(StdI)

    def namelist_entries(self, StdI: StdIntList) -> list:
        from ...writer.common_writer import _namelist_entries_hphi
        return _namelist_entries_hphi(StdI)

    @property
    def reset_scalars(self) -> list[tuple[str, object]]:
        return _RESET_SCALARS

    @property
    def reset_arrays(self) -> list[tuple[str, object]]:
        return _RESET_ARRAYS

    def post_lattice(self, StdI: StdIntList) -> None:
        """Compute LargeValue and optionally run Boost builder."""
        from ...lattice import get_lattice
        from ...writer.common_writer import unsupported_system

        large_value(StdI)

        if StdI.lBoost == 1:
            try:
                lattice_plugin = get_lattice(StdI.lattice)
            except KeyError:
                unsupported_system(StdI.model, StdI.lattice)
            else:
                lattice_plugin.boost(StdI)

    def write(self, StdI: StdIntList) -> None:
        """Write HPhi Expert-mode files (common files + excitation/pump/calcmod)."""
        from ...writer.common_writer import print_namelist
        self._write_common_files(StdI)
        print_excitation(StdI)
        if StdI.method == MethodType.TIME_EVOLUTION:
            print_pump(StdI)
        print_calc_mod(StdI)
        print_namelist(StdI)


# -----------------------------------------------------------------------
#  Keyword table
# -----------------------------------------------------------------------

_HPHI_KEYWORDS: dict[str, tuple] = {
    "calcspec":        (store_with_check_dup_sl, "CalcSpec"),
    "exct":            (store_with_check_dup_i,  "exct"),
    "eigenvecio":      (store_with_check_dup_sl, "EigenVecIO"),
    "expandcoef":      (store_with_check_dup_i,  "ExpandCoef"),
    "expecinterval":   (store_with_check_dup_i,  "ExpecInterval"),
    "cdatafilehead":   (store_with_check_dup_s,  "CDataFileHead"),
    "dt":              (store_with_check_dup_d,  "dt"),
    "flgtemp":         (store_with_check_dup_i,  "FlgTemp"),
    "freq":            (store_with_check_dup_d,  "freq"),
    "hamio":           (store_with_check_dup_sl, "HamIO"),
    "initialvectype":  (store_with_check_dup_sl, "InitialVecType"),
    "initial_iv":      (store_with_check_dup_i,  "initial_iv"),
    "lanczoseps":      (store_with_check_dup_i,  "LanczosEps"),
    "lanczostarget":   (store_with_check_dup_i,  "LanczosTarget"),
    "lanczos_max":     (store_with_check_dup_i,  "Lanczos_max"),
    "largevalue":      (store_with_check_dup_d,  "LargeValue"),
    "method":          (store_with_check_dup_sl, "method"),
    "nomega":          (store_with_check_dup_i,  "Nomega"),
    "numave":          (store_with_check_dup_i,  "NumAve"),
    "nvec":            (store_with_check_dup_i,  "nvec"),
    "omegamax":        (store_with_check_dup_d,  "OmegaMax"),
    "omegamin":        (store_with_check_dup_d,  "OmegaMin"),
    "omegaorg":        (store_with_check_dup_d,  "OmegaOrg"),
    "omegaim":         (store_with_check_dup_d,  "OmegaIm"),
    "outputexcitedvec": (store_with_check_dup_sl, "OutputExVec"),
    "pumptype":        (store_with_check_dup_sl, "PumpType"),
    "restart":         (store_with_check_dup_sl, "Restart"),
    "spectrumqh":      (store_with_check_dup_d, "SpectrumQ", 2, float),
    "spectrumql":      (store_with_check_dup_d, "SpectrumQ", 1, float),
    "spectrumqw":      (store_with_check_dup_d, "SpectrumQ", 0, float),
    "spectrumtype":    (store_with_check_dup_sl, "SpectrumType"),
    "tdump":           (store_with_check_dup_d,  "tdump"),
    "tshift":          (store_with_check_dup_d,  "tshift"),
    "uquench":         (store_with_check_dup_d,  "Uquench"),
    "vecpoth":         (store_with_check_dup_d, "VecPot", 2, float),
    "vecpotl":         (store_with_check_dup_d, "VecPot", 1, float),
    "vecpotw":         (store_with_check_dup_d, "VecPot", 0, float),
    "2s":              (store_with_check_dup_i,  "S2"),
    "ngpu":            (store_with_check_dup_i,  "NGPU"),
    "scalapack":       (store_with_check_dup_i,  "Scalapack"),
}

# -----------------------------------------------------------------------
#  Reset tables
# -----------------------------------------------------------------------

_RESET_SCALARS: list[tuple[str, object]] = [
    ("LargeValue", None),
    ("OmegaMax", None),
    ("OmegaMin", None),
    ("OmegaOrg", None),
    ("OmegaIm", None),
    ("Nomega", None),
    ("FlgTemp", 1),
    ("Lanczos_max", None),
    ("initial_iv", None),
    ("nvec", None),
    ("exct", None),
    ("LanczosEps", None),
    ("LanczosTarget", None),
    ("NumAve", None),
    ("ExpecInterval", None),
    ("dt", None),
    ("tdump", None),
    ("tshift", None),
    ("freq", None),
    ("Uquench", None),
    ("ExpandCoef", None),
    ("NGPU", None),
    ("Scalapack", None),
]

_RESET_ARRAYS: list[tuple[str, object]] = [
    ("SpectrumQ", np.nan),
    ("VecPot", np.nan),
]


# Auto-register on import
register(HPhiPlugin())
