"""
Read input file and write files for Expert mode; initialize variables; check parameters.

This module is the Python translation of ``StdFace_main.c``.  It provides the
top-level entry point for the Standard-mode input generator used by HPhi, mVMC,
UHF, and H-wave.

The following lattices are supported:

- 1D Chain
- 1D Ladder
- 2D Tetragonal
- 2D Triangular
- 2D Honeycomb
- 2D Kagome
- 3D Simple Orthorhombic
- 3D Face Centered Orthorhombic
- 3D Pyrochlore

License
-------
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""

from __future__ import annotations

import cmath
import math
import sys

import numpy as np

from stdface_vals import StdIntList
from stdface_model_util import (
    exit_program,
    not_used_i,
    print_val_d,
    print_val_i,
    required_val_i,
    generate_orb,
    proj,
    print_jastrow,
)
from export_wannier90 import export_geometry, export_interaction
import chain_lattice
import square_lattice
import ladder
import triangular_lattice
import honeycomb_lattice
import kagome
import orthorhombic
import fc_ortho
import pyrochlore
import wannier90 as wannier90_mod

# ---------------------------------------------------------------------------
#  Sentinel constants (matching the C code)
# ---------------------------------------------------------------------------
NaN_i: int = 2147483647
"""Sentinel for an unset integer parameter (same as ``INT_MAX`` in C)."""

NaN_d: float = float("nan")
"""Sentinel for an unset float parameter (IEEE NaN)."""

NaN_c: complex = complex(float("nan"), 0.0)
"""Sentinel for an unset complex parameter (real part is NaN)."""


# ===================================================================
#  _reset_vals
# ===================================================================


def _reset_vals(StdI: StdIntList) -> None:
    """Clear / initialize every field in *StdI* to its sentinel value.

    This is the Python translation of the C function
    ``StdFace_ResetVals()``.  Fields that have not been specified by
    the user are filled with NaN sentinels so that duplicate-input
    detection and default-value assignment work correctly later.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure whose fields are reset **in
        place**.
    """
    # NaN sentinel and mathematical constants
    StdI.NaN_i = NaN_i
    StdI.pi = math.acos(-1.0)

    # --- Lattice scalars / vectors -----------------------------------------
    StdI.a = NaN_d
    StdI.length[:] = NaN_d
    StdI.box[:, :] = NaN_i
    StdI.direct[:, :] = NaN_d

    # --- Magnetic field ----------------------------------------------------
    StdI.Gamma = NaN_d
    StdI.Gamma_y = NaN_d
    StdI.h = NaN_d

    # --- Lattice dimensions ------------------------------------------------
    StdI.Height = NaN_i
    StdI.L = NaN_i
    StdI.W = NaN_i

    # --- Isotropic / anisotropic scalar spin couplings ---------------------
    StdI.JAll = NaN_d
    StdI.JpAll = NaN_d
    StdI.JppAll = NaN_d
    StdI.J0All = NaN_d
    StdI.J0pAll = NaN_d
    StdI.J0ppAll = NaN_d
    StdI.J1All = NaN_d
    StdI.J1pAll = NaN_d
    StdI.J1ppAll = NaN_d
    StdI.J2All = NaN_d
    StdI.J2pAll = NaN_d
    StdI.J2ppAll = NaN_d

    # --- 3x3 spin coupling matrices ----------------------------------------
    StdI.J[:, :] = NaN_d
    StdI.Jp[:, :] = NaN_d
    StdI.Jpp[:, :] = NaN_d
    StdI.J0[:, :] = NaN_d
    StdI.J0p[:, :] = NaN_d
    StdI.J0pp[:, :] = NaN_d
    StdI.J1[:, :] = NaN_d
    StdI.J1p[:, :] = NaN_d
    StdI.J1pp[:, :] = NaN_d
    StdI.J2[:, :] = NaN_d
    StdI.J2p[:, :] = NaN_d
    StdI.J2pp[:, :] = NaN_d

    # D matrix: zero everywhere except D[2][2] = NaN_d
    StdI.D[:, :] = 0.0
    StdI.D[2, 2] = NaN_d

    StdI.K = NaN_d

    # --- Chemical potential / spin -----------------------------------------
    StdI.mu = NaN_d
    StdI.S2 = NaN_i

    # --- Hopping parameters (complex) --------------------------------------
    StdI.t = NaN_c
    StdI.tp = NaN_c
    StdI.tpp = NaN_c
    StdI.t0 = NaN_c
    StdI.t0p = NaN_c
    StdI.t0pp = NaN_c
    StdI.t1 = NaN_c
    StdI.t1p = NaN_c
    StdI.t1pp = NaN_c
    StdI.t2 = NaN_c
    StdI.t2p = NaN_c
    StdI.t2pp = NaN_c

    # --- Coulomb parameters (float) ----------------------------------------
    StdI.U = NaN_d
    StdI.V = NaN_d
    StdI.Vp = NaN_d
    StdI.Vpp = NaN_d
    StdI.V0 = NaN_d
    StdI.V0p = NaN_d
    StdI.V0pp = NaN_d
    StdI.V1 = NaN_d
    StdI.V1p = NaN_d
    StdI.V1pp = NaN_d
    StdI.V2 = NaN_d
    StdI.V2p = NaN_d
    StdI.V2pp = NaN_d

    # --- Phase / boundary --------------------------------------------------
    StdI.phase[:] = NaN_d
    StdI.pi180 = StdI.pi / 180.0

    # --- Calculation conditions (strings / ints) ---------------------------
    StdI.ncond = NaN_i
    StdI.Sz2 = NaN_i
    StdI.model = "****"
    StdI.lattice = "****"
    StdI.outputmode = "****"
    StdI.CDataFileHead = "****"
    StdI.double_counting_mode = "****"

    # --- Wannier90 cutoffs -------------------------------------------------
    StdI.cutoff_t = NaN_d
    StdI.cutoff_u = NaN_d
    StdI.cutoff_j = NaN_d
    StdI.cutoff_length_t = NaN_d
    StdI.cutoff_length_U = NaN_d
    StdI.cutoff_length_J = NaN_d
    StdI.lambda_ = NaN_d
    StdI.lambda_U = NaN_d
    StdI.lambda_J = NaN_d
    StdI.alpha = NaN_d
    StdI.cutoff_tR[:] = NaN_i
    StdI.cutoff_UR[:] = NaN_i
    StdI.cutoff_JR[:] = NaN_i
    StdI.cutoff_tVec[:, :] = NaN_d
    StdI.cutoff_UVec[:, :] = NaN_d
    StdI.cutoff_JVec[:, :] = NaN_d

    # --- HPhi-specific fields ----------------------------------------------
    if StdI.solver == "HPhi":
        StdI.LargeValue = NaN_d
        StdI.OmegaMax = NaN_d
        StdI.OmegaMin = NaN_d
        StdI.OmegaOrg = NaN_d
        StdI.OmegaIm = NaN_d
        StdI.Nomega = NaN_i
        StdI.SpectrumQ[:] = NaN_d
        StdI.method = "****"
        StdI.Restart = "****"
        StdI.EigenVecIO = "****"
        StdI.InitialVecType = "****"
        StdI.HamIO = "****"
        StdI.CalcSpec = "****"
        StdI.SpectrumType = "****"
        StdI.OutputExVec = "****"
        StdI.FlgTemp = 1
        StdI.Lanczos_max = NaN_i
        StdI.initial_iv = NaN_i
        StdI.nvec = NaN_i
        StdI.exct = NaN_i
        StdI.LanczosEps = NaN_i
        StdI.LanczosTarget = NaN_i
        StdI.NumAve = NaN_i
        StdI.ExpecInterval = NaN_i
        StdI.dt = NaN_d
        StdI.tdump = NaN_d
        StdI.tshift = NaN_d
        StdI.freq = NaN_d
        StdI.Uquench = NaN_d
        StdI.VecPot[:] = NaN_d
        StdI.PumpType = "****"
        StdI.ExpandCoef = NaN_i
        StdI.NGPU = NaN_i
        StdI.Scalapack = NaN_i

    # --- mVMC-specific fields ----------------------------------------------
    elif StdI.solver == "mVMC":
        StdI.CParaFileHead = "****"
        StdI.NVMCCalMode = NaN_i
        StdI.NLanczosMode = NaN_i
        StdI.NDataIdxStart = NaN_i
        StdI.NDataQtySmp = NaN_i
        StdI.NSPGaussLeg = NaN_i
        StdI.NSPStot = NaN_i
        StdI.NMPTrans = NaN_i
        StdI.NSROptItrStep = NaN_i
        StdI.NSROptItrSmp = NaN_i
        StdI.DSROptRedCut = NaN_d
        StdI.DSROptStaDel = NaN_d
        StdI.DSROptStepDt = NaN_d
        StdI.NVMCWarmUp = NaN_i
        StdI.NVMCInterval = NaN_i
        StdI.NVMCSample = NaN_i
        StdI.NExUpdatePath = NaN_i
        StdI.RndSeed = NaN_i
        StdI.NSplitSize = NaN_i
        StdI.NStore = NaN_i
        StdI.NSRCG = NaN_i
        StdI.ComplexType = NaN_i
        StdI.boxsub[:, :] = NaN_i
        StdI.Hsub = NaN_i
        StdI.Lsub = NaN_i
        StdI.Wsub = NaN_i

    # --- UHF-specific fields -----------------------------------------------
    elif StdI.solver == "UHF":
        StdI.NMPTrans = NaN_i
        StdI.RndSeed = NaN_i
        StdI.mix = NaN_d
        StdI.eps = NaN_i
        StdI.eps_slater = NaN_i
        StdI.Iteration_max = NaN_i
        StdI.boxsub[:, :] = NaN_i
        StdI.Hsub = NaN_i
        StdI.Lsub = NaN_i
        StdI.Wsub = NaN_i

    # --- HWAVE-specific fields ---------------------------------------------
    elif StdI.solver == "HWAVE":
        StdI.NMPTrans = NaN_i
        StdI.RndSeed = NaN_i
        StdI.mix = NaN_d
        StdI.eps = NaN_i
        StdI.eps_slater = NaN_i
        StdI.Iteration_max = NaN_i
        StdI.boxsub[:, :] = NaN_i
        StdI.Hsub = NaN_i
        StdI.Lsub = NaN_i
        StdI.Wsub = NaN_i
        StdI.calcmode = "****"
        StdI.fileprefix = "****"
        StdI.export_all = NaN_i
        StdI.lattice_gp = NaN_i

    # --- Boost (always zero, not NaN) --------------------------------------
    StdI.lBoost = 0


# ===================================================================
#  _text2lower
# ===================================================================


def _text2lower(text: str) -> str:
    """Convert *text* to lower case.

    This is the Python translation of the C helper ``Text2Lower()``.

    Parameters
    ----------
    text : str
        Keyword or value string.

    Returns
    -------
    str
        The input converted to lower case.
    """
    return text.lower()


# ===================================================================
#  _trim_space_quote
# ===================================================================


def _trim_space_quote(text: str) -> str:
    """Remove whitespace, colons, semicolons, quotes and backslashes from *text*.

    This is the Python translation of the C helper ``TrimSpaceQuote()``.
    The original C code strips the characters: space, colon, semicolon,
    double-quote, backspace (``\\b``), backslash, vertical-tab (``\\v``),
    newline, and null.

    Parameters
    ----------
    text : str
        Raw keyword or value string from an input file.

    Returns
    -------
    str
        The cleaned string with the above characters removed.
    """
    remove_chars = set(" :;\"\\\b\v\n\0")
    return "".join(ch for ch in text if ch not in remove_chars)


# ===================================================================
#  _store_with_check_dup_s
# ===================================================================


def _store_with_check_dup_s(keyword: str, value: str, current: str) -> str:
    """Store a string value after checking for duplicate assignment.

    If *current* has already been assigned (i.e. it is not the sentinel
    ``"****"``), the program prints an error and exits.

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file.
    current : str
        The current stored value (``"****"`` means unset).

    Returns
    -------
    str
        The accepted value (equal to *value*).

    Raises
    ------
    SystemExit
        If *current* is not the sentinel, indicating a duplicate keyword.
    """
    if current != "****":
        print(f"ERROR !  Keyword {keyword} is duplicated ! ")
        exit_program(-1)
    return value


# ===================================================================
#  _store_with_check_dup_sl
# ===================================================================


def _store_with_check_dup_sl(
    keyword: str, value: str, current: str, maxlen: int = 256
) -> str:
    """Store a string value (forced lower-case) after checking for duplicates.

    Behaves like :func:`_store_with_check_dup_s` but additionally
    converts *value* to lower case and truncates it to *maxlen*
    characters.

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file.
    current : str
        The current stored value (``"****"`` means unset).
    maxlen : int, optional
        Maximum number of characters to keep (default 256).

    Returns
    -------
    str
        The accepted value, lower-cased and truncated.

    Raises
    ------
    SystemExit
        If *current* is not the sentinel, indicating a duplicate keyword.
    """
    if current != "****":
        print(f"ERROR !  Keyword {keyword} is duplicated ! ")
        exit_program(-1)
    return _text2lower(value[:maxlen])


# ===================================================================
#  _store_with_check_dup_i
# ===================================================================


def _store_with_check_dup_i(keyword: str, value: str, current: int) -> int:
    """Store an integer value after checking for duplicate assignment.

    If *current* differs from the integer sentinel (``2147483647``),
    the program prints an error and exits.

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file (will be converted to int).
    current : int
        The current stored value (``NaN_i`` means unset).

    Returns
    -------
    int
        The parsed integer value.

    Raises
    ------
    SystemExit
        If *current* is not the sentinel, indicating a duplicate keyword.
    """
    if current != NaN_i:
        print(f"ERROR !  Keyword {keyword} is duplicated ! ")
        exit_program(-1)
    # C sscanf("%d") truncates floats like "2.0" → 2
    return int(float(value))


# ===================================================================
#  _store_with_check_dup_d
# ===================================================================


def _store_with_check_dup_d(keyword: str, value: str, current: float) -> float:
    """Store a float value after checking for duplicate assignment.

    If *current* is **not** NaN the program prints an error and exits.

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file (will be converted to float).
    current : float
        The current stored value (``NaN`` means unset).

    Returns
    -------
    float
        The parsed float value.

    Raises
    ------
    SystemExit
        If *current* is not NaN, indicating a duplicate keyword.
    """
    if not math.isnan(current):
        print(f"ERROR !  Keyword {keyword} is duplicated ! ")
        exit_program(-1)
    return float(value)


# ===================================================================
#  _store_with_check_dup_c
# ===================================================================


def _store_with_check_dup_c(keyword: str, value: str, current: complex) -> complex:
    """Store a complex value after checking for duplicate assignment.

    The input string *value* may be in one of the following forms:

    * ``"real,imag"`` -- both parts specified
    * ``"real"``      -- imaginary part defaults to 0
    * ``",imag"``     -- real part defaults to 0

    If the real part of *current* is **not** NaN the program prints an
    error and exits.

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file.
    current : complex
        The current stored value (real-part ``NaN`` means unset).

    Returns
    -------
    complex
        The parsed complex value.

    Raises
    ------
    SystemExit
        If *current* is already set, indicating a duplicate keyword.
    """
    if not cmath.isnan(current):
        print(f"ERROR !  Keyword {keyword} is duplicated ! ")
        exit_program(-1)

    # Split on comma, mirroring the C strtok(",") logic
    if "," in value:
        parts = value.split(",", 1)
        real_str = parts[0].strip()
        imag_str = parts[1].strip() if len(parts) > 1 else ""
    else:
        real_str = value.strip()
        imag_str = ""

    # Parse real part
    if real_str == "":
        real_part = 0.0
    else:
        try:
            real_part = float(real_str)
        except ValueError:
            real_part = 0.0

    # Parse imaginary part
    if imag_str == "":
        imag_part = 0.0
    else:
        try:
            imag_part = float(imag_str)
        except ValueError:
            imag_part = 0.0

    return complex(real_part, imag_part)


# ===================================================================
#  _print_loc_spin
# ===================================================================


def _print_loc_spin(StdI: StdIntList) -> None:
    """Write ``locspn.def`` listing local-spin flags for every site.

    This is the Python translation of the C function ``PrintLocSpin()``.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``nsite`` : int -- total number of sites.
        - ``locspinflag`` : list of int -- per-site flag (0 = itinerant
          electron, nonzero = local spin with :math:`S` given by the value).
    """
    nlocspin = 0
    for isite in range(StdI.nsite):
        if StdI.locspinflag[isite] != 0:
            nlocspin += 1

    with open("locspn.def", "w") as fp:
        fp.write("================================ \n")
        fp.write(f"NlocalSpin {nlocspin:5d}  \n")
        fp.write("================================ \n")
        fp.write("========i_1LocSpn_0IteElc ====== \n")
        fp.write("================================ \n")
        for isite in range(StdI.nsite):
            fp.write(f"{isite:5d}  {StdI.locspinflag[isite]:5d}\n")

    print("    locspn.def is written.")


# ===================================================================
#  _print_trans
# ===================================================================


def _print_trans(StdI: StdIntList) -> None:
    """Write ``trans.def`` listing one-body transfer integrals.

    This is the Python translation of the C function ``PrintTrans()``.
    Duplicate index quadruples are merged (their amplitudes summed) and
    entries whose absolute value is below ``1e-6`` are suppressed.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read
        and (for merging) modified **in place**:

        - ``ntrans`` : int -- number of registered transfer terms.
        - ``transindx`` : 2-D array of int, shape ``(ntrans, 4)`` --
          site/spin indices ``(i, s_i, j, s_j)`` for each term.
        - ``trans`` : 1-D array of complex -- transfer amplitudes.
    """
    # --- merge duplicate index quadruples ---
    for jtrans in range(StdI.ntrans):
        for ktrans in range(jtrans + 1, StdI.ntrans):
            if (StdI.transindx[jtrans][0] == StdI.transindx[ktrans][0]
                    and StdI.transindx[jtrans][1] == StdI.transindx[ktrans][1]
                    and StdI.transindx[jtrans][2] == StdI.transindx[ktrans][2]
                    and StdI.transindx[jtrans][3] == StdI.transindx[ktrans][3]):
                StdI.trans[jtrans] = StdI.trans[jtrans] + StdI.trans[ktrans]
                StdI.trans[ktrans] = 0.0

    # --- count non-negligible entries ---
    ntrans0 = 0
    for ktrans in range(StdI.ntrans):
        if abs(StdI.trans[ktrans]) > 0.000001:
            ntrans0 += 1

    # --- write file ---
    with open("trans.def", "w") as fp:
        fp.write("======================== \n")
        fp.write(f"NTransfer {ntrans0:7d}  \n")
        fp.write("======================== \n")
        fp.write("========i_j_s_tijs====== \n")
        fp.write("======================== \n")
        for ktrans in range(StdI.ntrans):
            if abs(StdI.trans[ktrans]) > 0.000001:
                fp.write(
                    f"{StdI.transindx[ktrans][0]:5d} "
                    f"{StdI.transindx[ktrans][1]:5d} "
                    f"{StdI.transindx[ktrans][2]:5d} "
                    f"{StdI.transindx[ktrans][3]:5d} "
                    f"{StdI.trans[ktrans].real:25.15f} "
                    f"{StdI.trans[ktrans].imag:25.15f}\n"
                )

    print("      trans.def is written.")


# ===================================================================
#  _print_namelist
# ===================================================================


def _print_namelist(StdI: StdIntList) -> None:
    """Write ``namelist.def`` that lists all definition files for the solver.

    This is the Python translation of the C function ``PrintNamelist()``.
    The content depends on which solver is active (``StdI.solver``).

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Relevant fields include:

        - ``solver`` : str -- ``"HPhi"``, ``"mVMC"``, ``"UHF"``, or
          ``"HWAVE"``.
        - ``LCintra``, ``LCinter``, ``LHund``, ``LEx``, ``LPairLift``,
          ``LPairHopp``, ``Lintr`` : int -- flags indicating whether the
          corresponding interaction file should be listed.
        - ``ioutputmode`` : int -- if nonzero, Green-function files are
          listed.
        - ``SpectrumBody`` : int -- (HPhi) 1 for single excitation, else
          pair excitation.
        - ``method`` : str -- (HPhi) calculation method name.
        - ``PumpBody`` : int -- (HPhi) 1 or 2 for time-evolution pump type.
        - ``CDataFileHead`` : str -- (HPhi) base name for eigenvector files.
        - ``lBoost`` : int -- (HPhi) 1 if Boost is used.
        - ``lGC`` : int -- (mVMC) 1 if grand-canonical.
        - ``Sz2`` : int -- twice the total :math:`S_z`.
        - ``NaN_i`` : int -- integer sentinel for "not set".
    """
    with open("namelist.def", "w") as fp:
        fp.write("         ModPara  modpara.def\n")
        fp.write("         LocSpin  locspn.def\n")
        fp.write("           Trans  trans.def\n")
        if StdI.LCintra == 1:
            fp.write("    CoulombIntra  coulombintra.def\n")
        if StdI.LCinter == 1:
            fp.write("    CoulombInter  coulombinter.def\n")
        if StdI.LHund == 1:
            fp.write("            Hund  hund.def\n")
        if StdI.LEx == 1:
            fp.write("        Exchange  exchange.def\n")
        if StdI.LPairLift == 1:
            fp.write("        PairLift  pairlift.def\n")
        if StdI.LPairHopp == 1:
            fp.write("         PairHop  pairhopp.def\n")
        if StdI.Lintr == 1:
            fp.write("        InterAll  interall.def\n")
        if StdI.ioutputmode != 0:
            fp.write("        OneBodyG  greenone.def\n")
            if StdI.solver in ("HPhi", "mVMC"):
                fp.write("        TwoBodyG  greentwo.def\n")

        # --- HPhi-specific entries ---
        if StdI.solver == "HPhi":
            fp.write("         CalcMod  calcmod.def\n")
            if StdI.SpectrumBody == 1:
                fp.write("SingleExcitation  single.def\n")
            else:
                fp.write("  PairExcitation  pair.def\n")
            if StdI.method == "timeevolution":
                if StdI.PumpBody == 1:
                    fp.write("       TEOneBody  teone.def\n")
                elif StdI.PumpBody == 2:
                    fp.write("       TETwoBody  tetwo.def\n")
            fp.write(f"     SpectrumVec  {StdI.CDataFileHead}_eigenvec_0\n")
            if StdI.lBoost == 1:
                fp.write("           Boost  boost.def\n")

        # --- mVMC-specific entries ---
        elif StdI.solver == "mVMC":
            fp.write("      Gutzwiller  gutzwilleridx.def\n")
            fp.write("         Jastrow  jastrowidx.def\n")
            fp.write("         Orbital  orbitalidx.def\n")
            if StdI.lGC == 1 or (StdI.Sz2 != 0 and StdI.Sz2 != StdI.NaN_i):
                fp.write(" OrbitalParallel  orbitalidxpara.def\n")
                fp.write("# OrbitalGeneral  orbitalidxgen.def\n")
            fp.write("        TransSym  qptransidx.def\n")

    print("    namelist.def is written.")


# ===================================================================
#  _print_mod_para
# ===================================================================


def _print_mod_para(StdI: StdIntList) -> None:
    """Write ``modpara.def`` containing model / calculation parameters.

    This is the Python translation of the C function ``PrintModPara()``.
    The file layout depends on the active solver (``StdI.solver``):

    - ``"HPhi"``  -- HPhi calculation parameters (Lanczos, spectrum, etc.).
    - ``"mVMC"``  -- VMC calculation parameters.
    - ``"UHF"``   -- UHF (unrestricted Hartree-Fock) parameters.
    - ``"HWAVE"`` -- H-wave parameters (same body as UHF, different header).

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  A large number of solver-specific
        fields are read; see the C source (``StdFace_main.c``, lines
        1441--1543) for the complete list.
    """
    with open("modpara.def", "w") as fp:
        fp.write("--------------------\n")
        fp.write("Model_Parameters   0\n")
        fp.write("--------------------\n")

        # ---------------------------------------------------------------
        #  HPhi
        # ---------------------------------------------------------------
        if StdI.solver == "HPhi":
            fp.write("HPhi_Cal_Parameters\n")
            fp.write("--------------------\n")
            fp.write(f"CDataFileHead  {StdI.CDataFileHead}\n")
            fp.write("CParaFileHead  zqp\n")
            fp.write("--------------------\n")
            fp.write(f"Nsite          {StdI.nsite:<5d}\n")
            if StdI.Sz2 != StdI.NaN_i:
                fp.write(f"2Sz            {StdI.Sz2:<5d}\n")
            if StdI.ncond != StdI.NaN_i:
                fp.write(f"Ncond          {StdI.ncond:<5d}\n")
            fp.write(f"Lanczos_max    {StdI.Lanczos_max:<5d}\n")
            fp.write(f"initial_iv     {StdI.initial_iv:<5d}\n")
            if StdI.nvec != StdI.NaN_i:
                fp.write(f"nvec           {StdI.nvec:<5d}\n")
            fp.write(f"exct           {StdI.exct:<5d}\n")
            fp.write(f"LanczosEps     {StdI.LanczosEps:<5d}\n")
            fp.write(f"LanczosTarget  {StdI.LanczosTarget:<5d}\n")
            fp.write(f"LargeValue     {StdI.LargeValue:<25.15e}\n")
            fp.write(f"NumAve         {StdI.NumAve:<5d}\n")
            fp.write(f"ExpecInterval  {StdI.ExpecInterval:<5d}\n")
            fp.write(f"NOmega         {StdI.Nomega:<5d}\n")
            fp.write(
                f"OmegaMax       {StdI.OmegaMax:<25.15e} {StdI.OmegaIm:<25.15e}\n"
            )
            fp.write(
                f"OmegaMin       {StdI.OmegaMin:<25.15e} {StdI.OmegaIm:<25.15e}\n"
            )
            fp.write(
                f"OmegaOrg       {StdI.OmegaOrg:<25.15e} {0.0:<25.15e}\n"
            )
            fp.write(f"PreCG          {0:<5d}\n")
            if StdI.method == "timeevolution":
                fp.write(f"ExpandCoef     {StdI.ExpandCoef:<5d}\n")

        # ---------------------------------------------------------------
        #  mVMC
        # ---------------------------------------------------------------
        elif StdI.solver == "mVMC":
            fp.write("VMC_Cal_Parameters\n")
            fp.write("--------------------\n")
            fp.write(f"CDataFileHead  {StdI.CDataFileHead}\n")
            fp.write(f"CParaFileHead  {StdI.CParaFileHead}\n")
            fp.write("--------------------\n")
            fp.write(f"NVMCCalMode    {StdI.NVMCCalMode}\n")
            fp.write(f"NLanczosMode   {StdI.NLanczosMode}\n")
            fp.write("--------------------\n")
            fp.write(f"NDataIdxStart  {StdI.NDataIdxStart}\n")
            fp.write(f"NDataQtySmp    {StdI.NDataQtySmp}\n")
            fp.write("--------------------\n")
            fp.write(f"Nsite          {StdI.nsite}\n")
            fp.write(f"Ncond          {StdI.ncond:<5d}\n")
            if StdI.Sz2 != StdI.NaN_i:
                fp.write(f"2Sz            {StdI.Sz2}\n")
            if StdI.NSPGaussLeg != StdI.NaN_i:
                fp.write(f"NSPGaussLeg    {StdI.NSPGaussLeg}\n")
            if StdI.NSPStot != StdI.NaN_i:
                fp.write(f"NSPStot        {StdI.NSPStot}\n")
            fp.write(f"NMPTrans       {StdI.NMPTrans}\n")
            fp.write(f"NSROptItrStep  {StdI.NSROptItrStep}\n")
            fp.write(f"NSROptItrSmp   {StdI.NSROptItrSmp}\n")
            fp.write(f"DSROptRedCut   {StdI.DSROptRedCut:.10f}\n")
            fp.write(f"DSROptStaDel   {StdI.DSROptStaDel:.10f}\n")
            fp.write(f"DSROptStepDt   {StdI.DSROptStepDt:.10f}\n")
            fp.write(f"NVMCWarmUp     {StdI.NVMCWarmUp}\n")
            fp.write(f"NVMCInterval   {StdI.NVMCInterval}\n")
            fp.write(f"NVMCSample     {StdI.NVMCSample}\n")
            fp.write(f"NExUpdatePath  {StdI.NExUpdatePath}\n")
            fp.write(f"RndSeed        {StdI.RndSeed}\n")
            fp.write(f"NSplitSize     {StdI.NSplitSize}\n")
            fp.write(f"NStore         {StdI.NStore}\n")
            fp.write(f"NSRCG          {StdI.NSRCG}\n")

        # ---------------------------------------------------------------
        #  UHF
        # ---------------------------------------------------------------
        elif StdI.solver == "UHF":
            fp.write("UHF_Cal_Parameters\n")
            fp.write("--------------------\n")
            fp.write(f"CDataFileHead  {StdI.CDataFileHead}\n")
            fp.write("CParaFileHead  zqp\n")
            fp.write("--------------------\n")
            fp.write(f"Nsite          {StdI.nsite}\n")
            if StdI.Sz2 != StdI.NaN_i:
                fp.write(f"2Sz            {StdI.Sz2:<5d}\n")
            fp.write(f"Ncond          {StdI.ncond:<5d}\n")
            fp.write(f"IterationMax   {StdI.Iteration_max}\n")
            fp.write(f"EPS            {StdI.eps}\n")
            fp.write(f"Mix            {StdI.mix:.10f}\n")
            fp.write(f"RndSeed        {StdI.RndSeed}\n")
            fp.write(f"EpsSlater      {StdI.eps_slater}\n")
            fp.write(f"NMPTrans       {StdI.NMPTrans}\n")

        # ---------------------------------------------------------------
        #  HWAVE
        # ---------------------------------------------------------------
        elif StdI.solver == "HWAVE":
            fp.write("HWAVE_Cal_Parameters\n")
            fp.write("--------------------\n")
            fp.write(f"CDataFileHead  {StdI.CDataFileHead}\n")
            fp.write("CParaFileHead  zqp\n")
            fp.write("--------------------\n")
            fp.write(f"Nsite          {StdI.nsite}\n")
            if StdI.Sz2 != StdI.NaN_i:
                fp.write(f"2Sz            {StdI.Sz2:<5d}\n")
            fp.write(f"Ncond          {StdI.ncond:<5d}\n")
            fp.write(f"IterationMax   {StdI.Iteration_max}\n")
            fp.write(f"EPS            {StdI.eps}\n")
            fp.write(f"Mix            {StdI.mix:.10f}\n")
            fp.write(f"RndSeed        {StdI.RndSeed}\n")
            fp.write(f"EpsSlater      {StdI.eps_slater}\n")
            fp.write(f"NMPTrans       {StdI.NMPTrans}\n")

    print("     modpara.def is written.")


# --- END OF PART 2 ---


# ===================================================================
#  _print_1_green
# ===================================================================


def _print_1_green(StdI: StdIntList) -> None:
    """Write ``greenone.def`` listing one-body Green-function indices.

    This is the Python translation of the C function ``Print1Green()``.
    A two-pass approach is used: the first pass counts the number of
    index tuples and the second pass collects them.  The output depends
    on the output mode stored in ``StdI.ioutputmode``:

    - **mode 1** (correlation): only same-spin pairs (``ispin == jspin``)
      are kept.  For Kondo models the unit-cell size is doubled
      (``NsiteUC * 2``).
    - **mode 2** (raw / full): all site-spin combinations are emitted,
      subject to the constraint that pairs of local-spin sites at
      different positions are skipped.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``ioutputmode`` : int -- 0 (none), 1 (correlation), or 2 (raw).
        - ``model`` : str -- model name (``"kondo"`` triggers doubled
          unit cell).
        - ``NsiteUC`` : int -- number of sites in the unit cell.
        - ``nsite`` : int -- total number of sites.
        - ``locspinflag`` : list of int -- per-site local-spin flag.
    """
    if StdI.ioutputmode == 0:
        return

    greenindx: list[tuple[int, int, int, int]] = []

    if StdI.model == "kondo":
        xkondo = 2
    else:
        xkondo = 1

    if StdI.ioutputmode == 1:
        for isite in range(StdI.NsiteUC * xkondo):
            if isite >= StdI.NsiteUC:
                isite2 = isite - StdI.NsiteUC + StdI.nsite // 2
            else:
                isite2 = isite

            if StdI.locspinflag[isite2] == 0:
                SiMax = 1
            else:
                SiMax = StdI.locspinflag[isite2]

            for ispin in range(SiMax + 1):
                for jsite in range(StdI.nsite):
                    if StdI.locspinflag[jsite] == 0:
                        SjMax = 1
                    else:
                        SjMax = StdI.locspinflag[jsite]

                    for jspin in range(SjMax + 1):
                        if (isite2 != jsite
                                and StdI.locspinflag[isite2] != 0
                                and StdI.locspinflag[jsite] != 0):
                            continue

                        if ispin == jspin:
                            greenindx.append((isite2, ispin, jsite, jspin))
    else:
        # ioutputmode == 2
        for isite in range(StdI.nsite):
            if StdI.locspinflag[isite] == 0:
                SiMax = 1
            else:
                SiMax = StdI.locspinflag[isite]

            for ispin in range(SiMax + 1):
                for jsite in range(StdI.nsite):
                    if StdI.locspinflag[jsite] == 0:
                        SjMax = 1
                    else:
                        SjMax = StdI.locspinflag[jsite]

                    for jspin in range(SjMax + 1):
                        if (isite != jsite
                                and StdI.locspinflag[isite] != 0
                                and StdI.locspinflag[jsite] != 0):
                            continue

                        greenindx.append((isite, ispin, jsite, jspin))

    ngreen = len(greenindx)

    with open("greenone.def", "w") as fp:
        fp.write("===============================\n")
        fp.write(f"NCisAjs {ngreen:10d}\n")
        fp.write("===============================\n")
        fp.write("======== Green functions ======\n")
        fp.write("===============================\n")
        for ig in range(ngreen):
            fp.write(
                f"{greenindx[ig][0]:5d} {greenindx[ig][1]:5d} "
                f"{greenindx[ig][2]:5d} {greenindx[ig][3]:5d}\n"
            )

    print("    greenone.def is written.")


# ===================================================================
#  _print_2_green
# ===================================================================


def _print_2_green(StdI: StdIntList) -> None:
    """Write ``greentwo.def`` listing two-body Green-function indices.

    This is the Python translation of the C function ``Print2Green()``.
    Like :func:`_print_1_green`, a two-pass approach is replaced by
    collecting tuples in a list.  Each entry is an 8-element tuple of
    ``(site1, spin1, site2, spin2, site3, spin3, site4, spin4)``.

    The behaviour depends on ``StdI.ioutputmode``:

    - **mode 1** (correlation): for each pair ``(site1, site3)`` only
      spin combinations satisfying
      ``spin1 - spin2 + spin3 - spin4 == 0`` are kept.  For mVMC, when
      ``spin1 != spin2`` or ``spin3 != spin4`` a special index order is
      used; otherwise the standard HPhi order is used.
    - **mode 2** (raw / full): all four-site combinations are emitted
      with local-spin pair filtering on ``(site1, site2)`` and
      ``(site3, site4)``.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``ioutputmode`` : int
        - ``model`` : str
        - ``solver`` : str
        - ``NsiteUC`` : int
        - ``nsite`` : int
        - ``locspinflag`` : list of int
    """
    ngreen = 0
    greenindx: list[tuple[int, int, int, int, int, int, int, int]] = []

    if StdI.ioutputmode == 1:
        if StdI.model == "kondo":
            xkondo = 2
        else:
            xkondo = 1

        for site1 in range(StdI.NsiteUC * xkondo):
            if site1 >= StdI.NsiteUC:
                site1k = site1 - StdI.NsiteUC + StdI.nsite // 2
            else:
                site1k = site1

            if StdI.locspinflag[site1k] == 0:
                S1Max = 1
            else:
                S1Max = StdI.locspinflag[site1k]

            for spin1 in range(S1Max + 1):
                for spin2 in range(S1Max + 1):
                    for site3 in range(StdI.nsite):
                        if StdI.locspinflag[site3] == 0:
                            S3Max = 1
                        else:
                            S3Max = StdI.locspinflag[site3]

                        for spin3 in range(S3Max + 1):
                            for spin4 in range(S3Max + 1):
                                if spin1 - spin2 + spin3 - spin4 == 0:
                                    if StdI.solver == "mVMC" and (
                                        spin1 != spin2 or spin3 != spin4
                                    ):
                                        greenindx.append((
                                            site1k, spin1,
                                            site3, spin4,
                                            site3, spin3,
                                            site1k, spin2,
                                        ))
                                    else:
                                        greenindx.append((
                                            site1k, spin1,
                                            site1k, spin2,
                                            site3, spin3,
                                            site3, spin4,
                                        ))

    elif StdI.ioutputmode == 2:
        for site1 in range(StdI.nsite):
            if StdI.locspinflag[site1] == 0:
                S1Max = 1
            else:
                S1Max = StdI.locspinflag[site1]

            for spin1 in range(S1Max + 1):
                for site2 in range(StdI.nsite):
                    if (StdI.locspinflag[site1] != 0
                            and StdI.locspinflag[site2] != 0
                            and site1 != site2):
                        continue

                    if StdI.locspinflag[site2] == 0:
                        S2Max = 1
                    else:
                        S2Max = StdI.locspinflag[site2]

                    for spin2 in range(S2Max + 1):
                        for site3 in range(StdI.nsite):
                            if StdI.locspinflag[site3] == 0:
                                S3Max = 1
                            else:
                                S3Max = StdI.locspinflag[site3]

                            for spin3 in range(S3Max + 1):
                                for site4 in range(StdI.nsite):
                                    if (StdI.locspinflag[site3] != 0
                                            and StdI.locspinflag[site4] != 0
                                            and site3 != site4):
                                        continue

                                    if StdI.locspinflag[site4] == 0:
                                        S4Max = 1
                                    else:
                                        S4Max = StdI.locspinflag[site4]

                                    for spin4 in range(S4Max + 1):
                                        greenindx.append((
                                            site1, spin1,
                                            site2, spin2,
                                            site3, spin3,
                                            site4, spin4,
                                        ))

    if StdI.ioutputmode != 0:
        ngreen = len(greenindx)
        with open("greentwo.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"NCisAjsCktAltDC {ngreen:10d}\n")
            fp.write("=============================================\n")
            fp.write("======== Green functions for Sq AND Nq ======\n")
            fp.write("=============================================\n")
            for ig in range(ngreen):
                fp.write(
                    f"{greenindx[ig][0]:5d} {greenindx[ig][1]:5d} "
                    f"{greenindx[ig][2]:5d} {greenindx[ig][3]:5d} "
                    f"{greenindx[ig][4]:5d} {greenindx[ig][5]:5d} "
                    f"{greenindx[ig][6]:5d} {greenindx[ig][7]:5d}\n"
                )

        print("    greentwo.def is written.")


# ===================================================================
#  _unsupported_system
# ===================================================================


def _unsupported_system(model: str, lattice: str) -> None:
    """Print an error message and abort for an unsupported model/lattice pair.

    This is the Python translation of the C function
    ``UnsupportedSystem()``.

    Parameters
    ----------
    model : str
        The model name specified by the user.
    lattice : str
        The lattice name specified by the user.

    Raises
    ------
    SystemExit
        Always raised after printing the error message.
    """
    print("\nSorry, specified combination, ")
    print(f"    MODEL : {model}  ")
    print(f"  LATTICE : {lattice}, ")
    print("is unsupported in the STANDARD MODE...")
    print("Please use the EXPART MODE, or write a NEW FUNCTION and post us.")
    exit_program(-1)


# ===================================================================
#  _check_output_mode
# ===================================================================


def _check_output_mode(StdI: StdIntList) -> None:
    """Verify and set the integer output-mode flag from the string keyword.

    This is the Python translation of the C function
    ``CheckOutputMode()``.  The mapping is:

    - ``"non"`` / ``"none"`` / ``"off"`` -> ``ioutputmode = 0``
    - ``"cor"`` / ``"corr"`` / ``"correlation"`` -> ``ioutputmode = 1``
    - ``"****"`` (default sentinel) -> ``ioutputmode = 1``
    - ``"raw"`` / ``"all"`` / ``"full"`` -> ``ioutputmode = 2``
    - anything else -> error and exit

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  ``StdI.outputmode`` is read and
        ``StdI.ioutputmode`` is set **in place**.

    Raises
    ------
    SystemExit
        If ``StdI.outputmode`` does not match any recognised keyword.
    """
    if StdI.outputmode in ("non", "none", "off"):
        StdI.ioutputmode = 0
        print(f"      ioutputmode = {StdI.ioutputmode:<10d}")
    elif StdI.outputmode in ("cor", "corr", "correlation"):
        StdI.ioutputmode = 1
        print(f"      ioutputmode = {StdI.ioutputmode:<10d}")
    elif StdI.outputmode == "****":
        StdI.ioutputmode = 1
        print(
            f"      ioutputmode = {StdI.ioutputmode:<10d}"
            "  ######  DEFAULT VALUE IS USED  ######"
        )
    elif StdI.outputmode in ("raw", "all", "full"):
        StdI.ioutputmode = 2
        print(f"      ioutputmode = {StdI.ioutputmode:<10d}")
    else:
        print(f"\n ERROR ! Unsupported OutPutMode : {StdI.outputmode}")
        exit_program(-1)


# ===================================================================
#  _check_mod_para
# ===================================================================


def _check_mod_para(StdI: StdIntList) -> None:
    """Validate and set default values for solver-specific model parameters.

    This is the Python translation of the C function ``CheckModPara()``.
    Depending on ``StdI.solver``, different parameter groups are
    validated:

    - **HPhi**: Lanczos parameters, spectrum frequency grid, large-value
      cutoff.
    - **mVMC**: VMC sampling parameters, exchange-update path count,
      SR-optimisation settings.
    - **UHF**: random seed, iteration limit, mixing parameters.
    - **HWAVE**: same as UHF.

    In all cases the conserved quantities (``ncond``, ``2Sz``) are
    checked for consistency with the model type and
    grand-canonical flag.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Many fields are read and
        modified **in place** via the helper functions
        :func:`print_val_i`, :func:`print_val_d`, :func:`not_used_i`,
        and :func:`required_val_i`.
    """
    # ------------------------------------------------------------------
    #  HPhi-specific defaults
    # ------------------------------------------------------------------
    if StdI.solver == "HPhi":
        StdI.Lanczos_max = print_val_i("Lanczos_max", StdI.Lanczos_max, 2000)
        StdI.initial_iv = print_val_i("initial_iv", StdI.initial_iv, -1)
        # nvec is not given a default here (commented out in C)
        StdI.exct = print_val_i("exct", StdI.exct, 1)
        StdI.LanczosEps = print_val_i("LanczosEps", StdI.LanczosEps, 14)
        StdI.LanczosTarget = print_val_i("LanczosTarget", StdI.LanczosTarget, 2)
        if StdI.LanczosTarget < StdI.exct:
            StdI.LanczosTarget = StdI.exct
        StdI.NumAve = print_val_i("NumAve", StdI.NumAve, 5)
        StdI.ExpecInterval = print_val_i("ExpecInterval", StdI.ExpecInterval, 20)
        StdI.Nomega = print_val_i("NOmega", StdI.Nomega, 200)
        StdI.OmegaMax = print_val_d(
            "OmegaMax", StdI.OmegaMax, StdI.LargeValue * StdI.nsite
        )
        StdI.OmegaMin = print_val_d(
            "OmegaMin", StdI.OmegaMin, -StdI.LargeValue * StdI.nsite
        )
        StdI.OmegaOrg = print_val_d("OmegaOrg", StdI.OmegaOrg, 0.0)
        StdI.OmegaIm = print_val_d(
            "OmegaIm", StdI.OmegaIm, 0.01 * int(StdI.LargeValue)
        )

    # ------------------------------------------------------------------
    #  mVMC-specific defaults
    # ------------------------------------------------------------------
    elif StdI.solver == "mVMC":
        if StdI.CParaFileHead == "****":
            StdI.CParaFileHead = "zqp"
            print(
                f"    CParaFileHead = {StdI.CParaFileHead:<12s}"
                "######  DEFAULT VALUE IS USED  ######"
            )
        else:
            print(f"    CParaFileHead = {StdI.CParaFileHead}")

        StdI.NVMCCalMode = print_val_i("NVMCCalMode", StdI.NVMCCalMode, 0)
        StdI.NLanczosMode = print_val_i("NLanczosMode", StdI.NLanczosMode, 0)
        StdI.NDataIdxStart = print_val_i("NDataIdxStart", StdI.NDataIdxStart, 1)

        if StdI.NVMCCalMode == 0:
            not_used_i("NDataQtySmp", StdI.NDataQtySmp)
        StdI.NDataQtySmp = print_val_i("NDataQtySmp", StdI.NDataQtySmp, 1)

        if StdI.lGC == 0 and (StdI.Sz2 == 0 or StdI.Sz2 == StdI.NaN_i):
            StdI.NSPGaussLeg = print_val_i("NSPGaussLeg", StdI.NSPGaussLeg, 8)
            StdI.NSPStot = print_val_i("NSPStot", StdI.NSPStot, 0)
        else:
            not_used_i("NSPGaussLeg", StdI.NSPGaussLeg)
            not_used_i("NSPStot", StdI.NSPStot)

        StdI.NMPTrans = print_val_i("NMPTrans", StdI.NMPTrans, -1)

        StdI.NSROptItrStep = print_val_i("NSROptItrStep", StdI.NSROptItrStep, 1000)

        if StdI.NVMCCalMode == 1:
            not_used_i("NSROptItrSmp", StdI.NSROptItrSmp)
        StdI.NSROptItrSmp = print_val_i(
            "NSROptItrSmp", StdI.NSROptItrSmp, StdI.NSROptItrStep // 10
        )

        StdI.NVMCWarmUp = print_val_i("NVMCWarmUp", StdI.NVMCWarmUp, 10)
        StdI.NVMCInterval = print_val_i("NVMCInterval", StdI.NVMCInterval, 1)
        StdI.NVMCSample = print_val_i("NVMCSample", StdI.NVMCSample, 1000)

        if StdI.model == "hubbard":
            StdI.NExUpdatePath = 0
        elif StdI.model == "spin":
            StdI.NExUpdatePath = 2
        elif StdI.model == "kondo":
            if StdI.lGC == 0:
                StdI.NExUpdatePath = 1
            else:
                StdI.NExUpdatePath = 3
        print(f"  {'NExUpdatePath':>15s} = {StdI.NExUpdatePath:<10d}")

        StdI.RndSeed = print_val_i("RndSeed", StdI.RndSeed, 123456789)
        StdI.NSplitSize = print_val_i("NSplitSize", StdI.NSplitSize, 1)
        StdI.NStore = print_val_i("NStore", StdI.NStore, 1)
        StdI.NSRCG = print_val_i("NSRCG", StdI.NSRCG, 0)

        StdI.DSROptRedCut = print_val_d("DSROptRedCut", StdI.DSROptRedCut, 0.001)
        StdI.DSROptStaDel = print_val_d("DSROptStaDel", StdI.DSROptStaDel, 0.02)
        StdI.DSROptStepDt = print_val_d("DSROptStepDt", StdI.DSROptStepDt, 0.02)

    # ------------------------------------------------------------------
    #  UHF-specific defaults
    # ------------------------------------------------------------------
    elif StdI.solver == "UHF":
        StdI.RndSeed = print_val_i("RndSeed", StdI.RndSeed, 123456789)
        StdI.Iteration_max = print_val_i("Iteration_max", StdI.Iteration_max, 1000)
        StdI.mix = print_val_d("Mix", StdI.mix, 0.5)
        StdI.eps = print_val_i("eps", StdI.eps, 8)
        StdI.eps_slater = print_val_i("EpsSlater", StdI.eps_slater, 6)
        StdI.NMPTrans = print_val_i("NMPTrans", StdI.NMPTrans, 0)

    # ------------------------------------------------------------------
    #  HWAVE-specific defaults
    # ------------------------------------------------------------------
    elif StdI.solver == "HWAVE":
        StdI.RndSeed = print_val_i("RndSeed", StdI.RndSeed, 123456789)
        StdI.Iteration_max = print_val_i("Iteration_max", StdI.Iteration_max, 1000)
        StdI.mix = print_val_d("Mix", StdI.mix, 0.5)
        StdI.eps = print_val_i("eps", StdI.eps, 8)
        StdI.eps_slater = print_val_i("EpsSlater", StdI.eps_slater, 6)
        StdI.NMPTrans = print_val_i("NMPTrans", StdI.NMPTrans, 0)

    # ------------------------------------------------------------------
    #  Conserved quantities: ncond and 2Sz
    # ------------------------------------------------------------------
    if StdI.model == "hubbard":
        if StdI.solver == "HPhi":
            if StdI.lGC == 0:
                required_val_i("nelec", StdI.ncond)
            else:
                not_used_i("nelec", StdI.ncond)
                not_used_i("2Sz", StdI.Sz2)
        else:
            required_val_i("ncond", StdI.ncond)
            if StdI.lGC == 0:
                StdI.Sz2 = print_val_i("2Sz", StdI.Sz2, 0)
            else:
                not_used_i("2Sz", StdI.Sz2)

    elif StdI.model == "spin":
        not_used_i("ncond", StdI.ncond)
        if StdI.solver == "mVMC":
            StdI.ncond = 0
        if StdI.lGC == 0:
            required_val_i("2Sz", StdI.Sz2)
        else:
            not_used_i("2Sz", StdI.Sz2)

    elif StdI.model == "kondo":
        if StdI.solver == "HPhi":
            if StdI.lGC == 0:
                required_val_i("ncond", StdI.ncond)
            else:
                not_used_i("nelec", StdI.ncond)
                not_used_i("2Sz", StdI.Sz2)
        else:
            required_val_i("ncond", StdI.ncond)
            if StdI.lGC == 0:
                StdI.Sz2 = print_val_i("2Sz", StdI.Sz2, 0)
            else:
                not_used_i("2Sz", StdI.Sz2)


# ===================================================================
#  _large_value  (HPhi only)
# ===================================================================


def _large_value(StdI: StdIntList) -> None:
    """Compute and set the ``LargeValue`` parameter for TPQ calculations.

    ``LargeValue`` is the sum of the absolute values of all one-body
    (transfer) and two-body (interaction / exchange / Hund / pair-lift /
    Coulomb-intra / Coulomb-inter) terms, divided by the number of
    sites.  The result is stored in ``StdI.LargeValue`` via
    :func:`print_val_d`.

    This is the Python translation of the C function
    ``StdFace_LargeValue()`` (lines 52--80 of ``StdFace_main.c``).

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``ntrans``, ``trans`` -- one-body transfer terms.
        - ``nintr``, ``intr`` -- general two-body interaction terms.
        - ``NCintra``, ``Cintra`` -- intra-site Coulomb terms.
        - ``NCinter``, ``Cinter`` -- inter-site Coulomb terms.
        - ``NEx``, ``Ex`` -- exchange terms.
        - ``NPairLift``, ``PairLift`` -- pair-lift terms.
        - ``NHund``, ``Hund`` -- Hund coupling terms.
        - ``nsite`` -- total number of sites.
        - ``LargeValue`` -- set **in place** to the computed value
          (unless the user already specified it).
    """
    large_value0 = 0.0

    for ktrans in range(StdI.ntrans):
        large_value0 += abs(StdI.trans[ktrans])

    for kintr in range(StdI.nintr):
        large_value0 += abs(StdI.intr[kintr])

    for kintr in range(StdI.NCintra):
        large_value0 += abs(StdI.Cintra[kintr])

    for kintr in range(StdI.NCinter):
        large_value0 += abs(StdI.Cinter[kintr])

    for kintr in range(StdI.NEx):
        large_value0 += 2.0 * abs(StdI.Ex[kintr])

    for kintr in range(StdI.NPairLift):
        large_value0 += 2.0 * abs(StdI.PairLift[kintr])

    for kintr in range(StdI.NHund):
        large_value0 += 2.0 * abs(StdI.Hund[kintr])

    large_value0 /= float(StdI.nsite)

    StdI.LargeValue = print_val_d("LargeValue", StdI.LargeValue, large_value0)


# ===================================================================
#  _print_calc_mod  (HPhi only)
# ===================================================================


def _print_calc_mod(StdI: StdIntList) -> None:
    """Write ``calcmod.def`` containing calculation-mode integers for HPhi.

    Maps the string-valued parameters ``method``, ``model``, ``Restart``,
    ``InitialVecType``, ``EigenVecIO``, ``HamIO``, ``CalcSpec``, and
    ``OutputExVec`` to the integer codes that HPhi expects.  Also
    validates ``NGPU`` and ``Scalapack``.

    This is the Python translation of the C function
    ``PrintCalcMod()`` (lines 85--306 of ``StdFace_main.c``).

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  String fields such as
        ``method``, ``model``, ``Restart``, etc. are read (and may be
        overwritten with default values).  ``NGPU`` and ``Scalapack``
        are validated if set.
    """
    print("\n  @ CalcMod\n")

    # ------------------------------------------------------------------
    #  Method
    # ------------------------------------------------------------------
    iCalcEigenvec = 0

    if StdI.method == "****":
        print("ERROR ! Method is NOT specified !")
        exit_program(-1)
    elif StdI.method == "lanczos":
        iCalcType = 0
    elif StdI.method == "lanczosenergy":
        iCalcType = 0
        iCalcEigenvec = 1
    elif StdI.method == "tpq":
        iCalcType = 1
    elif StdI.method == "fulldiag":
        iCalcType = 2
    elif StdI.method == "cg":
        iCalcType = 3
    elif StdI.method == "timeevolution":
        iCalcType = 4
    elif StdI.method == "ctpq":
        iCalcType = 5
    else:
        print(f"\n ERROR ! Unsupported Solver : {StdI.method}")
        exit_program(-1)

    if iCalcType != 4:
        StdI.PumpBody = 0

    # ------------------------------------------------------------------
    #  Model
    # ------------------------------------------------------------------
    if StdI.model == "hubbard":
        iCalcModel = 0 if StdI.lGC == 0 else 3
    elif StdI.model == "spin":
        iCalcModel = 1 if StdI.lGC == 0 else 4
    elif StdI.model == "kondo":
        iCalcModel = 2 if StdI.lGC == 0 else 5

    # ------------------------------------------------------------------
    #  Restart
    # ------------------------------------------------------------------
    if StdI.Restart == "****":
        StdI.Restart = "none"
        print("          Restart = none        ######  DEFAULT VALUE IS USED  ######")
        iRestart = 0
    else:
        print(f"          Restart = {StdI.Restart}")
        if StdI.Restart == "none":
            iRestart = 0
        elif StdI.Restart in ("restart_out", "save"):
            iRestart = 1
        elif StdI.Restart in ("restartsave", "restart"):
            iRestart = 2
        elif StdI.Restart == "restart_in":
            iRestart = 3
        else:
            print(f"\n ERROR ! Restart Mode : {StdI.Restart}")
            exit_program(-1)

    # ------------------------------------------------------------------
    #  InitialVecType
    # ------------------------------------------------------------------
    if StdI.InitialVecType == "****":
        StdI.InitialVecType = "c"
        print("   InitialVecType = c           ######  DEFAULT VALUE IS USED  ######")
        if StdI.method in ("tpq", "ctpq"):
            iInitialVecType = -1
        else:
            iInitialVecType = 0
    else:
        print(f"   InitialVecType = {StdI.InitialVecType}")
        if StdI.InitialVecType == "c":
            iInitialVecType = 0
        elif StdI.InitialVecType == "r":
            iInitialVecType = 1
        else:
            print(f"\n ERROR ! Restart Mode : {StdI.Restart}")
            exit_program(-1)

    # ------------------------------------------------------------------
    #  EigenVecIO
    # ------------------------------------------------------------------
    InputEigenVec = 0
    OutputEigenVec = 0

    if StdI.EigenVecIO == "****":
        StdI.EigenVecIO = "none"
        print("       EigenVecIO = none        ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"       EigenVecIO = {StdI.EigenVecIO}")
        if StdI.EigenVecIO == "none":
            InputEigenVec = 0
        elif StdI.EigenVecIO == "in":
            InputEigenVec = 1
        elif StdI.EigenVecIO == "out":
            OutputEigenVec = 1
        elif StdI.EigenVecIO == "inout":
            InputEigenVec = 1
            OutputEigenVec = 1
        else:
            print(f"\n ERROR ! EigenVecIO Mode : {StdI.Restart}")
            exit_program(-1)

    if StdI.method == "timeevolution":
        InputEigenVec = 1

    # ------------------------------------------------------------------
    #  HamIO
    # ------------------------------------------------------------------
    iOutputHam = 0
    iInputHam = 0

    if StdI.HamIO == "****":
        StdI.HamIO = "none"
        print("         HamIO = none        ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"         HamIO = {StdI.HamIO}")
        if StdI.HamIO == "none":
            iOutputHam = 0
            iInputHam = 0
        elif StdI.HamIO == "out":
            iOutputHam = 1
        elif StdI.HamIO == "in":
            iInputHam = 1
        else:
            print(f"\n ERROR ! HamIO mode : {StdI.HamIO}")
            exit_program(-1)

    # ------------------------------------------------------------------
    #  CalcSpec
    # ------------------------------------------------------------------
    if StdI.CalcSpec == "****":
        StdI.CalcSpec = "none"
        print("         CalcSpec = none        ######  DEFAULT VALUE IS USED  ######")
        iCalcSpec = 0
    else:
        print(f"         CalcSpec = {StdI.CalcSpec}")
        if StdI.CalcSpec == "none":
            iCalcSpec = 0
        elif StdI.CalcSpec == "normal":
            iCalcSpec = 1
        elif StdI.CalcSpec == "noiteration":
            iCalcSpec = 2
        elif StdI.CalcSpec == "restart_out":
            iCalcSpec = 3
        elif StdI.CalcSpec == "restart_in":
            iCalcSpec = 4
        elif StdI.CalcSpec in ("restartsave", "restart"):
            iCalcSpec = 5
        else:
            print(f"\n ERROR ! CalcSpec : {StdI.CalcSpec}")
            exit_program(-1)

    # ------------------------------------------------------------------
    #  OutputExcitedVec
    # ------------------------------------------------------------------
    iOutputExVec = 0

    if StdI.OutputExVec == "****":
        StdI.OutputExVec = "none"
        print("         OutputExcitedVec = none        ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"         OutputExcitedVec = {StdI.OutputExVec}")
        if StdI.OutputExVec == "none":
            iOutputExVec = 0
        elif StdI.OutputExVec == "out":
            iOutputExVec = 1
        else:
            print(f"\n ERROR ! OutputExcitedVec : {StdI.OutputExVec}")
            exit_program(-1)

    # ------------------------------------------------------------------
    #  NGPU validation
    # ------------------------------------------------------------------
    if StdI.NGPU != StdI.NaN_i:
        print(f"         NGPU = {StdI.NGPU}")
        if StdI.NGPU < 1:
            print(f"\n ERROR ! NGPU : {StdI.NGPU}")
            print("         NGPU should be a positive integer.")
            exit_program(-1)

    # ------------------------------------------------------------------
    #  Scalapack validation
    # ------------------------------------------------------------------
    if StdI.Scalapack != StdI.NaN_i:
        print(f"         Scalapack = {StdI.Scalapack}")
        if StdI.Scalapack < 0 or StdI.Scalapack > 1:
            print(f"\n ERROR ! Scalapack : {StdI.Scalapack}")
            print("         Scalapack should be 0 or 1.")
            exit_program(-1)

    # ------------------------------------------------------------------
    #  Write calcmod.def
    # ------------------------------------------------------------------
    with open("calcmod.def", "w") as fp:
        fp.write("#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution 5:cTPQ\n")
        fp.write("#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC\n")
        fp.write("#Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart\n")
        fp.write("#CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save\n")
        if StdI.NGPU != StdI.NaN_i:
            fp.write("#NGPU (for FullDiag): The number of GPU\n")
        if StdI.Scalapack != StdI.NaN_i:
            fp.write("#Scalapack (for FullDiag) = 0:w/o ScaLAPACK, 1:w/ ScaLAPACK\n")
        fp.write(f"CalcType {iCalcType:3d}\n")
        fp.write(f"CalcModel {iCalcModel:3d}\n")
        fp.write(f"ReStart {iRestart:3d}\n")
        fp.write(f"CalcSpec {iCalcSpec:3d}\n")
        fp.write(f"CalcEigenVec {iCalcEigenvec:3d}\n")
        fp.write(f"InitialVecType {iInitialVecType:3d}\n")
        fp.write(f"InputEigenVec {InputEigenVec:3d}\n")
        fp.write(f"OutputEigenVec {OutputEigenVec:3d}\n")
        fp.write(f"InputHam {iInputHam:3d}\n")
        fp.write(f"OutputHam {iOutputHam:3d}\n")
        fp.write(f"OutputExVec {iOutputExVec:3d}\n")
        if StdI.NGPU != StdI.NaN_i:
            fp.write(f"NGPU {StdI.NGPU:3d}\n")
        if StdI.Scalapack != StdI.NaN_i:
            fp.write(f"Scalapack {StdI.Scalapack:3d}\n")

    print("     calcmod.def is written.\n")


# ===================================================================
#  _print_excitation  (HPhi only)
# ===================================================================


def _print_excitation(StdI: StdIntList) -> None:
    """Write ``single.def`` or ``pair.def`` for spectrum calculations.

    Depending on ``SpectrumType``, this function generates excitation
    operators with the appropriate Fourier coefficients and spin
    structure:

    - ``"szsz"`` (default) or ``"****"``: pair excitation with Sz
      diagonal coefficients.
    - ``"s+s-"``: pair excitation with S+S- ladder coefficients.
    - ``"density"``: pair excitation with unit coefficients on both
      spins.
    - ``"up"``: single excitation with spin-up.
    - ``"down"``: single excitation with spin-down.

    For the Kondo model, the Fourier coefficients for the itinerant
    sites are duplicated to the local-spin sites.

    This is the Python translation of the C function
    ``PrintExcitation()`` (lines 311--502 of ``StdFace_main.c``).

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Key fields read:

        - ``model``, ``S2`` -- model type and twice the spin quantum
          number.
        - ``SpectrumType`` -- string selecting the excitation channel.
        - ``SpectrumQ`` -- wavevector components (W, L, H).
        - ``NCell``, ``NsiteUC``, ``Cell``, ``tau`` -- lattice cell and
          basis-site data.
        - ``nsite`` -- total number of sites.
        - ``pi`` -- the mathematical constant pi.

        Modified in place:

        - ``SpectrumBody`` -- set to 1 (single) or 2 (pair).
        - ``SpectrumQ`` -- defaults filled via :func:`print_val_d`.
    """
    # --- Allocate coefficient and spin-index arrays ---
    if StdI.model == "spin" and StdI.S2 > 1:
        n_alloc = StdI.S2 + 1
    else:
        n_alloc = 2

    coef = [0.0] * n_alloc
    spin = [[0, 0] for _ in range(n_alloc)]

    fourier_r = [0.0] * StdI.nsite
    fourier_i = [0.0] * StdI.nsite

    print("\n  @ Spectrum\n")

    StdI.SpectrumQ[0] = print_val_d("SpectrumQW", StdI.SpectrumQ[0], 0.0)
    StdI.SpectrumQ[1] = print_val_d("SpectrumQL", StdI.SpectrumQ[1], 0.0)
    StdI.SpectrumQ[2] = print_val_d("SpectrumQH", StdI.SpectrumQ[2], 0.0)

    # ------------------------------------------------------------------
    #  Determine NumOp, coef, spin based on SpectrumType
    # ------------------------------------------------------------------
    if StdI.SpectrumType == "****":
        StdI.SpectrumType = "szsz"
        print("     SpectrumType = szsz        ######  DEFAULT VALUE IS USED  ######")
        if StdI.model == "spin":
            NumOp = StdI.S2 + 1
            for ispin in range(StdI.S2 + 1):
                Sz = float(ispin) - float(StdI.S2) * 0.5
                coef[ispin] = Sz
                spin[ispin][0] = ispin
                spin[ispin][1] = ispin
        else:
            NumOp = 2
            coef[0] = 0.5
            coef[1] = -0.5
            spin[0][0] = 0
            spin[0][1] = 0
            spin[1][0] = 1
            spin[1][1] = 1
        StdI.SpectrumBody = 2
    else:
        print(f"     SpectrumType = {StdI.SpectrumType}")
        if StdI.SpectrumType == "szsz":
            if StdI.model == "spin":
                NumOp = StdI.S2 + 1
                for ispin in range(StdI.S2 + 1):
                    Sz = float(ispin) - float(StdI.S2) * 0.5
                    coef[ispin] = Sz
                    spin[ispin][0] = ispin
                    spin[ispin][1] = ispin
            else:
                NumOp = 2
                coef[0] = 0.5
                coef[1] = -0.5
                spin[0][0] = 0
                spin[0][1] = 0
                spin[1][0] = 1
                spin[1][1] = 1
            StdI.SpectrumBody = 2

        elif StdI.SpectrumType == "s+s-":
            if StdI.model == "spin" and StdI.S2 > 1:
                NumOp = StdI.S2
                S = float(StdI.S2) * 0.5
                for ispin in range(1, StdI.S2 + 1):
                    Sz = float(StdI.S2) * 0.5 - float(ispin)
                    coef[ispin - 1] = math.sqrt(S * (S + 1.0) - Sz * (Sz + 1.0))
                    spin[ispin - 1][0] = ispin
                    spin[ispin - 1][1] = ispin - 1
            else:
                NumOp = 1
                coef[0] = 1.0
                spin[0][0] = 1
                spin[0][1] = 0
            StdI.SpectrumBody = 2

        elif StdI.SpectrumType == "density":
            NumOp = 2
            coef[0] = 1.0
            coef[1] = 1.0
            spin[0][0] = 0
            spin[0][1] = 0
            spin[1][0] = 1
            spin[1][1] = 1
            StdI.SpectrumBody = 2

        elif StdI.SpectrumType == "up":
            NumOp = 1
            coef[0] = 1.0
            spin[0][0] = 0
            StdI.SpectrumBody = 1

        elif StdI.SpectrumType == "down":
            NumOp = 1
            coef[0] = 1.0
            spin[0][0] = 1
            StdI.SpectrumBody = 1

        else:
            print(f"\n ERROR ! SpectrumType : {StdI.SpectrumType}")
            exit_program(-1)

    # ------------------------------------------------------------------
    #  Compute Fourier coefficients
    # ------------------------------------------------------------------
    isite = 0
    for icell in range(StdI.NCell):
        for itau in range(StdI.NsiteUC):
            Cphase = ((StdI.Cell[icell][0] + StdI.tau[itau][0]) * StdI.SpectrumQ[0]
                      + (StdI.Cell[icell][1] + StdI.tau[itau][1]) * StdI.SpectrumQ[1]
                      + (StdI.Cell[icell][2] + StdI.tau[itau][2]) * StdI.SpectrumQ[2])
            fourier_r[isite] = math.cos(2.0 * StdI.pi * Cphase)
            fourier_i[isite] = math.sin(2.0 * StdI.pi * Cphase)
            isite += 1

    # For the Kondo model, duplicate Fourier coefficients
    if StdI.model == "kondo":
        half = StdI.nsite // 2
        for isite in range(half):
            fourier_r[isite + half] = fourier_r[isite]
            fourier_i[isite + half] = fourier_i[isite]

    # ------------------------------------------------------------------
    #  Write single.def or pair.def
    # ------------------------------------------------------------------
    if StdI.SpectrumBody == 1:
        with open("single.def", "w") as fp:
            fp.write("=============================================\n")
            if StdI.model == "kondo":
                fp.write(f"NSingle {StdI.nsite // 2 * NumOp}\n")
            else:
                fp.write(f"NSingle {StdI.nsite * NumOp}\n")
            fp.write("=============================================\n")
            fp.write("============== Single Excitation ============\n")
            fp.write("=============================================\n")
            if StdI.model == "kondo":
                for isite in range(StdI.nsite // 2, StdI.nsite):
                    fp.write(f"{isite} {spin[0][0]} 0 "
                             f"{fourier_r[isite] * coef[0]:25.15f} "
                             f"{fourier_i[isite] * coef[0]:25.15f}\n")
            else:
                for isite in range(StdI.nsite):
                    fp.write(f"{isite} {spin[0][0]} 0 "
                             f"{fourier_r[isite] * coef[0]:25.15f} "
                             f"{fourier_i[isite] * coef[0]:25.15f}\n")
        print("      single.def is written.\n")

    else:
        with open("pair.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"NPair {StdI.nsite * NumOp}\n")
            fp.write("=============================================\n")
            fp.write("=============== Pair Excitation =============\n")
            fp.write("=============================================\n")
            for isite in range(StdI.nsite):
                for ispin in range(NumOp):
                    fp.write(f"{isite} {spin[ispin][0]} {isite} {spin[ispin][1]} 1 "
                             f"{fourier_r[isite] * coef[ispin]:25.15f} "
                             f"{fourier_i[isite] * coef[ispin]:25.15f}\n")
        print("        pair.def is written.\n")


# ===================================================================
#  _vector_potential  (HPhi only)
# ===================================================================


def _vector_potential(StdI: StdIntList) -> None:
    """Compute the vector potential A(t) and electric field E(t) for time evolution.

    Depending on ``PumpType``, the time-dependent vector potential and
    electric field are computed as follows:

    - ``"quench"`` (default): no laser field; uses two-body pump
      (``PumpBody = 2``).
    - ``"pulselaser"``: Gaussian-enveloped cosine pulse
      (``PumpBody = 1``).
    - ``"aclaser"``: sinusoidal AC laser (``PumpBody = 1``).
    - ``"dclaser"``: linearly ramped DC laser (``PumpBody = 1``).

    For ``PumpBody == 1``, the file ``potential.dat`` is written with
    columns ``time, A_W, A_L, A_H, E_W, E_L, E_H``.

    This is the Python translation of the C function
    ``VectorPotential()`` (lines 506--597 of ``StdFace_main.c``).

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Key fields read / set:

        - ``VecPot`` -- initial vector potential components (W, L, H).
        - ``Lanczos_max`` -- number of time steps.
        - ``dt`` -- time step width.
        - ``freq`` -- laser frequency.
        - ``tshift`` -- time shift for the pulse center.
        - ``tdump`` -- Gaussian pulse width.
        - ``Uquench`` -- quench interaction strength.
        - ``ExpandCoef`` -- expansion coefficient order.
        - ``PumpType`` -- string selecting the pump type.
        - ``PumpBody`` -- set to 1 (one-body) or 2 (two-body).
        - ``At`` -- allocated as ``Lanczos_max x 3`` list of floats.
    """
    print("\n  @ Time-evolution\n")

    StdI.VecPot[0] = print_val_d("VecPotW", StdI.VecPot[0], 0.0)
    StdI.VecPot[1] = print_val_d("VecPotL", StdI.VecPot[1], 0.0)
    StdI.VecPot[2] = print_val_d("VecPotH", StdI.VecPot[2], 0.0)
    StdI.Lanczos_max = print_val_i("Lanczos_max", StdI.Lanczos_max, 1000)
    StdI.dt = print_val_d("dt", StdI.dt, 0.01)
    StdI.freq = print_val_d("freq", StdI.freq, 0.1)
    StdI.tshift = print_val_d("tshift", StdI.tshift, 0.0)
    StdI.tdump = print_val_d("tdump", StdI.tdump, 0.1)
    StdI.Uquench = print_val_d("Uquench", StdI.Uquench, 0.0)
    StdI.ExpandCoef = print_val_i("ExpandCoef", StdI.ExpandCoef, 10)

    # Allocate A(t) and E(t) arrays: Lanczos_max x 3
    StdI.At = [[0.0] * 3 for _ in range(StdI.Lanczos_max)]
    Et = [[0.0] * 3 for _ in range(StdI.Lanczos_max)]

    if StdI.PumpType == "****":
        StdI.PumpType = "quench"
        print("     PumpType = quench        ######  DEFAULT VALUE IS USED  ######")
        StdI.PumpBody = 2
    else:
        print(f"     PumpType = {StdI.PumpType}")
        if StdI.PumpType == "quench":
            StdI.PumpBody = 2

        elif StdI.PumpType == "pulselaser":
            for it in range(StdI.Lanczos_max):
                time = StdI.dt * float(it)
                for ii in range(3):
                    StdI.At[it][ii] = (StdI.VecPot[ii]
                                       * math.cos(StdI.freq * (time - StdI.tshift))
                                       * math.exp(-0.5 * (time - StdI.tshift) ** 2
                                                   / (StdI.tdump * StdI.tdump)))
                    Et[it][ii] = (-StdI.VecPot[ii]
                                  * ((StdI.tshift - time) / (StdI.tdump * StdI.tdump)
                                     * math.cos(StdI.freq * (time - StdI.tshift))
                                     - StdI.freq * math.sin(StdI.freq * (time - StdI.tshift)))
                                  * math.exp(-0.5 * (time - StdI.tshift) ** 2
                                              / (StdI.tdump * StdI.tdump)))
            StdI.PumpBody = 1

        elif StdI.PumpType == "aclaser":
            for it in range(StdI.Lanczos_max):
                time = StdI.dt * float(it)
                for ii in range(3):
                    StdI.At[it][ii] = StdI.VecPot[ii] * math.sin(StdI.freq * (time - StdI.tshift))
                    Et[it][ii] = StdI.VecPot[ii] * math.cos(StdI.freq * (time - StdI.tshift)) * StdI.freq
            StdI.PumpBody = 1

        elif StdI.PumpType == "dclaser":
            for it in range(StdI.Lanczos_max):
                time = StdI.dt * float(it)
                for ii in range(3):
                    StdI.At[it][ii] = StdI.VecPot[ii] * time
                    Et[it][ii] = -StdI.VecPot[ii]
            StdI.PumpBody = 1

        else:
            print(f"\n ERROR ! PumpType : {StdI.PumpType}")
            exit_program(-1)

    # ------------------------------------------------------------------
    #  Write potential.dat for one-body pump
    # ------------------------------------------------------------------
    if StdI.PumpBody == 1:
        with open("potential.dat", "w") as fp:
            fp.write("# Time A_W A_L A_H E_W E_L E_H\n")
            for it in range(StdI.Lanczos_max):
                time = StdI.dt * float(it)
                fp.write(f"{time:f} "
                         f"{StdI.At[it][0]:f} {StdI.At[it][1]:f} {StdI.At[it][2]:f} "
                         f"{Et[it][0]:f} {Et[it][1]:f} {Et[it][2]:f}\n")


# ===================================================================
#  _print_pump  (HPhi only)
# ===================================================================


def _print_pump(StdI: StdIntList) -> None:
    """Write ``teone.def`` or ``tetwo.def`` for time-evolution pump terms.

    - ``PumpBody == 1``: writes ``teone.def`` -- one-body pump terms.
      Equivalent pump terms (same index quadruples) are merged and
      entries below ``1e-6`` in magnitude are suppressed.
    - ``PumpBody == 2``: writes ``tetwo.def`` -- two-body pump terms
      using ``Uquench`` for on-site Hubbard quench interaction.

    This is the Python translation of the C function
    ``PrintPump()`` (lines 602--667 of ``StdFace_main.c``).

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Key fields read:

        - ``PumpBody`` -- 1 for one-body, 2 for two-body.
        - ``Lanczos_max`` -- number of time steps.
        - ``dt`` -- time step width.
        - ``npump`` -- list of int, per-timestep pump term counts.
        - ``pumpindx`` -- 3-D list of int, shape
          ``(Lanczos_max, npump[it], 4)`` -- site/spin indices.
        - ``pump`` -- 2-D list of complex, shape
          ``(Lanczos_max, npump[it])`` -- pump amplitudes.
        - ``nsite`` -- total number of sites.
        - ``Uquench`` -- quench interaction strength.
    """
    if StdI.PumpBody == 1:
        with open("teone.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"AllTimeStep {StdI.Lanczos_max}\n")
            fp.write("=============================================\n")
            fp.write("=========  OneBody Time Evolution  ==========\n")
            fp.write("=============================================\n")

            for it in range(StdI.Lanczos_max):
                # Sum equivalent pumping terms
                for ipump in range(StdI.npump[it]):
                    for jpump in range(ipump + 1, StdI.npump[it]):
                        if (StdI.pumpindx[it][ipump][0] == StdI.pumpindx[it][jpump][0]
                                and StdI.pumpindx[it][ipump][1] == StdI.pumpindx[it][jpump][1]
                                and StdI.pumpindx[it][ipump][2] == StdI.pumpindx[it][jpump][2]
                                and StdI.pumpindx[it][ipump][3] == StdI.pumpindx[it][jpump][3]):
                            StdI.pump[it][ipump] = StdI.pump[it][ipump] + StdI.pump[it][jpump]
                            StdI.pump[it][jpump] = 0.0

                # Count the number of finite pumping terms
                npump0 = 0
                for ipump in range(StdI.npump[it]):
                    if abs(StdI.pump[it][ipump]) > 0.000001:
                        npump0 += 1

                fp.write(f"{StdI.dt * float(it):f}  {npump0}\n")

                for ipump in range(StdI.npump[it]):
                    if abs(StdI.pump[it][ipump]) <= 0.000001:
                        continue
                    fp.write(f"{StdI.pumpindx[it][ipump][0]:5d} "
                             f"{StdI.pumpindx[it][ipump][1]:5d} "
                             f"{StdI.pumpindx[it][ipump][2]:5d} "
                             f"{StdI.pumpindx[it][ipump][3]:5d} "
                             f"{StdI.pump[it][ipump].real:25.15f} "
                             f"{StdI.pump[it][ipump].imag:25.15f}\n")

        print("      teone.def is written.\n")

    else:
        with open("tetwo.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"AllTimeStep {StdI.Lanczos_max}\n")
            fp.write("=============================================\n")
            fp.write("========== TwoBody Time Evolution ===========\n")
            fp.write("=============================================\n")

            for it in range(StdI.Lanczos_max):
                fp.write(f"{StdI.dt * float(it):f}  {StdI.nsite}\n")
                for isite in range(StdI.nsite):
                    fp.write(f"{isite:5d} {0:5d} {isite:5d} {0:5d} "
                             f"{isite:5d} {1:5d} {isite:5d} {1:5d} "
                             f"{StdI.Uquench:25.15f}  {0.0:25.15f}\n")

        print("        tetwo.def is written.\n")


# ===================================================================
#  Part 5 -- mVMC-specific output functions
# ===================================================================


def _print_orb(StdI: StdIntList) -> None:
    """Write the anti-parallel orbital index file ``orbitalidx.def``.

    This is the Python translation of the C function ``PrintOrb()``
    (lines 673-704 of ``StdFace_main.c``).  The file records the
    orbital pairing indices used by mVMC for the anti-parallel-spin
    part of the variational wave function.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``nsite`` -- total number of sites
        - ``NOrb`` -- number of orbital indices
        - ``ComplexType`` -- 0 for real, 1 for complex
        - ``Orb`` -- ``nsite x nsite`` orbital index matrix
        - ``AntiOrb`` -- ``nsite x nsite`` anti-periodic sign matrix
        - ``AntiPeriod`` -- length-3 array of anti-periodic boundary flags

    Notes
    -----
    The C version frees ``StdI->Orb`` at the end; in Python the
    garbage collector handles deallocation, so no explicit free is
    performed.
    """
    with open("orbitalidx.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NOrbitalIdx {StdI.NOrb:10d}\n")
        fp.write(f"ComplexType {StdI.ComplexType:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")

        has_anti = (StdI.AntiPeriod[0] == 1
                    or StdI.AntiPeriod[1] == 1
                    or StdI.AntiPeriod[2] == 1)

        for isite in range(StdI.nsite):
            for jsite in range(StdI.nsite):
                if has_anti:
                    fp.write(f"{isite:5d}  {jsite:5d}  "
                             f"{StdI.Orb[isite][jsite]:5d}  "
                             f"{StdI.AntiOrb[isite][jsite]:5d}\n")
                else:
                    fp.write(f"{isite:5d}  {jsite:5d}  "
                             f"{StdI.Orb[isite][jsite]:5d}\n")

        for iOrb in range(StdI.NOrb):
            fp.write(f"{iOrb:5d}  {1:5d}\n")

    print("    orbitalidx.def is written.")


def _print_orb_para(StdI: StdIntList) -> None:
    """Write parallel orbital index files ``orbitalidxpara.def`` and ``orbitalidxgen.def``.

    This is the Python translation of the C function ``PrintOrbPara()``
    (lines 709-839 of ``StdFace_main.c``).  The routine performs the
    following steps:

    1. Copy the anti-parallel orbital matrix (``StdI.Orb``) into a local
       ``OrbGC`` matrix and the sign matrix (``StdI.AntiOrb``) into
       ``reverse``.
    2. Symmetrise: for each known orbital index, set
       ``OrbGC[j][i] = OrbGC[i][j]`` and ``reverse[j][i] = -reverse[i][j]``.
    3. Renumber: walk the strict lower triangle (``isite > jsite``),
       assign negative indices to each newly seen orbital, then invert
       all indices so that they become non-negative.
    4. Write ``orbitalidxpara.def`` with the upper triangle
       (``isite < jsite``), four columns per line.
    5. Write ``orbitalidxgen.def`` combining anti-parallel (all pairs) and
       parallel (upper triangle) sections, offsetting parallel orbital
       indices by ``NOrb``.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``nsite`` -- total number of sites
        - ``NOrb`` -- number of anti-parallel orbital indices
        - ``ComplexType`` -- 0 for real, 1 for complex
        - ``Orb`` -- ``nsite x nsite`` orbital index matrix
        - ``AntiOrb`` -- ``nsite x nsite`` anti-periodic sign matrix
        - ``AntiPeriod`` -- length-3 array of anti-periodic boundary flags
    """
    nsite = StdI.nsite

    # -----------------------------------------------------------------
    # (1) Copy from anti-parallel orbital index
    # -----------------------------------------------------------------
    OrbGC = [[0] * nsite for _ in range(nsite)]
    reverse = [[0] * nsite for _ in range(nsite)]
    for isite in range(nsite):
        for jsite in range(nsite):
            OrbGC[isite][jsite] = int(StdI.Orb[isite][jsite])
            reverse[isite][jsite] = int(StdI.AntiOrb[isite][jsite])

    # -----------------------------------------------------------------
    # (2) Symmetrise
    # -----------------------------------------------------------------
    for iorb in range(StdI.NOrb):
        for isite in range(nsite):
            for jsite in range(nsite):
                if OrbGC[isite][jsite] == iorb:
                    OrbGC[jsite][isite] = OrbGC[isite][jsite]
                    reverse[jsite][isite] = -reverse[isite][jsite]

    # -----------------------------------------------------------------
    # (3) Renumber -- lower triangle (isite > jsite)
    # -----------------------------------------------------------------
    NOrbGC = 0
    for isite in range(nsite):
        for jsite in range(isite):
            if OrbGC[isite][jsite] >= 0:
                iOrbGC = OrbGC[isite][jsite]
                NOrbGC -= 1
                for isite1 in range(nsite):
                    for jsite1 in range(nsite):
                        if OrbGC[isite1][jsite1] == iOrbGC:
                            OrbGC[isite1][jsite1] = NOrbGC

    NOrbGC = -NOrbGC
    for isite in range(nsite):
        for jsite in range(nsite):
            OrbGC[isite][jsite] = -1 - OrbGC[isite][jsite]

    # -----------------------------------------------------------------
    # Write orbitalidxpara.def
    # -----------------------------------------------------------------
    with open("orbitalidxpara.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NOrbitalIdx {NOrbGC:10d}\n")
        fp.write(f"ComplexType {StdI.ComplexType:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")

        for isite in range(nsite):
            for jsite in range(nsite):
                if isite >= jsite:
                    continue
                fp.write(f"{isite:5d}  {jsite:5d}  "
                         f"{OrbGC[isite][jsite]:5d}  "
                         f"{reverse[isite][jsite]:5d}\n")

        for iOrbGC in range(NOrbGC):
            fp.write(f"{iOrbGC:5d}  {1:5d}\n")

    print("    orbitalidxpara.def is written.")

    # -----------------------------------------------------------------
    # Write orbitalidxgen.def
    # -----------------------------------------------------------------
    has_anti = (StdI.AntiPeriod[0] == 1
                or StdI.AntiPeriod[1] == 1
                or StdI.AntiPeriod[2] == 1)

    with open("orbitalidxgen.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NOrbitalIdx {StdI.NOrb + 2 * NOrbGC:10d}\n")
        fp.write(f"ComplexType {StdI.ComplexType:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")

        # -- anti-parallel section --
        for isite in range(nsite):
            for jsite in range(nsite):
                if has_anti:
                    fp.write(f"{isite:5d}  0  {jsite:5d}  1  "
                             f"{StdI.Orb[isite][jsite]:5d}  "
                             f"{StdI.AntiOrb[isite][jsite]:5d}\n")
                else:
                    fp.write(f"{isite:5d}  0  {jsite:5d}  1  "
                             f"{StdI.Orb[isite][jsite]:5d}  {1:5d}\n")

        # -- parallel section (upper triangle) --
        for isite in range(nsite):
            for jsite in range(nsite):
                if isite >= jsite:
                    continue
                if has_anti:
                    fp.write(f"{isite:5d}  0  {jsite:5d}  0  "
                             f"{OrbGC[isite][jsite] + StdI.NOrb:5d}  "
                             f"{reverse[isite][jsite]:5d}\n")
                    fp.write(f"{isite:5d}  1  {jsite:5d}  1  "
                             f"{OrbGC[isite][jsite] + StdI.NOrb + NOrbGC:5d}  "
                             f"{reverse[isite][jsite]:5d}\n")
                else:
                    fp.write(f"{isite:5d}  0  {jsite:5d}  0  "
                             f"{OrbGC[isite][jsite] + StdI.NOrb:5d}  "
                             f"{reverse[isite][jsite]:5d}\n")
                    fp.write(f"{isite:5d}  1  {jsite:5d}  1  "
                             f"{OrbGC[isite][jsite] + StdI.NOrb + NOrbGC:5d}  "
                             f"{reverse[isite][jsite]:5d}\n")

        for iOrbGC in range(StdI.NOrb):
            fp.write(f"{iOrbGC:5d}  {1:5d}\n")

        for iOrbGC in range(NOrbGC * 2):
            fp.write(f"{iOrbGC + StdI.NOrb:5d}  {1:5d}\n")

    print("    orbitalidxgen.def is written.")


def _print_gutzwiller(StdI: StdIntList) -> None:
    """Write the Gutzwiller variational-parameter file ``gutzwilleridx.def``.

    This is the Python translation of the C function
    ``PrintGutzwiller()`` (lines 843-921 of ``StdFace_main.c``).  The
    file defines the Gutzwiller projection indices used by mVMC.

    Two modes are supported depending on the value of ``NMPTrans``:

    **Momentum-projected mode** (``abs(NMPTrans) == 1`` or ``NMPTrans``
    is unset):

    - For the Hubbard model, ``NGutzwiller`` starts at 0.
    - For other models, ``NGutzwiller`` starts at -1.
    - Diagonal orbital indices (``Orb[i][i]``) are used; local-spin
      sites are excluded (set to -1).
    - Unique Gutzwiller indices are renumbered with negative temporaries
      and then inverted.

    **Global-optimisation mode** (all other ``NMPTrans``):

    - Hubbard: ``NGutzwiller = NsiteUC``, site index modulo ``NsiteUC``.
    - Spin: ``NGutzwiller = 1``, all sites map to index 0.
    - Kondo: ``NGutzwiller = NsiteUC + 1``, conduction sites map to 0,
      localised sites map to ``isite + 1``.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  The following fields are read:

        - ``nsite`` -- total number of sites
        - ``NMPTrans`` -- momentum-projection control
        - ``NaN_i`` -- integer sentinel value
        - ``model`` -- model name string
        - ``Orb`` -- ``nsite x nsite`` orbital index matrix
        - ``locspinflag`` -- per-site local-spin flag array
        - ``NsiteUC`` -- number of sites per unit cell
        - ``NCell`` -- number of unit cells
    """
    nsite = StdI.nsite
    Gutz = [0] * nsite

    if abs(StdI.NMPTrans) == 1 or StdI.NMPTrans == StdI.NaN_i:
        # ---- Momentum-projected mode ----
        if StdI.model == "hubbard":
            NGutzwiller = 0
        else:
            NGutzwiller = -1

        for isite in range(nsite):
            Gutz[isite] = int(StdI.Orb[isite][isite])

        for isite in range(nsite):
            # For local spin sites
            if StdI.locspinflag[isite] != 0:
                Gutz[isite] = -1
                continue
            # Renumber
            if Gutz[isite] >= 0:
                iGutz = Gutz[isite]
                NGutzwiller -= 1
                for jsite in range(nsite):
                    if Gutz[jsite] == iGutz:
                        Gutz[jsite] = NGutzwiller

        NGutzwiller = -NGutzwiller
        for isite in range(nsite):
            Gutz[isite] = -1 - Gutz[isite]

    else:
        # ---- Global-optimisation mode ----
        if StdI.model == "hubbard":
            NGutzwiller = StdI.NsiteUC
        elif StdI.model == "spin":
            NGutzwiller = 1
        else:
            NGutzwiller = StdI.NsiteUC + 1

        for iCell in range(StdI.NCell):
            for isite in range(StdI.NsiteUC):
                if StdI.model == "hubbard":
                    Gutz[isite + StdI.NsiteUC * iCell] = isite
                elif StdI.model == "spin":
                    Gutz[isite + StdI.NsiteUC * iCell] = 0
                else:
                    Gutz[isite + StdI.NsiteUC * iCell] = 0
                    Gutz[isite + StdI.NsiteUC * (iCell + StdI.NCell)] = isite + 1

    # -----------------------------------------------------------------
    # Write gutzwilleridx.def
    # -----------------------------------------------------------------
    with open("gutzwilleridx.def", "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NGutzwillerIdx {NGutzwiller:10d}\n")
        fp.write(f"ComplexType {0:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")

        for isite in range(nsite):
            fp.write(f"{isite:5d}  {Gutz[isite]:5d}\n")

        for iGutz in range(NGutzwiller):
            if StdI.model == "hubbard" or iGutz > 0:
                fp.write(f"{iGutz:5d}  {1:5d}\n")
            else:
                fp.write(f"{iGutz:5d}  {0:5d}\n")

    print("    gutzwilleridx.def is written.")


def _print_interactions(StdI: "StdIntList") -> None:
    """Process and write definition files for all interaction types.

    For each interaction type (CoulombIntra, CoulombInter, Hund, Exchange,
    PairLift, PairHopp, InterAll), this function:

    1. Merges duplicate terms by summing their coefficients and zeroing
       out the duplicate entry.
    2. Counts the number of non-zero terms.
    3. Sets a flag (e.g. ``LCintra``, ``LCinter``, ...) based on the count
       and whether boost mode is active.
    4. If the flag is set, writes the corresponding ``.def`` file.

    The InterAll section is more complex, involving three merge/reorder
    passes before file output.

    Parameters
    ----------
    StdI : StdIntList
        The central data structure holding all interaction arrays, index
        arrays, counts, and flags.  Modified in place.
    """
    # =================================================================
    #  Coulomb Intra
    # =================================================================
    for kintr in range(StdI.NCintra):
        for jintr in range(kintr + 1, StdI.NCintra):
            if StdI.CintraIndx[jintr][0] == StdI.CintraIndx[kintr][0]:
                StdI.Cintra[kintr] += StdI.Cintra[jintr]
                StdI.Cintra[jintr] = 0.0

    nintr0 = 0
    for kintr in range(StdI.NCintra):
        if abs(StdI.Cintra[kintr]) > 0.000001:
            nintr0 += 1

    if nintr0 == 0 or StdI.lBoost == 1:
        StdI.LCintra = 0
    else:
        StdI.LCintra = 1

    if StdI.LCintra == 1:
        with open("coulombintra.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"NCoulombIntra {nintr0:10d}\n")
            fp.write("=============================================\n")
            fp.write("================== CoulombIntra ================\n")
            fp.write("=============================================\n")
            for kintr in range(StdI.NCintra):
                if abs(StdI.Cintra[kintr]) > 0.000001:
                    fp.write(f"{StdI.CintraIndx[kintr][0]:5d} "
                             f"{StdI.Cintra[kintr]:25.15f}\n")
        print("    coulombintra.def is written.")

    # =================================================================
    #  Coulomb Inter
    # =================================================================
    for kintr in range(StdI.NCinter):
        for jintr in range(kintr + 1, StdI.NCinter):
            j0 = StdI.CinterIndx[jintr][0]
            j1 = StdI.CinterIndx[jintr][1]
            k0 = StdI.CinterIndx[kintr][0]
            k1 = StdI.CinterIndx[kintr][1]
            if (j0 == k0 and j1 == k1) or (j0 == k1 and j1 == k0):
                StdI.Cinter[kintr] += StdI.Cinter[jintr]
                StdI.Cinter[jintr] = 0.0

    nintr0 = 0
    for kintr in range(StdI.NCinter):
        if abs(StdI.Cinter[kintr]) > 0.000001:
            nintr0 += 1

    if nintr0 == 0 or StdI.lBoost == 1:
        StdI.LCinter = 0
    else:
        StdI.LCinter = 1

    if StdI.LCinter == 1:
        with open("coulombinter.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"NCoulombInter {nintr0:10d}\n")
            fp.write("=============================================\n")
            fp.write("================== CoulombInter ================\n")
            fp.write("=============================================\n")
            for kintr in range(StdI.NCinter):
                if abs(StdI.Cinter[kintr]) > 0.000001:
                    fp.write(f"{StdI.CinterIndx[kintr][0]:5d} "
                             f"{StdI.CinterIndx[kintr][1]:5d} "
                             f"{StdI.Cinter[kintr]:25.15f}\n")
        print("    coulombinter.def is written.")

    # =================================================================
    #  Hund
    # =================================================================
    for kintr in range(StdI.NHund):
        for jintr in range(kintr + 1, StdI.NHund):
            j0 = StdI.HundIndx[jintr][0]
            j1 = StdI.HundIndx[jintr][1]
            k0 = StdI.HundIndx[kintr][0]
            k1 = StdI.HundIndx[kintr][1]
            if (j0 == k0 and j1 == k1) or (j0 == k1 and j1 == k0):
                StdI.Hund[kintr] += StdI.Hund[jintr]
                StdI.Hund[jintr] = 0.0

    nintr0 = 0
    for kintr in range(StdI.NHund):
        if abs(StdI.Hund[kintr]) > 0.000001:
            nintr0 += 1

    if nintr0 == 0 or StdI.lBoost == 1:
        StdI.LHund = 0
    else:
        StdI.LHund = 1

    if StdI.LHund == 1:
        with open("hund.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"NHund {nintr0:10d}\n")
            fp.write("=============================================\n")
            fp.write("=============== Hund coupling ===============\n")
            fp.write("=============================================\n")
            for kintr in range(StdI.NHund):
                if abs(StdI.Hund[kintr]) > 0.000001:
                    fp.write(f"{StdI.HundIndx[kintr][0]:5d} "
                             f"{StdI.HundIndx[kintr][1]:5d} "
                             f"{StdI.Hund[kintr]:25.15f}\n")
        print("    hund.def is written.")

    # =================================================================
    #  Exchange
    # =================================================================
    for kintr in range(StdI.NEx):
        for jintr in range(kintr + 1, StdI.NEx):
            j0 = StdI.ExIndx[jintr][0]
            j1 = StdI.ExIndx[jintr][1]
            k0 = StdI.ExIndx[kintr][0]
            k1 = StdI.ExIndx[kintr][1]
            if (j0 == k0 and j1 == k1) or (j0 == k1 and j1 == k0):
                StdI.Ex[kintr] += StdI.Ex[jintr]
                StdI.Ex[jintr] = 0.0

    nintr0 = 0
    for kintr in range(StdI.NEx):
        if abs(StdI.Ex[kintr]) > 0.000001:
            nintr0 += 1

    if nintr0 == 0 or StdI.lBoost == 1:
        StdI.LEx = 0
    else:
        StdI.LEx = 1

    if StdI.LEx == 1:
        with open("exchange.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"NExchange {nintr0:10d}\n")
            fp.write("=============================================\n")
            fp.write("====== ExchangeCoupling coupling ============\n")
            fp.write("=============================================\n")
            for kintr in range(StdI.NEx):
                if abs(StdI.Ex[kintr]) > 0.000001:
                    fp.write(f"{StdI.ExIndx[kintr][0]:5d} "
                             f"{StdI.ExIndx[kintr][1]:5d} "
                             f"{StdI.Ex[kintr]:25.15f}\n")
        print("    exchange.def is written.")

    # =================================================================
    #  PairLift
    # =================================================================
    for kintr in range(StdI.NPairLift):
        for jintr in range(kintr + 1, StdI.NPairLift):
            j0 = StdI.PLIndx[jintr][0]
            j1 = StdI.PLIndx[jintr][1]
            k0 = StdI.PLIndx[kintr][0]
            k1 = StdI.PLIndx[kintr][1]
            if (j0 == k0 and j1 == k1) or (j0 == k1 and j1 == k0):
                StdI.PairLift[kintr] += StdI.PairLift[jintr]
                StdI.PairLift[jintr] = 0.0

    nintr0 = 0
    for kintr in range(StdI.NPairLift):
        if abs(StdI.PairLift[kintr]) > 0.000001:
            nintr0 += 1

    if nintr0 == 0 or StdI.lBoost == 1:
        StdI.LPairLift = 0
    else:
        StdI.LPairLift = 1

    if StdI.LPairLift == 1:
        with open("pairlift.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"NPairLift {nintr0:10d}\n")
            fp.write("=============================================\n")
            fp.write("====== Pair-Lift term ============\n")
            fp.write("=============================================\n")
            for kintr in range(StdI.NPairLift):
                if abs(StdI.PairLift[kintr]) > 0.000001:
                    fp.write(f"{StdI.PLIndx[kintr][0]:5d} "
                             f"{StdI.PLIndx[kintr][1]:5d} "
                             f"{StdI.PairLift[kintr]:25.15f}\n")
        print("    pairlift.def is written.")

    # =================================================================
    #  PairHopp
    # =================================================================
    for kintr in range(StdI.NPairHopp):
        for jintr in range(kintr + 1, StdI.NPairHopp):
            j0 = StdI.PHIndx[jintr][0]
            j1 = StdI.PHIndx[jintr][1]
            k0 = StdI.PHIndx[kintr][0]
            k1 = StdI.PHIndx[kintr][1]
            if (j0 == k0 and j1 == k1) or (j0 == k1 and j1 == k0):
                StdI.PairHopp[kintr] += StdI.PairHopp[jintr]
                StdI.PairHopp[jintr] = 0.0

    nintr0 = 0
    for kintr in range(StdI.NPairHopp):
        if abs(StdI.PairHopp[kintr]) > 0.000001:
            nintr0 += 1

    if nintr0 == 0 or StdI.lBoost == 1:
        StdI.LPairHopp = 0
    else:
        StdI.LPairHopp = 1

    if StdI.LPairHopp == 1:
        with open("pairhopp.def", "w") as fp:
            fp.write("=============================================\n")
            fp.write(f"NPairHopp {nintr0:10d}\n")
            fp.write("=============================================\n")
            fp.write("====== Pair-Hopping term ============\n")
            fp.write("=============================================\n")
            for kintr in range(StdI.NPairHopp):
                if abs(StdI.PairHopp[kintr]) > 0.000001:
                    fp.write(f"{StdI.PHIndx[kintr][0]:5d} "
                             f"{StdI.PHIndx[kintr][1]:5d} "
                             f"{StdI.PairHopp[kintr]:25.15f}\n")
        print("    pairhopp.def is written.")

    # =================================================================
    #  InterAll
    # =================================================================

    # -----------------------------------------------------------------
    #  Pass 1: Merge equivalent terms
    # -----------------------------------------------------------------
    for jintr in range(StdI.nintr):
        j = StdI.intrindx[jintr]  # shorthand for the 8-element index row
        for kintr in range(jintr + 1, StdI.nintr):
            k = StdI.intrindx[kintr]

            # Case A: all 8 indices match exactly
            # Case B: Hermitian conjugate match (indices 0-3 <-> 4-7)
            #         with two exclusion conditions
            if (
                (j[0] == k[0] and j[1] == k[1]
                 and j[2] == k[2] and j[3] == k[3]
                 and j[4] == k[4] and j[5] == k[5]
                 and j[6] == k[6] and j[7] == k[7])
                or
                (j[0] == k[4] and j[1] == k[5]
                 and j[2] == k[6] and j[3] == k[7]
                 and j[4] == k[0] and j[5] == k[1]
                 and j[6] == k[2] and j[7] == k[3]
                 and not (j[0] == j[6] and j[1] == j[7])
                 and not (j[2] == j[4] and j[3] == j[5]))
            ):
                StdI.intr[jintr] = StdI.intr[jintr] + StdI.intr[kintr]
                StdI.intr[kintr] = 0.0

            # Case C: Partial swap patterns -> SUBTRACT
            elif (
                (j[0] == k[4] and j[1] == k[5]
                 and j[2] == k[2] and j[3] == k[3]
                 and j[4] == k[0] and j[5] == k[1]
                 and j[6] == k[6] and j[7] == k[7]
                 and not (j[2] == j[0] and j[3] == j[1])
                 and not (j[2] == j[4] and j[3] == j[5]))
                or
                (j[0] == k[0] and j[1] == k[1]
                 and j[2] == k[6] and j[3] == k[7]
                 and j[4] == k[4] and j[5] == k[5]
                 and j[6] == k[2] and j[7] == k[3]
                 and not (j[4] == j[2] and j[5] == j[3])
                 and not (j[4] == j[6] and j[5] == j[7]))
            ):
                StdI.intr[jintr] = StdI.intr[jintr] - StdI.intr[kintr]
                StdI.intr[kintr] = 0.0

    # -----------------------------------------------------------------
    #  Pass 2: Force Hermitian ordering
    #  (c1+ c2 c3+ c4)+ = c4+ c3 c2+ c1
    # -----------------------------------------------------------------
    for jintr in range(StdI.nintr):
        j = StdI.intrindx[jintr]
        for kintr in range(jintr + 1, StdI.nintr):
            k = StdI.intrindx[kintr]

            # First reorder pattern (direct Hermitian conjugate match)
            if (j[6] == k[4] and j[7] == k[5]
                    and j[4] == k[6] and j[5] == k[7]
                    and j[2] == k[0] and j[3] == k[1]
                    and j[0] == k[2] and j[1] == k[3]
                    and not (k[0] == k[6] and k[1] == k[7])
                    and not (k[2] == k[4] and k[3] == k[5])):
                StdI.intrindx[kintr][0] = j[6]
                StdI.intrindx[kintr][1] = j[7]
                StdI.intrindx[kintr][2] = j[4]
                StdI.intrindx[kintr][3] = j[5]
                StdI.intrindx[kintr][4] = j[2]
                StdI.intrindx[kintr][5] = j[3]
                StdI.intrindx[kintr][6] = j[0]
                StdI.intrindx[kintr][7] = j[1]

            # Second reorder pattern (two sub-cases, with sign flip)
            elif (
                (j[6] == k[4] and j[7] == k[5]
                 and j[4] == k[2] and j[5] == k[3]
                 and j[2] == k[0] and j[3] == k[1]
                 and j[0] == k[6] and j[1] == k[7]
                 and not (k[2] == k[0] and k[3] == k[1])
                 and not (k[2] == k[4] and k[3] == k[5]))
                or
                (j[6] == k[0] and j[7] == k[1]
                 and j[4] == k[6] and j[5] == k[7]
                 and j[2] == k[4] and j[3] == k[5]
                 and j[0] == k[2] and j[1] == k[3]
                 and not (k[4] == k[2] and k[5] == k[3])
                 and not (k[4] == k[6] and k[5] == k[7]))
            ):
                StdI.intrindx[kintr][0] = j[6]
                StdI.intrindx[kintr][1] = j[7]
                StdI.intrindx[kintr][2] = j[4]
                StdI.intrindx[kintr][3] = j[5]
                StdI.intrindx[kintr][4] = j[2]
                StdI.intrindx[kintr][5] = j[3]
                StdI.intrindx[kintr][6] = j[0]
                StdI.intrindx[kintr][7] = j[1]

                StdI.intr[kintr] = -StdI.intr[kintr]

    # -----------------------------------------------------------------
    #  Pass 3: Remove diagonal terms
    #  If (site0,spin0)==(site4,spin4) or (site2,spin2)==(site6,spin6),
    #  zero the term unless any diagonal pair matches.
    # -----------------------------------------------------------------
    for jintr in range(StdI.nintr):
        idx = StdI.intrindx[jintr]

        has_diagonal_pair = (
            (idx[0] == idx[4] and idx[1] == idx[5])
            or (idx[2] == idx[6] and idx[3] == idx[7])
        )

        if has_diagonal_pair:
            any_matching_diagonal = (
                (idx[0] == idx[2] and idx[1] == idx[3])
                or (idx[0] == idx[6] and idx[1] == idx[7])
                or (idx[4] == idx[2] and idx[5] == idx[3])
                or (idx[4] == idx[6] and idx[5] == idx[7])
            )
            if not any_matching_diagonal:
                StdI.intr[jintr] = 0.0

    # -----------------------------------------------------------------
    #  Count non-zero InterAll terms and write file
    # -----------------------------------------------------------------
    nintr0 = 0
    for kintr in range(StdI.nintr):
        if abs(StdI.intr[kintr]) > 0.000001:
            nintr0 += 1

    if nintr0 == 0 or StdI.lBoost == 1:
        StdI.Lintr = 0
    else:
        StdI.Lintr = 1

    if StdI.Lintr == 1:
        with open("interall.def", "w") as fp:
            fp.write("====================== \n")
            fp.write(f"NInterAll {nintr0:7d}  \n")
            fp.write("====================== \n")
            fp.write("========zInterAll===== \n")
            fp.write("====================== \n")

            if StdI.lBoost == 0:
                for kintr in range(StdI.nintr):
                    if abs(StdI.intr[kintr]) > 0.000001:
                        idx = StdI.intrindx[kintr]
                        re_val = StdI.intr[kintr].real
                        im_val = StdI.intr[kintr].imag
                        fp.write(
                            f"{idx[0]:5d} {idx[1]:5d} "
                            f"{idx[2]:5d} {idx[3]:5d} "
                            f"{idx[4]:5d} {idx[5]:5d} "
                            f"{idx[6]:5d} {idx[7]:5d} "
                            f"{re_val:25.15f}  {im_val:25.15f}\n"
                        )

        print("    interall.def is written.")


# ===================================================================
#  Part 7 -- keyword parsing helpers and stdface_main() entry point
# ===================================================================


def _parse_common_keyword(keyword: str, value: str, StdI: StdIntList) -> bool:
    """Parse a keyword common to all solvers.

    Returns True if the keyword was recognised, False otherwise.
    """
    # --- scalars a -------------------------------------------------------
    if keyword == "a":
        StdI.a = _store_with_check_dup_d(keyword, value, StdI.a)
    # --- box (supercell) -------------------------------------------------
    elif keyword == "a0h":
        StdI.box[0, 2] = _store_with_check_dup_i(keyword, value, int(StdI.box[0, 2]))
    elif keyword == "a0l":
        StdI.box[0, 1] = _store_with_check_dup_i(keyword, value, int(StdI.box[0, 1]))
    elif keyword == "a0w":
        StdI.box[0, 0] = _store_with_check_dup_i(keyword, value, int(StdI.box[0, 0]))
    elif keyword == "a1h":
        StdI.box[1, 2] = _store_with_check_dup_i(keyword, value, int(StdI.box[1, 2]))
    elif keyword == "a1l":
        StdI.box[1, 1] = _store_with_check_dup_i(keyword, value, int(StdI.box[1, 1]))
    elif keyword == "a1w":
        StdI.box[1, 0] = _store_with_check_dup_i(keyword, value, int(StdI.box[1, 0]))
    elif keyword == "a2h":
        StdI.box[2, 2] = _store_with_check_dup_i(keyword, value, int(StdI.box[2, 2]))
    elif keyword == "a2l":
        StdI.box[2, 1] = _store_with_check_dup_i(keyword, value, int(StdI.box[2, 1]))
    elif keyword == "a2w":
        StdI.box[2, 0] = _store_with_check_dup_i(keyword, value, int(StdI.box[2, 0]))
    # --- cutoff J --------------------------------------------------------
    elif keyword == "cutoff_j":
        StdI.cutoff_j = _store_with_check_dup_d(keyword, value, StdI.cutoff_j)
    elif keyword == "cutoff_jh":
        StdI.cutoff_JR[2] = _store_with_check_dup_i(keyword, value, int(StdI.cutoff_JR[2]))
    elif keyword == "cutoff_jl":
        StdI.cutoff_JR[1] = _store_with_check_dup_i(keyword, value, int(StdI.cutoff_JR[1]))
    elif keyword == "cutoff_jw":
        StdI.cutoff_JR[0] = _store_with_check_dup_i(keyword, value, int(StdI.cutoff_JR[0]))
    elif keyword == "cutoff_j_a0w":
        StdI.cutoff_JVec[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_JVec[0, 0]))
    elif keyword == "cutoff_j_a0l":
        StdI.cutoff_JVec[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_JVec[0, 1]))
    elif keyword == "cutoff_j_a0h":
        StdI.cutoff_JVec[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_JVec[0, 2]))
    elif keyword == "cutoff_j_a1w":
        StdI.cutoff_JVec[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_JVec[1, 0]))
    elif keyword == "cutoff_j_a1l":
        StdI.cutoff_JVec[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_JVec[1, 1]))
    elif keyword == "cutoff_j_a1h":
        StdI.cutoff_JVec[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_JVec[1, 2]))
    elif keyword == "cutoff_j_a2w":
        StdI.cutoff_JVec[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_JVec[2, 0]))
    elif keyword == "cutoff_j_a2l":
        StdI.cutoff_JVec[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_JVec[2, 1]))
    elif keyword == "cutoff_j_a2h":
        StdI.cutoff_JVec[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_JVec[2, 2]))
    elif keyword == "cutoff_length_j":
        StdI.cutoff_length_J = _store_with_check_dup_d(keyword, value, StdI.cutoff_length_J)
    elif keyword == "cutoff_length_u":
        StdI.cutoff_length_U = _store_with_check_dup_d(keyword, value, StdI.cutoff_length_U)
    elif keyword == "cutoff_length_t":
        StdI.cutoff_length_t = _store_with_check_dup_d(keyword, value, StdI.cutoff_length_t)
    # --- cutoff t --------------------------------------------------------
    elif keyword == "cutoff_t":
        StdI.cutoff_t = _store_with_check_dup_d(keyword, value, StdI.cutoff_t)
    elif keyword == "cutoff_th":
        StdI.cutoff_tR[2] = _store_with_check_dup_i(keyword, value, int(StdI.cutoff_tR[2]))
    elif keyword == "cutoff_tl":
        StdI.cutoff_tR[1] = _store_with_check_dup_i(keyword, value, int(StdI.cutoff_tR[1]))
    elif keyword == "cutoff_tw":
        StdI.cutoff_tR[0] = _store_with_check_dup_i(keyword, value, int(StdI.cutoff_tR[0]))
    elif keyword == "cutoff_t_a0w":
        StdI.cutoff_tVec[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_tVec[0, 0]))
    elif keyword == "cutoff_t_a0l":
        StdI.cutoff_tVec[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_tVec[0, 1]))
    elif keyword == "cutoff_t_a0h":
        StdI.cutoff_tVec[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_tVec[0, 2]))
    elif keyword == "cutoff_t_a1w":
        StdI.cutoff_tVec[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_tVec[1, 0]))
    elif keyword == "cutoff_t_a1l":
        StdI.cutoff_tVec[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_tVec[1, 1]))
    elif keyword == "cutoff_t_a1h":
        StdI.cutoff_tVec[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_tVec[1, 2]))
    elif keyword == "cutoff_t_a2w":
        StdI.cutoff_tVec[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_tVec[2, 0]))
    elif keyword == "cutoff_t_a2l":
        StdI.cutoff_tVec[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_tVec[2, 1]))
    elif keyword == "cutoff_t_a2h":
        StdI.cutoff_tVec[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_tVec[2, 2]))
    # --- cutoff U --------------------------------------------------------
    elif keyword == "cutoff_u":
        StdI.cutoff_u = _store_with_check_dup_d(keyword, value, StdI.cutoff_u)
    elif keyword == "cutoff_uh":
        StdI.cutoff_UR[2] = _store_with_check_dup_i(keyword, value, int(StdI.cutoff_UR[2]))
    elif keyword == "cutoff_ul":
        StdI.cutoff_UR[1] = _store_with_check_dup_i(keyword, value, int(StdI.cutoff_UR[1]))
    elif keyword == "cutoff_uw":
        StdI.cutoff_UR[0] = _store_with_check_dup_i(keyword, value, int(StdI.cutoff_UR[0]))
    elif keyword == "cutoff_u_a0w":
        StdI.cutoff_UVec[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_UVec[0, 0]))
    elif keyword == "cutoff_u_a0l":
        StdI.cutoff_UVec[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_UVec[0, 1]))
    elif keyword == "cutoff_u_a0h":
        StdI.cutoff_UVec[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_UVec[0, 2]))
    elif keyword == "cutoff_u_a1w":
        StdI.cutoff_UVec[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_UVec[1, 0]))
    elif keyword == "cutoff_u_a1l":
        StdI.cutoff_UVec[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_UVec[1, 1]))
    elif keyword == "cutoff_u_a1h":
        StdI.cutoff_UVec[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_UVec[1, 2]))
    elif keyword == "cutoff_u_a2w":
        StdI.cutoff_UVec[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_UVec[2, 0]))
    elif keyword == "cutoff_u_a2l":
        StdI.cutoff_UVec[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_UVec[2, 1]))
    elif keyword == "cutoff_u_a2h":
        StdI.cutoff_UVec[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.cutoff_UVec[2, 2]))
    # --- lambda, alpha, D ------------------------------------------------
    elif keyword == "lambda":
        StdI.lambda_ = _store_with_check_dup_d(keyword, value, StdI.lambda_)
    elif keyword == "lambda_u":
        StdI.lambda_U = _store_with_check_dup_d(keyword, value, StdI.lambda_U)
    elif keyword == "lambda_j":
        StdI.lambda_J = _store_with_check_dup_d(keyword, value, StdI.lambda_J)
    elif keyword == "alpha":
        StdI.alpha = _store_with_check_dup_d(keyword, value, StdI.alpha)
    elif keyword == "d":
        StdI.D[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.D[2, 2]))
    elif keyword == "doublecounting":
        StdI.double_counting_mode = _store_with_check_dup_sl(keyword, value, StdI.double_counting_mode)
    # --- magnetic field and related ---------------------------------------
    elif keyword == "gamma":
        StdI.Gamma = _store_with_check_dup_d(keyword, value, StdI.Gamma)
    elif keyword == "h":
        StdI.h = _store_with_check_dup_d(keyword, value, StdI.h)
    elif keyword == "gamma_y":
        StdI.Gamma_y = _store_with_check_dup_d(keyword, value, StdI.Gamma_y)
    elif keyword == "height":
        StdI.Height = _store_with_check_dup_i(keyword, value, StdI.Height)
    elif keyword == "hlength":
        StdI.length[2] = _store_with_check_dup_d(keyword, value, float(StdI.length[2]))
    elif keyword == "hx":
        StdI.direct[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.direct[2, 0]))
    elif keyword == "hy":
        StdI.direct[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.direct[2, 1]))
    elif keyword == "hz":
        StdI.direct[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.direct[2, 2]))
    # --- J exchange couplings --------------------------------------------
    elif keyword == "j":
        StdI.JAll = _store_with_check_dup_d(keyword, value, StdI.JAll)
    elif keyword == "jx":
        StdI.J[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J[0, 0]))
    elif keyword == "jxy":
        StdI.J[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J[0, 1]))
    elif keyword == "jxz":
        StdI.J[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J[0, 2]))
    elif keyword == "jy":
        StdI.J[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J[1, 1]))
    elif keyword == "jyx":
        StdI.J[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J[1, 0]))
    elif keyword == "jyz":
        StdI.J[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J[1, 2]))
    elif keyword == "jz":
        StdI.J[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J[2, 2]))
    elif keyword == "jzx":
        StdI.J[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J[2, 0]))
    elif keyword == "jzy":
        StdI.J[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J[2, 1]))
    # --- J0 --------------------------------------------------------------
    elif keyword == "j0":
        StdI.J0All = _store_with_check_dup_d(keyword, value, StdI.J0All)
    elif keyword == "j0x":
        StdI.J0[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J0[0, 0]))
    elif keyword == "j0xy":
        StdI.J0[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J0[0, 1]))
    elif keyword == "j0xz":
        StdI.J0[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J0[0, 2]))
    elif keyword == "j0y":
        StdI.J0[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J0[1, 1]))
    elif keyword == "j0yx":
        StdI.J0[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J0[1, 0]))
    elif keyword == "j0yz":
        StdI.J0[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J0[1, 2]))
    elif keyword == "j0z":
        StdI.J0[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J0[2, 2]))
    elif keyword == "j0zx":
        StdI.J0[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J0[2, 0]))
    elif keyword == "j0zy":
        StdI.J0[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J0[2, 1]))
    # --- J0' -------------------------------------------------------------
    elif keyword == "j0'":
        StdI.J0pAll = _store_with_check_dup_d(keyword, value, StdI.J0pAll)
    elif keyword == "j0'x":
        StdI.J0p[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J0p[0, 0]))
    elif keyword == "j0'xy":
        StdI.J0p[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J0p[0, 1]))
    elif keyword == "j0'xz":
        StdI.J0p[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J0p[0, 2]))
    elif keyword == "j0'y":
        StdI.J0p[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J0p[1, 1]))
    elif keyword == "j0'yx":
        StdI.J0p[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J0p[1, 0]))
    elif keyword == "j0'yz":
        StdI.J0p[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J0p[1, 2]))
    elif keyword == "j0'z":
        StdI.J0p[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J0p[2, 2]))
    elif keyword == "j0'zx":
        StdI.J0p[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J0p[2, 0]))
    elif keyword == "j0'zy":
        StdI.J0p[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J0p[2, 1]))
    # --- J0'' ------------------------------------------------------------
    elif keyword == "j0''":
        StdI.J0ppAll = _store_with_check_dup_d(keyword, value, StdI.J0ppAll)
    elif keyword == "j0''x":
        StdI.J0pp[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J0pp[0, 0]))
    elif keyword == "j0''xy":
        StdI.J0pp[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J0pp[0, 1]))
    elif keyword == "j0''xz":
        StdI.J0pp[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J0pp[0, 2]))
    elif keyword == "j0''y":
        StdI.J0pp[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J0pp[1, 1]))
    elif keyword == "j0''yx":
        StdI.J0pp[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J0pp[1, 0]))
    elif keyword == "j0''yz":
        StdI.J0pp[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J0pp[1, 2]))
    elif keyword == "j0''z":
        StdI.J0pp[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J0pp[2, 2]))
    elif keyword == "j0''zx":
        StdI.J0pp[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J0pp[2, 0]))
    elif keyword == "j0''zy":
        StdI.J0pp[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J0pp[2, 1]))
    # --- J1 --------------------------------------------------------------
    elif keyword == "j1":
        StdI.J1All = _store_with_check_dup_d(keyword, value, StdI.J1All)
    elif keyword == "j1x":
        StdI.J1[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J1[0, 0]))
    elif keyword == "j1xy":
        StdI.J1[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J1[0, 1]))
    elif keyword == "j1xz":
        StdI.J1[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J1[0, 2]))
    elif keyword == "j1y":
        StdI.J1[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J1[1, 1]))
    elif keyword == "j1yx":
        StdI.J1[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J1[1, 0]))
    elif keyword == "j1yz":
        StdI.J1[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J1[1, 2]))
    elif keyword == "j1z":
        StdI.J1[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J1[2, 2]))
    elif keyword == "j1zx":
        StdI.J1[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J1[2, 0]))
    elif keyword == "j1zy":
        StdI.J1[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J1[2, 1]))
    # --- J1' -------------------------------------------------------------
    elif keyword == "j1'":
        StdI.J1pAll = _store_with_check_dup_d(keyword, value, StdI.J1pAll)
    elif keyword == "j1'x":
        StdI.J1p[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J1p[0, 0]))
    elif keyword == "j1'xy":
        StdI.J1p[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J1p[0, 1]))
    elif keyword == "j1'xz":
        StdI.J1p[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J1p[0, 2]))
    elif keyword == "j1'y":
        StdI.J1p[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J1p[1, 1]))
    elif keyword == "j1'yx":
        StdI.J1p[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J1p[1, 0]))
    elif keyword == "j1'yz":
        StdI.J1p[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J1p[1, 2]))
    elif keyword == "j1'z":
        StdI.J1p[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J1p[2, 2]))
    elif keyword == "j1'zx":
        StdI.J1p[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J1p[2, 0]))
    elif keyword == "j1'zy":
        StdI.J1p[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J1p[2, 1]))
    # --- J1'' ------------------------------------------------------------
    elif keyword == "j1''":
        StdI.J1ppAll = _store_with_check_dup_d(keyword, value, StdI.J1ppAll)
    elif keyword == "j1''x":
        StdI.J1pp[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J1pp[0, 0]))
    elif keyword == "j1''xy":
        StdI.J1pp[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J1pp[0, 1]))
    elif keyword == "j1''xz":
        StdI.J1pp[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J1pp[0, 2]))
    elif keyword == "j1''y":
        StdI.J1pp[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J1pp[1, 1]))
    elif keyword == "j1''yx":
        StdI.J1pp[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J1pp[1, 0]))
    elif keyword == "j1''yz":
        StdI.J1pp[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J1pp[1, 2]))
    elif keyword == "j1''z":
        StdI.J1pp[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J1pp[2, 2]))
    elif keyword == "j1''zx":
        StdI.J1pp[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J1pp[2, 0]))
    elif keyword == "j1''zy":
        StdI.J1pp[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J1pp[2, 1]))
    # --- J2 --------------------------------------------------------------
    elif keyword == "j2":
        StdI.J2All = _store_with_check_dup_d(keyword, value, StdI.J2All)
    elif keyword == "j2x":
        StdI.J2[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J2[0, 0]))
    elif keyword == "j2xy":
        StdI.J2[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J2[0, 1]))
    elif keyword == "j2xz":
        StdI.J2[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J2[0, 2]))
    elif keyword == "j2y":
        StdI.J2[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J2[1, 1]))
    elif keyword == "j2yx":
        StdI.J2[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J2[1, 0]))
    elif keyword == "j2yz":
        StdI.J2[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J2[1, 2]))
    elif keyword == "j2z":
        StdI.J2[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J2[2, 2]))
    elif keyword == "j2zx":
        StdI.J2[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J2[2, 0]))
    elif keyword == "j2zy":
        StdI.J2[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J2[2, 1]))
    # --- J2' -------------------------------------------------------------
    elif keyword == "j2'":
        StdI.J2pAll = _store_with_check_dup_d(keyword, value, StdI.J2pAll)
    elif keyword == "j2'x":
        StdI.J2p[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J2p[0, 0]))
    elif keyword == "j2'xy":
        StdI.J2p[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J2p[0, 1]))
    elif keyword == "j2'xz":
        StdI.J2p[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J2p[0, 2]))
    elif keyword == "j2'y":
        StdI.J2p[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J2p[1, 1]))
    elif keyword == "j2'yx":
        StdI.J2p[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J2p[1, 0]))
    elif keyword == "j2'yz":
        StdI.J2p[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J2p[1, 2]))
    elif keyword == "j2'z":
        StdI.J2p[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J2p[2, 2]))
    elif keyword == "j2'zx":
        StdI.J2p[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J2p[2, 0]))
    elif keyword == "j2'zy":
        StdI.J2p[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J2p[2, 1]))
    # --- J2'' ------------------------------------------------------------
    elif keyword == "j2''":
        StdI.J2ppAll = _store_with_check_dup_d(keyword, value, StdI.J2ppAll)
    elif keyword == "j2''x":
        StdI.J2pp[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J2pp[0, 0]))
    elif keyword == "j2''xy":
        StdI.J2pp[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J2pp[0, 1]))
    elif keyword == "j2''xz":
        StdI.J2pp[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J2pp[0, 2]))
    elif keyword == "j2''y":
        StdI.J2pp[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J2pp[1, 1]))
    elif keyword == "j2''yx":
        StdI.J2pp[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J2pp[1, 0]))
    elif keyword == "j2''yz":
        StdI.J2pp[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J2pp[1, 2]))
    elif keyword == "j2''z":
        StdI.J2pp[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.J2pp[2, 2]))
    elif keyword == "j2''zx":
        StdI.J2pp[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.J2pp[2, 0]))
    elif keyword == "j2''zy":
        StdI.J2pp[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.J2pp[2, 1]))
    # --- J' (nearest-neighbour prime) ------------------------------------
    elif keyword == "j'":
        StdI.JpAll = _store_with_check_dup_d(keyword, value, StdI.JpAll)
    elif keyword == "j'x":
        StdI.Jp[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.Jp[0, 0]))
    elif keyword == "j'xy":
        StdI.Jp[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.Jp[0, 1]))
    elif keyword == "j'xz":
        StdI.Jp[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.Jp[0, 2]))
    elif keyword == "j'y":
        StdI.Jp[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.Jp[1, 1]))
    elif keyword == "j'yx":
        StdI.Jp[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.Jp[1, 0]))
    elif keyword == "j'yz":
        StdI.Jp[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.Jp[1, 2]))
    elif keyword == "j'z":
        StdI.Jp[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.Jp[2, 2]))
    elif keyword == "j'zx":
        StdI.Jp[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.Jp[2, 0]))
    elif keyword == "j'zy":
        StdI.Jp[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.Jp[2, 1]))
    # --- J'' (next-nearest-neighbour double-prime) -----------------------
    elif keyword == "j''":
        StdI.JppAll = _store_with_check_dup_d(keyword, value, StdI.JppAll)
    elif keyword == "j''x":
        StdI.Jpp[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.Jpp[0, 0]))
    elif keyword == "j''xy":
        StdI.Jpp[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.Jpp[0, 1]))
    elif keyword == "j''xz":
        StdI.Jpp[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.Jpp[0, 2]))
    elif keyword == "j''y":
        StdI.Jpp[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.Jpp[1, 1]))
    elif keyword == "j''yx":
        StdI.Jpp[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.Jpp[1, 0]))
    elif keyword == "j''yz":
        StdI.Jpp[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.Jpp[1, 2]))
    elif keyword == "j''z":
        StdI.Jpp[2, 2] = _store_with_check_dup_d(keyword, value, float(StdI.Jpp[2, 2]))
    elif keyword == "j''zx":
        StdI.Jpp[2, 0] = _store_with_check_dup_d(keyword, value, float(StdI.Jpp[2, 0]))
    elif keyword == "j''zy":
        StdI.Jpp[2, 1] = _store_with_check_dup_d(keyword, value, float(StdI.Jpp[2, 1]))
    # --- K, L, lattice, length -------------------------------------------
    elif keyword == "k":
        StdI.K = _store_with_check_dup_d(keyword, value, StdI.K)
    elif keyword == "l":
        StdI.L = _store_with_check_dup_i(keyword, value, StdI.L)
    elif keyword == "lattice":
        StdI.lattice = _store_with_check_dup_sl(keyword, value, StdI.lattice)
    elif keyword == "llength":
        StdI.length[1] = _store_with_check_dup_d(keyword, value, float(StdI.length[1]))
    elif keyword == "lx":
        StdI.direct[1, 0] = _store_with_check_dup_d(keyword, value, float(StdI.direct[1, 0]))
    elif keyword == "ly":
        StdI.direct[1, 1] = _store_with_check_dup_d(keyword, value, float(StdI.direct[1, 1]))
    elif keyword == "lz":
        StdI.direct[1, 2] = _store_with_check_dup_d(keyword, value, float(StdI.direct[1, 2]))
    # --- model, mu, ncond, outputmode ------------------------------------
    elif keyword == "model":
        StdI.model = _store_with_check_dup_sl(keyword, value, StdI.model)
    elif keyword == "mu":
        StdI.mu = _store_with_check_dup_d(keyword, value, StdI.mu)
    elif keyword == "ncond":
        StdI.ncond = _store_with_check_dup_i(keyword, value, StdI.ncond)
    elif keyword == "nelec":
        StdI.ncond = _store_with_check_dup_i(keyword, value, StdI.ncond)
    elif keyword == "outputmode":
        StdI.outputmode = _store_with_check_dup_sl(keyword, value, StdI.outputmode)
    # --- phase -----------------------------------------------------------
    elif keyword == "phase0":
        StdI.phase[0] = _store_with_check_dup_d(keyword, value, float(StdI.phase[0]))
    elif keyword == "phase1":
        StdI.phase[1] = _store_with_check_dup_d(keyword, value, float(StdI.phase[1]))
    elif keyword == "phase2":
        StdI.phase[2] = _store_with_check_dup_d(keyword, value, float(StdI.phase[2]))
    # --- hopping t -------------------------------------------------------
    elif keyword == "t":
        StdI.t = _store_with_check_dup_c(keyword, value, StdI.t)
    elif keyword == "t0":
        StdI.t0 = _store_with_check_dup_c(keyword, value, StdI.t0)
    elif keyword == "t0'":
        StdI.t0p = _store_with_check_dup_c(keyword, value, StdI.t0p)
    elif keyword == "t0''":
        StdI.t0pp = _store_with_check_dup_c(keyword, value, StdI.t0pp)
    elif keyword == "t1":
        StdI.t1 = _store_with_check_dup_c(keyword, value, StdI.t1)
    elif keyword == "t1'":
        StdI.t1p = _store_with_check_dup_c(keyword, value, StdI.t1p)
    elif keyword == "t1''":
        StdI.t1pp = _store_with_check_dup_c(keyword, value, StdI.t1pp)
    elif keyword == "t2":
        StdI.t2 = _store_with_check_dup_c(keyword, value, StdI.t2)
    elif keyword == "t2'":
        StdI.t2p = _store_with_check_dup_c(keyword, value, StdI.t2p)
    elif keyword == "t2''":
        StdI.t2pp = _store_with_check_dup_c(keyword, value, StdI.t2pp)
    elif keyword == "t'":
        StdI.tp = _store_with_check_dup_c(keyword, value, StdI.tp)
    elif keyword == "t''":
        StdI.tpp = _store_with_check_dup_c(keyword, value, StdI.tpp)
    # --- U, V ------------------------------------------------------------
    elif keyword == "u":
        StdI.U = _store_with_check_dup_d(keyword, value, StdI.U)
    elif keyword == "v":
        StdI.V = _store_with_check_dup_d(keyword, value, StdI.V)
    elif keyword == "v0":
        StdI.V0 = _store_with_check_dup_d(keyword, value, StdI.V0)
    elif keyword == "v0'":
        StdI.V0p = _store_with_check_dup_d(keyword, value, StdI.V0p)
    elif keyword == "v0''":
        StdI.V0pp = _store_with_check_dup_d(keyword, value, StdI.V0pp)
    elif keyword == "v1":
        StdI.V1 = _store_with_check_dup_d(keyword, value, StdI.V1)
    elif keyword == "v1'":
        StdI.V1p = _store_with_check_dup_d(keyword, value, StdI.V1p)
    elif keyword == "v1''":
        StdI.V1pp = _store_with_check_dup_d(keyword, value, StdI.V1pp)
    elif keyword == "v2":
        StdI.V2 = _store_with_check_dup_d(keyword, value, StdI.V2)
    elif keyword == "v2'":
        StdI.V2p = _store_with_check_dup_d(keyword, value, StdI.V2p)
    elif keyword == "v2''":
        StdI.V2pp = _store_with_check_dup_d(keyword, value, StdI.V2pp)
    elif keyword == "v'":
        StdI.Vp = _store_with_check_dup_d(keyword, value, StdI.Vp)
    elif keyword == "v''":
        StdI.Vpp = _store_with_check_dup_d(keyword, value, StdI.Vpp)
    # --- W, wlength, wx/wy/wz -------------------------------------------
    elif keyword == "w":
        StdI.W = _store_with_check_dup_i(keyword, value, StdI.W)
    elif keyword == "wlength":
        StdI.length[0] = _store_with_check_dup_d(keyword, value, float(StdI.length[0]))
    elif keyword == "wx":
        StdI.direct[0, 0] = _store_with_check_dup_d(keyword, value, float(StdI.direct[0, 0]))
    elif keyword == "wy":
        StdI.direct[0, 1] = _store_with_check_dup_d(keyword, value, float(StdI.direct[0, 1]))
    elif keyword == "wz":
        StdI.direct[0, 2] = _store_with_check_dup_d(keyword, value, float(StdI.direct[0, 2]))
    # --- 2Sz -------------------------------------------------------------
    elif keyword == "2sz":
        StdI.Sz2 = _store_with_check_dup_i(keyword, value, StdI.Sz2)
    else:
        return False
    return True


def _parse_solver_keyword(keyword: str, value: str, StdI: StdIntList,
                          solver: str) -> bool:
    """Parse a solver-specific keyword.

    Returns True if the keyword was recognised, False otherwise.
    """
    if solver == "HPhi":
        return _parse_hphi_keyword(keyword, value, StdI)
    elif solver == "mVMC":
        return _parse_mvmc_keyword(keyword, value, StdI)
    elif solver == "UHF":
        return _parse_uhf_keyword(keyword, value, StdI)
    elif solver == "HWAVE":
        return _parse_hwave_keyword(keyword, value, StdI)
    return False


def _parse_hphi_keyword(keyword: str, value: str, StdI: StdIntList) -> bool:
    """Parse HPhi-specific keywords."""
    if keyword == "calcspec":
        StdI.CalcSpec = _store_with_check_dup_sl(keyword, value, StdI.CalcSpec)
    elif keyword == "exct":
        StdI.exct = _store_with_check_dup_i(keyword, value, StdI.exct)
    elif keyword == "eigenvecio":
        StdI.EigenVecIO = _store_with_check_dup_sl(keyword, value, StdI.EigenVecIO)
    elif keyword == "expandcoef":
        StdI.ExpandCoef = _store_with_check_dup_i(keyword, value, StdI.ExpandCoef)
    elif keyword == "expecinterval":
        StdI.ExpecInterval = _store_with_check_dup_i(keyword, value, StdI.ExpecInterval)
    elif keyword == "cdatafilehead":
        StdI.CDataFileHead = _store_with_check_dup_s(keyword, value, StdI.CDataFileHead)
    elif keyword == "dt":
        StdI.dt = _store_with_check_dup_d(keyword, value, StdI.dt)
    elif keyword == "flgtemp":
        StdI.FlgTemp = _store_with_check_dup_i(keyword, value, StdI.FlgTemp)
    elif keyword == "freq":
        StdI.freq = _store_with_check_dup_d(keyword, value, StdI.freq)
    elif keyword == "hamio":
        StdI.HamIO = _store_with_check_dup_sl(keyword, value, StdI.HamIO)
    elif keyword == "initialvectype":
        StdI.InitialVecType = _store_with_check_dup_sl(keyword, value, StdI.InitialVecType)
    elif keyword == "initial_iv":
        StdI.initial_iv = _store_with_check_dup_i(keyword, value, StdI.initial_iv)
    elif keyword == "lanczoseps":
        StdI.LanczosEps = _store_with_check_dup_i(keyword, value, StdI.LanczosEps)
    elif keyword == "lanczostarget":
        StdI.LanczosTarget = _store_with_check_dup_i(keyword, value, StdI.LanczosTarget)
    elif keyword == "lanczos_max":
        StdI.Lanczos_max = _store_with_check_dup_i(keyword, value, StdI.Lanczos_max)
    elif keyword == "largevalue":
        StdI.LargeValue = _store_with_check_dup_d(keyword, value, StdI.LargeValue)
    elif keyword == "method":
        StdI.method = _store_with_check_dup_sl(keyword, value, StdI.method)
    elif keyword == "nomega":
        StdI.Nomega = _store_with_check_dup_i(keyword, value, StdI.Nomega)
    elif keyword == "numave":
        StdI.NumAve = _store_with_check_dup_i(keyword, value, StdI.NumAve)
    elif keyword == "nvec":
        StdI.nvec = _store_with_check_dup_i(keyword, value, StdI.nvec)
    elif keyword == "omegamax":
        StdI.OmegaMax = _store_with_check_dup_d(keyword, value, StdI.OmegaMax)
    elif keyword == "omegamin":
        StdI.OmegaMin = _store_with_check_dup_d(keyword, value, StdI.OmegaMin)
    elif keyword == "omegaorg":
        StdI.OmegaOrg = _store_with_check_dup_d(keyword, value, StdI.OmegaOrg)
    elif keyword == "omegaim":
        StdI.OmegaIm = _store_with_check_dup_d(keyword, value, StdI.OmegaIm)
    elif keyword == "outputexcitedvec":
        StdI.OutputExVec = _store_with_check_dup_sl(keyword, value, StdI.OutputExVec)
    elif keyword == "pumptype":
        StdI.PumpType = _store_with_check_dup_sl(keyword, value, StdI.PumpType)
    elif keyword == "restart":
        StdI.Restart = _store_with_check_dup_sl(keyword, value, StdI.Restart)
    elif keyword == "spectrumqh":
        StdI.SpectrumQ[2] = _store_with_check_dup_d(keyword, value, float(StdI.SpectrumQ[2]))
    elif keyword == "spectrumql":
        StdI.SpectrumQ[1] = _store_with_check_dup_d(keyword, value, float(StdI.SpectrumQ[1]))
    elif keyword == "spectrumqw":
        StdI.SpectrumQ[0] = _store_with_check_dup_d(keyword, value, float(StdI.SpectrumQ[0]))
    elif keyword == "spectrumtype":
        StdI.SpectrumType = _store_with_check_dup_sl(keyword, value, StdI.SpectrumType)
    elif keyword == "tdump":
        StdI.tdump = _store_with_check_dup_d(keyword, value, StdI.tdump)
    elif keyword == "tshift":
        StdI.tshift = _store_with_check_dup_d(keyword, value, StdI.tshift)
    elif keyword == "uquench":
        StdI.Uquench = _store_with_check_dup_d(keyword, value, StdI.Uquench)
    elif keyword == "vecpoth":
        StdI.VecPot[2] = _store_with_check_dup_d(keyword, value, float(StdI.VecPot[2]))
    elif keyword == "vecpotl":
        StdI.VecPot[1] = _store_with_check_dup_d(keyword, value, float(StdI.VecPot[1]))
    elif keyword == "vecpotw":
        StdI.VecPot[0] = _store_with_check_dup_d(keyword, value, float(StdI.VecPot[0]))
    elif keyword == "2s":
        StdI.S2 = _store_with_check_dup_i(keyword, value, StdI.S2)
    elif keyword == "ngpu":
        StdI.NGPU = _store_with_check_dup_i(keyword, value, StdI.NGPU)
    elif keyword == "scalapack":
        StdI.Scalapack = _store_with_check_dup_i(keyword, value, StdI.Scalapack)
    else:
        return False
    return True


def _parse_mvmc_keyword(keyword: str, value: str, StdI: StdIntList) -> bool:
    """Parse mVMC-specific keywords."""
    if keyword == "a0hsub":
        StdI.boxsub[0, 2] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[0, 2]))
    elif keyword == "a0lsub":
        StdI.boxsub[0, 1] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[0, 1]))
    elif keyword == "a0wsub":
        StdI.boxsub[0, 0] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[0, 0]))
    elif keyword == "a1hsub":
        StdI.boxsub[1, 2] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[1, 2]))
    elif keyword == "a1lsub":
        StdI.boxsub[1, 1] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[1, 1]))
    elif keyword == "a1wsub":
        StdI.boxsub[1, 0] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[1, 0]))
    elif keyword == "a2hsub":
        StdI.boxsub[2, 2] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[2, 2]))
    elif keyword == "a2lsub":
        StdI.boxsub[2, 1] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[2, 1]))
    elif keyword == "a2wsub":
        StdI.boxsub[2, 0] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[2, 0]))
    elif keyword == "complextype":
        StdI.ComplexType = _store_with_check_dup_i(keyword, value, StdI.ComplexType)
    elif keyword == "cparafilehead":
        StdI.CParaFileHead = _store_with_check_dup_s(keyword, value, StdI.CParaFileHead)
    elif keyword == "dsroptredcut":
        StdI.DSROptRedCut = _store_with_check_dup_d(keyword, value, StdI.DSROptRedCut)
    elif keyword == "dsroptstadel":
        StdI.DSROptStaDel = _store_with_check_dup_d(keyword, value, StdI.DSROptStaDel)
    elif keyword == "dsroptstepdt":
        StdI.DSROptStepDt = _store_with_check_dup_d(keyword, value, StdI.DSROptStepDt)
    elif keyword == "hsub":
        StdI.Hsub = _store_with_check_dup_i(keyword, value, StdI.Hsub)
    elif keyword == "lsub":
        StdI.Lsub = _store_with_check_dup_i(keyword, value, StdI.Lsub)
    elif keyword == "nvmccalmode":
        StdI.NVMCCalMode = _store_with_check_dup_i(keyword, value, StdI.NVMCCalMode)
    elif keyword == "ndataidxstart":
        StdI.NDataIdxStart = _store_with_check_dup_i(keyword, value, StdI.NDataIdxStart)
    elif keyword == "ndataqtysmp":
        StdI.NDataQtySmp = _store_with_check_dup_i(keyword, value, StdI.NDataQtySmp)
    elif keyword == "nlanczosmode":
        StdI.NLanczosMode = _store_with_check_dup_i(keyword, value, StdI.NLanczosMode)
    elif keyword == "nmptrans":
        StdI.NMPTrans = _store_with_check_dup_i(keyword, value, StdI.NMPTrans)
    elif keyword == "nspgaussleg":
        StdI.NSPGaussLeg = _store_with_check_dup_i(keyword, value, StdI.NSPGaussLeg)
    elif keyword == "nsplitsize":
        StdI.NSplitSize = _store_with_check_dup_i(keyword, value, StdI.NSplitSize)
    elif keyword == "nspstot":
        StdI.NSPStot = _store_with_check_dup_i(keyword, value, StdI.NSPStot)
    elif keyword == "nsroptitrsmp":
        StdI.NSROptItrSmp = _store_with_check_dup_i(keyword, value, StdI.NSROptItrSmp)
    elif keyword == "nsroptitrstep":
        StdI.NSROptItrStep = _store_with_check_dup_i(keyword, value, StdI.NSROptItrStep)
    elif keyword == "nstore":
        StdI.NStore = _store_with_check_dup_i(keyword, value, StdI.NStore)
    elif keyword == "nsrcg":
        StdI.NSRCG = _store_with_check_dup_i(keyword, value, StdI.NSRCG)
    elif keyword == "nvmcinterval":
        StdI.NVMCInterval = _store_with_check_dup_i(keyword, value, StdI.NVMCInterval)
    elif keyword == "nvmcsample":
        StdI.NVMCSample = _store_with_check_dup_i(keyword, value, StdI.NVMCSample)
    elif keyword == "nvmcwarmup":
        StdI.NVMCWarmUp = _store_with_check_dup_i(keyword, value, StdI.NVMCWarmUp)
    elif keyword == "rndseed":
        StdI.RndSeed = _store_with_check_dup_i(keyword, value, StdI.RndSeed)
    elif keyword == "wsub":
        StdI.Wsub = _store_with_check_dup_i(keyword, value, StdI.Wsub)
    else:
        return False
    return True


def _parse_uhf_keyword(keyword: str, value: str, StdI: StdIntList) -> bool:
    """Parse UHF-specific keywords."""
    if keyword == "iteration_max":
        StdI.Iteration_max = _store_with_check_dup_i(keyword, value, StdI.Iteration_max)
    elif keyword == "rndseed":
        StdI.RndSeed = _store_with_check_dup_i(keyword, value, StdI.RndSeed)
    elif keyword == "nmptrans":
        StdI.NMPTrans = _store_with_check_dup_i(keyword, value, StdI.NMPTrans)
    elif keyword == "a0hsub":
        StdI.boxsub[0, 2] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[0, 2]))
    elif keyword == "a0lsub":
        StdI.boxsub[0, 1] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[0, 1]))
    elif keyword == "a0wsub":
        StdI.boxsub[0, 0] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[0, 0]))
    elif keyword == "a1hsub":
        StdI.boxsub[1, 2] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[1, 2]))
    elif keyword == "a1lsub":
        StdI.boxsub[1, 1] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[1, 1]))
    elif keyword == "a1wsub":
        StdI.boxsub[1, 0] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[1, 0]))
    elif keyword == "a2hsub":
        StdI.boxsub[2, 2] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[2, 2]))
    elif keyword == "a2lsub":
        StdI.boxsub[2, 1] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[2, 1]))
    elif keyword == "a2wsub":
        StdI.boxsub[2, 0] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[2, 0]))
    elif keyword == "hsub":
        StdI.Hsub = _store_with_check_dup_i(keyword, value, StdI.Hsub)
    elif keyword == "lsub":
        StdI.Lsub = _store_with_check_dup_i(keyword, value, StdI.Lsub)
    elif keyword == "wsub":
        StdI.Wsub = _store_with_check_dup_i(keyword, value, StdI.Wsub)
    elif keyword == "eps":
        StdI.eps = _store_with_check_dup_i(keyword, value, StdI.eps)
    elif keyword == "epsslater":
        StdI.eps_slater = _store_with_check_dup_i(keyword, value, StdI.eps_slater)
    elif keyword == "mix":
        StdI.mix = _store_with_check_dup_d(keyword, value, StdI.mix)
    else:
        return False
    return True


def _parse_hwave_keyword(keyword: str, value: str, StdI: StdIntList) -> bool:
    """Parse HWAVE-specific keywords."""
    if keyword == "iteration_max":
        StdI.Iteration_max = _store_with_check_dup_i(keyword, value, StdI.Iteration_max)
    elif keyword == "rndseed":
        StdI.RndSeed = _store_with_check_dup_i(keyword, value, StdI.RndSeed)
    elif keyword == "nmptrans":
        StdI.NMPTrans = _store_with_check_dup_i(keyword, value, StdI.NMPTrans)
    elif keyword == "a0hsub":
        StdI.boxsub[0, 2] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[0, 2]))
    elif keyword == "a0lsub":
        StdI.boxsub[0, 1] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[0, 1]))
    elif keyword == "a0wsub":
        StdI.boxsub[0, 0] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[0, 0]))
    elif keyword == "a1hsub":
        StdI.boxsub[1, 2] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[1, 2]))
    elif keyword == "a1lsub":
        StdI.boxsub[1, 1] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[1, 1]))
    elif keyword == "a1wsub":
        StdI.boxsub[1, 0] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[1, 0]))
    elif keyword == "a2hsub":
        StdI.boxsub[2, 2] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[2, 2]))
    elif keyword == "a2lsub":
        StdI.boxsub[2, 1] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[2, 1]))
    elif keyword == "a2wsub":
        StdI.boxsub[2, 0] = _store_with_check_dup_i(keyword, value, int(StdI.boxsub[2, 0]))
    elif keyword == "hsub":
        StdI.Hsub = _store_with_check_dup_i(keyword, value, StdI.Hsub)
    elif keyword == "lsub":
        StdI.Lsub = _store_with_check_dup_i(keyword, value, StdI.Lsub)
    elif keyword == "wsub":
        StdI.Wsub = _store_with_check_dup_i(keyword, value, StdI.Wsub)
    elif keyword == "eps":
        StdI.eps = _store_with_check_dup_i(keyword, value, StdI.eps)
    elif keyword == "epsslater":
        StdI.eps_slater = _store_with_check_dup_i(keyword, value, StdI.eps_slater)
    elif keyword == "mix":
        StdI.mix = _store_with_check_dup_d(keyword, value, StdI.mix)
    elif keyword == "calcmode":
        StdI.calcmode = _store_with_check_dup_sl(keyword, value, StdI.calcmode)
    elif keyword == "fileprefix":
        StdI.fileprefix = _store_with_check_dup_sl(keyword, value, StdI.fileprefix)
    elif keyword == "exportall":
        StdI.export_all = _store_with_check_dup_i(keyword, value, StdI.export_all)
    elif keyword == "lattice_gp":
        StdI.lattice_gp = _store_with_check_dup_i(keyword, value, StdI.lattice_gp)
    else:
        return False
    return True


# ===================================================================
#  stdface_main -- top-level entry point
# ===================================================================


def stdface_main(fname: str, solver: str = "HPhi") -> None:
    """Read a Standard-mode input file and generate Expert-mode definition files.

    This is the Python translation of the C function ``StdFace_main()``
    (``StdFace_main.c``, lines 2455--3067).  It performs the following
    sequence of operations:

    1. Create and initialise a :class:`StdIntList` parameter structure.
    2. Parse the input file, mapping keywords to structure fields.
    3. Validate and normalise the model name, lattice name and method.
    4. Dispatch to the appropriate lattice-generation function.
    5. Write all Expert-mode definition files required by the chosen
       solver.

    Parameters
    ----------
    fname : str
        Path to the Standard-mode input file.
    solver : str, optional
        Name of the solver that will consume the generated files.
        Must be one of ``"HPhi"``, ``"mVMC"``, ``"UHF"`` or ``"HWAVE"``.
        Defaults to ``"HPhi"``.  This replaces the compile-time
        ``#ifdef`` mechanism of the original C code.

    Raises
    ------
    SystemExit
        If the input file cannot be opened, a keyword is duplicated or
        unrecognised, or a required parameter is missing.
    """
    # ------------------------------------------------------------------
    #  Initialise
    # ------------------------------------------------------------------
    StdI = StdIntList()
    StdI.solver = solver

    print("\n######  Input Parameter of Standard Intarface  ######")

    try:
        fp = open(fname, "r")
    except OSError:
        print(f"\n  ERROR !  Cannot open input file {fname} !\n")
        exit_program(-1)
    else:
        print(f"\n  Open Standard-Mode Inputfile {fname} \n")

    _reset_vals(StdI)

    # ------------------------------------------------------------------
    #  Parse input file
    # ------------------------------------------------------------------
    for raw_line in fp:
        line = _trim_space_quote(raw_line)

        if line.startswith("//"):
            print("  Skipping a line.")
            continue
        if line == "":
            print("  Skipping a line.")
            continue

        parts = line.split("=", 1)
        if len(parts) < 2:
            print('\n  ERROR !  "=" is NOT found !\n')
            exit_program(-1)

        keyword = _text2lower(parts[0])
        value = parts[1]
        print(f"  KEYWORD : {keyword:<20s} | VALUE : {value} ")

        if not _parse_common_keyword(keyword, value, StdI):
            if not _parse_solver_keyword(keyword, value, StdI, solver):
                print("ERROR ! Unsupported Keyword in Standard mode!")
                exit_program(-1)

    fp.close()

    # ------------------------------------------------------------------
    #  Construct Model
    # ------------------------------------------------------------------
    print("")
    print("#######  Construct Model  #######")
    print("")

    # CDataFileHead default
    if StdI.CDataFileHead == "****":
        StdI.CDataFileHead = "zvo"
        print(f"    CDataFileHead = {'zvo':<12s}######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"    CDataFileHead = {StdI.CDataFileHead}")

    # ------------------------------------------------------------------
    #  Check / normalise the model name
    # ------------------------------------------------------------------
    StdI.lGC = 0
    StdI.lBoost = 0

    if StdI.model in ("fermionhubbard", "hubbard"):
        StdI.model = "hubbard"
    elif StdI.model in ("fermionhubbardgc", "hubbardgc"):
        StdI.model = "hubbard"
        StdI.lGC = 1
    elif StdI.model == "spin":
        StdI.model = "spin"
    elif StdI.model == "spingc":
        StdI.model = "spin"
        StdI.lGC = 1
    elif solver == "HPhi" and StdI.model in ("spingcboost", "spingccma"):
        StdI.model = "spin"
        StdI.lGC = 1
        StdI.lBoost = 1
    elif StdI.model in ("kondolattice", "kondo"):
        StdI.model = "kondo"
    elif StdI.model in ("kondolatticegc", "kondogc"):
        StdI.model = "kondo"
        StdI.lGC = 1
    else:
        _unsupported_system(StdI.model, StdI.lattice)

    # ------------------------------------------------------------------
    #  Check / normalise the method (HPhi only)
    # ------------------------------------------------------------------
    if solver == "HPhi":
        if StdI.method in ("direct", "alldiag"):
            StdI.method = "fulldiag"
        elif StdI.method in ("te", "time-evolution"):
            StdI.method = "timeevolution"

        # Compute vector potential and electrical field
        if StdI.method == "timeevolution":
            _vector_potential(StdI)

    # ------------------------------------------------------------------
    #  Generate Hamiltonian definition files -- lattice dispatch
    # ------------------------------------------------------------------
    lattice = StdI.lattice
    if lattice in ("chain", "chainlattice"):
        chain_lattice.chain(StdI)
    elif lattice in ("face-centeredorthorhombic", "fcorthorhombic", "fco",
                     "face-centeredcubic", "fccubic", "fcc"):
        fc_ortho.fc_ortho(StdI)
    elif lattice in ("honeycomb", "honeycomblattice"):
        honeycomb_lattice.honeycomb(StdI)
    elif lattice in ("kagome", "kagomelattice"):
        kagome.kagome(StdI)
    elif lattice in ("ladder", "ladderlattice"):
        ladder.ladder(StdI)
    elif lattice in ("orthorhombic", "simpleorthorhombic",
                     "cubic", "simplecubic"):
        orthorhombic.orthorhombic(StdI)
    elif lattice == "pyrochlore":
        pyrochlore.pyrochlore(StdI)
    elif lattice in ("tetragonal", "tetragonallattice",
                     "square", "squarelattice"):
        square_lattice.tetragonal(StdI)
    elif lattice in ("triangular", "triangularlattice"):
        triangular_lattice.triangular(StdI)
    elif lattice == "wannier90":
        wannier90_mod.wannier90(StdI)
    else:
        _unsupported_system(StdI.model, StdI.lattice)

    # ------------------------------------------------------------------
    #  HPhi extras: LargeValue and Boost
    # ------------------------------------------------------------------
    if solver == "HPhi":
        _large_value(StdI)

        if StdI.lBoost == 1:
            if lattice in ("chain", "chainlattice"):
                chain_lattice.chain_boost(StdI)
            elif lattice in ("honeycomb", "honeycomblattice"):
                honeycomb_lattice.honeycomb_boost(StdI)
            elif lattice in ("kagome", "kagomelattice"):
                kagome.kagome_boost(StdI)
            elif lattice in ("ladder", "ladderlattice"):
                ladder.ladder_boost(StdI)
            else:
                _unsupported_system(StdI.model, StdI.lattice)

    # ------------------------------------------------------------------
    #  Print Expert input files
    # ------------------------------------------------------------------
    print("")
    print("######  Print Expert input files  ######")
    print("")

    if solver == "HPhi":
        _print_loc_spin(StdI)
        _print_trans(StdI)
        _print_interactions(StdI)
        _check_mod_para(StdI)
        _print_mod_para(StdI)
        _print_excitation(StdI)
        if StdI.method == "timeevolution":
            _print_pump(StdI)
        _print_calc_mod(StdI)
        _check_output_mode(StdI)
        _print_1_green(StdI)
        _print_2_green(StdI)
        _print_namelist(StdI)

    elif solver == "mVMC":
        _print_loc_spin(StdI)
        _print_trans(StdI)
        _print_interactions(StdI)
        _check_mod_para(StdI)
        _print_mod_para(StdI)

        if StdI.lGC == 0 and (StdI.Sz2 == 0 or StdI.Sz2 == NaN_i):
            StdI.ComplexType = print_val_i("ComplexType", StdI.ComplexType, 0)
        else:
            StdI.ComplexType = print_val_i("ComplexType", StdI.ComplexType, 1)

        generate_orb(StdI)
        proj(StdI)
        print_jastrow(StdI)
        if StdI.lGC == 1 or (StdI.Sz2 != 0 and StdI.Sz2 != NaN_i):
            _print_orb_para(StdI)
        _print_gutzwiller(StdI)
        _print_orb(StdI)
        _check_output_mode(StdI)
        _print_1_green(StdI)
        _print_2_green(StdI)
        _print_namelist(StdI)

    elif solver == "UHF":
        _print_loc_spin(StdI)
        _print_trans(StdI)
        _print_interactions(StdI)
        _check_mod_para(StdI)
        _print_mod_para(StdI)
        _check_output_mode(StdI)
        _print_1_green(StdI)
        _print_namelist(StdI)

    elif solver == "HWAVE":
        if StdI.calcmode == "uhfr":
            _print_trans(StdI)
            _print_interactions(StdI)
            _check_mod_para(StdI)
            _check_output_mode(StdI)
            _print_1_green(StdI)
        else:
            export_geometry(StdI)
            export_interaction(StdI)

    # ------------------------------------------------------------------
    #  Finalise
    # ------------------------------------------------------------------
    print("\n######  Input files are generated.  ######\n")
