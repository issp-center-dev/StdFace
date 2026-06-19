"""HPhi solver-specific output functions.

This module contains functions that generate definition files specific to the
HPhi (Exact Diagonalization) solver.  These were extracted from
``stdface_main.py`` during refactoring.

Functions
---------
large_value
    Compute the ``LargeValue`` parameter for TPQ calculations.
print_calc_mod
    Write ``calcmod.def`` with calculation-mode integers.
print_excitation
    Write ``single.def`` or ``pair.def`` for spectrum calculations.
vector_potential
    Compute vector potential A(t) and electric field E(t) for time evolution.
print_pump
    Write ``teone.def`` or ``tetwo.def`` for time-evolution pump terms.

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

import logging
import math
from collections.abc import Callable
from typing import NamedTuple

import numpy as np

from ...core.stdface_vals import StdIntList, ModelType, MethodType, AMPLITUDE_EPS
from ...core.param_check import print_val_d, print_val_i
from ...writer.common_writer import _merge_duplicate_terms

# ---------------------------------------------------------------------------
#  String → integer dispatch tables for calcmod.def
# ---------------------------------------------------------------------------

METHOD_TO_CALC_TYPE: dict[MethodType, int] = {
    MethodType.LANCZOS:        0,
    MethodType.LANCZOS_ENERGY: 0,
    MethodType.TPQ:            1,
    MethodType.FULLDIAG:       2,
    MethodType.CG:             3,
    MethodType.TIME_EVOLUTION: 4,
    MethodType.CTPQ:           5,
}
"""Maps each HPhi calculation method to its ``CalcType`` integer code."""

RESTART_TO_INT: dict[str, int] = {
    "none":         0,
    "restart_out":  1,
    "save":         1,
    "restartsave":  2,
    "restart":      2,
    "restart_in":   3,
}
"""Maps each HPhi ``Restart`` string to its integer code."""

CALC_SPEC_TO_INT: dict[str, int] = {
    "none":         0,
    "normal":       1,
    "noiteration":  2,
    "restart_out":  3,
    "restart_in":   4,
    "restartsave":  5,
    "restart":      5,
}
"""Maps each HPhi ``CalcSpec`` string to its integer code."""

MODEL_GC_TO_CALC_MODEL: dict[tuple, int] = {
    (ModelType.HUBBARD, 0): 0,
    (ModelType.HUBBARD, 1): 3,
    (ModelType.SPIN,    0): 1,
    (ModelType.SPIN,    1): 4,
    (ModelType.KONDO,   0): 2,
    (ModelType.KONDO,   1): 5,
}
"""Maps ``(ModelType, lGC)`` to the ``CalcModel`` integer code."""

INITIAL_VEC_TYPE_TO_INT: dict[str, int] = {
    "c": 0,
    "r": 1,
}
"""Maps the ``InitialVecType`` string to its integer code."""

logger = logging.getLogger(__name__)


class _IOFlags(NamedTuple):
    """Pair of input/output flag integers for I/O dispatch tables.

    Attributes
    ----------
    flag0 : int
        First flag (InputEigenVec for EigenVecIO, iOutputHam for HamIO).
    flag1 : int
        Second flag (OutputEigenVec for EigenVecIO, iInputHam for HamIO).
    """

    flag0: int
    flag1: int


EIGENVEC_IO_TO_FLAGS: dict[str, _IOFlags] = {
    "none":  _IOFlags(0, 0),
    "in":    _IOFlags(1, 0),
    "out":   _IOFlags(0, 1),
    "inout": _IOFlags(1, 1),
}
"""Maps ``EigenVecIO`` string to ``(InputEigenVec, OutputEigenVec)`` flags."""

HAM_IO_TO_FLAGS: dict[str, _IOFlags] = {
    "none": _IOFlags(0, 0),
    "out":  _IOFlags(1, 0),
    "in":   _IOFlags(0, 1),
}
"""Maps ``HamIO`` string to ``(iOutputHam, iInputHam)`` flags."""

OUTPUT_EX_VEC_TO_INT: dict[str, int] = {
    "none": 0,
    "out":  1,
}
"""Maps ``OutputExVec`` string to its integer code."""


def _resolve_string_param(
    StdI: StdIntList,
    field: str,
    label: str,
    default: str,
    default_value: int | tuple[int, int],
    dispatch: dict,
) -> int | tuple[int, int]:
    """Resolve a string-valued parameter to its integer code(s).

    If the field on *StdI* is ``None`` (unset), it is set to *default*
    and *default_value* is returned.  Otherwise the field value is looked
    up in *dispatch*; a missing key raises ``ValueError``.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure (field is read/written in-place).
    field : str
        Attribute name on *StdI* (e.g. ``"Restart"``).
    label : str
        Display label for log messages (e.g. ``"Restart"``).
    default : str
        Default string value when unset (e.g. ``"none"``).
    default_value : int or tuple of int
        Value(s) to return when the field is unset.
    dispatch : dict
        Mapping from string values to integer code(s).

    Returns
    -------
    int or tuple of int
        The resolved integer code(s).
    """
    field_val = getattr(StdI, field)
    if field_val is None:
        setattr(StdI, field, default)
        logger.info(f"  {label:>20s} = {default:<12s}######  DEFAULT VALUE IS USED  ######")
        return default_value
    logger.info(f"  {label:>20s} = {field_val}")
    result = dispatch.get(field_val)
    if result is None:
        msg = f"\n ERROR ! {label} : {field_val}"
        logger.error(msg)
        raise ValueError(msg)
    return result


def large_value(StdI: StdIntList) -> None:
    """Compute and set the ``LargeValue`` parameter for TPQ calculations.

    ``LargeValue`` is the sum of the absolute values of all one-body
    (transfer) and two-body (interaction / exchange / Hund / pair-lift /
    Coulomb-intra / Coulomb-inter) terms, divided by the number of
    sites.  The result is stored in ``StdI.LargeValue`` via
    :func:`print_val_d`.

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
    large_value0 = (
        sum(abs(t[0]) for t in StdI.trans_list)
        + sum(abs(t[0]) for t in StdI.intr_list)
        + sum(abs(t[0]) for t in StdI.Cintra_list)
        + sum(abs(t[0]) for t in StdI.Cinter_list)
        + 2.0 * sum(abs(t[0]) for t in StdI.Ex_list)
        + 2.0 * sum(abs(t[0]) for t in StdI.PairLift_list)
        + 2.0 * sum(abs(t[0]) for t in StdI.Hund_list)
    )

    large_value0 /= float(StdI.nsite)

    StdI.LargeValue = print_val_d("LargeValue", StdI.LargeValue, large_value0)


def _validate_ngpu_scalapack(StdI: StdIntList) -> None:
    """Validate ``NGPU`` and ``Scalapack`` parameters if set.

    Logs the parameter values and raises ``ValueError`` if
    they are out of range.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Reads ``NGPU`` and ``Scalapack``.
    """
    if StdI.NGPU is not None:
        logger.info(f"         NGPU = {StdI.NGPU}")
        if StdI.NGPU < 1:
            logger.info(f"\n ERROR ! NGPU : {StdI.NGPU}")
            msg = "         NGPU should be a positive integer."
            logger.error(msg)
            raise ValueError(msg)

    if StdI.Scalapack is not None:
        logger.info(f"         Scalapack = {StdI.Scalapack}")
        if StdI.Scalapack < 0 or StdI.Scalapack > 1:
            logger.info(f"\n ERROR ! Scalapack : {StdI.Scalapack}")
            msg = "         Scalapack should be 0 or 1."
            logger.error(msg)
            raise ValueError(msg)


class _CalcModParams(NamedTuple):
    """Resolved integer parameters for ``calcmod.def``.

    Attributes
    ----------
    iCalcType : int
        Calculation type code (0–5).
    iCalcModel : int
        Calculation model code (0–5).
    iCalcEigenvec : int
        Whether to calculate eigenvectors (0 or 1).
    iRestart : int
        Restart flag (0–3).
    iCalcSpec : int
        Spectrum calculation flag (0–5).
    iInitialVecType : int
        Initial vector type code.
    InputEigenVec : int
        Whether to read eigenvectors from file.
    OutputEigenVec : int
        Whether to write eigenvectors to file.
    iInputHam : int
        Whether to read Hamiltonian from file.
    iOutputHam : int
        Whether to write Hamiltonian to file.
    iOutputExVec : int
        Whether to output excited-state vectors.
    """

    iCalcType: int
    iCalcModel: int
    iCalcEigenvec: int
    iRestart: int
    iCalcSpec: int
    iInitialVecType: int
    InputEigenVec: int
    OutputEigenVec: int
    iInputHam: int
    iOutputHam: int
    iOutputExVec: int


def _write_calcmod_file(StdI: StdIntList, params: _CalcModParams) -> None:
    """Write ``calcmod.def`` with the resolved integer parameters.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  Reads ``NGPU`` and ``Scalapack``
        for conditional output lines.
    params : _CalcModParams
        Resolved calculation-mode parameters.
    """
    with open("calcmod.def", "w") as fp:
        fp.write("#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution 5:cTPQ\n")
        fp.write("#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC\n")
        fp.write("#Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart\n")
        fp.write("#CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save\n")
        if StdI.NGPU is not None:
            fp.write("#NGPU (for FullDiag): The number of GPU\n")
        if StdI.Scalapack is not None:
            fp.write("#Scalapack (for FullDiag) = 0:w/o ScaLAPACK, 1:w/ ScaLAPACK\n")
        fp.write(f"CalcType {params.iCalcType:3d}\n")
        fp.write(f"CalcModel {params.iCalcModel:3d}\n")
        fp.write(f"ReStart {params.iRestart:3d}\n")
        fp.write(f"CalcSpec {params.iCalcSpec:3d}\n")
        fp.write(f"CalcEigenVec {params.iCalcEigenvec:3d}\n")
        fp.write(f"InitialVecType {params.iInitialVecType:3d}\n")
        fp.write(f"InputEigenVec {params.InputEigenVec:3d}\n")
        fp.write(f"OutputEigenVec {params.OutputEigenVec:3d}\n")
        fp.write(f"InputHam {params.iInputHam:3d}\n")
        fp.write(f"OutputHam {params.iOutputHam:3d}\n")
        fp.write(f"OutputExVec {params.iOutputExVec:3d}\n")
        if StdI.NGPU is not None:
            fp.write(f"NGPU {StdI.NGPU:3d}\n")
        if StdI.Scalapack is not None:
            fp.write(f"Scalapack {StdI.Scalapack:3d}\n")


def print_calc_mod(StdI: StdIntList) -> None:
    """Write ``calcmod.def`` containing calculation-mode integers for HPhi.

    Maps the string-valued parameters ``method``, ``model``, ``Restart``,
    ``InitialVecType``, ``EigenVecIO``, ``HamIO``, ``CalcSpec``, and
    ``OutputExVec`` to the integer codes that HPhi expects.  Validation
    is delegated to :func:`_validate_ngpu_scalapack` and file writing to
    :func:`_write_calcmod_file`.

    Parameters
    ----------
    StdI : StdIntList
        The global parameter structure.  String fields such as
        ``method``, ``model``, ``Restart``, etc. are read (and may be
        overwritten with default values).  ``NGPU`` and ``Scalapack``
        are validated if set.
    """
    logger.info("\n  @ CalcMod\n")

    # ------------------------------------------------------------------
    #  Method  →  CalcType integer
    # ------------------------------------------------------------------
    iCalcEigenvec = 0

    if StdI.method is None:
        msg = "ERROR ! Method is NOT specified !"
        logger.error(msg)
        raise ValueError(msg)

    iCalcType = METHOD_TO_CALC_TYPE.get(StdI.method)
    if iCalcType is None:
        msg = f"\n ERROR ! Unsupported Solver : {StdI.method}"
        logger.error(msg)
        raise ValueError(msg)

    if StdI.method == MethodType.LANCZOS_ENERGY:
        iCalcEigenvec = 1

    if iCalcType != 4:
        StdI.PumpBody = 0

    # ------------------------------------------------------------------
    #  Model
    # ------------------------------------------------------------------
    iCalcModel = MODEL_GC_TO_CALC_MODEL.get((StdI.model, StdI.lGC))
    if iCalcModel is None:
        msg = (f"Unsupported Model / GC combination : "
               f"{StdI.model}, lGC={StdI.lGC}")
        logger.error(msg)
        raise ValueError(msg)

    # ------------------------------------------------------------------
    #  Resolve string parameters to integer codes
    # ------------------------------------------------------------------
    iRestart = _resolve_string_param(
        StdI, "Restart", "Restart", "none", 0, RESTART_TO_INT)

    # InitialVecType: default depends on method (TPQ/CTPQ → -1, else → 0)
    ivt_default = -1 if StdI.method in (MethodType.TPQ, MethodType.CTPQ) else 0
    iInitialVecType = _resolve_string_param(
        StdI, "InitialVecType", "InitialVecType", "c",
        ivt_default, INITIAL_VEC_TYPE_TO_INT)

    InputEigenVec, OutputEigenVec = _resolve_string_param(
        StdI, "EigenVecIO", "EigenVecIO", "none", (0, 0), EIGENVEC_IO_TO_FLAGS)
    if StdI.method == MethodType.TIME_EVOLUTION:
        InputEigenVec = 1

    iOutputHam, iInputHam = _resolve_string_param(
        StdI, "HamIO", "HamIO", "none", (0, 0), HAM_IO_TO_FLAGS)

    iCalcSpec = _resolve_string_param(
        StdI, "CalcSpec", "CalcSpec", "none", 0, CALC_SPEC_TO_INT)

    iOutputExVec = _resolve_string_param(
        StdI, "OutputExVec", "OutputExcitedVec", "none",
        0, OUTPUT_EX_VEC_TO_INT)

    # ------------------------------------------------------------------
    #  Validate and write
    # ------------------------------------------------------------------
    _validate_ngpu_scalapack(StdI)
    _write_calcmod_file(
        StdI,
        _CalcModParams(
            iCalcType, iCalcModel, iCalcEigenvec,
            iRestart, iCalcSpec, iInitialVecType,
            InputEigenVec, OutputEigenVec, iInputHam, iOutputHam, iOutputExVec,
        ),
    )

    logger.info("     calcmod.def is written.\n")


# ---------------------------------------------------------------------------
#  Spectrum-type operator configuration
# ---------------------------------------------------------------------------

def _spectrum_szsz(
    model: ModelType, S2: int,
    coef: list[float], spin: list[list[int]],
) -> tuple[int, int]:
    """Configure SzSz spectrum operators.

    For SPIN models, produces ``S2 + 1`` operators with coefficients equal
    to :math:`S_z` values.  For Hubbard/Kondo models, produces two
    operators with coefficients :math:`\\pm 0.5`.

    Parameters
    ----------
    model : ModelType
        The model type.
    S2 : int
        Twice the spin quantum number.
    coef : list of float
        Coefficient array (modified in place).
    spin : list of list of int
        Spin-index array (modified in place).

    Returns
    -------
    tuple of (int, int)
        ``(NumOp, SpectrumBody)`` — always ``SpectrumBody = 2``.
    """
    if model == ModelType.SPIN:
        NumOp = S2 + 1
        for ispin in range(S2 + 1):
            Sz = float(ispin) - float(S2) * 0.5
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
    return NumOp, 2


def _spectrum_spsm(
    model: ModelType, S2: int,
    coef: list[float], spin: list[list[int]],
) -> tuple[int, int]:
    """Configure S+S- spectrum operators.

    For high-spin SPIN models (``S2 > 1``), produces ``S2`` operators with
    ladder coefficients :math:`\\sqrt{S(S+1) - S_z(S_z+1)}`.  Otherwise
    produces a single operator with coefficient 1.

    Parameters
    ----------
    model : ModelType
        The model type.
    S2 : int
        Twice the spin quantum number.
    coef : list of float
        Coefficient array (modified in place).
    spin : list of list of int
        Spin-index array (modified in place).

    Returns
    -------
    tuple of (int, int)
        ``(NumOp, SpectrumBody)`` — always ``SpectrumBody = 2``.
    """
    if model == ModelType.SPIN and S2 > 1:
        NumOp = S2
        S = float(S2) * 0.5
        for ispin in range(1, S2 + 1):
            Sz = float(S2) * 0.5 - float(ispin)
            coef[ispin - 1] = math.sqrt(S * (S + 1.0) - Sz * (Sz + 1.0))
            spin[ispin - 1][0] = ispin
            spin[ispin - 1][1] = ispin - 1
    else:
        NumOp = 1
        coef[0] = 1.0
        spin[0][0] = 1
        spin[0][1] = 0
    return NumOp, 2


def _spectrum_density(
    model: ModelType, S2: int,
    coef: list[float], spin: list[list[int]],
) -> tuple[int, int]:
    """Configure density spectrum operators.

    Always produces two operators (up and down) with coefficient 1.

    Parameters
    ----------
    model : ModelType
        The model type (unused, included for dispatch signature).
    S2 : int
        Twice the spin quantum number (unused).
    coef : list of float
        Coefficient array (modified in place).
    spin : list of list of int
        Spin-index array (modified in place).

    Returns
    -------
    tuple of (int, int)
        ``(2, 2)`` — two operators, pair body.
    """
    coef[0] = 1.0
    coef[1] = 1.0
    spin[0][0] = 0
    spin[0][1] = 0
    spin[1][0] = 1
    spin[1][1] = 1
    return 2, 2


def _spectrum_up(
    model: ModelType, S2: int,
    coef: list[float], spin: list[list[int]],
) -> tuple[int, int]:
    """Configure spin-up single-particle spectrum operator.

    Parameters
    ----------
    model : ModelType
        The model type (unused, included for dispatch signature).
    S2 : int
        Twice the spin quantum number (unused).
    coef : list of float
        Coefficient array (modified in place).
    spin : list of list of int
        Spin-index array (modified in place).

    Returns
    -------
    tuple of (int, int)
        ``(1, 1)`` — one operator, single body.
    """
    coef[0] = 1.0
    spin[0][0] = 0
    return 1, 1


def _spectrum_down(
    model: ModelType, S2: int,
    coef: list[float], spin: list[list[int]],
) -> tuple[int, int]:
    """Configure spin-down single-particle spectrum operator.

    Parameters
    ----------
    model : ModelType
        The model type (unused, included for dispatch signature).
    S2 : int
        Twice the spin quantum number (unused).
    coef : list of float
        Coefficient array (modified in place).
    spin : list of list of int
        Spin-index array (modified in place).

    Returns
    -------
    tuple of (int, int)
        ``(1, 1)`` — one operator, single body.
    """
    coef[0] = 1.0
    spin[0][0] = 1
    return 1, 1


_SPECTRUM_HANDLERS: dict[str, Callable] = {
    "szsz": _spectrum_szsz,
    "s+s-": _spectrum_spsm,
    "density": _spectrum_density,
    "up": _spectrum_up,
    "down": _spectrum_down,
}
"""Dispatch table mapping spectrum-type strings to handler functions."""


def _configure_spectrum_ops(
    spectrum_type: str,
    model: ModelType,
    S2: int,
    coef: list[float],
    spin: list[list[int]],
) -> tuple[int, int]:
    """Configure operator coefficients and spins for a spectrum type.

    Dispatches to a per-type handler function via ``_SPECTRUM_HANDLERS``.
    Populates *coef* and *spin* in place and returns ``(NumOp, SpectrumBody)``.

    Parameters
    ----------
    spectrum_type : str
        One of ``"szsz"``, ``"s+s-"``, ``"density"``, ``"up"``, ``"down"``.
    model : ModelType
        The model type (SPIN, HUBBARD, or KONDO).
    S2 : int
        Twice the spin quantum number.
    coef : list of float
        Coefficient array to populate (modified in place).
    spin : list of list of int
        Spin-index array to populate (modified in place, each element is
        ``[spin_out, spin_in]``).

    Returns
    -------
    tuple of (int, int)
        ``(NumOp, SpectrumBody)`` where *NumOp* is the number of operators
        and *SpectrumBody* is 1 (single) or 2 (pair).

    Raises
    ------
    ValueError
        If *spectrum_type* is not recognized.
    """
    handler = _SPECTRUM_HANDLERS.get(spectrum_type)
    if handler is not None:
        return handler(model, S2, coef, spin)

    msg = f"\n ERROR ! SpectrumType : {spectrum_type}"
    logger.error(msg)
    raise ValueError(msg)


def _compute_fourier_coefficients(StdI: StdIntList) -> tuple[list[float], list[float]]:
    """Compute site Fourier coefficients for spectrum excitation.

    For each site in the super-cell, compute the real and imaginary parts
    of :math:`\\exp(2\\pi i \\mathbf{q} \\cdot \\mathbf{r})`, where
    :math:`\\mathbf{q}` is ``SpectrumQ`` and :math:`\\mathbf{r}` is
    the cell + basis position.  For the Kondo model, the itinerant-site
    coefficients are duplicated to the local-spin sites.

    Parameters
    ----------
    StdI : StdIntList
        Global parameter structure.

    Returns
    -------
    fourier_r : list of float
        Real parts, length ``nsite``.
    fourier_i : list of float
        Imaginary parts, length ``nsite``.
    """
    fourier_r = np.zeros(StdI.nsite)
    fourier_i = np.zeros(StdI.nsite)
    n_computed = StdI.NCell * StdI.NsiteUC

    if n_computed > 0:
        # Compute position vectors: Cell[icell] + tau[itau] for all (icell, itau) pairs
        # Shape: (NCell, 1, 3) + (1, NsiteUC, 3) -> (NCell, NsiteUC, 3)
        positions = StdI.Cell[:StdI.NCell, :].astype(float)[:, np.newaxis, :] + \
                    StdI.tau[:StdI.NsiteUC, :][np.newaxis, :, :]
        # Dot with SpectrumQ: shape (NCell, NsiteUC), then flatten to site order
        Cphase_flat = (2.0 * math.pi * (positions @ StdI.SpectrumQ)).ravel()
        fourier_r[:n_computed] = np.cos(Cphase_flat)
        fourier_i[:n_computed] = np.sin(Cphase_flat)

    if StdI.model == ModelType.KONDO:
        half = StdI.nsite // 2
        fourier_r[half:] = fourier_r[:half]
        fourier_i[half:] = fourier_i[:half]

    return fourier_r, fourier_i


def _write_excitation_file(
    StdI: StdIntList,
    NumOp: int,
    coef: list[float],
    spin: list[list[int]],
    fourier_r: list[float],
    fourier_i: list[float],
) -> None:
    """Write ``single.def`` or ``pair.def`` excitation file.

    Depending on ``StdI.SpectrumBody``, writes either a single-body
    excitation file (``single.def``) or a pair-body file (``pair.def``).

    Parameters
    ----------
    StdI : StdIntList
        Global parameter structure.
    NumOp : int
        Number of operator components per site.
    coef : list of float
        Coefficients for each operator component.
    spin : list of list of int
        Spin indices ``[s1, s2]`` for each operator component.
    fourier_r : list of float
        Real parts of Fourier coefficients, length ``nsite``.
    fourier_i : list of float
        Imaginary parts of Fourier coefficients, length ``nsite``.
    """
    if StdI.SpectrumBody == 1:
        lines = ["=============================================\n"]
        if StdI.model == ModelType.KONDO:
            lines.append(f"NSingle {StdI.nsite // 2 * NumOp}\n")
        else:
            lines.append(f"NSingle {StdI.nsite * NumOp}\n")
        lines.append("=============================================\n")
        lines.append("============== Single Excitation ============\n")
        lines.append("=============================================\n")
        if StdI.model == ModelType.KONDO:
            for isite in range(StdI.nsite // 2, StdI.nsite):
                lines.append(f"{isite} {spin[0][0]} 0 "
                             f"{fourier_r[isite] * coef[0]:25.15f} "
                             f"{fourier_i[isite] * coef[0]:25.15f}\n")
        else:
            for isite in range(StdI.nsite):
                lines.append(f"{isite} {spin[0][0]} 0 "
                             f"{fourier_r[isite] * coef[0]:25.15f} "
                             f"{fourier_i[isite] * coef[0]:25.15f}\n")
        with open("single.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("      single.def is written.\n")
    else:
        lines = ["=============================================\n",
                 f"NPair {StdI.nsite * NumOp}\n",
                 "=============================================\n",
                 "=============== Pair Excitation =============\n",
                 "=============================================\n"]
        for isite in range(StdI.nsite):
            for ispin in range(NumOp):
                lines.append(f"{isite} {spin[ispin][0]} {isite} {spin[ispin][1]} 1 "
                             f"{fourier_r[isite] * coef[ispin]:25.15f} "
                             f"{fourier_i[isite] * coef[ispin]:25.15f}\n")
        with open("pair.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("        pair.def is written.\n")


def print_excitation(StdI: StdIntList) -> None:
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
    if StdI.model == ModelType.SPIN and StdI.S2 > 1:
        n_alloc = StdI.S2 + 1
    else:
        n_alloc = 2

    coef = [0.0] * n_alloc
    spin = [[0, 0] for _ in range(n_alloc)]

    logger.info("\n  @ Spectrum\n")

    StdI.SpectrumQ[0] = print_val_d("SpectrumQW", StdI.SpectrumQ[0], 0.0)
    StdI.SpectrumQ[1] = print_val_d("SpectrumQL", StdI.SpectrumQ[1], 0.0)
    StdI.SpectrumQ[2] = print_val_d("SpectrumQH", StdI.SpectrumQ[2], 0.0)

    # ------------------------------------------------------------------
    #  Resolve SpectrumType default and determine NumOp, coef, spin
    # ------------------------------------------------------------------
    if StdI.SpectrumType is None:
        StdI.SpectrumType = "szsz"
        logger.info("     SpectrumType = szsz        ######  DEFAULT VALUE IS USED  ######")
    else:
        logger.info(f"     SpectrumType = {StdI.SpectrumType}")

    NumOp, StdI.SpectrumBody = _configure_spectrum_ops(
        StdI.SpectrumType, StdI.model, StdI.S2, coef, spin)

    # Compute Fourier coefficients
    fourier_r, fourier_i = _compute_fourier_coefficients(StdI)

    # Write single.def or pair.def
    _write_excitation_file(StdI, NumOp, coef, spin, fourier_r, fourier_i)


# ---------------------------------------------------------------------------
#  Pump-type A(t) / E(t) computation functions
# ---------------------------------------------------------------------------

def _pulselaser_At_Et(
    time: float, V: float, freq: float, tshift: float, tdump: float,
) -> tuple[float, float]:
    """Compute A(t) and E(t) for Gaussian-enveloped cosine pulse laser.

    Parameters
    ----------
    time : float
        Current time.
    V : float
        Vector potential amplitude for one component.
    freq : float
        Laser frequency.
    tshift : float
        Time shift for the pulse center.
    tdump : float
        Gaussian pulse width parameter.

    Returns
    -------
    tuple of (float, float)
        ``(A_component, E_component)`` for this time and component.
    """
    dt_shift = time - tshift
    gauss = math.exp(-0.5 * dt_shift ** 2 / (tdump * tdump))
    cos_val = math.cos(freq * dt_shift)
    sin_val = math.sin(freq * dt_shift)
    At = V * cos_val * gauss
    Et = -V * ((-dt_shift) / (tdump * tdump) * cos_val - freq * sin_val) * gauss
    return At, Et


def _aclaser_At_Et(
    time: float, V: float, freq: float, tshift: float, tdump: float,
) -> tuple[float, float]:
    """Compute A(t) and E(t) for sinusoidal AC laser.

    Parameters
    ----------
    time : float
        Current time.
    V : float
        Vector potential amplitude for one component.
    freq : float
        Laser frequency.
    tshift : float
        Time shift.
    tdump : float
        Unused (kept for uniform signature).

    Returns
    -------
    tuple of (float, float)
        ``(A_component, E_component)`` for this time and component.
    """
    dt_shift = time - tshift
    At = V * math.sin(freq * dt_shift)
    Et = V * math.cos(freq * dt_shift) * freq
    return At, Et


def _dclaser_At_Et(
    time: float, V: float, freq: float, tshift: float, tdump: float,
) -> tuple[float, float]:
    """Compute A(t) and E(t) for linearly ramped DC laser.

    Parameters
    ----------
    time : float
        Current time.
    V : float
        Vector potential amplitude for one component.
    freq : float
        Unused (kept for uniform signature).
    tshift : float
        Unused (kept for uniform signature).
    tdump : float
        Unused (kept for uniform signature).

    Returns
    -------
    tuple of (float, float)
        ``(A_component, E_component)`` for this time and component.
    """
    return V * time, -V


class _PumpTypeConfig(NamedTuple):
    """Configuration for a pump-type entry in the dispatch table.

    Attributes
    ----------
    pump_body : int
        Pump body type (1 = one-body laser, 2 = two-body quench).
    handler_fn : Callable or None
        Field computation function with signature
        ``(time, V, freq, tshift, tdump) -> (At, Et)``, or ``None``
        for quench (no field computation).
    """

    pump_body: int
    handler_fn: Callable | None


_PUMP_TYPE_HANDLERS: dict[str, _PumpTypeConfig] = {
    "quench":     _PumpTypeConfig(2, None),
    "pulselaser": _PumpTypeConfig(1, _pulselaser_At_Et),
    "aclaser":    _PumpTypeConfig(1, _aclaser_At_Et),
    "dclaser":    _PumpTypeConfig(1, _dclaser_At_Et),
}
"""Maps each PumpType to a :class:`_PumpTypeConfig`.

For ``PumpBody == 2`` (quench), *handler_fn* is ``None`` — no field is
computed.  For ``PumpBody == 1``, *handler_fn* is called once per
(timestep, component) pair with signature
``(time, V, freq, tshift, tdump) -> (At, Et)``.
"""


def vector_potential(StdI: StdIntList) -> None:
    """Compute vector potential A(t) and electric field E(t) for time evolution.

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
    logger.info("\n  @ Time-evolution\n")

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
    StdI.At = np.zeros((StdI.Lanczos_max, 3))
    Et = np.zeros((StdI.Lanczos_max, 3))

    # Resolve PumpType default
    if StdI.PumpType is None:
        StdI.PumpType = "quench"
        logger.info("     PumpType = quench        ######  DEFAULT VALUE IS USED  ######")
    else:
        logger.info(f"     PumpType = {StdI.PumpType}")

    # Dispatch on PumpType
    handler_entry = _PUMP_TYPE_HANDLERS.get(StdI.PumpType)
    if handler_entry is None:
        msg = f"\n ERROR ! PumpType : {StdI.PumpType}"
        logger.error(msg)
        raise ValueError(msg)

    StdI.PumpBody, handler_fn = handler_entry

    # Compute A(t) and E(t) for one-body pump types
    if handler_fn is not None:
        for it in range(StdI.Lanczos_max):
            time = StdI.dt * float(it)
            results = [handler_fn(time, StdI.VecPot[ii], StdI.freq,
                                  StdI.tshift, StdI.tdump) for ii in range(3)]
            StdI.At[it, :] = [r[0] for r in results]
            Et[it, :] = [r[1] for r in results]

    # ------------------------------------------------------------------
    #  Write potential.dat for one-body pump
    # ------------------------------------------------------------------
    if StdI.PumpBody == 1:
        lines = ["# Time A_W A_L A_H E_W E_L E_H\n"]
        for it in range(StdI.Lanczos_max):
            time = StdI.dt * float(it)
            lines.append(f"{time:f} "
                         f"{StdI.At[it][0]:f} {StdI.At[it][1]:f} {StdI.At[it][2]:f} "
                         f"{Et[it][0]:f} {Et[it][1]:f} {Et[it][2]:f}\n")
        with open("potential.dat", "w") as fp:
            fp.write("".join(lines))


def print_pump(StdI: StdIntList) -> None:
    """Write ``teone.def`` or ``tetwo.def`` for time-evolution pump terms.

    - ``PumpBody == 1``: writes ``teone.def`` -- one-body pump terms.
      Equivalent pump terms (same index quadruples) are merged and
      entries below ``AMPLITUDE_EPS`` in magnitude are suppressed.
    - ``PumpBody == 2``: writes ``tetwo.def`` -- two-body pump terms
      using ``Uquench`` for on-site Hubbard quench interaction.

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
        lines = ["=============================================\n",
                 f"AllTimeStep {StdI.Lanczos_max}\n",
                 "=============================================\n",
                 "=========  OneBody Time Evolution  ==========\n",
                 "=============================================\n"]

        for it in range(StdI.Lanczos_max):
            npump0 = _merge_duplicate_terms(
                StdI.pumpindx[it], StdI.pump[it], StdI.npump[it])

            lines.append(f"{StdI.dt * float(it):f}  {npump0}\n")

            for ipump in range(StdI.npump[it]):
                val = StdI.pump[it][ipump]
                if abs(val) <= AMPLITUDE_EPS:
                    continue
                i0, s0, i1, s1 = StdI.pumpindx[it][ipump]
                lines.append(
                    f"{i0:5d} {s0:5d} {i1:5d} {s1:5d} "
                    f"{val.real:25.15f} {val.imag:25.15f}\n"
                )

        with open("teone.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("      teone.def is written.\n")

    else:
        lines = ["=============================================\n",
                 f"AllTimeStep {StdI.Lanczos_max}\n",
                 "=============================================\n",
                 "========== TwoBody Time Evolution ===========\n",
                 "=============================================\n"]

        for it in range(StdI.Lanczos_max):
            lines.append(f"{StdI.dt * float(it):f}  {StdI.nsite}\n")
            for isite in range(StdI.nsite):
                lines.append(f"{isite:5d} {0:5d} {isite:5d} {0:5d} "
                             f"{isite:5d} {1:5d} {isite:5d} {1:5d} "
                             f"{StdI.Uquench:25.15f}  {0.0:25.15f}\n")

        with open("tetwo.def", "w") as fp:
            fp.write("".join(lines))
        logger.info("        tetwo.def is written.\n")
