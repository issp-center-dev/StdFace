"""Parameter validation, printing, and program-exit utilities.

This module provides helper functions for printing parameter values with
default-value handling, aborting on unused or missing parameters, and
terminating the program.  These utilities are shared across all lattice
implementations, solver writers, and the main entry point.

Functions
---------
print_val_d
    Print / default a real-valued parameter.
print_val_dd
    Print / default a real-valued parameter with two fallback defaults.
print_val_c
    Print / default a complex-valued parameter.
print_val_i
    Print / default an integer parameter.
not_used_d
    Abort if a real or complex parameter is specified but unused.
not_used_j
    Abort if any component of a J-type interaction is specified but unused.
not_used_i
    Abort if an integer parameter is specified but unused.
required_val_i
    Abort if a required integer parameter is missing.

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

import itertools
import logging
import math

import numpy as np

from .stdface_vals import NaN_i


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
#  PrintVal / NotUsed / Required helpers
# ---------------------------------------------------------------------------


def print_val_d(valname: str, val: float, val0: float) -> float:
    """Print and optionally set a real-valued parameter.

    If *val* is NaN, set it to the default *val0* and print with a
    ``DEFAULT VALUE`` tag.

    Parameters
    ----------
    valname : str
        Name of the variable (for display).
    val : float
        Current value (may be NaN if not specified).
    val0 : float
        Default value to use when *val* is NaN.

    Returns
    -------
    float
        The (possibly updated) value.
    """
    if math.isnan(val):
        val = val0
        logger.info("  %15s = %-10.5f  ######  DEFAULT VALUE IS USED  ######", valname, val)
    else:
        logger.info("  %15s = %-10.5f", valname, val)
    return val


def print_val_dd(valname: str, val: float, val0: float, val1: float) -> float:
    """Print and optionally set a real-valued parameter with two defaults.

    If *val* is NaN, use the primary default *val0* if specified,
    otherwise fall back to the secondary default *val1*.  Delegates to
    :func:`print_val_d` after resolving the effective default.

    Parameters
    ----------
    valname : str
        Name of the variable (for display).
    val : float
        Current value (may be NaN).
    val0 : float
        Primary default (may itself be NaN).
    val1 : float
        Secondary default.

    Returns
    -------
    float
        The (possibly updated) value.
    """
    default = val1 if math.isnan(val0) else val0
    return print_val_d(valname, val, default)


def print_val_c(valname: str, val: complex, val0: complex) -> complex:
    """Print and optionally set a complex-valued parameter.

    Parameters
    ----------
    valname : str
        Name of the variable (for display).
    val : complex
        Current value (real part NaN if not specified).
    val0 : complex
        Default value.

    Returns
    -------
    complex
        The (possibly updated) value.
    """
    if math.isnan(val.real):
        val = val0
        logger.info(
            "  %15s = %-10.5f %-10.5f  ######  DEFAULT VALUE IS USED  ######",
            valname, val.real, val.imag,
        )
    else:
        logger.info("  %15s = %-10.5f %-10.5f", valname, val.real, val.imag)
    return val


def print_val_i(valname: str, val: int, val0: int) -> int:
    """Print and optionally set an integer parameter.

    If *val* equals the sentinel ``2147483647`` (NaN_i), set it to *val0*.

    Parameters
    ----------
    valname : str
        Name of the variable (for display).
    val : int
        Current value (sentinel if not specified).
    val0 : int
        Default value.

    Returns
    -------
    int
        The (possibly updated) value.
    """
    if val is None or val == NaN_i:
        val = val0
        logger.info("  %15s = %-10d  ######  DEFAULT VALUE IS USED  ######", valname, val)
    else:
        logger.info("  %15s = %-10d", valname, val)
    return val


def _fail_not_used(valname: str) -> None:
    """Log a "specified but not used" error and raise.

    Parameters
    ----------
    valname : str
        Name of the unused parameter.

    Raises
    ------
    ValueError
        Always raised after logging the error message.
    """
    msg = (
        f"{valname} is SPECIFIED but will NOT be USED. "
        "Please COMMENT-OUT this line, "
        "or check this input is REALLY APPROPRIATE for your purpose."
    )
    logger.error(msg)
    raise ValueError(msg)


def not_used_d(valname: str, val: float | complex) -> None:
    """Abort if a real parameter is specified but will not be used.

    Parameters
    ----------
    valname : str
        Name of the variable.
    val : float or complex
        Value to check (abort if not NaN).  If complex, the real part
        is tested (matching the C behaviour of implicitly casting
        ``double complex`` to ``double``).
    """
    check = val.real if isinstance(val, complex) else val
    if not math.isnan(check):
        _fail_not_used(valname)


# Spin-interaction suffix matrix for 3x3 J-coupling tensors.
# Shared by not_used_j and input_params module.
SPIN_SUFFIXES: list[list[str]] = [
    ["x", "xy", "xz"],
    ["yx", "y", "yz"],
    ["zx", "zy", "z"],
]
"""3x3 suffix matrix for spin-interaction tensor components (Jx, Jxy, etc.)."""


def not_used_j(valname: str, JAll: float, J: np.ndarray) -> None:
    """Abort if any component of a J-type interaction is specified but unused.

    Parameters
    ----------
    valname : str
        Base name of the variable (e.g. ``"J0"``).
    JAll : float
        Scalar (isotropic) part.
    J : numpy.ndarray
        3x3 matrix of anisotropic components.
    """
    not_used_d(valname, JAll)
    for i1, i2 in itertools.product(range(3), repeat=2):
        not_used_d(f"{valname}{SPIN_SUFFIXES[i1][i2]}", J[i1, i2])


def not_used_i(valname: str, val: int) -> None:
    """Abort if an integer parameter is specified but will not be used.

    Parameters
    ----------
    valname : str
        Name of the variable.
    val : int
        Value to check (abort if not the sentinel 2147483647).
    """
    if val is not None and val != NaN_i:
        _fail_not_used(valname)


def required_val_i(valname: str, val: int) -> None:
    """Abort if a required integer parameter is missing.

    Parameters
    ----------
    valname : str
        Name of the variable.
    val : int
        Value to check (abort if equals the sentinel 2147483647).

    Raises
    ------
    ValueError
        If ``val`` equals the sentinel (i.e. the parameter is unset).
    """
    if val is None or val == NaN_i:
        msg = f"{valname} is NOT specified."
        logger.error(msg)
        raise ValueError(msg)
    logger.info("  %s = %d", valname, val)
