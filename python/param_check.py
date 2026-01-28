"""Parameter validation, printing, and program-exit utilities.

This module provides helper functions for printing parameter values with
default-value handling, aborting on unused or missing parameters, and
terminating the program.  These utilities are shared across all lattice
implementations, solver writers, and the main entry point.

Functions
---------
exit_program
    Terminate the program with a given error code.
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

import math
import sys

import numpy as np

from stdface_vals import NaN_i


# ---------------------------------------------------------------------------
#  Exit wrapper
# ---------------------------------------------------------------------------


def exit_program(errorcode: int) -> None:
    """Terminate the program with the given error code.

    Parameters
    ----------
    errorcode : int
        Exit code passed to ``sys.exit``.

    Notes
    -----
    In the original C code this was ``StdFace_exit`` which wrapped
    ``MPI_Abort`` / ``MPI_Finalize``.  The Python version simply calls
    ``sys.exit``.
    """
    sys.stdout.flush()
    sys.stderr.flush()
    sys.exit(errorcode)


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
        print(f"  {valname:>15s} = {val:<10.5f}  ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"  {valname:>15s} = {val:<10.5f}")
    return val


def print_val_dd(valname: str, val: float, val0: float, val1: float) -> float:
    """Print and optionally set a real-valued parameter with two defaults.

    If *val* is NaN, use the primary default *val0* if specified,
    otherwise fall back to the secondary default *val1*.

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
    if math.isnan(val):
        if math.isnan(val0):
            val = val1
        else:
            val = val0
        print(f"  {valname:>15s} = {val:<10.5f}  ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"  {valname:>15s} = {val:<10.5f}")
    return val


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
        print(f"  {valname:>15s} = {val.real:<10.5f} {val.imag:<10.5f}"
              f"  ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"  {valname:>15s} = {val.real:<10.5f} {val.imag:<10.5f}")
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
    if val == NaN_i:
        val = val0
        print(f"  {valname:>15s} = {val:<10d}  ######  DEFAULT VALUE IS USED  ######")
    else:
        print(f"  {valname:>15s} = {val:<10d}")
    return val


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
        print(f"\n Check !  {valname} is SPECIFIED but will NOT be USED. ")
        print("            Please COMMENT-OUT this line ")
        print("            or check this input is REALLY APPROPRIATE for your purpose !\n")
        exit_program(-1)


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
    suffixes = [["x", "xy", "xz"],
                ["yx", "y", "yz"],
                ["zx", "zy", "z"]]
    not_used_d(valname, JAll)
    for i1 in range(3):
        for i2 in range(3):
            not_used_d(f"{valname}{suffixes[i1][i2]}", J[i1, i2])


def not_used_i(valname: str, val: int) -> None:
    """Abort if an integer parameter is specified but will not be used.

    Parameters
    ----------
    valname : str
        Name of the variable.
    val : int
        Value to check (abort if not the sentinel 2147483647).
    """
    if val != NaN_i:
        print(f"\n Check !  {valname} is SPECIFIED but will NOT be USED. ")
        print("            Please COMMENT-OUT this line ")
        print("            or check this input is REALLY APPROPRIATE for your purpose !\n")
        exit_program(-1)


def required_val_i(valname: str, val: int) -> None:
    """Abort if a required integer parameter is missing.

    Parameters
    ----------
    valname : str
        Name of the variable.
    val : int
        Value to check (abort if equals the sentinel 2147483647).
    """
    if val == NaN_i:
        print(f"ERROR ! {valname} is NOT specified !")
        exit_program(-1)
    else:
        print(f"  {valname:>15s} = {val:<3d}")
