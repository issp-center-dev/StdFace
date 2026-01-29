"""Input parameter resolution helpers for spin, Coulomb, and hopping.

This module provides functions that resolve user-specified input parameters
into concrete values, handling conflicts between isotropic and anisotropic
specifications.  These are used by all lattice modules when setting up
Hamiltonian parameters.

Functions
---------
input_spin_nn
    Resolve nearest-neighbour spin-spin interaction (handles J/J0 conflicts).
input_spin
    Resolve spin-spin interaction for beyond-nearest-neighbour bonds.
input_coulomb_v
    Resolve off-site Coulomb interaction.
input_hopp
    Resolve hopping integral.

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
import math

import numpy as np

from param_check import exit_program, SPIN_SUFFIXES



# ---------------------------------------------------------------------------
#  Vectorised conflict-detection helpers
# ---------------------------------------------------------------------------


def _has_set_elements(mat: np.ndarray) -> bool:
    """Return True if *mat* contains any non-NaN elements.

    Parameters
    ----------
    mat : numpy.ndarray
        A 3×3 float array.

    Returns
    -------
    bool
        True if at least one element is finite (not NaN).
    """
    return bool(np.any(~np.isnan(mat)))


def _first_set_index(mat: np.ndarray) -> tuple[int, int] | None:
    """Return the ``(row, col)`` of the first non-NaN element, or None.

    Parameters
    ----------
    mat : numpy.ndarray
        A 3×3 float array.

    Returns
    -------
    tuple[int, int] or None
        Index of the first set element, or None if all NaN.
    """
    mask = ~np.isnan(mat)
    if not np.any(mask):
        return None
    idx = np.argmax(mask)  # index of first True in flattened array
    return int(idx // 3), int(idx % 3)


def _check_scalar_vs_matrix(scalar: float, mat: np.ndarray,
                             scalar_name: str, mat_name: str) -> None:
    """Abort if both a scalar and any matrix element are set.

    Parameters
    ----------
    scalar : float
        Isotropic scalar value (NaN means unset).
    mat : numpy.ndarray
        3×3 anisotropic matrix (NaN elements are unset).
    scalar_name : str
        Display name for the scalar in error messages.
    mat_name : str
        Base name for the matrix in error messages.
    """
    if math.isnan(scalar):
        return
    idx = _first_set_index(mat)
    if idx is not None:
        i1, i2 = idx
        print(f"\n ERROR! {scalar_name} and {mat_name}{SPIN_SUFFIXES[i1][i2]} conflict !\n")
        exit_program(-1)


def _check_matrix_vs_matrix(mat_a: np.ndarray, mat_b: np.ndarray,
                              name_a: str, name_b: str) -> None:
    """Abort if any element in *mat_a* and any element in *mat_b* are both set.

    Parameters
    ----------
    mat_a : numpy.ndarray
        First 3×3 matrix (NaN elements are unset).
    mat_b : numpy.ndarray
        Second 3×3 matrix (NaN elements are unset).
    name_a : str
        Base name for *mat_a* in error messages.
    name_b : str
        Base name for *mat_b* in error messages.
    """
    idx_a = _first_set_index(mat_a)
    idx_b = _first_set_index(mat_b)
    if idx_a is not None and idx_b is not None:
        i1, i2 = idx_a
        i3, i4 = idx_b
        print(f"\n ERROR! {name_a}{SPIN_SUFFIXES[i1][i2]} "
              f"and {name_b}{SPIN_SUFFIXES[i3][i4]} conflict !\n")
        exit_program(-1)


def _resolve_spin_matrix(
    J0: np.ndarray,
    J0name: str,
    *,
    J: np.ndarray | None = None,
    J0All: float = float("nan"),
    JAll: float = float("nan"),
) -> None:
    """Fill a 3×3 spin interaction matrix in-place using the resolution cascade.

    Priority order for each element ``(i, j)``:
    1. Explicit ``J0[i, j]`` (already set)
    2. Fallback ``J[i, j]`` (if provided and set)
    3. Diagonal ``J0All`` (if ``i == j``)
    4. Diagonal ``JAll`` (if ``i == j``)
    5. Default ``0.0``

    Parameters
    ----------
    J0 : numpy.ndarray
        Output 3×3 matrix (modified in-place).
    J0name : str
        Display name for the interaction.
    J : numpy.ndarray or None, optional
        Fallback 3×3 matrix. Default is None (no fallback).
    J0All : float, optional
        Isotropic value for this bond. Default is NaN (unset).
    JAll : float, optional
        Global isotropic value. Default is NaN (unset).
    """
    for i1, i2 in itertools.product(range(3), repeat=2):
        resolved = False
        if not math.isnan(J0[i1, i2]):
            resolved = True
        elif J is not None and not math.isnan(J[i1, i2]):
            J0[i1, i2] = J[i1, i2]
            resolved = True
        elif i1 == i2 and not math.isnan(J0All):
            J0[i1, i2] = J0All
            resolved = True
        elif i1 == i2 and not math.isnan(JAll):
            J0[i1, i2] = JAll
            resolved = True
        else:
            J0[i1, i2] = 0.0

        if resolved:
            label = J0name + SPIN_SUFFIXES[i1][i2]
            print(f"  {label:>14s} = {J0[i1, i2]:<10.5f}")


def input_spin_nn(
    J: np.ndarray,
    JAll: float,
    J0: np.ndarray,
    J0All: float,
    J0name: str,
) -> None:
    """Input nearest-neighbour spin-spin interaction.

    Resolves conflicts between isotropic and anisotropic specifications
    and fills in the output matrix *J0* in-place.

    Parameters
    ----------
    J : numpy.ndarray
        Anisotropic spin interaction (3x3).
    JAll : float
        Isotropic interaction.
    J0 : numpy.ndarray
        Output anisotropic spin interaction (3x3, modified in-place).
    J0All : float
        Isotropic interaction for this bond.
    J0name : str
        Name of this spin interaction (e.g. ``"J1"``).
    """
    # Scalar-scalar conflict: JAll vs J0All
    if not math.isnan(JAll) and not math.isnan(J0All):
        print(f"\n ERROR! {J0name} conflict !\n")
        exit_program(-1)

    # Scalar-matrix conflicts
    _check_scalar_vs_matrix(JAll, J, "J", "J")
    _check_scalar_vs_matrix(J0All, J, J0name, "J")
    _check_scalar_vs_matrix(J0All, J0, J0name, J0name)
    _check_scalar_vs_matrix(JAll, J0, J0name, J0name)

    # Cross-matrix conflict: any J0 element vs any J element
    _check_matrix_vs_matrix(J0, J, J0name, "J")

    # Resolve values using the cascade
    _resolve_spin_matrix(J0, J0name, J=J, J0All=J0All, JAll=JAll)


def input_spin(Jp: np.ndarray, JpAll: float, Jpname: str) -> None:
    """Input spin-spin interaction other than nearest-neighbour.

    Parameters
    ----------
    Jp : numpy.ndarray
        Fully anisotropic spin interaction (3x3, modified in-place).
    JpAll : float
        Isotropic interaction value.
    Jpname : str
        Name of this spin interaction (e.g. ``"J'"``).
    """
    # Scalar-matrix conflict
    _check_scalar_vs_matrix(JpAll, Jp, Jpname, Jpname)

    # Resolve values
    _resolve_spin_matrix(Jp, Jpname, J0All=JpAll)


def input_coulomb_v(V: float, V0: float, V0name: str) -> float:
    """Input off-site Coulomb interaction from the input file.

    Parameters
    ----------
    V : float
        Isotropic Coulomb interaction (may be NaN).
    V0 : float
        Specific Coulomb parameter (may be NaN).
    V0name : str
        Name of the parameter (e.g. ``"V1"``).

    Returns
    -------
    float
        The resolved value of *V0*.
    """
    if not math.isnan(V) and not math.isnan(V0):
        print(f"\n ERROR! {V0name} conflicts !\n")
        exit_program(-1)
    elif not math.isnan(V0):
        print(f"  {V0name:>15s} = {V0:<10.5f}")
    elif not math.isnan(V):
        V0 = V
        print(f"  {V0name:>15s} = {V0:<10.5f}")
    else:
        V0 = 0.0
    return V0


def input_hopp(t: complex, t0: complex, t0name: str) -> complex:
    """Input hopping integral from the input file.

    Parameters
    ----------
    t : complex
        Isotropic hopping (may have NaN real part).
    t0 : complex
        Specific hopping parameter (may have NaN real part).
    t0name : str
        Name of the parameter (e.g. ``"t1"``).

    Returns
    -------
    complex
        The resolved value of *t0*.
    """
    if not math.isnan(t.real) and not math.isnan(t0.real):
        print(f"\n ERROR! {t0name} conflicts !\n")
        exit_program(-1)
    elif not math.isnan(t0.real):
        print(f"  {t0name:>15s} = {t0.real:<10.5f} {t0.imag:<10.5f}")
    elif not math.isnan(t.real):
        t0 = t
        print(f"  {t0name:>15s} = {t0.real:<10.5f} {t0.imag:<10.5f}")
    else:
        t0 = 0.0 + 0j
    return t0
