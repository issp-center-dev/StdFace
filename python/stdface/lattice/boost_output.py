"""Boost-output helpers for writing ``boost.def`` sections.

This module provides functions used by lattice builders that support HPhi's
Boost method (chain, honeycomb, kagome, ladder).  Each function writes a
specific section of the ``boost.def`` file.

Functions
---------
write_boost_mag_field
    Write the magnetic-field line to ``boost.def``.
write_boost_j_full
    Write a full 3×3 J-coupling matrix to ``boost.def``.
write_boost_j_symmetric
    Write an upper-triangle symmetric 3×3 J-coupling matrix to ``boost.def``.
write_boost_6spin_star
    Write the ``list_6spin_star`` array to ``boost.def``.
write_boost_6spin_pair
    Write the ``list_6spin_pair`` array to ``boost.def``.

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

from typing import TextIO

import numpy as np

from ..core.stdface_vals import StdIntList


def write_boost_mag_field(fp: TextIO, StdI: StdIntList) -> None:
    """Write the magnetic-field header line to ``boost.def``.

    Outputs ``# Magnetic field`` followed by a line containing the three
    scaled field components ``-0.5 * Gamma``, ``-0.5 * Gamma_y``, and
    ``-0.5 * h``.

    Parameters
    ----------
    fp : TextIO
        Open file handle for ``boost.def``.
    StdI : StdIntList
        Model parameter structure (read-only access to ``Gamma``,
        ``Gamma_y``, ``h``).
    """
    fp.write("# Magnetic field\n")
    fp.write(
        f"{-0.5 * StdI.Gamma:25.15e} "
        f"{-0.5 * StdI.Gamma_y:25.15e} "
        f"{-0.5 * StdI.h:25.15e}\n"
    )


def write_boost_j_full(fp: TextIO, J: np.ndarray, scale: float = 0.25) -> None:
    """Write a full 3×3 J-coupling matrix to ``boost.def``.

    Outputs three rows of three values each, where each element is
    ``scale * J[i, j]``.  Used by the chain lattice.

    Parameters
    ----------
    fp : TextIO
        Open file handle for ``boost.def``.
    J : numpy.ndarray
        3×3 coupling matrix.
    scale : float, optional
        Multiplicative scaling factor.  Default is ``0.25``.
    """
    for row in J:
        scaled = scale * row
        fp.write(f"{scaled[0]:25.15e} {scaled[1]:25.15e} {scaled[2]:25.15e}\n")


def write_boost_j_symmetric(fp: TextIO, J: np.ndarray, scale: float = 0.25) -> None:
    """Write an upper-triangle symmetric 3×3 J-coupling matrix to ``boost.def``.

    Outputs three rows using only upper-triangle elements::

        J[0,0]  J[0,1]  J[0,2]
        J[0,1]  J[1,1]  J[1,2]
        J[0,2]  J[1,2]  J[2,2]

    Each element is multiplied by *scale*.  Used by honeycomb, ladder,
    and kagome lattices.

    Parameters
    ----------
    fp : TextIO
        Open file handle for ``boost.def``.
    J : numpy.ndarray
        3×3 coupling matrix (only upper-triangle elements are read).
    scale : float, optional
        Multiplicative scaling factor.  Default is ``0.25``.
    """
    Jsym = np.triu(J) + np.triu(J, 1).T
    write_boost_j_full(fp, Jsym, scale)


def write_boost_6spin_star(fp: TextIO, StdI: StdIntList) -> None:
    """Write the ``list_6spin_star`` array to ``boost.def``.

    Outputs each pivot's 7-element row from ``StdI.list_6spin_star``.

    Parameters
    ----------
    fp : TextIO
        Open file handle for ``boost.def``.
    StdI : StdIntList
        Structure containing ``num_pivot`` and ``list_6spin_star``.
    """
    fp.write("# StdI->list_6spin_star\n")
    for ipivot in range(StdI.num_pivot):
        fp.write(f"# pivot {ipivot}\n")
        row = StdI.list_6spin_star[ipivot, :7]
        fp.write(" ".join(str(x) for x in row) + " \n")


def write_boost_6spin_pair(fp: TextIO, StdI: StdIntList) -> None:
    """Write the ``list_6spin_pair`` array to ``boost.def``.

    Outputs each pivot's interaction rows from ``StdI.list_6spin_pair``,
    with the number of interactions per pivot given by
    ``StdI.list_6spin_star[ipivot, 0]``.

    Parameters
    ----------
    fp : TextIO
        Open file handle for ``boost.def``.
    StdI : StdIntList
        Structure containing ``num_pivot``, ``list_6spin_star``,
        and ``list_6spin_pair``.
    """
    fp.write("# StdI->list_6spin_pair\n")
    for ipivot in range(StdI.num_pivot):
        fp.write(f"# pivot {ipivot}\n")
        for kintr in range(StdI.list_6spin_star[ipivot, 0]):
            row = StdI.list_6spin_pair[ipivot, :7, kintr]
            fp.write(" ".join(str(x) for x in row) + " \n")
