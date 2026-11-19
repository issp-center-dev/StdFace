"""
Version information for StdFace.

mVMC - A numerical solver package for a wide range of quantum lattice models
based on many-variable Variational Monte Carlo method Copyright (C) 2016 The
University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.
"""

from __future__ import annotations

# Version of the Python package (PEP 440).  It is independent of the C
# generator's version in src/version.h, and this is the only place it is
# written: pyproject.toml reads it from here.
__version__ = "0.6.0-dev"


def print_version() -> None:
    """Print the version of the Python implementation to standard output.

    The ``StdFace(python)`` label tells this apart from the C dry-run
    binaries, which print ``StdFace version <version>``.
    """
    print(f"StdFace(python) version {__version__}")
