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

# Semantic Versioning http://semver.org
# <major>.<minor>.<patch>-<prerelease>
VERSION_MAJOR = 0
VERSION_MINOR = 5
VERSION_PATCH = 0
VERSION_PRERELEASE = ""  # "alpha", "beta.1", etc.


def print_version():
    """Print the StdFace version string to standard output.

    Prints a version string in the format ``major.minor.patch``.
    If ``VERSION_PRERELEASE`` is non-empty, ``-<prerelease>`` is appended.
    """
    version_str = f"StdFace version {VERSION_MAJOR}.{VERSION_MINOR}.{VERSION_PATCH}"
    if VERSION_PRERELEASE:
        version_str += f"-{VERSION_PRERELEASE}"
    print(version_str)
