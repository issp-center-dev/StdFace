# StdFace

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://issp-center-dev.github.io/StdFace/)
[![Build Documentation](https://github.com/issp-center-dev/StdFace/actions/workflows/docs.yml/badge.svg)](https://github.com/issp-center-dev/StdFace/actions/workflows/docs.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

An input file generator for quantum lattice model solvers developed at ISSP, University of Tokyo.

## Overview

StdFace reads a simple configuration file specifying the physical model and lattice geometry, then generates solver-specific input files for:

- [HPhi](https://github.com/issp-center-dev/HPhi) - Exact Diagonalization
- [mVMC](https://github.com/issp-center-dev/mVMC) - Variational Monte Carlo
- [UHF](https://github.com/issp-center-dev/UHF-dev) - Unrestricted Hartree-Fock
- [H-wave](https://github.com/issp-center-dev/H-wave) - Mean-field solver

## Features

- Simple `key = value` input format
- Support for various lattice geometries and physical models
- Wannier90 format support for ab-initio calculations
- Single codebase for multiple solver backends

## Supported Lattices

| Lattice Type | Description |
|--------------|-------------|
| `chain` | 1D chain |
| `ladder` | 2-leg ladder |
| `square` | 2D square lattice |
| `triangular` | 2D triangular lattice |
| `honeycomb` | 2D honeycomb lattice |
| `kagome` | 2D kagome lattice |
| `tetragonal` | 3D tetragonal lattice |
| `orthorhombic` | 3D orthorhombic lattice |
| `fcortho` | Face-centered orthorhombic |
| `pyrochlore` | 3D pyrochlore lattice |
| `wannier90` | Import from Wannier90 files |

## Supported Models

- Hubbard model
- Spin models
- Kondo lattice model

## Requirements

- CMake 2.8.12 or later
- C99-compatible compiler (GCC, Clang, Intel, Fujitsu, etc.)

## Installation

```bash
git clone https://github.com/issp-center-dev/StdFace
cd StdFace
cmake -B build -DHPHI=ON   # Enable HPhi mode
cmake --build build
cmake --install build --prefix /path/to/install
```

### Build Options

Enable one or more solver modes:

| Option | Executable | Target Solver |
|--------|------------|---------------|
| `-DHPHI=ON` | `hphi_dry.out` | HPhi (Exact Diagonalization) |
| `-DMVMC=ON` | `mvmc_dry.out` | mVMC (Variational Monte Carlo) |
| `-DUHF=ON` | `uhf_dry.out` | UHF (Unrestricted Hartree-Fock) |
| `-DHWAVE=ON` | `hwave_dry.out` | H-wave |

Build all solvers:

```bash
cmake -B build -DHPHI=ON -DMVMC=ON -DUHF=ON -DHWAVE=ON
cmake --build build
```

## Quick Start

1. Create an input file `stan.in`:

```text
model = "Hubbard"
lattice = square
W = 2
L = 2
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
```

2. Run StdFace:

```bash
./hphi_dry.out stan.in
```

3. Input files for the target solver are generated in the current directory.

## Documentation

Full documentation is available at: https://issp-center-dev.github.io/StdFace/

- [Quickstart Guide](https://issp-center-dev.github.io/StdFace/user/quickstart.html)
- [Input/Output Reference](https://issp-center-dev.github.io/StdFace/user/input_output.html)
- [Examples](https://issp-center-dev.github.io/StdFace/user/examples.html)

## License

StdFace is distributed under the [GNU General Public License version 3 (GPL v3)](http://www.gnu.org/licenses/gpl-3.0.en.html).

## Authors

Kazuyoshi Yoshimi, Mitsuaki Kawamura, Kota Ido, Yuichi Motoyama, and Tatsumi Aoyama.

## Related Projects

- [HPhi](https://github.com/issp-center-dev/HPhi) - Exact Diagonalization package
- [mVMC](https://github.com/issp-center-dev/mVMC) - Variational Monte Carlo package
- [UHF](https://github.com/issp-center-dev/UHF-dev) - Unrestricted Hartree-Fock package
- [H-wave](https://github.com/issp-center-dev/H-wave) - Mean-field solver
