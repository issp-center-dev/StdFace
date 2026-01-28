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

### C Implementation

- CMake 2.8.12 or later
- C99-compatible compiler (GCC, Clang, Intel, Fujitsu, etc.)

### Python Implementation

- Python 3.10 or later
- NumPy

## Installation

### C Implementation

```bash
git clone https://github.com/issp-center-dev/StdFace
cd StdFace
cmake -B build -DHPHI=ON   # Enable HPhi mode
cmake --build build
cmake --install build --prefix /path/to/install
```

#### Build Options

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

### Python Implementation

The Python implementation requires no installation. Simply ensure Python 3.10+ and NumPy are available:

```bash
# Check Python version
python3 --version  # Should be 3.10 or later

# Install NumPy if needed
pip install numpy
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

**C Implementation:**
```bash
./hphi_dry.out stan.in
```

**Python Implementation:**
```bash
PYTHONPATH=python python3 python/__main__.py stan.in
```

To select a solver (Python):
```bash
PYTHONPATH=python python3 python/__main__.py stan.in --solver mVMC
PYTHONPATH=python python3 python/__main__.py stan.in --solver UHF
PYTHONPATH=python python3 python/__main__.py stan.in --solver HWAVE
```

3. Input files for the target solver are generated in the current directory.

Both implementations produce identical output files.

## Python Implementation

The `python/` directory contains a fully-featured Python port of StdFace that produces byte-identical output to the C implementation. The Python codebase has been refactored into idiomatic Python with:

- **Modular architecture**: Organized into `lattice/` and `writer/` subpackages
- **Comprehensive testing**: 1,252 unit tests and 83 integration tests
- **Python idioms**: Enums, dict dispatch, context managers, and helper functions
- **Full feature parity**: Supports all lattices, models, and solvers

### Python Project Structure

```
python/
  __main__.py              # CLI entry point
  stdface_main.py          # Main logic
  stdface_vals.py          # Data structures
  stdface_model_util.py    # Shared utilities
  keyword_parser.py        # Keyword parsing
  param_check.py           # Parameter validation
  lattice/                 # Lattice implementations
    chain_lattice.py
    square_lattice.py
    honeycomb_lattice.py
    kagome.py
    wannier90.py
    ...
  writer/                  # Solver-specific writers
    common_writer.py
    hphi_writer.py
    mvmc_writer.py
    ...
```

### Running Python Tests

**Unit Tests:**
```bash
python3 -m pytest test/unit/ -v
```

**Integration Tests:**
```bash
# Run all integration tests
bash test/run_all_integration.sh
```

For more details, see [python/README.md](python/README.md).

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
