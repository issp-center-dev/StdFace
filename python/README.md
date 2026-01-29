# StdFace (Python)

The `python/` directory contains a Python port of the C-based StdFace input generator.
It reads the same input files as the C version and produces solver-specific definition files
for HPhi, mVMC, UHF, and H-wave.

## Requirements

- Python 3.10 or later
- NumPy

## Installation

```bash
cd python
pip install -e .          # install with runtime dependencies (numpy)
pip install -e ".[dev]"   # also install dev dependencies (pytest, pytest-cov)
```

After installation, the `stdface` command becomes available:

```bash
stdface stan.in
stdface stan.in --solver mVMC
stdface -v
```

The editable install (`-e`) is recommended for development. Source changes take effect immediately without reinstalling.

## Architecture

StdFace uses a **plugin architecture** for both solvers and lattices.
New solvers and lattices can be added without modifying any existing dispatch logic.

### Plugin System Overview

```
stdface/
  plugin.py              # SolverPlugin ABC + solver registry
  lattice/__init__.py    # LatticePlugin ABC + lattice registry
```

- **SolverPlugin** (`plugin.py`): Defines how a solver parses keywords, resets fields, and writes output files. Each solver (HPhi, mVMC, UHF, H-wave) is a plugin registered at import time.
- **LatticePlugin** (`lattice/__init__.py`): Defines lattice geometry, aliases, and the setup/boost methods. Each lattice (chain, square, kagome, etc.) is a plugin registered at import time.

See [docs/tutorial_plugin.md](../docs/tutorial_plugin.md) for a step-by-step guide on adding new solvers and lattices.

## Directory Structure

```
python/
  __main__.py              # CLI entry point (port of dry.c)
  stdface/
    plugin.py              # SolverPlugin ABC + solver registry
    core/
      stdface_main.py      # Main logic (port of StdFace_main.c)
      stdface_vals.py      # StdIntList dataclass (port of StdFace_vals.h)
      keyword_parser.py    # Keyword parsing subsystem
      param_check.py       # Parameter validation utilities
    lattice/               # Lattice plugins
      __init__.py          # LatticePlugin ABC + lattice registry
      chain_lattice.py     # 1D chain (ChainPlugin)
      square_lattice.py    # 2D square (SquarePlugin)
      ladder.py            # 2-leg ladder (LadderPlugin)
      triangular_lattice.py # 2D triangular (TriangularPlugin)
      honeycomb_lattice.py # 2D honeycomb (HoneycombPlugin)
      kagome.py            # 2D kagome (KagomePlugin)
      orthorhombic.py      # 3D orthorhombic (OrthorhombicPlugin)
      fc_ortho.py          # 3D face-centered orthorhombic (FCOrthoPlugin)
      pyrochlore.py        # 3D pyrochlore (PyrochlorePlugin)
      wannier90.py         # Wannier90 interface (Wannier90Plugin)
      boost_output.py      # Boost output utilities
      geometry_output.py   # Geometry output functions
      input_params.py      # Input parameter resolution
      interaction_builder.py # Interaction building utilities
      site_util.py         # Site utility functions
    solvers/               # Solver plugins
      __init__.py          # Auto-imports all solver plugins
      hphi/_plugin.py      # HPhi plugin (HPhiPlugin)
      mvmc/_plugin.py      # mVMC plugin (MVMCPlugin)
      uhf/_plugin.py       # UHF plugin (UHFPlugin)
      hwave/_plugin.py     # H-wave plugin (HWavePlugin)
    writer/                # Output writers (shared)
      common_writer.py     # Common output functions
      interaction_writer.py # Interaction file writer
      export_wannier90.py  # Wannier90 format export
  history/
    refactoring_log.md     # Refactoring change log
```

## Usage

### Command Line

From the project root:

```bash
./stdface stan.in
```

To select a solver:

```bash
./stdface stan.in --solver mVMC
./stdface stan.in --solver UHF
./stdface stan.in --solver HWAVE
```

The default solver is HPhi.

### Print Version

```bash
./stdface -v
```

### Calling from Python

```python
import sys
sys.path.insert(0, "python")

from stdface.core.stdface_main import stdface_main

stdface_main("stan.in", solver="HPhi")
```

## Input File Examples

Parameters are specified in `stan.in` using the `key = value` format.

### Hubbard Model (Square Lattice)

```text
model = Hubbard
lattice = square
W = 2
L = 2
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
```

### Spin Model (Kagome Lattice)

```text
model = Spin
lattice = kagome
W = 2
L = 2
J = 1.0
2S = 1
2Sz = 0
```

### Kondo Lattice Model (1D Chain)

```text
model = Kondo
lattice = chain
L = 4
t = 1.0
U = 4.0
J = 0.8
nelec = 4
2Sz = 0
```

## Running Tests

### Unit Tests

```bash
# Run all unit tests
python3 -m pytest test/unit/ -v

# Run with coverage
python3 -m pytest test/unit/ --cov=python --cov-report=html
```

**Test Coverage:**
- 1,268 unit tests covering all modules
- Tests for lattice implementations, writers, parsers, plugin registries, and utilities
- All tests maintain byte-identical output with C version

### Integration Tests

**Run all integration tests:**
```bash
bash test/run_all_integration.sh
```

**Single test case:**
```bash
# Copy test input to a working directory and run
cp test/hphi/lanczos_hubbard_square/stan.in /tmp/work/
cd /tmp/work
/path/to/StdFace/stdface stan.in

# Compare outputs against reference
diff /path/to/StdFace/test/hphi/lanczos_hubbard_square/ref/modpara.def modpara.def
```

**Integration Test Results:**
- 83/83 integration tests passing
- All solvers (HPhi, mVMC, UHF, H-wave) verified
- Byte-identical output confirmed for all test cases

## Code Quality

The Python implementation has been refactored from a direct C translation into idiomatic Python:

- **Plugin architecture**: Solvers and lattices are self-registering plugins
- **Modular design**: Code organized into logical subpackages (`lattice/`, `solvers/`, `writer/`)
- **Python idioms**: Enums, ABC, context managers, type hints
- **Reduced duplication**: Helper functions extracted to eliminate code duplication
- **Comprehensive testing**: 1,268 unit tests with full coverage
- **Documentation**: NumPy-style docstrings throughout

See `history/refactoring_log.md` for detailed refactoring history.

## Mapping to C Sources

| C Source File | Python Module |
|---|---|
| `dry.c` | `__main__.py` |
| `StdFace_main.c` | `stdface/core/stdface_main.py` |
| `StdFace_vals.h` | `stdface/core/stdface_vals.py` |
| `StdFace_ModelUtil.c/h` | `stdface/core/stdface_model_util.py` |
| `version.h` | `stdface/version.py` |
| `ChainLattice.c` | `stdface/lattice/chain_lattice.py` |
| `SquareLattice.c` | `stdface/lattice/square_lattice.py` |
| `Ladder.c` | `stdface/lattice/ladder.py` |
| `TriangularLattice.c` | `stdface/lattice/triangular_lattice.py` |
| `HoneycombLattice.c` | `stdface/lattice/honeycomb_lattice.py` |
| `Kagome.c` | `stdface/lattice/kagome.py` |
| `Orthorhombic.c` | `stdface/lattice/orthorhombic.py` |
| `FCOrtho.c` | `stdface/lattice/fc_ortho.py` |
| `Pyrochlore.c` | `stdface/lattice/pyrochlore.py` |
| `Wannier90.c` | `stdface/lattice/wannier90.py` |
| `export_wannier90.c` | `stdface/writer/export_wannier90.py` |

Note: The Python implementation has been significantly refactored beyond the original C structure for better maintainability and Pythonic code style.

## License

GNU General Public License v3 (GPLv3)
