# StdFace (Python)

The `python/` directory contains a Python port of the C-based StdFace input generator.
It reads the same input files as the C version and produces solver-specific definition files
for HPhi, mVMC, UHF, and H-wave.

## Requirements

- Python 3.10 or later
- NumPy

## Directory Structure

```
python/
  __main__.py              # CLI entry point (port of dry.c)
  stdface_main.py          # Main logic (port of StdFace_main.c)
  stdface_vals.py          # StdIntList dataclass (port of StdFace_vals.h)
  stdface_model_util.py    # Shared utilities (port of StdFace_ModelUtil.c)
  version.py               # Version information
  keyword_parser.py        # Keyword parsing subsystem
  param_check.py           # Parameter validation utilities
  history/
    refactoring_log.md     # Refactoring change log
  lattice/                 # Lattice implementations
    __init__.py
    chain_lattice.py       # 1D chain lattice
    square_lattice.py      # 2D square lattice
    ladder.py              # 2-leg ladder lattice
    triangular_lattice.py  # 2D triangular lattice
    honeycomb_lattice.py   # 2D honeycomb lattice
    kagome.py              # 2D kagome lattice
    orthorhombic.py        # 3D orthorhombic lattice
    fc_ortho.py            # Face-centered orthorhombic lattice
    pyrochlore.py          # Pyrochlore lattice
    wannier90.py           # Wannier90 input reader
    boost_output.py        # Boost output utilities
    geometry_output.py     # Geometry output functions
    input_params.py        # Input parameter resolution
    interaction_builder.py # Interaction building utilities
    site_util.py           # Site utility functions
  writer/                  # Solver-specific writers
    __init__.py
    common_writer.py       # Common output functions
    hphi_writer.py         # HPhi-specific writer
    mvmc_writer.py         # mVMC-specific writer
    mvmc_variational.py    # mVMC variational functions
    interaction_writer.py  # Interaction file writer
    solver_writer.py       # Solver writer base classes
    export_wannier90.py    # Wannier90 format export
```

## Usage

### Command Line

From the project root:

```bash
PYTHONPATH=python python3 python/__main__.py stan.in
```

To select a solver:

```bash
PYTHONPATH=python python3 python/__main__.py stan.in --solver mVMC
PYTHONPATH=python python3 python/__main__.py stan.in --solver UHF
PYTHONPATH=python python3 python/__main__.py stan.in --solver HWAVE
```

The default solver is HPhi.

### Print Version

```bash
PYTHONPATH=python python3 python/__main__.py -v
```

### Calling from Python

```python
import sys
sys.path.insert(0, "python")

from stdface_main import stdface_main

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
- 1,252 unit tests covering all modules
- Tests for lattice implementations, writers, parsers, and utilities
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
PYTHONPATH=/path/to/StdFace/python python3 -c \
  "from stdface_main import stdface_main; stdface_main('stan.in', solver='HPhi')"

# Compare outputs against reference
diff /path/to/StdFace/test/hphi/lanczos_hubbard_square/ref/modpara.def modpara.def
```

**Integration Test Results:**
- 83/83 integration tests passing
- All solvers (HPhi, mVMC, UHF, H-wave) verified
- Byte-identical output confirmed for all test cases

## Code Quality

The Python implementation has been refactored from a direct C translation into idiomatic Python:

- **Modular design**: Code organized into logical subpackages (`lattice/`, `writer/`)
- **Python idioms**: Enums, dict dispatch, context managers, type hints
- **Reduced duplication**: Helper functions extracted to eliminate code duplication
- **Comprehensive testing**: 1,252 unit tests with full coverage
- **Documentation**: NumPy-style docstrings throughout

See `history/refactoring_log.md` for detailed refactoring history (Steps 1-77).

## Mapping to C Sources

| C Source File | Python File |
|---|---|
| `dry.c` | `__main__.py` |
| `StdFace_main.c` | `stdface_main.py` (refactored into multiple modules) |
| `StdFace_vals.h` | `stdface_vals.py` |
| `StdFace_ModelUtil.c/h` | `stdface_model_util.py` |
| `version.h` | `version.py` |
| `ChainLattice.c` | `lattice/chain_lattice.py` |
| `SquareLattice.c` | `lattice/square_lattice.py` |
| `Ladder.c` | `lattice/ladder.py` |
| `TriangularLattice.c` | `lattice/triangular_lattice.py` |
| `HoneycombLattice.c` | `lattice/honeycomb_lattice.py` |
| `Kagome.c` | `lattice/kagome.py` |
| `Orthorhombic.c` | `lattice/orthorhombic.py` |
| `FCOrtho.c` | `lattice/fc_ortho.py` |
| `Pyrochlore.c` | `lattice/pyrochlore.py` |
| `Wannier90.c` | `lattice/wannier90.py` |
| `export_wannier90.c` | `writer/export_wannier90.py` |

Note: The Python implementation has been significantly refactored beyond the original C structure for better maintainability and Pythonic code style.

## License

GNU General Public License v3 (GPLv3)
