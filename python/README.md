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
  __main__.py            # CLI entry point (port of dry.c)
  stdface_main.py        # Main logic (port of StdFace_main.c)
  stdface_vals.py        # StdIntList dataclass (port of StdFace_vals.h)
  stdface_model_util.py  # Shared utilities (port of StdFace_ModelUtil.c)
  version.py             # Version information
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
python3 -m pytest test/unit/ -v
```

### Integration Tests (Single Test Case)

```bash
# Copy test input to a working directory and run
cp test/hphi/lanczos_hubbard_square/stan.in /tmp/work/
cd /tmp/work
PYTHONPATH=/path/to/StdFace/python python3 -c \
  "from stdface_main import stdface_main; stdface_main('stan.in', solver='HPhi')"

# Compare outputs against reference
diff /path/to/StdFace/test/hphi/lanczos_hubbard_square/ref/modpara.def modpara.def
```

## Mapping to C Sources

| C Source File | Python File |
|---|---|
| `dry.c` | `__main__.py` |
| `StdFace_main.c` | `stdface_main.py` |
| `StdFace_vals.h` | `stdface_vals.py` |
| `StdFace_ModelUtil.c/h` | `stdface_model_util.py` |
| `version.h` | `version.py` |
| `ChainLattice.c` | `chain_lattice.py` |
| `SquareLattice.c` | `square_lattice.py` |
| `Ladder.c` | `ladder.py` |
| `TriangularLattice.c` | `triangular_lattice.py` |
| `HoneycombLattice.c` | `honeycomb_lattice.py` |
| `Kagome.c` | `kagome.py` |
| `Orthorhombic.c` | `orthorhombic.py` |
| `FCOrtho.c` | `fc_ortho.py` |
| `Pyrochlore.c` | `pyrochlore.py` |
| `Wannier90.c` | `wannier90.py` |
| `export_wannier90.c` | `export_wannier90.py` |

## License

GNU General Public License v3 (GPLv3)
