#!/usr/bin/env python3
"""Generate matrix test cases covering all (solver, lattice, model) combinations.

Usage:
    python3 gen_matrix_tests.py [--gen-ref] [--update-cmake]
                                [--only REGEX] [--build-dir DIR]

    --gen-ref       Run dry executables to populate ref/ directories.
    --update-cmake  Rewrite the matrix block of each solver's CMakeLists.txt.
    --only REGEX    With --gen-ref, regenerate only the cases whose name matches.
    --build-dir DIR Directory holding the *_dry.out executables
                    (default: ../build/src relative to this script).

References are the C generators' output, except for the H-wave uhfk cases:
their Wannier-format ``*.dat`` files come from the Python CLI (``stdface``
must be on PATH), which normalizes the small-supercell fold-back that the C
generator leaves in.  Use --only when regenerating: some references outside
the cases touched are kept from the Python implementation on purpose.
"""

import argparse
import re
import shutil
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()
BUILD_DIR = SCRIPT_DIR.parent / "build" / "src"

# ---------------------------------------------------------------------------
# Lattice definitions
# params: geometry keywords written verbatim into stan.in
# nsite_uc: number of sites per unit cell (used to compute nelec/ncond)
# ncell: product of L*W*Height
# ---------------------------------------------------------------------------
LATTICES = {
    "chain":        {"params": {"L": 4},                  "nsite_uc": 1},
    "square":       {"params": {"L": 2, "W": 2},          "nsite_uc": 1},
    "triangular":   {"params": {"L": 2, "W": 2},          "nsite_uc": 1},
    "honeycomb":    {"params": {"L": 2, "W": 2},          "nsite_uc": 2},
    "kagome":       {"params": {"L": 2, "W": 2},          "nsite_uc": 3},
    # ladder: NsiteUC = W (set by Ladder.c), then W = 1 → nsite = W * L
    # W must be explicitly specified; use nsite_uc=1 because ncell = L*W.
    "ladder":       {"params": {"L": 4, "W": 2},          "nsite_uc": 1},
    # The cell height is "Height": "H" is the magnetic field keyword.
    "orthorhombic": {"params": {"L": 2, "W": 2, "Height": 2},  "nsite_uc": 1},
    "pyrochlore":   {"params": {"L": 1, "W": 1, "Height": 1},  "nsite_uc": 4},
    "fco":          {"params": {"L": 2, "W": 2, "Height": 2},  "nsite_uc": 1},
}

# Magnetic-field variants: the 3D lattices with a longitudinal and a
# transverse field. Grand-canonical models only (Gamma does not conserve Sz).
FIELD_LATTICES = ()
FIELD_PARAMS = {"h": 0.5, "Gamma": 0.3}


def nsite_lattice(lat_name: str) -> int:
    lat = LATTICES[lat_name]
    p = lat["params"]
    ncell = p.get("L", 1) * p.get("W", 1) * p.get("Height", 1)
    return ncell * lat["nsite_uc"]


# Combinations that are known to be unsupported at the source level.
# Each entry is (solver, calcmode, base_model); None matches any value.
UNSUPPORTED = [
    # HWave uhfk + kondo: segfaults (not implemented in the uhfk path)
    ("hwave", "uhfk", "kondo"),
]


def is_unsupported(solver: str, calcmode, base: str) -> bool:
    for s, c, b in UNSUPPORTED:
        if (s is None or s == solver) and (c is None or c == calcmode) and (b is None or b == base):
            return True
    return False


# Models: (user-facing name, is_gc, base_model)
MODELS = [
    ("hubbard",   False, "hubbard"),
    ("hubbardgc", True,  "hubbard"),
    ("spin",      False, "spin"),
    ("spingc",    True,  "spin"),
    ("kondo",     False, "kondo"),
    ("kondogc",   True,  "kondo"),
]

# Solver configuration.
# "calcmodes": list of calcmode values to generate; None means no calcmode parameter.
SOLVERS = {
    "hphi": {
        "exe":       "hphi_dry.out",
        "infile":    "stan.in",
        "calcmodes": [None],
    },
    "mvmc": {
        "exe":       "mvmc_dry.out",
        "infile":    "StdFace.def",
        "calcmodes": [None],
    },
    "uhf": {
        "exe":       "uhf_dry.out",
        "infile":    "stan.in",
        "calcmodes": [None],
    },
    "hwave": {
        "exe":       "hwave_dry.out",
        "infile":    "stan.in",
        "calcmodes": ["uhfr", "uhfk"],
    },
}

# Solver names of the Python CLI (stdface --solver <name>)
PYTHON_SOLVER = {"hphi": "HPhi", "mvmc": "mVMC", "uhf": "UHF", "hwave": "HWAVE"}


def model_input_name(model_name: str) -> str:
    """Return the string written to the model = field in stan.in."""
    table = {
        "hubbard":   "FermionHubbard",
        "hubbardgc": "FermionHubbardGC",
        "spin":      "Spin",
        "spingc":    "SpinGC",
        "kondo":     "Kondo",
        "kondogc":   "KondoGC",
    }
    return table[model_name]


def build_stan_in(solver: str, lat_name: str, model_name: str, is_gc: bool, base: str,
                  calcmode=None, field: bool = False) -> str:
    """Return the content of stan.in (or StdFace.def) for the given combination."""
    lat = LATTICES[lat_name]
    ns = nsite_lattice(lat_name)

    lines = []

    # Geometry
    for key, val in lat["params"].items():
        lines.append(f"{key} = {val}")

    lines.append(f'model = "{model_input_name(model_name)}"')
    lines.append(f'lattice = "{lat_name}"')

    # Physics parameters
    # Ladder uses bond-specific parameter names (t0/t1, J0/J1) instead of t/J.
    if lat_name == "ladder":
        if base == "hubbard":
            lines.append("t0 = 1.0")
            lines.append("t1 = 1.0")
            lines.append("U = 4.0")
        elif base == "spin":
            lines.append("J0 = 1.0")
            lines.append("J1 = 1.0")
        elif base == "kondo":
            lines.append("t0 = 1.0")
            lines.append("t1 = 1.0")
            # J is marked NotUsed for ladder before model dispatch → cannot specify
            lines.append("U = 4.0")
    else:
        if base == "hubbard":
            lines.append("t = 1.0")
            lines.append("U = 4.0")
        elif base == "spin":
            lines.append("J = 1.0")
        elif base == "kondo":
            lines.append("t = 1.0")
            lines.append("J = 1.0")
            lines.append("U = 4.0")

    if field:
        for key, val in FIELD_PARAMS.items():
            lines.append(f"{key} = {val}")

    # Method / calcmode
    if solver == "hphi":
        lines.append('method = "FullDiag"')
    elif calcmode is not None:
        lines.append(f'calcmode = "{calcmode}"')

    # Conserved-quantity parameters.
    #
    # HPhi hubbard canonical:  nelec required
    # HPhi hubbardgc:          nelec/2Sz are NotUsed → must NOT specify
    # mVMC/UHF/HWave hubbard: ncond always required (canonical and GC)
    #
    # HPhi kondo canonical:    ncond required
    # HPhi kondogc:            ncond is NotUsed → must NOT specify
    # mVMC/UHF/HWave kondo:   ncond always required (canonical and GC)
    #
    # spin (all solvers):      2Sz required for canonical, NotUsed for GC
    if base == "hubbard":
        if solver == "hphi":
            if not is_gc:
                lines.append(f"nelec = {ns}")
                lines.append("2Sz = 0")
        else:
            lines.append(f"nelec = {ns}")
            if not is_gc:
                lines.append("2Sz = 0")
    elif base == "spin":
        if not is_gc:
            lines.append("2Sz = 0")
    elif base == "kondo":
        if solver == "hphi":
            if not is_gc:
                lines.append(f"ncond = {ns}")
                lines.append("2Sz = 0")
        else:
            lines.append(f"ncond = {ns}")
            if not is_gc:
                lines.append("2Sz = 0")

    # Solver-specific extras for mVMC reproducibility
    if solver == "mvmc":
        lines.append("RndSeed = 1")

    return "\n".join(lines) + "\n"


def test_name(lat_name: str, model_name: str, calcmode=None, field: bool = False) -> str:
    name = f"matrix_{lat_name}_{model_name}"
    if field:
        name += "_field"
    if calcmode is not None:
        name += f"_{calcmode}"
    return name


def cases(solver_name: str):
    """Yield (name, lattice, model, is_gc, base, calcmode, field) for a solver.

    The field variants come after all plain cases, so adding them leaves the
    order of the existing CMake entries unchanged.
    """
    calcmodes = SOLVERS[solver_name]["calcmodes"]
    for lat_name in LATTICES:
        for model_name, is_gc, base in MODELS:
            for calcmode in calcmodes:
                if not is_unsupported(solver_name, calcmode, base):
                    yield (test_name(lat_name, model_name, calcmode),
                           lat_name, model_name, is_gc, base, calcmode, False)
    for lat_name in FIELD_LATTICES:
        for model_name, is_gc, base in MODELS:
            # uhfk is left out: its small-cell references are Python-only.
            for calcmode in [c for c in calcmodes if c != "uhfk"]:
                if is_gc and not is_unsupported(solver_name, calcmode, base):
                    yield (test_name(lat_name, model_name, calcmode, field=True),
                           lat_name, model_name, is_gc, base, calcmode, True)


def generate_test_dirs(dry_run: bool = False) -> None:
    """Create all test directories with input files."""
    for solver_name, solver_cfg in SOLVERS.items():
        solver_dir = SCRIPT_DIR / solver_name
        infile = solver_cfg["infile"]
        for name, lat_name, model_name, is_gc, base, calcmode, field in cases(solver_name):
            if dry_run:
                print(f"  [dry] {solver_name}/{name}")
                continue
            test_dir = solver_dir / name
            test_dir.mkdir(parents=True, exist_ok=True)
            content = build_stan_in(solver_name, lat_name, model_name, is_gc, base, calcmode, field)
            (test_dir / infile).write_text(content)
            print(f"  wrote {solver_name}/{name}/{infile}")


def generate_refs(build_dir: Path, only=None) -> None:
    """Run each dry executable inside ref/ to produce reference output files."""
    for solver_name, solver_cfg in SOLVERS.items():
        exe_path = (build_dir / solver_cfg["exe"]).resolve()
        if not exe_path.exists():
            print(f"WARNING: {exe_path} not found, skipping {solver_name}")
            continue

        solver_dir = SCRIPT_DIR / solver_name
        infile = solver_cfg["infile"]
        python_cli = ["stdface", "--solver", PYTHON_SOLVER[solver_name]]

        for name, _lat, _model, _is_gc, _base, calcmode, _field in cases(solver_name):
            if only is not None and not re.search(only, name):
                continue
            test_dir = solver_dir / name

            if not (test_dir / infile).exists():
                print(f"  SKIP {solver_name}/{name} (no input file)")
                continue

            ref_dir = test_dir / "ref"
            if ref_dir.exists():
                shutil.rmtree(ref_dir)
            ref_dir.mkdir()

            result = subprocess.run(
                [str(exe_path), f"../{infile}"],
                cwd=str(ref_dir),
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                print(f"  FAIL {solver_name}/{name}: {result.stderr.strip()}")
                continue
            if calcmode == "uhfk":
                # Wannier-format files from the Python CLI (see module docstring)
                py_dir = test_dir / "ref.python"
                if py_dir.exists():
                    shutil.rmtree(py_dir)
                py_dir.mkdir()
                result = subprocess.run(python_cli + [f"../{infile}"], cwd=str(py_dir),
                                        capture_output=True, text=True)
                if result.returncode != 0:
                    print(f"  FAIL {solver_name}/{name} (python): {result.stderr.strip()}")
                    shutil.rmtree(py_dir)
                    continue
                for dat in py_dir.glob("*.dat"):
                    shutil.copy(dat, ref_dir / dat.name)
                shutil.rmtree(py_dir)
            print(f"  OK   {solver_name}/{name}")


def cmake_entries(solver_name: str) -> str:
    """Return CMakeLists.txt add_test() entries for all matrix tests of a solver."""
    lines = []
    for name, *_ in cases(solver_name):
        lines.append(
            f"add_test(\n"
            f"  NAME {name}\n"
            f"  COMMAND ${{CMAKE_SOURCE_DIR}}/test/{solver_name}/check.sh"
            f" {name} ${{CMAKE_SOURCE_DIR}}/test/{solver_name} ${{TEST_MODE}}\n"
            f")"
        )
    return "\n".join(lines) + "\n"


_MATRIX_ENTRY = re.compile(r"add_test\(\n  NAME matrix_[^\n]*\n  COMMAND [^\n]*\n\)\n")


def update_cmake() -> None:
    """Rewrite the matrix block of each solver's CMakeLists.txt.

    The block is the marker line and the matrix_* entries that follow it;
    whatever comes after them (other test groups) is kept.
    """
    marker = "# --- matrix tests (auto-generated) ---\n"
    for solver_name in SOLVERS:
        cmake_path = SCRIPT_DIR / solver_name / "CMakeLists.txt"
        original = cmake_path.read_text()

        if marker in original:
            start = original.index(marker)
            end = start + len(marker)
            while (m := _MATRIX_ENTRY.match(original, end)):
                end = m.end()
            head, tail = original[:start], original[end:]
        else:
            head, tail = original.rstrip("\n") + "\n\n", ""

        cmake_path.write_text(head + marker + cmake_entries(solver_name) + tail)
        print(f"  updated {solver_name}/CMakeLists.txt")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gen-ref", action="store_true")
    parser.add_argument("--update-cmake", action="store_true")
    parser.add_argument("--only", metavar="REGEX")
    parser.add_argument("--build-dir", type=Path, default=BUILD_DIR)
    args = parser.parse_args()

    print("==> Generating test directories...")
    generate_test_dirs()

    if args.gen_ref:
        print("==> Generating ref/ files...")
        generate_refs(args.build_dir, args.only)

    if args.update_cmake:
        print("==> Updating CMakeLists.txt...")
        update_cmake()

    print("Done.")
