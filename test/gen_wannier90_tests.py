#!/usr/bin/env python3
"""Generate wannier90 regression tests for all (solver, model) combinations.

Test system: 2D square lattice, a=3 Å, 1 Wannier orbital per unit cell,
nearest-neighbor hopping t=-1.0 eV.  Supercell W=2, L=2, Height=1 → 4 sites.

Usage:
    python3 gen_wannier90_tests.py [--gen-ref] [--update-cmake]

    --gen-ref       Run dry executables to populate ref/ directories.
    --update-cmake  Append wannier90 test entries to each solver's CMakeLists.txt.

Executables are expected in ../build/src/ relative to this script.

Key difference from gen_matrix_tests.py:
  The wannier90 dry-run reads zvo_geom.dat / zvo_hr.dat from the WORKING directory.
  check_case.sh runs the executable from the test directory, so the data files must
  live alongside stan.in.  Reference generation also runs from the test directory
  (not from ref/) so that the same data files are found; generated output files are
  then moved into ref/.
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()
BUILD_DIR  = SCRIPT_DIR.parent / "build" / "src"
DATA_DIR   = SCRIPT_DIR / "wannier90_data"

# Data file names (prefix = CDataFileHead = "zvo")
DATA_FILES = ["zvo_geom.dat", "zvo_hr.dat", "zvo_ur.dat"]

# Supercell: W=2, L=2, Height=1  → nsite = W*L*Height*NsiteUC = 4
NSITE = 4

# Valid models for wannier90 (kondo excluded: StdFace_exit in source)
MODELS = [
    ("hubbard",   False, "hubbard"),
    ("hubbardgc", True,  "hubbard"),
    ("spin",      False, "spin"),
    ("spingc",    True,  "spin"),
]

# Solver configuration  (same structure as gen_matrix_tests.py)
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

MODEL_INPUT = {
    "hubbard":   "FermionHubbard",
    "hubbardgc": "FermionHubbardGC",
    "spin":      "Spin",
    "spingc":    "SpinGC",
}


def build_stan_in(solver, model_name, is_gc, base, calcmode=None):
    lines = []
    # CDataFileHead is HPhi-only; other solvers default to "zvo" automatically.
    if solver == "hphi":
        lines.append('CDataFileHead = "zvo"')
    lines += [
        f'model = "{MODEL_INPUT[model_name]}"',
        'lattice = "wannier90"',
        "W = 2",
        "L = 2",
        "Height = 1",
    ]

    if solver == "hphi":
        lines.append('method = "FullDiag"')
    elif calcmode is not None:
        lines.append(f'calcmode = "{calcmode}"')

    # Conserved-quantity parameters (same rules as matrix tests)
    if base == "hubbard":
        if solver == "hphi":
            if not is_gc:
                lines.append(f"nelec = {NSITE}")
                lines.append("2Sz = 0")
        else:
            lines.append(f"nelec = {NSITE}")
            if not is_gc:
                lines.append("2Sz = 0")
    elif base == "spin":
        if not is_gc:
            lines.append("2Sz = 0")

    if solver == "mvmc":
        lines.append("RndSeed = 1")

    return "\n".join(lines) + "\n"


def test_name(model_name, calcmode=None):
    suffix = f"_{calcmode}" if calcmode is not None else ""
    return f"wannier90_{model_name}{suffix}"


def generate_test_dirs():
    for solver_name, solver_cfg in SOLVERS.items():
        solver_dir = SCRIPT_DIR / solver_name
        infile = solver_cfg["infile"]
        for model_name, is_gc, base in MODELS:
            for calcmode in solver_cfg["calcmodes"]:
                name = test_name(model_name, calcmode)
                test_dir = solver_dir / name
                test_dir.mkdir(parents=True, exist_ok=True)

                # Input file
                content = build_stan_in(solver_name, model_name, is_gc, base, calcmode)
                (test_dir / infile).write_text(content)

                # Data files (must be co-located with stan.in for check_case.sh)
                for df in DATA_FILES:
                    shutil.copy(DATA_DIR / df, test_dir / df)

                print(f"  wrote {solver_name}/{name}/")


def generate_refs():
    """Run each dry executable from the TEST directory (not from ref/).

    wannier90 reads zvo_geom.dat relative to cwd, so we must run from the
    directory that contains the data files.  Generated output files are then
    collected into ref/.
    """
    # Files that belong to the test fixture, not to the generated output.
    fixture_files = set(DATA_FILES)

    for solver_name, solver_cfg in SOLVERS.items():
        exe_path = (BUILD_DIR / solver_cfg["exe"]).resolve()
        if not exe_path.exists():
            print(f"WARNING: {exe_path} not found, skipping {solver_name}")
            continue

        solver_dir  = SCRIPT_DIR / solver_name
        infile      = solver_cfg["infile"]

        for model_name, is_gc, base in MODELS:
            for calcmode in solver_cfg["calcmodes"]:
                name     = test_name(model_name, calcmode)
                test_dir = solver_dir / name

                if not (test_dir / infile).exists():
                    print(f"  SKIP {solver_name}/{name} (no input file)")
                    continue

                ref_dir = test_dir / "ref"
                if ref_dir.exists():
                    shutil.rmtree(ref_dir)
                ref_dir.mkdir()

                # Snapshot of files BEFORE running
                before = set(p.name for p in test_dir.iterdir()
                             if p.is_file() and p.name != infile)

                result = subprocess.run(
                    [str(exe_path), infile],
                    cwd=str(test_dir),
                    capture_output=True,
                    text=True,
                )

                if result.returncode != 0:
                    print(f"  FAIL {solver_name}/{name}: {result.stdout.strip()[-200:]}")
                    continue

                # Collect newly generated files → ref/
                after = set(p.name for p in test_dir.iterdir()
                            if p.is_file() and p.name != infile)
                new_files = after - before - fixture_files
                for fname in sorted(new_files):
                    shutil.move(str(test_dir / fname), str(ref_dir / fname))

                print(f"  OK   {solver_name}/{name}  ({len(new_files)} files)")


def cmake_entries(solver_name):
    solver_cfg = SOLVERS[solver_name]
    lines = []
    for model_name, _is_gc, _base in MODELS:
        for calcmode in solver_cfg["calcmodes"]:
            name = test_name(model_name, calcmode)
            lines.append(
                f"add_test(\n"
                f"  NAME {name}\n"
                f"  COMMAND ${{CMAKE_SOURCE_DIR}}/test/{solver_name}/check.sh"
                f" {name} ${{CMAKE_SOURCE_DIR}}/test/{solver_name}\n"
                f")"
            )
    return "\n".join(lines) + "\n"


def update_cmake():
    marker = "# --- wannier90 tests (auto-generated) ---"
    for solver_name in SOLVERS:
        cmake_path = SCRIPT_DIR / solver_name / "CMakeLists.txt"
        original = cmake_path.read_text()

        if marker in original:
            original = original[: original.index(marker)]

        new_content = original.rstrip("\n") + "\n\n" + marker + "\n" + cmake_entries(solver_name)
        cmake_path.write_text(new_content)
        print(f"  updated {solver_name}/CMakeLists.txt")


if __name__ == "__main__":
    args = set(sys.argv[1:])

    print("==> Generating test directories...")
    generate_test_dirs()

    if "--gen-ref" in args:
        print("==> Generating ref/ files...")
        generate_refs()

    if "--update-cmake" in args:
        print("==> Updating CMakeLists.txt...")
        update_cmake()

    print("Done.")
