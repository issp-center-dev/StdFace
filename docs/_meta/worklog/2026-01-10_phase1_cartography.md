# Worklog: Phase 1 (Repo Cartography)

Date: 2026-01-10

## Scope
Identify public API surface, entry points, CMake options, and seed terminology for StdFace.

## What was done

### 1. Public Headers Identification
- **Finding**: No `include/` directory exists. All headers are in `src/`.
- Identified 5 header files:
  - `src/StdFace_ModelUtil.h` - 37 functions (core lattice/model utilities)
  - `src/StdFace_vals.h` - StdIntList struct with ~125 fields
  - `src/setmemory.h` - 24 memory allocation functions
  - `src/export_wannier90.h` - 2 functions (HWAVE-only, conditional)
  - `src/version.h` - 1 function + version macros

### 2. API Surface Extraction
- Total: 64 public functions, 1 struct type, 7 macros (including guards)
- Used Serena `get_symbols_overview` for each header
- Used `search_for_pattern` to find `#define` macros

### 3. Entry Points Identification
- **Main entry**: `src/dry.c::main(int argc, char *argv[])`
- **Core function**: `src/StdFace_main.c::StdFace_main(char *fname)`
- 4 binary variants: uhf_dry.out, mvmc_dry.out, hphi_dry.out, hwave_dry.out

### 4. CMake Build Options
- `UHF` (OFF) - builds uhf_dry.out with `_UHF` define
- `MVMC` (OFF) - builds mvmc_dry.out with `_mVMC` define
- `HPHI` (OFF) - builds hphi_dry.out with `_HPhi` define
- `HWAVE` (OFF) - builds hwave_dry.out with `_HWAVE` define
- `TestStdFace` (ON) - enables testing

### 5. Terminology Seeding
- Documented 25 terms from README, CMake, and samples
- Categories: project/tools, models, lattice types, parameters

## Tools used
- **serena**: `find_file`, `list_dir`, `get_symbols_overview`, `search_for_pattern`
- **filesystem**: Read CMakeLists.txt, README.md, samples/
- **git**: Status check

## Review loop summary
- **Readability issues**: N/A (meta files, not documentation)
- **Correctness issues**:
  - All symbols verified via Serena
  - All file paths confirmed to exist

## Changes summary

Files created/updated:
- `docs/_meta/header_map.yaml` (166 lines) - complete API surface
- `docs/_meta/entry_points.yaml` (102 lines) - entry points and binaries
- `docs/_meta/terminology.yaml` (154 lines) - seed glossary
- `docs/_meta/progress.json` - updated phase1_cartography to done

## Key findings

1. **No include/ directory**: Headers are in `src/`, not following typical C project layout.
   This is documented but does not affect API extraction.

2. **Conditional compilation**: `export_wannier90.h` functions are only available when
   `_HWAVE` is defined.

3. **Single entry point**: All 4 binaries share the same `main()` in `dry.c`, differentiated
   by compile-time defines.

4. **Large parameter structure**: `StdIntList` has ~125 fields covering lattice geometry,
   model parameters, and interaction coefficients.

## Next actions
- Phase 2: Architecture documentation (use header_map.yaml and entry_points.yaml as input)
- Do NOT proceed to Phase 2 per user instruction

## Suggested git commit message
```
docs: complete Phase 1 Repo Cartography

- Generate header_map.yaml with 5 headers, 64 functions, 1 type
- Generate entry_points.yaml with main(), binaries, and CMake options
- Generate terminology.yaml with 25 seed terms
- Document finding: headers in src/, not include/
- Update progress.json: phase1_cartography -> done

All symbols verified via Serena. Ready for Phase 2 (Architecture).
```
