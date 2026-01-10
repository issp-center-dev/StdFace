# Phase 2 Worklog: Developer Architecture Documentation

**Date**: 2026-01-10
**Phase**: 2 (Developer Architecture Documentation)
**Status**: Complete

## Objective

Generate the 4 architecture documentation files under `docs/dev/architecture/`:
- `overview.rst`
- `directory_structure.rst`
- `execution_flow.rst`
- `dependency_graph.rst`

## Inputs Used

- `docs/_meta/header_map.yaml` - Header file inventory and function counts
- `docs/_meta/entry_points.yaml` - Entry point documentation
- `docs/_meta/terminology.yaml` - Seed terminology
- Source code analysis via Serena MCP tools

## Artifacts Produced

### docs/dev/architecture/overview.rst
- Purpose and target solvers (HPhi, mVMC, UHF, H-wave)
- Design philosophy: single-source, multi-target
- Header location note (headers in `src/`, not `include/`)
- Core components: entry point, main processing, lattices
- Build-time mode selection table
- Thread/MPI safety noted as "unspecified in the current code"

### docs/dev/architecture/directory_structure.rst
- Repository layout diagram
- Header files table (5 headers with purposes)
- Source files tables (core and lattice-specific)
- samples/ directory structure
- Build output structure

### docs/dev/architecture/execution_flow.rst
- High-level flow diagram
- CLI entry point code example (src/dry.c:34-47)
- StdFace_main() 6-step processing:
  1. Initialization
  2. Input file parsing
  3. Model validation
  4. Lattice construction
  5. Mode-specific output generation
  6. Cleanup
- Mode-specific behavior table
- Error handling documentation

### docs/dev/architecture/dependency_graph.rst
- Header inclusion graph (ASCII diagram)
- Header dependency diagram
- Main call path diagram
- Lattice constructor dependencies
- Memory allocation dependencies
- Component responsibilities table
- CMake dependency structure
- Compile-time dependencies table
- Data flow diagram

## Key Findings

1. **Headers in src/**: All 5 public headers reside in `src/`, not `include/`. This is explicitly documented with an `.. important::` directive in both `overview.rst` and `directory_structure.rst`.

2. **StdIntList structure**: Central data container with ~125 fields, defined in `src/StdFace_vals.h`.

3. **Build-time mode selection**: Four mutually exclusive modes (_HPhi, _mVMC, _UHF, _HWAVE) controlled by CMake options.

4. **Lattice abstraction**: 10 lattice constructors following a common pattern, all using utility functions from `StdFace_ModelUtil.c`.

5. **Conditional compilation**: `export_wannier90.h` only available when `_HWAVE` is defined.

## Source References

All claims in the architecture documents include source references in the format:
- `(file_path:line_range)` for code locations
- `(file_path)` for general file references

## Verification

- All 4 RST files follow Sphinx reStructuredText syntax
- Source references verified against actual source code
- ASCII diagrams render correctly in monospace
- Tables use proper RST grid or simple table format

## Notes

- Large file `src/StdFace_main.c` (3066 lines) required partial reads and Serena symbol extraction
- No issues discovered during Phase 2 that require tracking
