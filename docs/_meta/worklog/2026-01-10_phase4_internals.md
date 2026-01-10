# Phase 4 Worklog: Developer Internals Documentation

**Date**: 2026-01-10
**Phase**: 4 (Developer Internals)
**Status**: Complete

## Objective

Create internal documentation pages covering error handling, memory management,
and build system options.

## Inputs Used

- `docs/_meta/entry_points.yaml` - CLI entry points and CMake options
- `CMakeLists.txt` (root) - Build option definitions
- `src/CMakeLists.txt` - Build configuration for each mode
- `src/StdFace_ModelUtil.c` - Error handling implementations
- `src/setmemory.h`, `src/setmemory.c` - Memory allocation functions

## Artifacts Produced

### docs/dev/internals/error_handling.rst

- Documents `StdFace_exit()` function with source reference
- Covers error categories: required parameter missing, unused parameter
- Documents `StdFace_RequiredVal_i()` and `StdFace_NotUsed_*` functions
- Exit codes and error output format documented

### docs/dev/internals/memory_management.rst

- Documents setmemory.h module (24 functions)
- Explains naming convention (`<type>_<dim>d_allocate`)
- Describes allocation patterns for 1D, 2D, 3D arrays
- Documents contiguous memory layout for multi-dimensional arrays
- Function summary table with allocation/deallocation pairs

### docs/dev/internals/build_and_options.rst

- Documents CMake options (UHF, MVMC, HPHI, HWAVE)
- Explains preprocessor defines (`_UHF`, `_mVMC`, `_HPhi`, `_HWAVE`)
- Documents build artifacts (static libraries and executables)
- Covers conditional compilation patterns
- Details export_wannier90.h availability under `_HWAVE`
- Lists all source files compiled into variants

## Key Findings

1. **Error handling pattern**: All errors are fatal with no recovery mechanism.
   Messages are written to stdout (not stderr).

2. **MPI awareness**: `StdFace_exit()` properly handles MPI cleanup when
   compiled with MPI support.

3. **Memory allocation**: Uses `calloc()` for zero-initialization. No NULL
   checks after allocation (observed pattern).

4. **2D/3D array layout**: Multi-dimensional arrays use contiguous memory with
   pointer indexing for cache efficiency.

5. **Build modes**: Four distinct modes (UHF, MVMC, HPHI, HWAVE) each produce
   separate static library and executable.

## Review Loop Summary

**Readability Review**:
- Organized content with clear section headings
- Code examples included for key implementations
- Tables used for summarizing options and functions

**Correctness Review**:
- All code snippets verified against source using Serena
- Source file:line references included for key implementations
- CMake option names verified against CMakeLists.txt

## Post-Draft Fixes

1. **RST table alignment**: Fixed table column widths in error_handling.rst
   and build_and_options.rst to prevent "Malformed table" errors.

## Verification

- Sphinx build: 0 warnings, 0 errors
- All 3 internals pages render correctly

## Tools Used

- **Serena**: Symbol search and body extraction
- **filesystem**: Reading CMake files, writing RST files
- **Bash**: Sphinx build verification

## Metrics

| Metric | Value |
|--------|-------|
| Pages created | 3 |
| Functions documented | 24 (setmemory) + 6 (error handling) |
| Build options documented | 4 |
| Preprocessor defines documented | 4 |
