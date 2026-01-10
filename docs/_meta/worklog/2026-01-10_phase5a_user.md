# Phase 5 Worklog: User Documentation (Part A)

**Date**: 2026-01-10
**Phase**: 5 (User Documentation) - Part A
**Status**: Complete (Part A only)

## Objective

Create user-facing quickstart guide and examples documentation.

## Inputs Used

- `docs/_meta/entry_points.yaml` - CLI entry points and CMake options
- `docs/_meta/terminology.yaml` - Glossary of terms
- `docs/dev/internals/build_and_options.rst` - Build system documentation
- `samples/hubbard/default_model/stan.in` - Default model example
- `samples/hubbard/wannier/stan.in` - Wannier90 example
- `samples/hubbard/wannier/zvo_geom.dat` - Geometry data
- `samples/hubbard/wannier/zvo_hr.dat` - Hopping data
- `samples/hubbard/wannier/zvo_ur.dat` - Interaction data

## Artifacts Produced

### docs/user/quickstart.rst

- Prerequisites (CMake, compiler, libm)
- Build instructions with CMake options table
- Run command pattern (`xxx_dry.out <input_file>`)
- Minimal input file example (2x2 Hubbard square)
- Troubleshooting section (build errors, runtime errors)
- Cross-references to examples and input_output

### docs/user/examples.rst

- Hubbard Model (Default): Complete stan.in example with parameter table
- Hubbard Model (Wannier90-based): Advanced example with multiple files
- File relationships diagram (stan.in -> zvo_*.dat)
- Lattice types reference table

## Key Decisions

1. **Output files unspecified**: The dev docs do not explicitly document output
   file names/formats. Instructed users to "consult generated outputs" rather
   than guessing.

2. **Wannier90 example as advanced**: Clearly labeled as advanced use case,
   with explicit file relationship documentation.

3. **No mode-specific labeling for Wannier90**: The `lattice = "wannier90"`
   option works with all solver modes (not HWAVE-specific). HWAVE-specific
   functions (ExportGeometry, ExportInteraction) are for export, not import.

4. **Primary "first success" path**: Used the `samples/hubbard/default_model`
   example as the recommended starting point.

## Review Loop Summary

**Readability Review**:
- Clear section headings and progressive complexity
- Code blocks with proper syntax highlighting
- Tables for parameter explanations and references

**Correctness Review**:
- All CMake options verified against entry_points.yaml
- Sample files quoted exactly from source
- No invented output filenames or formats

## Post-Draft Fixes

1. **RST table alignment**: Fixed column separator width in quickstart.rst
   build options table (``-DHWAVE=ON`` required 14-character column).

## Verification

- Sphinx build: 0 warnings, 0 errors
- Both user pages render correctly

## Tools Used

- **Serena**: Reading sample files and dev docs
- **filesystem**: Writing RST files
- **Bash**: Sphinx build verification

## Metrics

| Metric | Value |
|--------|-------|
| User pages created | 2 |
| Examples documented | 2 (default, wannier) |
| Lattice types listed | 11 |

## Revision: quickstart.rst (Post-Draft)

**Issue**: Original quickstart.rst contained inferred input keys and claimed
parameters were "required" without code evidence.

**Changes made**:
1. Replaced "Minimal example" and "Your First Run" blocks with exact contents
   of `samples/hubbard/default_model/stan.in`
2. Changed "Key parameters" to "Keys used in the sample"
3. Removed claim about "Common required parameters" in troubleshooting
4. Added explicit source path for sample file

**Consistency pass**:
- Readability: Clear structure, properly labeled sample
- Correctness: All content verified against entry_points.yaml and sample file

## Remaining Work

- Part B: `docs/user/input_output.rst` (detailed parameter reference)
