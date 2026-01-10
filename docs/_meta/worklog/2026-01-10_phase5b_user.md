# Phase 5 Worklog: User Documentation (Part B)

**Date**: 2026-01-10
**Phase**: 5 (User Documentation) - Part B
**Status**: Complete

## Objective

Create comprehensive input/output reference documentation for StdFace users.

## Inputs Used

- `samples/hubbard/default_model/stan.in` - Default model sample
- `samples/hubbard/wannier/stan.in` - Wannier90 sample
- `samples/hubbard/wannier/zvo_geom.dat` - Geometry file
- `samples/hubbard/wannier/zvo_hr.dat` - Hopping file
- `samples/hubbard/wannier/zvo_ur.dat` - Interaction file
- `docs/dev/internals/error_handling.rst` - Error message patterns
- `docs/_meta/terminology.yaml` - Key definitions
- `docs/_meta/entry_points.yaml` - Build options

## Artifact Produced

### docs/user/input_output.rst

Seven-section structure as required:

1. **Overview**: Purpose and scope of StdFace
2. **Input File Format**: key=value syntax, quoting, ordering
3. **Common Keys (Observed in Samples)**: 10 keys from default_model sample
4. **Mode-Specific Inputs**: Wannier90 file format and relationships
5. **Output Files**: Solver-dependent, unspecified
6. **Validation and Errors**: Error messages and causes
7. **Limitations and Non-Goals**: What StdFace does not do

## Key Decisions

1. **No invented keys**: Only documented keys observed in sample files
2. **Labeling convention**: All keys labeled as "observed in sample" or
   "solver-specific", never "required" or "mandatory"
3. **"Unspecified" stance**: Used throughout for behaviors not explicitly
   documented in dev docs
4. **Error message accuracy**: Mapped directly from error_handling.rst

## Review Loop Summary

**Readability Review**:
- Clear 7-section structure
- Tables for key reference
- Code blocks from actual samples
- Consistent terminology

**Correctness Review**:
- All keys verified against sample files
- Error messages match dev docs
- No unverified claims about requirements
- Wannier90 format from actual file content

## Verification

- Sphinx build: 0 warnings, 0 errors
- All 7 sections render correctly

## Tools Used

- **Serena**: Reading sample and dev doc files
- **Bash**: Finding all stan.in files in samples
- **filesystem**: Writing RST file

## Metrics

| Metric | Value |
|--------|-------|
| Sections created | 7 |
| Keys documented | 10 (from default_model) |
| Wannier90 files documented | 3 |
| Error messages documented | 2 |
