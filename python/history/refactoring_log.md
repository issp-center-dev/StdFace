# StdFace Python Refactoring Log

This file tracks all refactoring changes applied to the Python codebase.
Each entry records the date, what was changed, why, and the test results.

---

## Refactoring Plan Overview

### Goal
Transform the mechanically-translated C→Python code into idiomatic, maintainable Python.

### Current state (before refactoring)
- 16 Python files, ~11,400 lines total
- `stdface_main.py` alone is ~3,900 lines (monolithic)
- No classes (all functions taking `StdI` dataclass as first argument)
- C-style patterns: sentinel values, manual dispatch, flat structure
- 252 unit tests, 83 integration tests — all passing

### Planned phases
1. Extract modules from `stdface_main.py` (keyword parser, solver writers, interaction writer)
2. Introduce classes (SolverWriter hierarchy, KeywordParser, InteractionBuilder)
3. Leverage Python idioms (enums, dict dispatch, dataclass improvements)

---

## Change Log

### 2026-01-27 — Phase 1, Step 1: Extract HPhi writer functions

**Target**: `stdface_main.py` → new `hphi_writer.py`

**What changed**:
- Created `python/hphi_writer.py` (~530 lines) containing 5 HPhi-specific output functions:
  - `large_value(StdI)` — compute LargeValue for TPQ
  - `print_calc_mod(StdI)` — write `calcmod.def`
  - `print_excitation(StdI)` — write `single.def` / `pair.def`
  - `vector_potential(StdI)` — compute A(t)/E(t), write `potential.dat`
  - `print_pump(StdI)` — write `teone.def` / `tetwo.def`
- Updated `stdface_main.py` to import these functions from `hphi_writer`
- Removed ~659 lines from `stdface_main.py` (reduced from ~3,908 to ~3,195 lines)

**Why**: First step in breaking up the monolithic `stdface_main.py`. HPhi solver
output functions are a self-contained group with clear boundaries.

**New tests**: `test/unit/test_hphi_writer.py` — 9 tests covering:
- `TestLargeValue` (3 tests): basic computation, interactions, user-specified value
- `TestPrintCalcMod` (4 tests): fulldiag/hubbard, spinGC, EigenVecIO, kondoGC
- `TestVectorPotential` (2 tests): quench default, aclaser potential output

**Test results**:
- Unit tests: 261 passed (252 existing + 9 new)
- Integration check: HPhi lanczos_hubbard_square — no diffs against reference

**Files created**: `python/hphi_writer.py`, `test/unit/test_hphi_writer.py`
**Files modified**: `python/stdface_main.py`

---

### 2026-01-27 — Phase 1, Step 2: Extract mVMC writer functions

**Target**: `stdface_main.py` → new `mvmc_writer.py`

**What changed**:
- Created `python/mvmc_writer.py` (~346 lines) containing 3 mVMC-specific output functions:
  - `print_orb(StdI)` — write `orbitalidx.def` (anti-parallel orbital indices)
  - `print_orb_para(StdI)` — write `orbitalidxpara.def` and `orbitalidxgen.def` (parallel orbital indices)
  - `print_gutzwiller(StdI)` — write `gutzwilleridx.def` (Gutzwiller variational parameters)
- Updated `stdface_main.py` to import these functions from `mvmc_writer`
- Removed ~311 lines from `stdface_main.py` (reduced from ~3,195 to ~2,884 lines)

**Why**: Second step in breaking up the monolithic `stdface_main.py`. mVMC solver
output functions are a self-contained group that write variational wave function
parameter files.

**New tests**: `test/unit/test_mvmc_writer.py` — 14 tests covering:
- `TestPrintOrb` (4 tests): file creation, 3-column format (no anti-periodic),
  4-column format (anti-periodic), optimization line count
- `TestPrintOrbPara` (3 tests): both files created, header content, NOrbitalIdx count
- `TestPrintGutzwiller` (7 tests): file creation, momentum-projected hubbard/spin,
  global-optimisation hubbard/spin, optimization flags

**Test results**:
- Unit tests: 275 passed (261 existing + 14 new)
- Integration check: HPhi lanczos_hubbard_square — no diffs; mVMC — OK

**Files created**: `python/mvmc_writer.py`, `test/unit/test_mvmc_writer.py`
**Files modified**: `python/stdface_main.py`

---

### 2026-01-27 — Phase 1, Step 3: Extract common output writer functions

**Target**: `stdface_main.py` → new `common_writer.py`

**What changed**:
- Created `python/common_writer.py` (~1,248 lines) containing 10 shared output functions:
  - `print_loc_spin(StdI)` — write `locspn.def`
  - `print_trans(StdI)` — write `trans.def`
  - `print_namelist(StdI)` — write `namelist.def`
  - `print_mod_para(StdI)` — write `modpara.def`
  - `print_1_green(StdI)` — write `greenone.def`
  - `print_2_green(StdI)` — write `greentwo.def`
  - `unsupported_system(model, lattice)` — error and exit
  - `check_output_mode(StdI)` — validate/set ioutputmode
  - `check_mod_para(StdI)` — validate solver-specific parameters
  - `print_interactions(StdI)` — write all interaction `.def` files
- Updated `stdface_main.py` to import these functions from `common_writer`
- Removed ~1,248 lines from `stdface_main.py` (reduced from ~2,884 to ~1,648 lines)

**Why**: Third and largest step in breaking up the monolithic `stdface_main.py`.
These functions are shared across all 4 solvers (HPhi, mVMC, UHF, HWAVE) and
write the common definition files (locspn, trans, modpara, namelist, green
functions, interactions). This extraction completes Phase 1 of the module
extraction plan.

**New tests**: `test/unit/test_common_writer.py` — 23 tests covering:
- `TestPrintLocSpin` (2 tests): file creation, zero local spins
- `TestPrintTrans` (3 tests): file creation, duplicate merging, small-value suppression
- `TestCheckOutputMode` (6 tests): none/correlation/default/raw/full modes, invalid exit
- `TestUnsupportedSystem` (1 test): exit behavior
- `TestPrintNamelist` (2 tests): HPhi entries, mVMC entries
- `TestPrint1Green` (3 tests): no-output mode, correlation, raw
- `TestPrint2Green` (1 test): correlation mode
- `TestPrintInteractions` (2 tests): no interactions, coulomb intra
- `TestPrintModPara` (3 tests): HPhi, UHF, HWAVE solver headers

**Test results**:
- Unit tests: 298 passed (275 existing + 23 new)
- Integration check: HPhi — no diffs; mVMC, UHF, HWAVE — all OK

**Files created**: `python/common_writer.py`, `test/unit/test_common_writer.py`
**Files modified**: `python/stdface_main.py`, `python/history/refactoring_log.md`

**Summary of Phase 1 progress** (module extraction from `stdface_main.py`):

| Step | Module | Functions | Lines extracted |
|------|--------|-----------|----------------|
| 1 | `hphi_writer.py` | 5 HPhi functions | ~530 |
| 2 | `mvmc_writer.py` | 3 mVMC functions | ~346 |
| 3 | `common_writer.py` | 10 shared functions | ~1,248 |
| **Total** | | **18 functions** | **~2,124** |

`stdface_main.py` reduced from ~3,908 lines (original) to ~1,648 lines (**58% reduction**).

---

### 2026-01-27 — Phase 1, Step 4: Extract keyword parser functions

**Target**: `stdface_main.py` → new `keyword_parser.py`

**What changed**:
- Created `python/keyword_parser.py` (~1,020 lines) containing 13 keyword parsing functions:
  - `text2lower(text)` — convert text to lower case
  - `trim_space_quote(text)` — remove whitespace, colons, semicolons, quotes, backslashes
  - `store_with_check_dup_s(keyword, value, current)` — store string with duplicate check
  - `store_with_check_dup_sl(keyword, value, current, maxlen)` — store string (lowercased) with duplicate check
  - `store_with_check_dup_i(keyword, value, current)` — store integer with duplicate check
  - `store_with_check_dup_d(keyword, value, current)` — store float with duplicate check
  - `store_with_check_dup_c(keyword, value, current)` — store complex with duplicate check
  - `parse_common_keyword(keyword, value, StdI)` — parse keywords common to all solvers
  - `parse_solver_keyword(keyword, value, StdI, solver)` — dispatch to solver-specific parser
  - `parse_hphi_keyword(keyword, value, StdI)` — parse HPhi-specific keywords
  - `parse_mvmc_keyword(keyword, value, StdI)` — parse mVMC-specific keywords
  - `parse_uhf_keyword(keyword, value, StdI)` — parse UHF-specific keywords
  - `parse_hwave_keyword(keyword, value, StdI)` — parse HWAVE-specific keywords
- Updated `stdface_main.py` to import these functions from `keyword_parser`
- Removed ~1,039 lines from `stdface_main.py` (reduced from ~1,648 to ~608 lines)
- Removed unused `import cmath` from `stdface_main.py`

**Why**: Fourth and final step of Phase 1, extracting the keyword parsing subsystem.
The parser helpers and solver-specific parsers form a cohesive, self-contained group
with clear boundaries. This completes the decomposition of the monolithic
`stdface_main.py`.

**New tests**: `test/unit/test_keyword_parser.py` — 54 tests covering:
- `TestText2Lower` (3 tests): lowercase, mixed case, already lower
- `TestTrimSpaceQuote` (5 tests): spaces, colons/semicolons, quotes, backslashes, preserves equals
- `TestStoreWithCheckDupS` (2 tests): new value, duplicate exit
- `TestStoreWithCheckDupSl` (3 tests): lowercase, truncation, duplicate exit
- `TestStoreWithCheckDupI` (3 tests): integer, float truncation, duplicate exit
- `TestStoreWithCheckDupD` (2 tests): float, duplicate exit
- `TestStoreWithCheckDupC` (4 tests): real+imag, real only, imag only, duplicate exit
- `TestParseCommonKeyword` (12 tests): model, lattice, L, W, U, t, J, 2Sz, phase0, unrecognised, box, J0 matrix
- `TestParseSolverKeyword` (5 tests): HPhi/mVMC/UHF/HWAVE dispatch, unknown solver
- `TestParseHPhiKeyword` (4 tests): method, exct, dt, unrecognised
- `TestParseMVMCKeyword` (4 tests): nvmcsample, complextype, boxsub, unrecognised
- `TestParseUHFKeyword` (3 tests): iteration_max, mix, unrecognised
- `TestParseHWAVEKeyword` (4 tests): calcmode, fileprefix, exportall, unrecognised

**Test results**:
- Unit tests: 352 passed (298 existing + 54 new)
- Integration check: HPhi, mVMC, UHF, HWAVE — all OK

**Files created**: `python/keyword_parser.py`, `test/unit/test_keyword_parser.py`
**Files modified**: `python/stdface_main.py`, `python/history/refactoring_log.md`

**Summary of Phase 1 (complete)** — module extraction from `stdface_main.py`:

| Step | Module | Functions | Lines extracted |
|------|--------|-----------|----------------|
| 1 | `hphi_writer.py` | 5 HPhi functions | ~530 |
| 2 | `mvmc_writer.py` | 3 mVMC functions | ~346 |
| 3 | `common_writer.py` | 10 shared functions | ~1,248 |
| 4 | `keyword_parser.py` | 13 parser functions | ~1,020 |
| **Total** | | **31 functions** | **~3,144** |

`stdface_main.py` reduced from ~3,908 lines (original) to ~608 lines (**84% reduction**).
Remaining in `stdface_main.py`: `_reset_vals()` and `stdface_main()` entry point.

---

### 2026-01-27 — Phase 2, Step 1: Introduce SolverWriter class hierarchy

**Target**: `stdface_main.py` solver dispatch → new `solver_writer.py`

**What changed**:
- Created `python/solver_writer.py` (~255 lines) containing:
  - `SolverWriter` — abstract base class with `write(StdI)` method
  - `HPhiWriter` — HPhi solver output (locspn, trans, interactions, modpara,
    excitation, pump, calcmod, green, namelist)
  - `MVMCWriter` — mVMC solver output (locspn, trans, interactions, modpara,
    ComplexType, orbital, Gutzwiller, Jastrow, projection, green, namelist)
  - `UHFWriter` — UHF solver output (locspn, trans, interactions, modpara,
    green, namelist)
  - `HWaveWriter` — H-wave solver output (uhfr mode: trans/interactions/green;
    otherwise: Wannier90 export)
  - `get_solver_writer(solver)` — factory function returning the correct writer
  - `_SOLVER_WRITERS` — registry dict mapping solver names to writer classes
- Updated `stdface_main.py`:
  - Replaced 60-line `if solver == "HPhi" ... elif ... elif ... elif ...` block
    with `writer = get_solver_writer(solver); writer.write(StdI)` (2 lines)
  - Removed unused imports: `sys`, `numpy`, `print_val_i`, `generate_orb`,
    `proj`, `print_jastrow`, `export_geometry`, `export_interaction`,
    `_store_with_check_dup_*`, all `_print_*` / `_check_*` aliases from
    common_writer, hphi_writer, mvmc_writer
  - `stdface_main.py` reduced from ~608 to ~517 lines

**Why**: First step of Phase 2 (introduce classes). The solver dispatch block
was a classic "switch on type" code smell — each solver's file-writing logic
is now encapsulated in its own class with a clean `write()` interface. The
factory function provides a single point of lookup. This makes it easy to
add new solvers and simplifies the main entry point.

**New tests**: `test/unit/test_solver_writer.py` — 11 tests covering:
- `TestGetSolverWriter` (5 tests): HPhi/mVMC/UHF/HWAVE factory, unknown solver error
- `TestSolverWriterIsAbstract` (1 test): cannot instantiate base class
- `TestHPhiWriter` (3 tests): namelist+calcmod+modpara+locspn, trans, green files
- `TestUHFWriter` (1 test): expected file set
- `TestHWaveWriter` (1 test): uhfr mode output

**Test results**:
- Unit tests: 363 passed (352 existing + 11 new)
- Integration check: HPhi, mVMC, UHF, HWAVE — all OK

**Files created**: `python/solver_writer.py`, `test/unit/test_solver_writer.py`
**Files modified**: `python/stdface_main.py`, `python/history/refactoring_log.md`

---

### 2026-01-27 — Phase 2, Step 2: Replace lattice dispatch if/elif with dict tables

**Target**: `stdface_main.py` lattice dispatch block

**What changed**:
- Added two module-level dispatch dicts to `stdface_main.py`:
  - `LATTICE_DISPATCH` — maps 25 lattice name aliases to 10 builder functions
  - `BOOST_DISPATCH` — maps 8 lattice name aliases to 4 boost builder functions
- Replaced the 25-line lattice dispatch `if/elif` chain (lines 460-484) with a
  4-line dict lookup: `lattice_builder = LATTICE_DISPATCH.get(lattice)`
- Replaced the 11-line boost dispatch `if/elif` chain (lines 492-502) with a
  4-line dict lookup: `boost_builder = BOOST_DISPATCH.get(lattice)`
- The `stdface_main()` function body shrank by ~27 lines net

**Why**: Second step of Phase 2. The lattice dispatch was a textbook case for
dict-based dispatch — many aliases mapping to a small number of functions, with
no per-branch logic beyond calling the function. The dict approach is more
readable, easier to extend with new lattice types, and eliminates the risk of
duplicating alias strings across the main dispatch and boost dispatch.

**New tests**: `test/unit/test_lattice_dispatch.py` — 18 tests covering:
- `TestLatticeDispatch` (12 tests): chain, square, ladder, triangular,
  honeycomb, kagome, orthorhombic, fco, pyrochlore, wannier90 aliases;
  unknown returns None; all entries callable
- `TestBoostDispatch` (6 tests): chain, honeycomb, kagome, ladder boost
  aliases; unsupported lattices absent; all entries callable

**Test results**:
- Unit tests: 381 passed (363 existing + 18 new)
- Integration check: HPhi, mVMC, UHF, HWAVE — all OK

**Files created**: `test/unit/test_lattice_dispatch.py`
**Files modified**: `python/stdface_main.py`, `python/history/refactoring_log.md`

---

### 2026-01-27 — Phase 2, Step 3: Replace model/method normalisation if/elif with dict tables

**Target**: `stdface_main.py` model name and method name normalisation blocks

**What changed**:
- Added three module-level normalisation dicts to `stdface_main.py`:
  - `MODEL_ALIASES` — maps 10 model name aliases to `(canonical_name, lGC, lBoost)` tuples
  - `MODEL_ALIASES_HPHI_BOOST` — maps 2 HPhi-only model aliases (spingcboost, spingccma) to `(canonical_name, lGC=1, lBoost=1)` tuples
  - `METHOD_ALIASES` — maps 4 HPhi method name aliases to canonical forms (fulldiag, timeevolution)
- Replaced the 21-line model normalisation `if/elif` chain with a 6-line dict lookup:
  first checks `MODEL_ALIASES`, then (for HPhi) checks `MODEL_ALIASES_HPHI_BOOST`,
  and unpacks `(model, lGC, lBoost)` from the result
- Replaced the 4-line method normalisation `if/elif` chain with a 1-line
  `METHOD_ALIASES.get(StdI.method, StdI.method)` lookup
- The `stdface_main()` function body shrank by ~18 lines net

**Why**: Third step of Phase 2. The model and method normalisation blocks were
switch-on-string patterns that mapped user-facing aliases to canonical names plus
flag values. The dict approach is more readable, groups all alias→canonical
mappings in one visible place, and makes it trivial to add new model or method
aliases without touching control flow.

**New tests**: `test/unit/test_model_method_dispatch.py` — 22 tests covering:
- `TestModelAliases` (10 tests): hubbard, hubbardGC, spin, spinGC, kondo, kondoGC
  aliases; unknown returns None; all entries are valid tuples; GC/Boost flags 0/1;
  boost models absent from common table
- `TestModelAliasesHPhiBoost` (5 tests): spingcboost, spingccma; all have boost=1
  and gc=1; no overlap with common aliases
- `TestMethodAliases` (7 tests): direct→fulldiag, alldiag→fulldiag, te→timeevolution,
  time-evolution→timeevolution; canonical names not aliased; unknown returns None;
  all values are strings

**Test results**:
- Unit tests: 403 passed (381 existing + 22 new)
- Integration check: HPhi lanczos_hubbard_square, fulldiag_spingcboost_chain,
  fulldiag_hubbard_chain_longrange — all no diffs against reference

**Files created**: `test/unit/test_model_method_dispatch.py`
**Files modified**: `python/stdface_main.py`, `python/history/refactoring_log.md`

---

### 2026-01-27 — Phase 2, Step 4: Extract solver-specific reset helpers with dict dispatch

**Target**: `stdface_main.py` `_reset_vals()` solver-specific if/elif chain

**What changed**:
- Extracted 4 helper functions from the 95-line solver-specific if/elif block in
  `_reset_vals()`:
  - `_reset_hphi_fields(StdI)` — reset HPhi-specific fields (35 assignments)
  - `_reset_mvmc_fields(StdI)` — reset mVMC-specific fields (25 assignments)
  - `_reset_uhf_fields(StdI)` — reset UHF-specific fields (10 assignments)
  - `_reset_hwave_fields(StdI)` — reset HWAVE-specific fields (14 assignments)
- Added `_SOLVER_RESET_DISPATCH` dict mapping solver names to their reset functions
- Replaced the 95-line if/elif chain in `_reset_vals()` with a 3-line dict lookup:
  ```python
  _solver_reset = _SOLVER_RESET_DISPATCH.get(StdI.solver)
  if _solver_reset is not None:
      _solver_reset(StdI)
  ```
- `_reset_vals()` reduced from ~230 lines to ~140 lines; the extracted helpers
  total ~130 lines with docstrings

**Why**: Fourth step of Phase 2. The solver-specific field reset block was the last
major if/elif chain in `stdface_main.py`. Each branch was an independent block of
assignments with no shared logic, making it ideal for extraction into named helper
functions with dict-based dispatch. This improves readability (each solver's fields
are now in a clearly labelled function) and makes it straightforward to add new
solvers.

**New tests**: `test/unit/test_reset_vals_dispatch.py` — 23 tests covering:
- `TestSolverResetDispatch` (7 tests): HPhi/mVMC/UHF/HWAVE entries, unknown
  returns None, all entries callable, exactly 4 entries
- `TestResetHPhiFields` (5 tests): method sentinel, NaN_d fields, NaN_i fields,
  FlgTemp=1, string sentinels
- `TestResetMVMCFields` (4 tests): CParaFileHead sentinel, NaN_i fields, NaN_d
  fields, boxsub sentinels
- `TestResetUHFFields` (3 tests): NaN_i fields, NaN_d fields, boxsub sentinels
- `TestResetHWaveFields` (4 tests): NaN_i fields, NaN_d fields, string sentinels,
  extra fields vs UHF

**Test results**:
- Unit tests: 426 passed (403 existing + 23 new)
- Integration check: HPhi lanczos_hubbard_square, fulldiag_spingcboost_chain — all
  no diffs against reference

**Files created**: `test/unit/test_reset_vals_dispatch.py`
**Files modified**: `python/stdface_main.py`, `python/history/refactoring_log.md`

---

### 2026-01-27 — Phase 2, Step 5: Extract mVMC variational functions into new module

**Target**: `stdface_model_util.py` → new `mvmc_variational.py`

**What changed**:
- Created `python/mvmc_variational.py` (~388 lines) containing 5 mVMC-specific
  variational parameter functions extracted from `stdface_model_util.py`:
  - `_fold_site_sub(StdI, iCellV)` — fold site coordinates into a sub-lattice cell
  - `_init_site_sub(StdI)` — initialize sub-cell parameters (Lsub/Wsub/Hsub or boxsub)
  - `proj(StdI)` — generate and write `qptransidx.def` (momentum projection)
  - `generate_orb(StdI)` — compute orbital indices (Orb, AntiOrb arrays) for mVMC
  - `print_jastrow(StdI)` — generate and write `jastrowidx.def` (Jastrow variational parameters)
- Removed these 5 functions (~386 lines) from `stdface_model_util.py`, replaced with
  a comment noting the new location
- Updated `python/solver_writer.py` imports: `generate_orb`, `proj`, `print_jastrow`
  now imported from `mvmc_variational` instead of `stdface_model_util`
- `stdface_model_util.py` reduced from ~1,721 to ~1,335 lines

**Why**: Fifth step of Phase 2. These 5 functions are exclusively used by the mVMC
solver (via `solver_writer.py`'s `MVMCWriter`). They form a cohesive group dealing
with variational parameter generation — sub-lattice folding, orbital index computation,
Jastrow indices, and momentum projection. Extracting them into a dedicated module
reduces `stdface_model_util.py`'s size and improves separation of concerns.

**New tests**: `test/unit/test_mvmc_variational.py` — 13 tests covering:
- `TestFoldSiteSub` (2 tests): identity fold (site inside sub-cell), fold outside sub-cell
- `TestInitSiteSub` (3 tests): Wsub mode, boxsub mode, conflicting Lsub+boxsub raises SystemExit
- `TestGenerateOrb` (2 tests): NOrb set > 0, Orb/AntiOrb arrays have correct shape
- `TestProj` (3 tests): qptransidx.def created, NSym > 0, header contains NQPTrans
- `TestPrintJastrow` (3 tests): jastrowidx.def created, header contains NJastrowIdx,
  spin model with NMPTrans != 1

**Test results**:
- Unit tests: 439 passed (426 existing + 13 new)
- Integration check: HPhi fulldiag_hubbard_chain_longrange — all 12 output files match reference

**Files created**: `python/mvmc_variational.py`, `test/unit/test_mvmc_variational.py`
**Files modified**: `python/stdface_model_util.py`, `python/solver_writer.py`, `python/history/refactoring_log.md`

---

### 2026-01-27 — Phase 2, Step 6: Extract parameter validation utilities into new module

**Target**: `stdface_model_util.py` → new `param_check.py`

**What changed**:
- Created `python/param_check.py` (~260 lines) containing 10 parameter validation,
  printing, and program-exit utility functions extracted from `stdface_model_util.py`:
  - `exit_program(errorcode)` — terminate with error code (wraps ``sys.exit``)
  - `print_val_d(valname, val, val0)` — print/default a real-valued parameter
  - `print_val_dd(valname, val, val0, val1)` — print/default with two fallback defaults
  - `print_val_c(valname, val, val0)` — print/default a complex-valued parameter
  - `print_val_i(valname, val, val0)` — print/default an integer parameter
  - `not_used_d(valname, val)` — abort if a real parameter is specified but unused
  - `not_used_c(valname, val)` — abort if a complex parameter is specified but unused
  - `not_used_j(valname, JAll, J)` — abort if any J-type interaction component is unused
  - `not_used_i(valname, val)` — abort if an integer parameter is specified but unused
  - `required_val_i(valname, val)` — abort if a required integer parameter is missing
- Removed these 10 functions (~210 lines) and the `import sys` from `stdface_model_util.py`,
  replaced with a re-export import block for backward compatibility:
  ```python
  from param_check import (  # noqa: F401 – re-exported for backward compatibility
      exit_program, print_val_d, print_val_dd, print_val_c, print_val_i,
      not_used_d, not_used_c, not_used_j, not_used_i, required_val_i,
  )
  ```
- Updated 7 modules that only imported Category A (param_check) functions to import
  directly from `param_check` instead of `stdface_model_util`:
  - `stdface_main.py`: `exit_program`
  - `keyword_parser.py`: `exit_program`
  - `export_wannier90.py`: `exit_program`
  - `solver_writer.py`: `print_val_i`
  - `common_writer.py`: `exit_program`, `print_val_i`, `print_val_d`, `required_val_i`, `not_used_i`
  - `hphi_writer.py`: `exit_program`, `print_val_d`, `print_val_i`
  - `mvmc_variational.py`: `exit_program`, `print_val_i` (split into two import lines)
- 9 lattice modules (chain, square, ladder, triangular, honeycomb, kagome,
  orthorhombic, fc_ortho, pyrochlore) and `wannier90.py` continue importing from
  `stdface_model_util` via re-exports — no import changes needed for these
- `stdface_model_util.py` reduced from ~1,335 to ~1,125 lines

**Why**: Sixth step of Phase 2. The 10 parameter validation/printing utilities are
the most widely imported functions in the codebase (used by all 9 lattice modules,
all 4 writer modules, the main entry point, keyword parser, and the Wannier90
exporter). They have zero dependency on lattice or interaction logic — only `math`,
`sys`, and `numpy`. Extracting them into a dedicated module creates a clean
dependency leaf that other modules can import without pulling in the heavier
`stdface_model_util`. The re-export in `stdface_model_util` ensures full backward
compatibility with no import changes needed in the lattice modules.

**New tests**: `test/unit/test_param_check.py` — 29 tests covering:
- `TestExitProgram` (2 tests): exit with code -1, exit with code 0
- `TestPrintValD` (4 tests): NaN→default, set value, default tag printed, no tag when set
- `TestPrintValDD` (3 tests): primary default, secondary default, set value
- `TestPrintValC` (2 tests): NaN→default, set value
- `TestPrintValI` (3 tests): sentinel→default, set value, default tag
- `TestNotUsedD` (4 tests): NaN ok, set exits, complex NaN ok, complex set exits
- `TestNotUsedC` (2 tests): NaN ok, set exits
- `TestNotUsedJ` (3 tests): all NaN ok, scalar set exits, matrix element set exits
- `TestNotUsedI` (2 tests): sentinel ok, set exits
- `TestRequiredValI` (3 tests): missing exits, present ok, prints value
- `TestBackwardCompatibility` (1 test): all 10 functions importable from both
  `param_check` and `stdface_model_util`, verified as identical objects

**Test results**:
- Unit tests: 468 passed (439 existing + 29 new)
- Integration check: HPhi fulldiag_hubbard_chain_longrange — all 12 output files match reference

**Files created**: `python/param_check.py`, `test/unit/test_param_check.py`
**Files modified**: `python/stdface_model_util.py`, `python/stdface_main.py`,
`python/keyword_parser.py`, `python/export_wannier90.py`, `python/solver_writer.py`,
`python/common_writer.py`, `python/hphi_writer.py`, `python/mvmc_variational.py`,
`python/history/refactoring_log.md`

---

### 2026-01-27 — Phase 2, Step 7: Extract input parameter resolution helpers into new module

**Target**: `stdface_model_util.py` → new `input_params.py`

**What changed**:
- Created `python/input_params.py` (~210 lines) containing 4 input parameter
  resolution functions extracted from `stdface_model_util.py`:
  - `input_spin_nn(J, JAll, J0, J0All, J0name)` — resolve nearest-neighbour
    spin-spin interaction, handling J/J0 anisotropic/isotropic conflicts
  - `input_spin(Jp, JpAll, Jpname)` — resolve beyond-nearest-neighbour spin
    interactions
  - `input_coulomb_v(V, V0, V0name)` — resolve off-site Coulomb interaction
  - `input_hopp(t, t0, t0name)` — resolve hopping integral
- Also extracted a shared module-level constant `_SUFFIXES` (the 3×3 spin-component
  suffix matrix) that was previously duplicated as a local variable in both
  `input_spin_nn` and `input_spin`
- Removed these 4 functions (~167 lines) from `stdface_model_util.py`, replaced
  with a comment noting the new location
- Added re-export in `stdface_model_util.py` for backward compatibility:
  ```python
  from input_params import (  # noqa: F401
      input_spin_nn, input_spin, input_coulomb_v, input_hopp,
  )
  ```
- All 9 lattice modules that import these functions continue to use the
  re-exports from `stdface_model_util` — no import changes needed
- `stdface_model_util.py` reduced from ~1,125 to ~960 lines

**Why**: Seventh step of Phase 2. These 4 functions form a cohesive group that
resolves user-specified input parameters into concrete values, handling conflicts
between isotropic (scalar) and anisotropic (matrix) specifications. They are used
by all 9 lattice modules but not by any writer, solver, or variational module.
They depend only on `math` and `param_check.exit_program` — no lattice or
interaction logic. Extracting them reduces `stdface_model_util.py`'s scope to
its core responsibilities: transfer/interaction term building, site initialization,
and geometry output.

**New tests**: `test/unit/test_input_params.py` — 21 tests covering:
- `TestInputSpinNN` (8 tests): all-NaN→zeros, J0All sets diagonal, JAll sets
  diagonal, J matrix copies to J0, explicit J0 takes priority, JAll/J0All conflict,
  JAll/J element conflict, J0/J element conflict
- `TestInputSpin` (4 tests): all-NaN→zeros, JpAll sets diagonal, explicit
  elements kept, JpAll/Jp conflict
- `TestInputCoulombV` (4 tests): both NaN→0, V0 explicit, V fallback, conflict
- `TestInputHopp` (4 tests): both NaN→0+0j, t0 explicit, t fallback, conflict
- `TestBackwardCompatibility` (1 test): all 4 functions importable from both
  `input_params` and `stdface_model_util`, verified as identical objects

**Test results**:
- Unit tests: 489 passed (468 existing + 21 new)
- Integration check: HPhi fulldiag_hubbard_chain_longrange — all 12 output files match reference

**Files created**: `python/input_params.py`, `test/unit/test_input_params.py`
**Files modified**: `python/stdface_model_util.py`, `python/history/refactoring_log.md`

---

### 2026-01-27 — Phase 2, Step 8: Extract geometry output functions into new module

**Target**: `stdface_model_util.py` → new `geometry_output.py`

**What changed**:
- Created `python/geometry_output.py` (~140 lines) containing 2 geometry output
  functions extracted from `stdface_model_util.py`:
  - `print_xsf(StdI)` — write ``lattice.xsf`` in XCrysDen format (primitive vectors,
    optional conventional vectors for orthorhombic/fco/pyrochlore, atomic coordinates)
  - `print_geometry(StdI)` — write ``geometry.dat`` for correlation-function
    post-processing (direct vectors, phases, box matrix, site coordinates; skips for
    HWAVE uhfk/rpa modes; doubles sites for kondo model)
- Removed these 2 functions (~93 lines) from `stdface_model_util.py`, replaced with
  a comment noting the new location
- Added re-export in `stdface_model_util.py` for backward compatibility:
  ```python
  from geometry_output import (  # noqa: F401
      print_xsf, print_geometry,
  )
  ```
- All 9 lattice modules and `wannier90.py` continue to use the re-exports from
  `stdface_model_util` — no import changes needed
- `stdface_model_util.py` reduced from ~960 to ~870 lines

**Why**: Eighth step of Phase 2. The 2 geometry output functions are self-contained
file-writing operations with no dependency on interaction building, site initialization,
or parameter validation logic. They depend only on `StdIntList`. They are used by
all 9 lattice modules and `wannier90.py`. Extracting them continues the decomposition
of `stdface_model_util.py`, which now contains only: transfer/interaction term
builders, site initialization/folding, set_label, and memory allocation.

**New tests**: `test/unit/test_geometry_output.py` — 12 tests covering:
- `TestPrintGeometry` (6 tests): file creation, direct vectors content, phase line,
  cell count (L=4), kondo doubles cell lines, HWAVE uhfk skips output
- `TestPrintXsf` (5 tests): file creation, CRYSTAL/PRIMVEC/PRIMCOORD headers,
  no CONVVEC for chain, CONVVEC present for orthorhombic, PRIMCOORD atom count
- `TestBackwardCompatibility` (1 test): both functions importable from both
  `geometry_output` and `stdface_model_util`, verified as identical objects

**Test results**:
- Unit tests: 501 passed (489 existing + 12 new)
- Integration check: HPhi fulldiag_hubbard_chain_longrange — all 12 output files
  match reference (including geometry.dat itself)

**Files created**: `python/geometry_output.py`, `test/unit/test_geometry_output.py`
**Files modified**: `python/stdface_model_util.py`, `python/history/refactoring_log.md`

---

### 2026-01-27 — Phase 2, Step 9: Extract interaction builder functions into new module

**Target**: `stdface_model_util.py` → new `interaction_builder.py`

**What changed**:
- Created `python/interaction_builder.py` (~370 lines) containing 8 interaction term
  builder functions extracted from `stdface_model_util.py`:
  - `trans(StdI, trans0, isite, ispin, jsite, jspin)` — add a single transfer
    (one-body) term; skips values below 1e-12 threshold
  - `hopping(StdI, trans0, isite, jsite, dR)` — add hopping for both spin channels;
    supports HPhi time-evolution pump mode
  - `hubbard_local(StdI, mu0, h0, Gamma0, Gamma0_y, U0, isite)` — add chemical
    potential, longitudinal/transverse magnetic field, and intra-site Coulomb U
  - `mag_field(StdI, S2, h, Gamma, Gamma_y, isite)` — add magnetic field terms for
    arbitrary spin S using Bogoliubov representation
  - `intr(StdI, intr0, site1, spin1, ..., site4, spin4)` — add a general two-body
    (InterAll) interaction term; skips values below 1e-12
  - `general_j(StdI, J, Si2, Sj2, isite, jsite)` — treat J as a 3×3 matrix for
    general spin interactions; uses optimised Hund/Exchange/PairLift shortcuts for
    S=1/2 and falls back to general InterAll for higher spin
  - `coulomb(StdI, V, isite, jsite)` — add off-site Coulomb interaction
  - `malloc_interactions(StdI, ntransMax, nintrMax)` — allocate all interaction arrays
    (transfer, InterAll, Coulomb intra/inter, Hund, Exchange, PairLift, PairHopp,
    and HPhi pump arrays for time-evolution mode)
- Removed these 8 functions (~360 lines) from `stdface_model_util.py`, replaced with
  a comment noting the new location
- Added re-export in `stdface_model_util.py` for backward compatibility:
  ```python
  from interaction_builder import (  # noqa: F401
      trans, hopping, hubbard_local, mag_field,
      intr, general_j, coulomb, malloc_interactions,
  )
  ```
- All 9 lattice modules and `wannier90.py` continue to use the re-exports from
  `stdface_model_util` — no import changes needed
- `stdface_model_util.py` reduced from ~870 to ~445 lines

**Why**: Ninth step of Phase 2. These 8 functions form a tightly cohesive group:
`malloc_interactions` allocates the arrays that `trans`, `intr`, and `coulomb` populate;
`hopping`, `hubbard_local`, and `mag_field` are higher-level wrappers that call `trans`;
`general_j` is a higher-level wrapper that calls `intr`. Together they constitute the
core Hamiltonian term builder — the central purpose of the StdFace library. Extracting
them into a dedicated module with a clear name makes their role explicit and reduces
`stdface_model_util.py` to only site/cell utility functions (`_fold_site`, `init_site`,
`find_site`, `set_label`).

**New tests**: `test/unit/test_interaction_builder.py` — 28 tests covering:
- `TestMallocInteractions` (5 tests): transfer arrays shape, interaction arrays shape,
  Coulomb arrays shape, exchange/pair arrays shape, pump arrays for time-evolution
- `TestTrans` (4 tests): adds transfer term, skips small value, increments counter,
  complex value stored
- `TestHopping` (2 tests): adds both spins (4 terms), hermitian conjugate pairs
- `TestHubbardLocal` (4 tests): chemical potential, magnetic field splitting,
  intra-Coulomb U, transverse field Gamma
- `TestMagField` (3 tests): S=1/2 longitudinal, S=1/2 transverse, S=1 (3 states)
- `TestIntr` (2 tests): adds interaction term, skips small value
- `TestGeneralJ` (5 tests): isotropic S=1/2 (Hund/Cinter/Ex/PairLift), Jzz-only,
  S=1 general path, kondo exchange sign, non-kondo exchange sign
- `TestCoulomb` (2 tests): adds Coulomb term, increments counter
- `TestBackwardCompatibility` (1 test): all 8 functions importable from both
  `interaction_builder` and `stdface_model_util`, verified as identical objects

**Test results**:
- Unit tests: 529 passed (501 existing + 28 new)
- Integration check: HPhi fulldiag_hubbard_chain_longrange — all 12 output files
  match reference

**Files created**: `python/interaction_builder.py`, `test/unit/test_interaction_builder.py`
**Files modified**: `python/stdface_model_util.py`, `python/history/refactoring_log.md`

**Suggested next step**: Extract the remaining site utility functions (`_fold_site`,
`init_site`, `find_site`, `set_label`) from `stdface_model_util.py` into a dedicated
`site_util.py` module. This would complete the decomposition of `stdface_model_util.py`,
leaving it as a pure re-export hub for backward compatibility.

---

### 2026-01-27 — Phase 2, Step 10: Extract site utilities into new module

**Target**: `stdface_model_util.py` → new `site_util.py`

**What changed**:
- Created `python/site_util.py` (~340 lines) containing 4 site utility functions
  extracted from `stdface_model_util.py`:
  - `_fold_site(StdI, iCellV)` — fold a fractional coordinate into the original
    super-cell; returns the super-cell index and folded coordinate
  - `init_site(StdI, fp, dim)` — initialise the super-cell: validate L/W/Height vs
    box parameters, compute reciprocal lattice, find cells, write gnuplot header
  - `find_site(StdI, iW, iL, iH, diW, diL, diH, isiteUC, jsiteUC)` — find global
    site indices and boundary phase factor for a pair of sites
  - `set_label(StdI, fp, iW, iL, diW, diL, isiteUC, jsiteUC, connect)` — write
    gnuplot labels and arrows for 2-D lattice visualisation
- Rewrote `python/stdface_model_util.py` as a **pure re-export hub** (~77 lines):
  - Module docstring updated to describe it as a backward-compatibility re-export hub
  - Contains only `from X import ...  # noqa: F401` statements re-exporting all
    30 public symbols from 5 sub-modules: `param_check` (10), `input_params` (4),
    `geometry_output` (2), `interaction_builder` (8), `site_util` (4)
  - Removed all direct `import math`, `from typing import TextIO`, `import numpy`
    — these are no longer needed since no code remains in the module
  - All function definitions, section comments, and "moved to" comments removed
- `stdface_model_util.py` reduced from ~447 lines to **~77 lines** (83% reduction
  from Step 9; 96% reduction from original ~1,700 lines)

**Why**: Tenth and final step of Phase 2 decomposition. The 4 site utility functions
were the last remaining code in `stdface_model_util.py`. They form a cohesive group
dealing with super-cell geometry: `_fold_site` is called by `init_site` during cell
enumeration and by `find_site` when locating neighbour sites; `set_label` wraps
`find_site` with gnuplot output. Extracting them completes the full decomposition of
the formerly monolithic utility module into 5 focused, single-responsibility modules.

**New tests**: `test/unit/test_site_util.py` — 21 tests covering:
- `TestFoldSite` (4 tests): site inside cell stays, boundary wraps, negative wraps,
  interior site stays
- `TestInitSite` (6 tests): L/W/Height sets box, NCell computed, Cell array shape,
  tau array shape, anti-periodic phase flag, gnuplot header written
- `TestFindSite` (6 tests): same site, nearest neighbour, boundary wraps, kondo
  offset, dR vector, anti-periodic phase gives Cphase=-1
- `TestSetLabel` (4 tests): returns site indices, writes gnuplot labels, suppresses
  arrow for connect>=3, fp=None no error
- `TestBackwardCompatibility` (1 test): all 4 functions importable from both
  `site_util` and `stdface_model_util`, verified as identical objects

**Test results**:
- Unit tests: 550 passed (529 existing + 21 new)
- Integration check: HPhi fulldiag_hubbard_chain_longrange — all 12 output files
  match reference

**Files created**: `python/site_util.py`, `test/unit/test_site_util.py`
**Files modified**: `python/stdface_model_util.py`, `python/history/refactoring_log.md`

**Summary of Phase 2 decomposition of `stdface_model_util.py`** (Steps 5–10):

| Step | Module | Functions | Lines extracted |
|------|--------|-----------|----------------|
| 5 | `mvmc_variational.py` | 5 mVMC functions | ~386 |
| 6 | `param_check.py` | 10 validation utilities | ~210 |
| 7 | `input_params.py` | 4 input resolvers | ~167 |
| 8 | `geometry_output.py` | 2 geometry writers | ~93 |
| 9 | `interaction_builder.py` | 8 term builders | ~360 |
| 10 | `site_util.py` | 4 site utilities | ~370 |
| **Total** | **5 new modules** | **33 functions** | **~1,586** |

`stdface_model_util.py` reduced from ~1,700 lines (original) to ~77 lines
(**pure re-export hub**, 96% reduction). All lattice modules continue to import from
`stdface_model_util` with zero changes — backward compatibility fully preserved.

**Suggested next step**: With Phase 2 module decomposition complete, consider
Phase 3 improvements: introduce `enum.Enum` for model/lattice/solver types,
replace sentinel values with `None`, or migrate lattice modules to import directly
from sub-modules instead of through the re-export hub.

---

### 2026-01-27 — Phase 2, Step 11: Migrate lattice module imports to direct sub-module imports

**Target**: All 11 files that import from `stdface_model_util` re-export hub

**What changed**:
- Replaced `from stdface_model_util import (...)` blocks in all 11 consumer files
  with focused imports directly from the 5 sub-modules:
  - `param_check` — validation/printing utilities
  - `input_params` — input parameter resolution helpers
  - `geometry_output` — geometry file writers
  - `interaction_builder` — Hamiltonian term builders
  - `site_util` — super-cell initialisation and site finding

- Files migrated (11 total):

  | File | param_check | input_params | geometry_output | interaction_builder | site_util |
  |------|-------------|--------------|-----------------|---------------------|-----------|
  | `chain_lattice.py` | exit_program, print_val_d, print_val_i, not_used_d/c/i/j, required_val_i | input_spin_nn, input_spin, input_hopp, input_coulomb_v | print_geometry | malloc_interactions, mag_field, general_j, hopping, coulomb, hubbard_local | init_site, set_label |
  | `square_lattice.py` | print_val_d, print_val_i, not_used_j/c/d/i | input_spin_nn, input_spin, input_hopp, input_coulomb_v | print_geometry | malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb | init_site, set_label |
  | `ladder.py` | exit_program, print_val_d, print_val_i, not_used_j/c/d/i, required_val_i | input_spin, input_hopp, input_coulomb_v | print_geometry | malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb | init_site, set_label |
  | `triangular_lattice.py` | print_val_d, print_val_i, not_used_j/c/d/i | input_spin_nn, input_spin, input_hopp, input_coulomb_v | print_geometry | malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb | init_site, set_label |
  | `honeycomb_lattice.py` | exit_program, print_val_d, print_val_i, not_used_j/c/d/i | input_spin_nn, input_spin, input_hopp, input_coulomb_v | print_geometry | malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb | init_site, set_label |
  | `kagome.py` | exit_program, print_val_d, print_val_i, not_used_j/c/d/i | input_spin_nn, input_spin, input_hopp, input_coulomb_v | print_geometry | malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb | init_site, set_label |
  | `orthorhombic.py` | print_val_d, print_val_c, print_val_i, not_used_j/c/d/i | input_spin_nn, input_spin, input_hopp, input_coulomb_v | print_geometry, print_xsf | malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb | init_site, find_site |
  | `fc_ortho.py` | print_val_d, print_val_i, not_used_j/c/d/i | input_spin_nn, input_spin, input_hopp, input_coulomb_v | print_geometry, print_xsf | malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb | init_site, find_site |
  | `pyrochlore.py` | print_val_d, print_val_i, not_used_j/c/d/i | input_spin_nn, input_spin, input_hopp, input_coulomb_v | print_geometry, print_xsf | malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb | init_site, find_site |
  | `wannier90.py` | exit_program, print_val_d, print_val_i, not_used_d | (none) | print_geometry, print_xsf | malloc_interactions, mag_field, general_j, hubbard_local, hopping, coulomb | init_site, find_site |
  | `mvmc_variational.py` | (already direct) | (none) | (none) | (none) | _fold_site, find_site |

- No logic changes — only import source paths changed
- `stdface_model_util.py` remains as a re-export hub for any future consumers

**Why**: Eleventh step of Phase 2. Now that all 5 sub-modules are stable and tested,
the lattice modules can import directly from them rather than going through the
re-export hub. This makes the dependency graph explicit: each lattice file's imports
clearly show which sub-modules it depends on, aiding navigation and reducing
coupling to the legacy monolithic module. The re-export hub remains available for
any third-party code that still uses the old import paths.

**New tests**: None — this is a pure import path change; existing tests cover all
function behaviour. The backward compatibility tests in each sub-module's test file
continue to verify that the re-export hub works.

**Test results**:
- Unit tests: 550 passed (no change)
- Integration check: all 10 lattice tests pass (chain, square, triangular, honeycomb,
  kagome, kondo_chain, kondo_honey, kondo_kagome, kondo_square, kondo_tri);
  Wannier90 tests have a pre-existing trans.def diff unrelated to this change

**Files modified**: `python/chain_lattice.py`, `python/square_lattice.py`,
`python/ladder.py`, `python/triangular_lattice.py`, `python/honeycomb_lattice.py`,
`python/kagome.py`, `python/orthorhombic.py`, `python/fc_ortho.py`,
`python/pyrochlore.py`, `python/wannier90.py`, `python/mvmc_variational.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Phase 2 is now complete. The codebase has been decomposed
from 2 monolithic modules into 13 focused modules with explicit dependency graphs.
Consider Phase 3 improvements: introduce `enum.Enum` for model/lattice/solver types,
replace sentinel values with `None`, or consolidate the `_reset_vals()` field
initialisations into `StdIntList.__init__`.

---

### 2026-01-27 — Phase 3, Step 1: Introduce ModelType enum for model name constants

**Target**: All 17 files with `StdI.model == "spin"/"hubbard"/"kondo"` comparisons

**What changed**:
- Created `ModelType(str, Enum)` class in `python/stdface_vals.py` with three members:
  - `ModelType.SPIN = "spin"` — pure spin model
  - `ModelType.HUBBARD = "hubbard"` — fermion Hubbard model
  - `ModelType.KONDO = "kondo"` — Kondo lattice model
- Inherits from both `str` and `Enum`, so `ModelType.SPIN == "spin"` is `True` —
  full backward compatibility with existing string comparisons
- Updated `MODEL_ALIASES` and `MODEL_ALIASES_HPHI_BOOST` in `stdface_main.py` to
  produce `ModelType` values instead of plain strings
- Replaced all 219 `StdI.model == "spin"/"hubbard"/"kondo"` string comparisons across
  17 source files with `ModelType.SPIN`/`ModelType.HUBBARD`/`ModelType.KONDO`
- Added `from stdface_vals import ModelType` to all 17 files

**Files updated** (17 source + 1 test):
- `python/stdface_vals.py` — added `ModelType` enum class
- `python/stdface_main.py` — updated `MODEL_ALIASES`/`MODEL_ALIASES_HPHI_BOOST` values
- `python/chain_lattice.py`, `python/square_lattice.py`, `python/ladder.py`,
  `python/triangular_lattice.py`, `python/honeycomb_lattice.py`, `python/kagome.py`,
  `python/orthorhombic.py`, `python/fc_ortho.py`, `python/pyrochlore.py` — lattice modules
- `python/wannier90.py`, `python/site_util.py`, `python/interaction_builder.py`,
  `python/geometry_output.py`, `python/hphi_writer.py`, `python/common_writer.py`,
  `python/mvmc_writer.py`, `python/mvmc_variational.py` — utility/writer modules
- `test/unit/test_stdface_vals.py` — added `TestModelType` class

**Why**: First step of Phase 3 (Python idioms). The model type string `"spin"`,
`"hubbard"`, or `"kondo"` is compared 219 times across 17 files — the single most
frequent string comparison in the codebase. Using a `str` enum:
- Prevents typos (e.g. `"spim"` would fail at import time, not silently)
- Provides IDE autocompletion and type checking
- Documents all valid model types in one place
- Retains full backward compatibility with existing string comparisons via `str` base class

**New tests**: 9 tests in `TestModelType`:
- `test_has_three_members` — exactly 3 enum members
- `test_values` — SPIN="spin", HUBBARD="hubbard", KONDO="kondo"
- `test_is_str_subclass` — all members are `str` instances
- `test_string_equality` — enum == string in both directions
- `test_string_inequality` — cross-type comparison
- `test_construction_from_string` — `ModelType("spin")` returns `ModelType.SPIN`
- `test_invalid_string_raises` — `ModelType("invalid")` raises `ValueError`
- `test_can_be_used_as_dict_key` — interchangeable with string dict keys
- `test_model_field_assignment` — `StdIntList.model = ModelType.HUBBARD` works

**Test results**:
- Unit tests: 559 passed (550 existing + 9 new)
- Integration check: all 6 representative cases pass (hubbard chain/honeycomb/square,
  kondo chain/honeycomb/square)

**Suggested next step**: Introduce `SolverType(str, Enum)` for the 35 solver string
comparisons, or `MethodType(str, Enum)` for the 15 method string comparisons.
Alternatively, consolidate sentinel constants into `stdface_vals.py`.

---

### 2026-01-27 — Phase 3, Step 2: Introduce SolverType enum for solver name constants

**Target**: All 12 files with `StdI.solver == "HPhi"/"mVMC"/"UHF"/"HWAVE"` comparisons

**What changed**:
- Created `SolverType(str, Enum)` class in `python/stdface_vals.py` with four members:
  - `SolverType.HPhi = "HPhi"` — exact-diagonalisation / Lanczos solver
  - `SolverType.mVMC = "mVMC"` — variational Monte Carlo solver
  - `SolverType.UHF = "UHF"` — unrestricted Hartree-Fock solver
  - `SolverType.HWAVE = "HWAVE"` — H-wave solver
- Inherits from both `str` and `Enum`, so `SolverType.HPhi == "HPhi"` is `True` —
  full backward compatibility with existing string comparisons
- Replaced all ~35 `StdI.solver == "HPhi"/"mVMC"/"UHF"/"HWAVE"` and
  `StdI.solver != "HWAVE"/"HPhi"` string comparisons across 12 source files with
  `SolverType` enum constants
- Updated `_SOLVER_RESET_DISPATCH` dict keys in `stdface_main.py` to use `SolverType`
- Updated `_SOLVER_WRITERS` dict keys and `super().__init__()` calls in `solver_writer.py`
- Updated `parse_solver_keyword()` dispatch in `keyword_parser.py`
- Added `from stdface_vals import SolverType` to all 12 files

**Files updated** (12 source + 1 test):
- `python/stdface_vals.py` — added `SolverType` enum class
- `python/stdface_main.py` — updated `_SOLVER_RESET_DISPATCH` keys, 3 comparisons
- `python/common_writer.py` — 14 comparisons (`==` HPhi/mVMC/UHF/HWAVE, `in` tuple)
- `python/interaction_builder.py` — 4 comparisons (`==` HPhi, `==` mVMC)
- `python/geometry_output.py` — 1 comparison (`==` HWAVE)
- `python/wannier90.py` — 1 comparison (`==` mVMC)
- `python/chain_lattice.py` — 2 `!=` comparisons (HWAVE, HPhi)
- `python/square_lattice.py` — 2 `!=` comparisons (HWAVE)
- `python/ladder.py` — 2 `!=` comparisons (HWAVE)
- `python/triangular_lattice.py` — 2 `!=` comparisons (HWAVE)
- `python/honeycomb_lattice.py` — 2 `!=` comparisons (HWAVE)
- `python/kagome.py` — 2 `!=` comparisons (HWAVE)
- `python/solver_writer.py` — import + registry dict keys + 4 `super().__init__()` calls
- `python/keyword_parser.py` — import + 4 comparisons in `parse_solver_keyword()`
- `test/unit/test_stdface_vals.py` — added `TestSolverType` class

**Why**: Second step of Phase 3 (Python idioms). The solver type string `"HPhi"`,
`"mVMC"`, `"UHF"`, or `"HWAVE"` is compared ~35 times across 12 files. Using a
`str` enum:
- Prevents typos (e.g. `"Hphi"` would be caught by IDE/type-checker, not silently
  fall through to the wrong branch)
- Provides IDE autocompletion and type checking
- Documents all valid solver types in one place
- Retains full backward compatibility with existing string comparisons via `str` base class

**New tests**: 9 tests in `TestSolverType`:
- `test_has_four_members` — exactly 4 enum members
- `test_values` — HPhi="HPhi", mVMC="mVMC", UHF="UHF", HWAVE="HWAVE"
- `test_is_str_subclass` — all members are `str` instances
- `test_string_equality` — enum == string in both directions
- `test_string_inequality` — cross-type comparison
- `test_construction_from_string` — `SolverType("HPhi")` returns `SolverType.HPhi`
- `test_invalid_string_raises` — `SolverType("invalid")` and `SolverType("hphi")`
  raise `ValueError`
- `test_can_be_used_as_dict_key` — interchangeable with string dict keys
- `test_solver_field_assignment` — `StdIntList.solver = SolverType.HPhi` works

**Test results**:
- Unit tests: 568 passed (559 existing + 9 new)
- Integration check: all 6 representative cases pass (hubbard chain/square/tri/kagome,
  kondo chain, spinaniso honeycomb)

**Suggested next step**: Introduce `MethodType(str, Enum)` for the ~15 method string
comparisons (`"fulldiag"`, `"timeevolution"`, etc.), or consolidate sentinel constants
into `stdface_vals.py`.

---

### 2026-01-27 — Phase 3, Step 3: Introduce MethodType enum for method name constants

**Target**: All 5 files with `StdI.method == "lanczos"/"timeevolution"/...` comparisons

**What changed**:
- Created `MethodType(str, Enum)` class in `python/stdface_vals.py` with seven members:
  - `MethodType.LANCZOS = "lanczos"` — Lanczos diagonalisation
  - `MethodType.LANCZOS_ENERGY = "lanczosenergy"` — Lanczos (energy-only)
  - `MethodType.TPQ = "tpq"` — Thermal Pure Quantum state
  - `MethodType.FULLDIAG = "fulldiag"` — Full diagonalisation
  - `MethodType.CG = "cg"` — Conjugate Gradient
  - `MethodType.TIME_EVOLUTION = "timeevolution"` — Real-time evolution
  - `MethodType.CTPQ = "ctpq"` — Canonical TPQ
- Inherits from both `str` and `Enum`, so `MethodType.LANCZOS == "lanczos"` is `True` —
  full backward compatibility with existing string comparisons
- Updated `METHOD_ALIASES` dict in `stdface_main.py` to produce `MethodType` values
  instead of plain strings
- Replaced all 16 `StdI.method == "lanczos"/"timeevolution"/...` string comparisons
  across 5 source files with `MethodType` enum constants
- Added `from stdface_vals import MethodType` to all 5 files
- The sentinel check `StdI.method == "****"` remains a plain string (sentinels are a
  separate concern)

**Files updated** (5 source + 1 test):
- `python/stdface_vals.py` — added `MethodType` enum class
- `python/stdface_main.py` — updated import, `METHOD_ALIASES` values, 1 comparison
- `python/hphi_writer.py` — updated import, 10 comparisons (7 method dispatch, 1 `in`
  tuple, 1 time-evolution check, plus 1 sentinel that stays as-is)
- `python/interaction_builder.py` — updated import, 2 comparisons
- `python/common_writer.py` — updated import, 2 comparisons
- `python/solver_writer.py` — updated import, 1 comparison
- `test/unit/test_stdface_vals.py` — added `TestMethodType` class

**Why**: Third step of Phase 3 (Python idioms). The method type string `"lanczos"`,
`"timeevolution"`, `"fulldiag"`, etc. is compared 16 times across 5 files. Using a
`str` enum:
- Prevents typos (e.g. `"timeevalution"` would be caught by IDE/type-checker)
- Provides IDE autocompletion and type checking
- Documents all valid method types in one place
- Retains full backward compatibility with existing string comparisons via `str` base class

**New tests**: 9 tests in `TestMethodType`:
- `test_has_seven_members` — exactly 7 enum members
- `test_values` — all 7 members have expected string values
- `test_is_str_subclass` — all members are `str` instances
- `test_string_equality` — enum == string in both directions
- `test_string_inequality` — cross-type comparison
- `test_construction_from_string` — all 7 members constructible from string
- `test_invalid_string_raises` — invalid/case-wrong strings raise `ValueError`
- `test_can_be_used_as_dict_key` — interchangeable with string dict keys
- `test_method_field_assignment` — `StdIntList.method = MethodType.FULLDIAG` works

**Test results**:
- Unit tests: 577 passed (568 existing + 9 new)
- Integration check: all 5 representative cases pass:
  - `fulldiag_hubbard_chain_longrange` (10 files match)
  - `fulldiag_spingcboost_chain_success` (9 files match)
  - `te_hubbard_chain_phase90` (9 files match — exercises timeevolution path)
  - `fulldiag_kondo_chain_longrange` (12 files match)
  - `fulldiag_spinaniso_honeycomb` (12 files match)

**Suggested next step**: Consolidate sentinel constants (`NaN_i`, `NaN_d`, `NaN_c`,
`"****"`) into `stdface_vals.py`, or introduce enums for remaining string parameters
(`CalcSpec`, `Restart`, `PumpType`, `SpectrumType`).

---

### 2026-01-27 — Phase 3, Step 4: Consolidate sentinel constants into `stdface_vals.py`

**Target**: 7 redundant sentinel definitions across 4 files + 15 `"****"` comparisons across 6 files

**What changed**:
- Added four module-level sentinel constants to `python/stdface_vals.py`:
  - `NaN_i: int = 2147483647` — sentinel for unset integer parameter (INT_MAX in C)
  - `NaN_d: float = float("nan")` — sentinel for unset float parameter (IEEE NaN)
  - `NaN_c: complex = complex(float("nan"), 0.0)` — sentinel for unset complex parameter
  - `UNSET_STRING: str = "****"` — sentinel for unset string parameter (C convention)
- Removed redundant `NaN_i` definitions from 3 files:
  - `stdface_main.py` (line 148) — module-level definition removed; now imports from `stdface_vals`
  - `solver_writer.py` (line 60) — module-level definition removed; now imports from `stdface_vals`
  - `keyword_parser.py` (line 15) — module-level definition removed; now imports from `stdface_vals`
- Removed redundant `NaN_d` and `NaN_c` definitions from `stdface_main.py` (lines 151, 154) —
  now imports from `stdface_vals`
- Removed 3 function-local `NaN_i = 2147483647` definitions in `param_check.py`
  (`print_val_i`, `not_used_i`, `required_val_i`) — now uses module-level import from `stdface_vals`
- Replaced all 15 `"****"` comparisons in executable code across 6 files with `UNSET_STRING`:
  - `stdface_main.py` — 18 assignments + 1 comparison
  - `hphi_writer.py` — 9 comparisons
  - `common_writer.py` — 2 comparisons
  - `keyword_parser.py` — 2 comparisons
  - `wannier90.py` — 1 comparison
  - `export_wannier90.py` — 1 comparison
- `"****"` literals in docstrings were left unchanged (documentation references)
- Updated all import statements to source sentinel constants from `stdface_vals`

**Files updated** (8 source + 2 test):
- `python/stdface_vals.py` — added `NaN_i`, `NaN_d`, `NaN_c`, `UNSET_STRING` module constants
- `python/stdface_main.py` — removed local definitions; imports from `stdface_vals`; `"****"` → `UNSET_STRING`
- `python/param_check.py` — added `from stdface_vals import NaN_i`; removed 3 function-local definitions
- `python/solver_writer.py` — removed local `NaN_i`; imports from `stdface_vals`
- `python/keyword_parser.py` — removed local `NaN_i`; imports `NaN_i`, `UNSET_STRING` from `stdface_vals`
- `python/hphi_writer.py` — imports `UNSET_STRING`; 9 comparisons updated
- `python/common_writer.py` — imports `UNSET_STRING`; 2 comparisons updated
- `python/wannier90.py` — imports `UNSET_STRING`; 1 comparison updated
- `python/export_wannier90.py` — imports `UNSET_STRING`; 1 comparison updated
- `test/unit/test_stdface_vals.py` — added `TestSentinelConstants` class
- `test/unit/test_param_check.py` — updated to import `NaN_i`, `NaN_d` from `stdface_vals`

**Why**: Fourth step of Phase 3 (Python idioms). The sentinel constants were defined
redundantly in 7 places across 4 files (plus the `StdIntList.NaN_i` dataclass field).
The `"****"` string literal was scattered across 6 files with no named constant.
Consolidating all sentinels in `stdface_vals.py`:
- Single source of truth — changing a sentinel value requires editing one line
- Named constant `UNSET_STRING` is self-documenting (vs magic string `"****"`)
- Eliminates risk of divergent definitions (e.g. if one file used a different value)
- `param_check.py` no longer shadows the module-level constant with function-local copies

**New tests**: 11 tests in `TestSentinelConstants`:
- `test_nan_i_value` — equals 2147483647
- `test_nan_i_type` — is `int`
- `test_nan_d_is_nan` — is IEEE NaN
- `test_nan_d_type` — is `float`
- `test_nan_c_real_is_nan` — real part is NaN
- `test_nan_c_imag_is_zero` — imaginary part is 0.0
- `test_nan_c_type` — is `complex`
- `test_unset_string_value` — equals `"****"`
- `test_unset_string_type` — is `str`
- `test_stdintlist_nan_i_matches` — `StdIntList().NaN_i == NaN_i`
- `test_sentinels_are_distinct` — different types for the three numeric sentinels

**Test results**:
- Unit tests: 588 passed (577 existing + 11 new)
- Integration check: 5/5 representative cases pass:
  - `fulldiag_hubbard_chain_longrange` (10 files match)
  - `fulldiag_spingcboost_chain_success` (9 files match)
  - `te_hubbard_chain_phase90` (9 files match)
  - `fulldiag_kondo_chain_longrange` (12 files match)
  - `fulldiag_spinaniso_honeycomb` (12 files match)
  - `fulldiag_wannier90_alphaBETS` (8/9 match; trans.def has pre-existing
    `-0.0` vs `0.0` sign diff, not a regression)

**Suggested next step**: Introduce enums for remaining string parameters that are
compared multiple times (`CalcSpec`, `Restart`, `PumpType`, `SpectrumType`), or
replace `StdI.NaN_i` field accesses with direct `NaN_i` constant imports across
lattice/writer modules.

---

### 2026-01-27 — Phase 3, Step 5: Replace `StdI.NaN_i` field accesses with direct `NaN_i` imports

**Target**: 9 files accessing the sentinel via `StdI.NaN_i` field indirection

**What changed**:
- Replaced all `StdI.NaN_i` comparisons (23 occurrences across 7 files) with the
  module-level `NaN_i` constant imported from `stdface_vals`
- Removed 4 local aliases `NaN_i = StdI.NaN_i` that existed as function-local caching
  of the field value:
  - `site_util.py` `init_site()` (line 98)
  - `wannier90.py` `_read_wannier90()` (line 201)
  - `wannier90.py` `wannier90()` (line 567)
  - `mvmc_variational.py` `_init_site_sub()` (line 144)
- Added `NaN_i` to the `from stdface_vals import` line in 7 files:
  - `site_util.py`, `wannier90.py`, `mvmc_variational.py`, `common_writer.py`,
    `hphi_writer.py`, `export_wannier90.py`, `mvmc_writer.py`
- Retained the single `StdI.NaN_i = NaN_i` assignment in `stdface_main.py` `_reset_vals()`
  (sets the dataclass field; removing it would be a semantic change)

**Why**: Fifth step of Phase 3 (Python idioms). After Step 4 consolidated the sentinel
constants into `stdface_vals.py`, many files still accessed `NaN_i` through the
`StdI` dataclass field (`StdI.NaN_i`) rather than importing it directly. This
indirection:
- Obscured the fact that `NaN_i` is a compile-time constant, not a per-instance value
- Required 4 functions to create local aliases `NaN_i = StdI.NaN_i` as an optimisation
- Made the code harder to understand (readers had to check whether `StdI.NaN_i`
  could vary between instances)

Using the module-level constant directly makes the code clearer and avoids the
unnecessary field access.

**New tests**: None — pure mechanical replacement; existing 588 tests cover all
affected code paths.

**Test results**:
- Unit tests: 588 passed (no change)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange` (10 files match)
  - `fulldiag_spingcboost_chain_success` (9 files match)
  - `te_hubbard_chain_phase90` (9 files match)
  - `fulldiag_kondo_chain_longrange` (12 files match)
  - `fulldiag_spinaniso_honeycomb` (12 files match)

**Files modified**: `python/site_util.py`, `python/wannier90.py`,
`python/mvmc_variational.py`, `python/common_writer.py`, `python/hphi_writer.py`,
`python/export_wannier90.py`, `python/mvmc_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Introduce enums for remaining string parameters that are
compared multiple times (`CalcSpec`, `Restart`, `PumpType`, `SpectrumType`), or
remove the now-redundant `NaN_i` field from the `StdIntList` dataclass (requires
auditing all code that reads `StdI.NaN_i`).

---

### 2026-01-27 — Phase 3, Step 6: Dict dispatch tables in `print_calc_mod()`

**Target**: `hphi_writer.py` — replace if/elif dispatch chains with dict lookups

**What changed**:
- Added 3 module-level dispatch dictionaries in `hphi_writer.py`:
  - `METHOD_TO_CALC_TYPE: dict[MethodType, int]` — maps 7 `MethodType` enum members
    to their `CalcType` integer codes
  - `RESTART_TO_INT: dict[str, int]` — maps 6 Restart string values (incl. aliases)
    to integer codes
  - `CALC_SPEC_TO_INT: dict[str, int]` — maps 7 CalcSpec string values (incl. aliases)
    to integer codes
- Rewrote 3 if/elif chains in `print_calc_mod()` to use `dict.get()` lookups:
  - Method → CalcType (was 7-branch if/elif, now `METHOD_TO_CALC_TYPE.get()`)
  - Restart → iRestart (was 5-branch if/elif, now `RESTART_TO_INT.get()`)
  - CalcSpec → iCalcSpec (was 7-branch if/elif, now `CALC_SPEC_TO_INT.get()`)
- Error handling preserved: `.get()` returns `None` on unknown keys, triggering the
  same error messages and `exit_program(-1)` as before.

**Why**: Dict dispatch is idiomatic Python for mapping string/enum keys to values.
The if/elif chains were direct transliterations from C `strcmp()` cascades. Dict
lookups are more concise, easier to maintain (add a new entry = add one dict entry),
and their mapping is visible at module level rather than buried inside a function.

**New tests**: 23 tests in `TestDispatchDicts` class (`test/unit/test_hphi_writer.py`):
- `METHOD_TO_CALC_TYPE`: 8 tests (7 known methods + unknown-returns-None)
- `RESTART_TO_INT`: 7 tests (6 known strings incl. aliases + unknown-returns-None)
- `CALC_SPEC_TO_INT`: 8 tests (7 known strings incl. aliases + unknown-returns-None)

**Test results**:
- Unit tests: 611 passed (588 → 611, +23 new dispatch dict tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_fcortho_smoke`

**Files modified**: `python/hphi_writer.py`, `test/unit/test_hphi_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Apply the same dict-dispatch pattern to `print_calc_mod()`'s
remaining if/elif chains (Model, InitialVecType, EigenVecIO, HamIO, OutputExVec) or
to other functions with string→int dispatch (e.g., `SpectrumType` in
`print_excitation()`, `PumpType` in `vector_potential()`).

---

### 2026-01-27 — Phase 3, Step 7: Complete dict dispatch in `print_calc_mod()`

**Target**: `hphi_writer.py` — convert 5 remaining if/elif chains to dict lookups

**What changed**:
- Added 5 module-level dispatch dictionaries in `hphi_writer.py`:
  - `MODEL_GC_TO_CALC_MODEL: dict[tuple, int]` — maps `(ModelType, lGC)` pairs to
    `CalcModel` integer codes (6 entries)
  - `INITIAL_VEC_TYPE_TO_INT: dict[str, int]` — maps `"c"→0`, `"r"→1`
  - `EIGENVEC_IO_TO_FLAGS: dict[str, tuple[int, int]]` — maps EigenVecIO string to
    `(InputEigenVec, OutputEigenVec)` flag pairs (4 entries)
  - `HAM_IO_TO_FLAGS: dict[str, tuple[int, int]]` — maps HamIO string to
    `(iOutputHam, iInputHam)` flag pairs (3 entries)
  - `OUTPUT_EX_VEC_TO_INT: dict[str, int]` — maps `"none"→0`, `"out"→1`
- Rewrote 5 if/elif chains in `print_calc_mod()` to use `dict.get()`:
  - Model: was 3-branch if/elif with ternary, now `MODEL_GC_TO_CALC_MODEL.get()`
  - InitialVecType: was 2-branch, now `INITIAL_VEC_TYPE_TO_INT.get()`
  - EigenVecIO: was 4-branch with dual assignments, now `EIGENVEC_IO_TO_FLAGS.get()`
    returning `(InputEigenVec, OutputEigenVec)` tuple
  - HamIO: was 3-branch with dual assignments, now `HAM_IO_TO_FLAGS.get()`
    returning `(iOutputHam, iInputHam)` tuple
  - OutputExVec: was 2-branch, now `OUTPUT_EX_VEC_TO_INT.get()`
- All 8 if/elif chains in `print_calc_mod()` now use dict dispatch (3 from Step 6 +
  5 from this step). The function body is more concise and all mappings are visible
  as module-level constants.

**Why**: Completes the dict-dispatch conversion of `print_calc_mod()` started in
Step 6. For multi-value outputs (EigenVecIO, HamIO), tuple-valued dicts keep
the mapping compact and avoid scattered variable assignments. The `(Model, lGC)`
composite-key dict replaces a pattern that was hard to extend (3 if/elif branches
each with a ternary).

**New tests**: 22 tests added to `TestDispatchDicts` in `test/unit/test_hphi_writer.py`:
- `MODEL_GC_TO_CALC_MODEL`: 7 tests (6 known combos + unknown-returns-None)
- `INITIAL_VEC_TYPE_TO_INT`: 3 tests (2 values + unknown-returns-None)
- `EIGENVEC_IO_TO_FLAGS`: 5 tests (4 values + unknown-returns-None)
- `HAM_IO_TO_FLAGS`: 4 tests (3 values + unknown-returns-None)
- `OUTPUT_EX_VEC_TO_INT`: 3 tests (2 values + unknown-returns-None)

**Test results**:
- Unit tests: 633 passed (611 → 633, +22 new dispatch dict tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_fcortho_smoke`

**Files modified**: `python/hphi_writer.py`, `test/unit/test_hphi_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Apply dict-dispatch to other functions with string→int
dispatch chains: `SpectrumType` in `print_excitation()`, `PumpType` in
`vector_potential()`. Or move on to a different Phase 3 task such as removing the
redundant `NaN_i` field from the `StdIntList` dataclass.

---

### 2026-01-27 — Phase 3, Step 8: Remove redundant `NaN_i` field from `StdIntList`

**Target**: `stdface_vals.py` dataclass and all code that referenced `StdI.NaN_i`

**What changed**:
- Removed the `NaN_i: int = 2147483647` field from the `StdIntList` dataclass in
  `stdface_vals.py` (line 600). The sentinel is now only the module-level constant
  `NaN_i` in `stdface_vals`, which all code already imports directly since Step 5.
- Removed the `StdI.NaN_i = NaN_i` assignment from `_reset_vals()` in
  `stdface_main.py` (line 310).
- Removed 15 `StdI.NaN_i = NaN_i` (or literal `= 2147483647`) setup lines from
  10 test files:
  - `test_mvmc_writer.py` (1)
  - `test_keyword_parser.py` (1)
  - `test_mvmc_variational.py` (1)
  - `test_solver_writer.py` (1)
  - `test_wannier90.py` (4)
  - `test_reset_vals_dispatch.py` (1)
  - `test_common_writer.py` (1)
  - `test_site_util.py` (2)
  - `test_geometry_output.py` (1)
  - `test_stdface_vals.py` (2 tests removed: `test_sentinel_value`,
    `test_stdintlist_nan_i_matches`)
- Updated `StdIntList` docstring to remove the `NaN_i` attribute entry.
- Renamed the dataclass section header from "Initial (undefined) / sentinel" to
  "Mathematical constants" since `pi` is now the only field in that block.

**Why**: The `NaN_i` field on the dataclass was a C-ism — C used a struct field
because there was no module-level constant facility. In Python, the sentinel is
properly a module-level constant (added in Step 4). After Step 5 replaced all
23 `StdI.NaN_i` read accesses with direct `NaN_i` imports, the field had zero
readers — only the one write in `_reset_vals()`. Removing it eliminates dead
state and prevents accidental reintroduction of field-access patterns.

**Tests removed**: 2 tests that explicitly asserted `StdIntList().NaN_i`:
- `TestStdIntListDefaults.test_sentinel_value`
- `TestSentinelConstants.test_stdintlist_nan_i_matches`

**Test results**:
- Unit tests: 631 passed (633 − 2 removed = 631)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_fcortho_smoke`

**Files modified**: `python/stdface_vals.py`, `python/stdface_main.py`,
`test/unit/test_mvmc_writer.py`, `test/unit/test_keyword_parser.py`,
`test/unit/test_mvmc_variational.py`, `test/unit/test_solver_writer.py`,
`test/unit/test_wannier90.py`, `test/unit/test_reset_vals_dispatch.py`,
`test/unit/test_common_writer.py`, `test/unit/test_site_util.py`,
`test/unit/test_geometry_output.py`, `test/unit/test_stdface_vals.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Replace the remaining `"****"` string sentinels in parameter
default values with `UNSET_STRING` (many `StdIntList` fields still default to `""`
rather than `UNSET_STRING`), or begin Phase 3 work on replacing `None`-compatible
sentinels (e.g., using `Optional[int]` with `None` instead of `NaN_i = 2147483647`
for unset integer parameters).

---

### 2026-01-27 — Phase 3, Step 9: Change string field defaults to `UNSET_STRING`

**Target**: `stdface_vals.py` — 17 string fields in `StdIntList` dataclass

**What changed**:
- Changed the default value of 17 string fields in the `StdIntList` dataclass from
  `""` to `UNSET_STRING` (`"****"`). These are all the string fields that
  `_reset_vals()` explicitly sets to `UNSET_STRING` — the dataclass defaults now
  match the runtime-initialized state:

  | Field | Section |
  |-------|---------|
  | `lattice` | LATTICE |
  | `model` | MODEL |
  | `outputmode` | Common |
  | `CDataFileHead` | Common |
  | `double_counting_mode` | Wannier90 |
  | `method` | HPhi |
  | `Restart` | HPhi |
  | `InitialVecType` | HPhi |
  | `EigenVecIO` | HPhi |
  | `HamIO` | HPhi |
  | `CalcSpec` | HPhi |
  | `SpectrumType` | HPhi |
  | `OutputExVec` | HPhi |
  | `PumpType` | HPhi |
  | `CParaFileHead` | mVMC |
  | `calcmode` | HWAVE |
  | `fileprefix` | HWAVE |

- Left `solver: str = ""` unchanged — it is not reset to `UNSET_STRING` (it
  receives its value from the function argument).
- Updated `test_string_defaults` in `test_stdface_vals.py` to assert the new defaults.

**Why**: In the original C code, `StdFace_ResetVals()` sets every string parameter
to `"****"` before parsing. The Python dataclass defaults of `""` didn't match this
convention — a freshly constructed `StdIntList()` had different string values from
one that had been through `_reset_vals()`. This made it easy to introduce bugs: code
that checked `== UNSET_STRING` would silently fail on a `StdIntList` that hadn't been
reset. Now the dataclass defaults match the intended "unset" sentinel, making the
`_reset_vals()` string assignments idempotent and the object self-consistent from
construction.

**New tests**: None — updated 1 existing test (`test_string_defaults`).

**Test results**:
- Unit tests: 631 passed (unchanged)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_fcortho_smoke`

**Files modified**: `python/stdface_vals.py`, `test/unit/test_stdface_vals.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Remove the now-redundant 17 `StdI.xxx = UNSET_STRING`
assignments from `_reset_vals()` / `_reset_hphi_fields()` / `_reset_mvmc_fields()`
/ `_reset_hwave_fields()` (they are idempotent since the dataclass defaults already
match). Or begin a larger Phase 3 task: replacing integer/float sentinel values
(`NaN_i`, `NaN_d`) with `Optional[int]`/`Optional[float]` using `None`.

---

### 2026-01-27 — Phase 3, Step 10: Remove redundant `UNSET_STRING` assignments from reset functions

**Target**: `stdface_main.py` — 4 reset functions

**What changed**:
- Removed 17 now-idempotent `StdI.xxx = UNSET_STRING` assignments across 4 functions:
  - `_reset_hphi_fields()`: removed 9 assignments (`method`, `Restart`,
    `EigenVecIO`, `InitialVecType`, `HamIO`, `CalcSpec`, `SpectrumType`,
    `OutputExVec`, `PumpType`)
  - `_reset_mvmc_fields()`: removed 1 assignment (`CParaFileHead`)
  - `_reset_hwave_fields()`: removed 2 assignments (`calcmode`, `fileprefix`)
  - `_reset_vals()`: removed 5 assignments (`model`, `lattice`, `outputmode`,
    `CDataFileHead`, `double_counting_mode`)
- Updated section comment from "Calculation conditions (strings / ints)" to
  "Calculation conditions" since the string fields are no longer listed there.
- Net reduction: 17 lines removed from `stdface_main.py`.
- `UNSET_STRING` import retained — still used at line 515 for a comparison.

**Why**: Step 9 changed the 17 `StdIntList` string field defaults from `""` to
`UNSET_STRING`. Since the dataclass constructor now initialises these fields to the
correct sentinel value, the explicit `StdI.xxx = UNSET_STRING` assignments in the
reset functions were performing a no-op (setting a field to the value it already
has). Removing them eliminates dead code and makes the reset functions shorter and
easier to read — they now only contain assignments that actually change a field's
value from the dataclass default.

**New tests**: None — pure dead-code removal; no behavior change.

**Test results**:
- Unit tests: 631 passed (unchanged)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_fcortho_smoke`

**Files modified**: `python/stdface_main.py`, `python/history/refactoring_log.md`

**Suggested next step**: The `UNSET_STRING` sentinel cleanup is now complete (Steps
4, 5, 8, 9, 10). The Phase 3 Python-idioms track can continue with:
- Apply dict-dispatch to `print_excitation()` (`SpectrumType` branches — though
  these are algorithmic, not simple mappings) or `vector_potential()` (`PumpType`)
- Replace integer/float sentinel values (`NaN_i`, `NaN_d`) with `Optional` types
  using `None` — a larger architectural change affecting many comparisons
- Extract remaining if/elif dispatch chains in `common_writer.py` or
  `solver_writer.py` into dict tables

---

### 2026-01-27 — Phase 3, Step 11: Dict dispatch for `check_output_mode()` in `common_writer.py`

**Target**: `common_writer.py` — `check_output_mode()` function

**What changed**:
- Added `OUTPUT_MODE_TO_INT: dict[str, int]` module-level dispatch dictionary
  mapping 9 `outputmode` keyword aliases to 3 integer codes:
  - `"non"`, `"none"`, `"off"` → 0 (no correlation output)
  - `"cor"`, `"corr"`, `"correlation"` → 1 (correlation functions)
  - `"raw"`, `"all"`, `"full"` → 2 (all output)
- Rewrote `check_output_mode()` to use `OUTPUT_MODE_TO_INT.get()`:
  - Default (`UNSET_STRING`) case handled first (sets `ioutputmode = 1`)
  - Other cases use dict lookup; `None` result triggers the existing error path
- Reduced the function body from a 5-branch if/elif chain (16 lines) to a clean
  2-branch structure (10 lines)

**Why**: Follows the same dict-dispatch pattern established in `hphi_writer.py`
(Steps 6-7). The `outputmode` mapping had 9 string aliases spread across 4 if/elif
branches with duplicated `print()` calls. The dict centralises the mapping and the
print is now a single line in the else branch.

**New tests**: 10 tests in `TestOutputModeToInt` class (`test/unit/test_common_writer.py`):
- 3 tests for "off" aliases (non, none, off → 0)
- 3 tests for "correlation" aliases (cor, corr, correlation → 1)
- 3 tests for "raw" aliases (raw, all, full → 2)
- 1 test for unknown-returns-None

**Test results**:
- Unit tests: 641 passed (631 → 641, +10 new dispatch dict tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_fcortho_smoke`

**Files modified**: `python/common_writer.py`, `test/unit/test_common_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Convert the `NExUpdatePath` mapping in `check_mod_para()`
to a `MODEL_GC_TO_EX_UPDATE_PATH` dict (composite-key `(ModelType, lGC) → int`,
similar to `MODEL_GC_TO_CALC_MODEL` in `hphi_writer.py`). Or extract the
solver-specific default-parameter blocks in `check_mod_para()` into separate
functions dispatched via `_SOLVER_DEFAULTS_DISPATCH` dict.

### 2026-01-27 — Phase 3, Step 12: Dict dispatch for `NExUpdatePath` in `common_writer.py`

**Target**: `common_writer.py` — `check_mod_para()` function, NExUpdatePath mapping

**What changed**:
- Added `MODEL_GC_TO_EX_UPDATE_PATH: dict[tuple, int]` module-level dispatch dictionary
  mapping `(ModelType, lGC)` composite keys to `NExUpdatePath` integer values:
  - `(HUBBARD, 0)` → 0, `(HUBBARD, 1)` → 0
  - `(SPIN, 0)` → 2, `(SPIN, 1)` → 2
  - `(KONDO, 0)` → 1, `(KONDO, 1)` → 3
- Rewrote the 9-line if/elif/nested-if chain (lines 783–791) to a 3-line dict lookup:
  ```python
  key = (StdI.model, StdI.lGC)
  ex_path = MODEL_GC_TO_EX_UPDATE_PATH.get(key)
  if ex_path is not None:
      StdI.NExUpdatePath = ex_path
  ```

**Why**: The NExUpdatePath mapping is a classic composite-key dispatch: the result
depends on both `model` and `lGC`. This is the same pattern used by
`MODEL_GC_TO_CALC_MODEL` in `hphi_writer.py` (Step 7). The dict makes the mapping
explicit and tabular rather than nested-conditional.

**New tests**: 8 tests in `TestModelGcToExUpdatePath` class (`test/unit/test_common_writer.py`):
- 2 tests for Hubbard (canonical + GC → 0)
- 2 tests for Spin (canonical + GC → 2)
- 2 tests for Kondo (canonical → 1, GC → 3)
- 1 test for dict completeness (6 entries)
- 1 test for unknown key returns None

**Test results**:
- Unit tests: 649 passed (641 → 649, +8 new dispatch dict tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spinaniso_honeycomb`
  - `fulldiag_hubbard_kagome_longrange`

**Files modified**: `python/common_writer.py`, `test/unit/test_common_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Extract the solver-specific default-parameter blocks in
`check_mod_para()` into separate helper functions dispatched via a solver-type dict.
The function currently has large `if solver == "mVMC"` / `elif solver == "HPhi"`
blocks that each set dozens of defaults — these could become `_check_mod_para_mvmc()`,
`_check_mod_para_hphi()`, etc., dispatched via a `_SOLVER_DEFAULTS` dict.

### 2026-01-27 — Phase 3, Step 13: Extract solver-specific blocks from `check_mod_para()` into dispatched helpers

**Target**: `common_writer.py` — `check_mod_para()` function

**What changed**:
- Extracted 4 solver-specific if/elif blocks into private helper functions:
  - `_check_mod_para_hphi(StdI)` — HPhi: Lanczos params, spectrum frequency grid
  - `_check_mod_para_mvmc(StdI)` — mVMC: VMC sampling, SR-optimisation, exchange-update
  - `_check_mod_para_uhf(StdI)` — UHF/H-wave: random seed, iteration, mixing, convergence
- UHF and H-wave share the identical handler (`_check_mod_para_uhf`)
- Added `_SOLVER_DEFAULTS_DISPATCH: dict[str, callable]` mapping `SolverType` → handler
- Rewrote the 4-branch if/elif chain in `check_mod_para()` to a 3-line dict dispatch:
  ```python
  handler = _SOLVER_DEFAULTS_DISPATCH.get(StdI.solver)
  if handler is not None:
      handler(StdI)
  ```
- The conserved-quantities section (ncond / 2Sz) remains inline in `check_mod_para()`
  since it cross-cuts solver and model types and doesn't benefit from extraction

**Why**: The 4-branch if/elif was ~100 lines of solver-specific defaults crammed into
a single function. Extracting each branch into its own function:
- Makes `check_mod_para()` shorter and clearer (dispatch + conserved quantities only)
- Makes individual solver-specific logic independently testable
- Eliminates the duplicate UHF/HWAVE code block (both now call the same function)
- Follows the dict-dispatch pattern established in earlier steps

**New tests**: 15 tests in `TestSolverDefaultsDispatch` class (`test/unit/test_common_writer.py`):
- 6 tests for dispatch dict structure (4 entries, correct function mappings, UHF=HWAVE, unknown→None)
- 3 tests for HPhi handler (Lanczos_max, exct, Omega defaults)
- 3 tests for mVMC handler (CParaFileHead, NVMCCalMode, RndSeed)
- 3 tests for UHF handler (RndSeed, Iteration_max, mix)

**Test results**:
- Unit tests: 664 passed (649 → 664, +15 new helper/dispatch tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spinaniso_honeycomb`
  - `fulldiag_hubbard_kagome_longrange`

**Files modified**: `python/common_writer.py`, `test/unit/test_common_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Apply dict dispatch to the conserved-quantities section of
`check_mod_para()`, or move to a different refactoring target. The conserved-quantities
block has complex nested `model × solver × lGC` logic that may benefit from a strategy
table. Alternatively, extract the `print_interactions()` function (currently ~400 lines)
into a dedicated `interaction_writer.py` module as suggested by Phase 1 of the
refactoring strategy.

### 2026-01-27 — Phase 3, Step 14: Extract `_check_conserved_quantities()` from `check_mod_para()`

**Target**: `common_writer.py` — `check_mod_para()` function, conserved-quantities section

**What changed**:
- Extracted the 36-line conserved-quantities block (ncond / 2Sz validation) from
  `check_mod_para()` into a new private function `_check_conserved_quantities(StdI)`.
- `check_mod_para()` is now just 4 lines: dispatch solver-specific defaults, then
  call `_check_conserved_quantities(StdI)`.
- Added full NumPy-style docstring to `_check_conserved_quantities()` documenting all
  11 model × solver × lGC combinations.
- Updated `check_mod_para()` docstring to reference the new helper.

**Why**: The conserved-quantities logic is a distinct concern from the solver-specific
defaults: it's cross-cutting (depends on model × solver × lGC) and applies to all
solvers. Extracting it makes the main function trivially readable and the validation
logic independently testable. A dict-dispatch approach was considered but rejected:
the branches contain different function calls (`required_val_i`, `not_used_i`,
`print_val_i`) and conditional assignments, not simple value mappings.

**New tests**: 11 tests in `TestCheckConservedQuantities` class (`test/unit/test_common_writer.py`):
- Hubbard: mVMC canonical (Sz2 defaults to 0, ncond preserved), mVMC GC (Sz2 unchanged),
  HPhi canonical (passes with ncond set), HPhi GC (passes)
- Spin: mVMC (ncond set to 0), HPhi canonical (passes with Sz2), HPhi GC (passes)
- Kondo: mVMC canonical (Sz2 defaults to 0), HPhi canonical (passes with ncond),
  HPhi GC (passes)

**Test results**:
- Unit tests: 675 passed (664 → 675, +11 new conserved-quantities tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spinaniso_honeycomb`
  - `fulldiag_hubbard_kagome_longrange`

**Files modified**: `python/common_writer.py`, `test/unit/test_common_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Extract the `print_interactions()` function (~400 lines) into
a dedicated `interaction_writer.py` module as suggested by Phase 1 of the refactoring
strategy. This is the last large function remaining in `common_writer.py` and would
complete the decomposition of the original monolithic code.

### 2026-01-27 — Phase 3, Step 15: Extract `print_interactions()` into `interaction_writer.py`

**Target**: `common_writer.py` → new `interaction_writer.py` module

**What changed**:
- Created `python/interaction_writer.py` containing the `print_interactions()` function
  (429 lines), extracted verbatim from `common_writer.py`.
- `common_writer.py` now re-exports `print_interactions` via:
  ```python
  from interaction_writer import print_interactions  # re-exported for backward compat
  ```
- Updated `solver_writer.py` to import `print_interactions` directly from
  `interaction_writer` instead of `common_writer`.
- Updated module docstring in `common_writer.py` to note the re-export.
- `common_writer.py` shrank from 1324 → 924 lines (−400 lines, −30%).

**Why**: `print_interactions()` was a self-contained 400-line function handling
interaction-term merging and `.def` file output. It has no dependency on any other
function in `common_writer.py` — it only needs `StdIntList`. Extracting it into its
own module:
- Completes the Phase 1 goal of decomposing `common_writer.py`
- Makes the interaction logic independently importable and testable
- Reduces `common_writer.py` to a manageable size focused on parameter validation
  and common output (locspn, trans, namelist, modpara, green functions)

**New files**: `python/interaction_writer.py`, `test/unit/test_interaction_writer.py`

**New tests**: 7 tests across 6 test classes in `test_interaction_writer.py`:
- `TestImportFromInteractionWriter` — verifies direct import
- `TestMergeDuplicateCoulombIntra` — duplicate site indices merged correctly
- `TestMergeDuplicateCoulombInter` — symmetric pair indices merged correctly
- `TestBoostSuppressesOutput` — `lBoost=1` sets flags to 0, no files written
- `TestExchangeWritten` — exchange.def created with correct header
- `TestInterAllWritten` — interall.def created for single term
- `TestReExportFromCommonWriter` — backward-compatible re-export verified

**Test results**:
- Unit tests: 682 passed (675 → 682, +7 new tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spinaniso_honeycomb`
  - `fulldiag_hubbard_kagome_longrange`

**Files created**: `python/interaction_writer.py`, `test/unit/test_interaction_writer.py`
**Files modified**: `python/common_writer.py`, `python/solver_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: The `print_interactions()` function has 6 nearly-identical
code blocks (CoulombIntra, CoulombInter, Hund, Exchange, PairLift, PairHopp) that each
perform: merge duplicates → count non-zero → set flag → write file. These could be
consolidated into a generic helper function parameterized by the interaction type
metadata (array name, index array, flag name, filename, header strings). This would
reduce ~200 lines of repetitive code to ~40 lines + a data table.

### 2026-01-27 — Phase 3, Step 16: Consolidate 6 interaction blocks into data-driven loop

**Target**: `interaction_writer.py` — the 6 repetitive interaction blocks in
`print_interactions()`

**What changed**:
- Added 4 private helper functions:
  - `_merge_1idx(nterms, indx, coeff)` — merge duplicates for 1-index interactions
  - `_merge_2idx(nterms, indx, coeff)` — merge duplicates for 2-index interactions
    (symmetric pair matching)
  - `_count_nonzero(nterms, coeff)` — count terms with |coeff| > 1e-6
  - `_process_interaction(StdI, ...)` — generic driver: merge → count → set flag → write
- Added 2 private file-writing functions:
  - `_write_1idx_file(filename, count_label, banner, nterms, indx, coeff)`
  - `_write_2idx_file(filename, count_label, banner, nterms, indx, coeff)`
- Added `_INTERACTION_TYPES` metadata table: a list of 6 dicts, one per interaction type,
  containing attribute names, filenames, header strings, and index count
- Replaced the 6 copy-pasted blocks (~200 lines) with a 2-line data-driven loop:
  ```python
  for spec in _INTERACTION_TYPES:
      _process_interaction(StdI, **spec)
  ```

**Why**: The 6 blocks (CoulombIntra, CoulombInter, Hund, Exchange, PairLift, PairHopp)
were structurally identical — each performed the same merge → count → flag → write
pipeline with only attribute names and header strings varying. The data-driven approach:
- Eliminates all code duplication (DRY)
- Makes the per-type metadata explicit and declarative
- Makes adding a new interaction type a 1-dict addition (no code changes)
- Makes merge/count/write logic independently testable

**New tests**: 17 tests across 5 classes in `test/unit/test_interaction_writer.py`:
- `TestMerge1Idx` (4 tests): no duplicates, duplicate merged, triple, empty
- `TestMerge2Idx` (4 tests): no duplicates, same order, reversed order, empty
- `TestCountNonzero` (5 tests): all nonzero, some zero, all below threshold,
  at threshold, empty
- `TestInteractionTypesTable` (4 tests): 6 entries, CoulombIntra is 1-idx,
  remaining are 2-idx, all have required keys

**Test results**:
- Unit tests: 699 passed (682 → 699, +17 new tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spinaniso_honeycomb`
  - `fulldiag_hubbard_kagome_longrange`

**Files modified**: `python/interaction_writer.py`, `test/unit/test_interaction_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: The InterAll section in `print_interactions()` (~160 lines
of complex merge/reorder logic) could potentially benefit from having its 3 passes
extracted into named helper functions (e.g. `_merge_interall_equivalent`,
`_reorder_interall_hermitian`, `_remove_interall_diagonal`). Alternatively, review
the overall Phase 3 progress and consider whether to advance to Phase 2 items
(introducing classes like `LatticeBuilder` or `SolverWriter` subclass hierarchy).

### 2026-01-27 — Phase 3, Step 17: Extract InterAll passes into named helpers

**Target**: `interaction_writer.py` — the ~160-line inline InterAll section in
`print_interactions()`

**What changed**:
- Extracted 4 private helper functions from the InterAll section:
  - `_merge_interall_equivalent(nintr, intrindx, intr)` — Pass 1: merge equivalent
    terms (exact match → add, Hermitian conjugate → add, partial swap → subtract)
  - `_reorder_interall_hermitian(nintr, intrindx, intr)` — Pass 2: force Hermitian
    ordering on index pairs (direct conjugate → reorder, two sub-cases → reorder
    with sign flip)
  - `_remove_interall_diagonal(nintr, intrindx, intr)` — Pass 3: zero spurious
    diagonal terms (diagonal pair without matching diagonal)
  - `_write_interall(StdI)` — count non-zero, set `Lintr` flag, write `interall.def`
- Replaced ~160 lines of inline InterAll code with 4 function calls:
  ```python
  _merge_interall_equivalent(StdI.nintr, StdI.intrindx, StdI.intr)
  _reorder_interall_hermitian(StdI.nintr, StdI.intrindx, StdI.intr)
  _remove_interall_diagonal(StdI.nintr, StdI.intrindx, StdI.intr)
  _write_interall(StdI)
  ```

**Why**: The InterAll section was a long block of inline code with 3 distinct
algorithmic passes that were difficult to understand as a monolith. Extracting
into named functions:
- Makes each pass independently testable
- Provides clear documentation of each pass's purpose via docstrings
- Reduces `print_interactions()` to a clean, readable sequence of function calls
- Each helper is self-contained with clear parameters

**New tests**: 20 tests across 4 classes in `test/unit/test_interaction_writer.py`:
- `TestMergeInterallEquivalent` (6 tests): empty, no match, Case A exact match,
  Case B Hermitian conjugate add, Case C partial swap subtract, three terms
- `TestReorderInterallHermitian` (4 tests): empty, no reorder needed, Pattern 1
  direct Hermitian, Pattern 2 sign flip
- `TestRemoveInterallDiagonal` (6 tests): empty, no diagonal pair kept, diagonal
  pair (0,4) removed, diagonal pair (2,6) removed, matching diagonal kept, mixed
- `TestWriteInterall` (4 tests): zero terms no file, boost suppresses, nonzero
  writes file, all zero no file

**Test results**:
- Unit tests: 719 passed (699 → 719, +20 new tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_square_longrange`

**Files modified**: `python/interaction_writer.py`, `test/unit/test_interaction_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Consider consolidating the common `tempfile.TemporaryDirectory`
+ `os.chdir` pattern in `test_interaction_writer.py` into a pytest fixture, or
review the overall Phase 3 progress and plan Phase 4 items. The `interaction_writer.py`
module is now well-decomposed: 6 standard interactions handled by data-driven loop,
InterAll handled by 4 named helpers. Potential next targets include:
- Extracting `print_trans()` from `common_writer.py` into its own module (~200 lines)
- Consolidating the lattice file writer pattern (all 8 lattice modules share similar
  structure)

### 2026-01-27 — Phase 3, Step 18: Extract `print_mod_para` body writers + dict dispatch

**Target**: `common_writer.py` — the 4-branch if/elif chain in `print_mod_para()`

**What changed**:
- Extracted 3 private body-writer functions:
  - `_write_modpara_hphi(fp, StdI)` — HPhi-specific fields (Lanczos, spectrum, etc.)
  - `_write_modpara_mvmc(fp, StdI)` — mVMC-specific fields (VMC sampling, SR opt, etc.)
  - `_write_modpara_uhf_hwave(fp, StdI)` — shared UHF/H-wave fields (iteration, mixing,
    convergence); selects the correct banner from `_MODPARA_BANNER`
- Added `_MODPARA_BANNER` dict mapping UHF/HWAVE solver types to their header strings
- Added `_MODPARA_BODY_DISPATCH` dict mapping all 4 solver types to their body-writer
  function (UHF and HWAVE share the same writer)
- Rewrote `print_mod_para()` from a 4-branch if/elif chain (~140 lines) to a 3-line
  dict dispatch:
  ```python
  writer = _MODPARA_BODY_DISPATCH.get(StdI.solver)
  if writer is not None:
      writer(fp, StdI)
  ```

**Why**: The UHF and HWAVE blocks were 100% identical except for the header banner
string — clear code duplication. The 4-branch if/elif pattern is the same dispatch
pattern already eliminated elsewhere (check_mod_para, NExUpdatePath). Benefits:
- Eliminates the UHF/HWAVE duplication entirely
- Makes `print_mod_para()` a clean 6-line dispatcher (open file + dispatch)
- Each solver's body is independently testable via `io.StringIO`
- Adding a new solver requires only a new function + one dict entry

**New tests**: 11 tests in `TestModparaBodyDispatch` in `test/unit/test_common_writer.py`:
- Dispatch table: has 4 solvers, maps all solver types, UHF/HWAVE share writer
- Banner dict: has UHF and HWAVE, correct values
- Body writers: HPhi content, mVMC content, UHF banner, HWAVE banner
- Conditional logic: ExpandCoef omitted for non-TE, included for TE

**Test results**:
- Unit tests: 730 passed (719 → 730, +11 new tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_square_longrange`

**Files modified**: `python/common_writer.py`, `test/unit/test_common_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Apply the same extract-and-dispatch pattern to
`print_namelist()`, which also has a 2-branch solver-specific section (HPhi vs mVMC).
Alternatively, extract `print_trans()` into its own module, or begin consolidating
the lattice modules which share similar structure.

---

### 2026-01-27 — Phase 3, Step 19: Package restructuring — `lattice/` and `writer/` subpackages

**Motivation**: The flat `python/` directory had 30 `.py` files with no organizational
structure. Related modules (lattice construction, output writers) were mixed together,
making navigation difficult. Reorganizing into subpackages makes the dependency graph
explicit and improves code discoverability.

**What changed**:

1. Created `python/lattice/` subpackage with `__init__.py` (re-exports all lattice modules)
2. Created `python/writer/` subpackage with `__init__.py` (re-exports `get_solver_writer`, `print_interactions`)
3. Moved 14 files into `lattice/`:
   - `chain_lattice.py`, `square_lattice.py`, `ladder.py`, `triangular_lattice.py`
   - `honeycomb_lattice.py`, `kagome.py`, `orthorhombic.py`, `fc_ortho.py`, `pyrochlore.py`
   - `wannier90.py`
   - `site_util.py`, `geometry_output.py`, `interaction_builder.py`, `input_params.py`
4. Moved 7 files into `writer/`:
   - `common_writer.py`, `interaction_writer.py`, `hphi_writer.py`
   - `mvmc_writer.py`, `mvmc_variational.py`, `solver_writer.py`, `export_wannier90.py`
5. 7 files remain at top level:
   - `__main__.py`, `stdface_main.py`, `version.py`, `stdface_vals.py`
   - `param_check.py`, `keyword_parser.py`, `stdface_model_util.py`

**Import strategy**:
- Within `lattice/`: relative imports for siblings (`from .site_util import ...`)
- Within `writer/`: relative imports for siblings (`from .common_writer import ...`)
- Cross-package: absolute (`from lattice.site_util import ...` in `writer/mvmc_variational.py`)
- Top-level files: absolute (`from lattice import chain_lattice`, `from writer.hphi_writer import ...`)
- `stdface_model_util.py`: updated re-export paths (`from lattice.input_params import ...`)

**Test updates**: 13 of 19 test files updated with new import paths. 6 test files
unchanged (import only top-level modules: `stdface_vals`, `param_check`, `stdface_main`,
`stdface_model_util`).

**Test results**:
- Unit tests: 730 passed (no change in count — import paths only)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_square_longrange`

**Files modified**: 14 lattice files (moved + relative imports), 7 writer files
(moved + relative imports), `stdface_main.py`, `stdface_model_util.py`,
13 test files (import paths), 2 new `__init__.py` files.

**Final directory structure**:
```
python/
├── __main__.py
├── stdface_main.py
├── version.py
├── stdface_vals.py
├── param_check.py
├── keyword_parser.py
├── stdface_model_util.py
├── lattice/
│   ├── __init__.py
│   ├── site_util.py
│   ├── geometry_output.py
│   ├── interaction_builder.py
│   ├── input_params.py
│   ├── chain_lattice.py
│   ├── square_lattice.py
│   ├── ladder.py
│   ├── triangular_lattice.py
│   ├── honeycomb_lattice.py
│   ├── kagome.py
│   ├── orthorhombic.py
│   ├── fc_ortho.py
│   ├── pyrochlore.py
│   └── wannier90.py
└── writer/
    ├── __init__.py
    ├── common_writer.py
    ├── interaction_writer.py
    ├── hphi_writer.py
    ├── mvmc_writer.py
    ├── mvmc_variational.py
    ├── solver_writer.py
    └── export_wannier90.py
```

**Suggested next step**: Continue Phase 3 idiom improvements (e.g., extract-and-dispatch
for `print_namelist()`) or begin Phase 2 class extraction (SolverWriter hierarchy).

---

### 2026-01-27 — Phase 3, Step 20: Extract `print_namelist()` body writers into dict dispatch

**Motivation**: `print_namelist()` had a 2-branch `if/elif` chain (HPhi vs mVMC),
with UHF/HWAVE having no solver-specific entries. The 7-flag interaction check was
also a repetitive `if` ladder. Both patterns are better expressed as data-driven dispatch.

**What changed** (in `python/writer/common_writer.py`):

1. Extracted `_write_namelist_hphi(fp, StdI)` — writes CalcMod, excitation type,
   time-evolution pump entries, SpectrumVec, and Boost entries.
2. Extracted `_write_namelist_mvmc(fp, StdI)` — writes Gutzwiller, Jastrow,
   Orbital, OrbitalParallel, and TransSym entries.
3. Added `_NAMELIST_BODY_DISPATCH` dict mapping `SolverType` → body-writer function.
   Only HPhi and mVMC have entries; UHF/HWAVE have no solver-specific namelist lines.
4. Added `_INTERACTION_FLAGS` list of `(attr_name, line)` tuples to replace the 7
   sequential `if StdI.LXxx == 1:` checks with a data-driven loop.
5. Rewrote `print_namelist()` to: write common prefix, loop over `_INTERACTION_FLAGS`,
   write Green-function entries, then dispatch to solver-specific body writer.

**Design**: Same pattern as Step 18 (`print_mod_para` dispatch). Each solver's
body is independently testable via `io.StringIO`. The `_INTERACTION_FLAGS` table
makes adding new interaction types a single-line table addition.

**New tests**: 14 tests in `TestNamelistBodyDispatch` in `test/unit/test_common_writer.py`:
- Dispatch table: has HPhi/mVMC, UHF/HWAVE absent, maps to correct functions
- Interaction flags table: 7 entries, valid attribute names
- HPhi body: CalcMod, SingleExcitation, PairExcitation, TEOneBody, TETwoBody, Boost
- mVMC body: basic entries, OrbitalParallel with lGC=1, OrbitalParallel with nonzero Sz2
- UHF: no solver-specific entries in output
- Interaction flags: conditional gating of CoulombInter, Exchange, InterAll

**Test results**:
- Unit tests: 744 passed (730 → 744, +14 new tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_square_longrange`

**Files modified**: `python/writer/common_writer.py`, `test/unit/test_common_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Apply the same extract-and-dispatch pattern to
`check_mod_para()` solver-defaults section, or extract `print_trans()` into its
own module, or begin Phase 2 class extraction (SolverWriter hierarchy).

---

### 2026-01-27 — Phase 3, Step 21: Extract `_spin_max()` and `_skip_local_spin_pair()` helpers

**Motivation**: The pattern `if locspinflag[site] == 0: SMax = 1; else: SMax = locspinflag[site]`
appeared 8 times across `print_1_green()` and `print_2_green()`. The local-spin pair
filtering pattern `if site_a != site_b and locspinflag[site_a] != 0 and locspinflag[site_b] != 0: continue`
appeared 4 times. Both are error-prone to repeat and obscure the core loop logic.

**What changed** (in `python/writer/common_writer.py`):

1. Extracted `_spin_max(locspinflag, site)` → returns `1` for itinerant sites
   (flag=0), or the flag value for local-spin sites. Replaces 8 inline if/else blocks.
2. Extracted `_skip_local_spin_pair(locspinflag, site_a, site_b)` → returns `True`
   when two distinct sites are both local-spin. Replaces 4 inline filter conditions.
3. Simplified `print_1_green()`: 102 → 68 lines. The 4-deep nested loops now read
   as `for ispin in range(_spin_max(...) + 1)` with `if _skip_local_spin_pair(...)`.
4. Simplified `print_2_green()`: 142 → 80 lines. Same helper usage across 8 call sites.

**Design**: Both helpers are pure functions (no side effects, no `StdI` dependency)
that take only the `locspinflag` array and site indices. This makes them independently
testable and reusable in any future Green-function code.

**New tests**: 9 tests in `TestSpinMax` and `TestSkipLocalSpinPair` in
`test/unit/test_common_writer.py`:
- `_spin_max`: itinerant (returns 1), S=1/2 (returns 1), S=1 (returns 2), S=3/2 (returns 3)
- `_skip_local_spin_pair`: both itinerant (False), same site (False), mixed (False),
  both local different sites (True), both local same site (False)

**Test results**:
- Unit tests: 753 passed (744 → 753, +9 new tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_square_longrange`

**Files modified**: `python/writer/common_writer.py`, `test/unit/test_common_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Split `print_2_green()` into mode-1 (correlation) and mode-2
(raw) index generators for further clarity, or begin Phase 2 class extraction
(SolverWriter hierarchy).

---

### 2026-01-27 — Phase 3, Step 22: Extract Green-function index generators from file-I/O functions

**Motivation**: `print_1_green()` and `print_2_green()` each mixed index-generation
algorithms with file I/O in a single function. This made the core logic hard to test
without temp directories, and the Kondo site-mapping appeared duplicated in both.
Separating index generation from file writing enables direct unit testing of the
algorithms as pure functions.

**What changed** (in `python/writer/common_writer.py`):

1. Extracted `_kondo_site(isite, NsiteUC, nsite)` — maps a Kondo unit-cell index
   to the physical site index. Used by both one-body and two-body correlation generators.
2. Extracted `_green1_indices_corr(NsiteUC, nsite, locspinflag, is_kondo)` — generates
   one-body correlation-mode index tuples (same-spin pairs only).
3. Extracted `_green1_indices_raw(nsite, locspinflag)` — generates one-body raw-mode
   index tuples (all combinations with local-spin filtering).
4. Extracted `_green2_indices_corr(NsiteUC, nsite, locspinflag, is_kondo, is_mvmc)` —
   generates two-body correlation-mode index tuples with spin conservation and
   mVMC alternative ordering.
5. Extracted `_green2_indices_raw(nsite, locspinflag)` — generates two-body raw-mode
   index tuples with local-spin pair filtering.
6. Simplified `print_1_green()` to a 2-branch dispatcher + file writer (~25 lines).
7. Simplified `print_2_green()` to a 2-branch dispatcher + file writer (~25 lines).

**Design**: All five generators are pure functions (no file I/O, no `StdI` dependency)
that take primitive parameters and return lists of tuples. This makes them independently
testable, composable, and easy to reason about. The `_kondo_site` helper eliminates
the duplicated if/else Kondo mapping from both `print_1_green` and `print_2_green`.

**New tests**: 15 tests in 5 test classes in `test/unit/test_common_writer.py`:
- `TestKondoSite` (3): within UC, beyond UC, single-site UC
- `TestGreen1IndicesCorr` (4): basic 2-site, same-spin only, Kondo doubles, local-spin skip
- `TestGreen1IndicesRaw` (3): 2-site count, 1-site all combos, local-spin skip
- `TestGreen2IndicesCorr` (3): spin conservation, mVMC reorder, Kondo doubles
- `TestGreen2IndicesRaw` (2): 1-site count, local-spin filtering

**Test results**:
- Unit tests: 768 passed (753 → 768, +15 new tests)
- Integration check: 5/5 pass:
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`
  - `fulldiag_hubbard_square_longrange`

**Files modified**: `python/writer/common_writer.py`, `test/unit/test_common_writer.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Begin Phase 2 class extraction (SolverWriter hierarchy),
or continue Phase 3 with `_check_conserved_quantities()` dict-dispatch refactoring.

---

## Step 23 — Data-driven solver-specific field resets (stdface_main.py)

**Date**: 2026-01-27
**Phase**: 3 — Idiomatic Python (data-driven tables)
**Scope**: `python/stdface_main.py`, `test/unit/test_reset_vals_dispatch.py`

### Motivation

The four solver-specific reset functions (`_reset_hphi_fields`,
`_reset_mvmc_fields`, `_reset_uhf_fields`, `_reset_hwave_fields`) each
performed a flat list of `StdI.field = sentinel` assignments, totalling ~72
scalar/array assignments across ~106 lines of code (including docstrings).
The UHF and HWAVE functions were near-identical — HWAVE differed only by
adding `export_all` and `lattice_gp`.  The `_SOLVER_RESET_DISPATCH` dict
mapped solver types to these functions.

### Changes

1. **Replaced 4 functions + 1 dispatch dict with 2 data tables + 1 generic function**:
   - `_SOLVER_RESET_SCALARS: dict[SolverType, list[tuple[str, object]]]` — scalar
     field resets (applied via `setattr`)
   - `_SOLVER_RESET_ARRAYS: dict[SolverType, list[tuple[str, object]]]` — array-fill
     field resets (applied via `getattr(StdI, name)[...] = value`)
   - `_apply_field_resets(StdI, solver)` — single 8-line function that iterates
     both tables

2. **Eliminated UHF/HWAVE duplication** with shared base lists:
   - `_UHF_BASE_SCALARS` — 9 shared scalar entries
   - `_UHF_BASE_ARRAYS` — 1 shared array entry (`boxsub`)
   - UHF uses the base lists directly; HWAVE extends scalars with `+ [("export_all", ...), ("lattice_gp", ...)]`

3. **Updated `_reset_vals()`** call site: replaced 3-line dispatch lookup
   (`_SOLVER_RESET_DISPATCH.get(...)` / `if ... is not None` / call) with
   a single call to `_apply_field_resets(StdI, StdI.solver)`.

4. **Rewrote tests** in `test/unit/test_reset_vals_dispatch.py`:
   - Old: 23 tests across 5 classes testing 4 individual functions + dispatch dict
   - New: 30 tests across 7 classes testing data tables, base sharing, and
     `_apply_field_resets` behaviour for all 4 solvers + unknown solver

### Design rationale

- **Data tables over functions**: Each reset function was pure boilerplate —
  a sequence of `field = value` with no branching or computation.  Data tables
  make the field lists scannable, eliminate function-call overhead, and make
  the UHF/HWAVE relationship explicit.
- **Separate scalar/array tables**: Array-fill assignments (`arr[...] = val`)
  require `getattr` + indexing, unlike scalar assignments which use `setattr`.
  Keeping them separate avoids needing sentinel markers or type-dispatch in
  the loop body.
- **`arr[...]` (Ellipsis)**: Covers both 1D (`SpectrumQ`, `VecPot`) and 2D
  (`boxsub`) arrays uniformly without needing to know the array shape.

### Test results

- **Baseline**: 768 passed
- **After**: 775 passed (768 − 23 old + 30 new = 775, net +7)
- **Integration**: 5/5 pass
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`

**Files modified**: `python/stdface_main.py`, `test/unit/test_reset_vals_dispatch.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Continue Phase 3 with `_check_conserved_quantities()`
dict-dispatch refactoring, or begin Phase 2 class extraction (SolverWriter hierarchy).

---

## Step 24 — Dict-dispatch solver keyword parsers (keyword_parser.py)

**Date**: 2026-01-27
**Phase**: 3 — Idiomatic Python (dict-dispatch tables)
**Scope**: `python/keyword_parser.py`, `test/unit/test_keyword_parser.py`

### Motivation

Four solver-specific keyword parser functions (`parse_hphi_keyword`,
`parse_mvmc_keyword`, `parse_uhf_keyword`, `parse_hwave_keyword`) each
consisted of long if/elif chains: 39, 32, 18, and 21 branches respectively
(~246 lines total).  Every branch followed one of two mechanical patterns:

1. **Scalar**: `StdI.field = store_func(keyword, value, StdI.field)`
2. **Array element**: `StdI.arr[idx] = store_func(keyword, value, cast(StdI.arr[idx]))`

The UHF and HWAVE parsers were near-identical (HWAVE added 4 extra keywords),
and the 9 `boxsub` keywords appeared identically in mVMC, UHF, and HWAVE.

### Changes

1. **Created keyword dispatch tables**:
   - `_BOXSUB_KEYWORDS` — 9 shared boxsub entries (mVMC, UHF, HWAVE)
   - `_UHF_BASE_KEYWORDS` — 18 shared UHF/HWAVE entries (includes boxsub)
   - `_HPHI_KEYWORDS` — 39 entries (includes 6 array-element entries for
     SpectrumQ and VecPot)
   - `_MVMC_KEYWORDS` — 32 entries (includes boxsub via spread)
   - `_UHF_KEYWORDS` — alias for `_UHF_BASE_KEYWORDS`
   - `_HWAVE_KEYWORDS` — extends `_UHF_BASE_KEYWORDS` with 4 extra entries
   - `_SOLVER_KEYWORD_TABLES` — maps `SolverType` to keyword table

2. **Two entry formats**:
   - Scalar: `(store_func, "field_name")` — 2-tuple
   - Array element: `(store_func, "array_name", index, cast)` — 4-tuple

3. **Generic dispatch function** `_apply_keyword_table(table, keyword, value, StdI)`:
   - Looks up keyword in table (O(1) dict lookup vs O(n) if/elif scan)
   - Applies scalar or array-element assignment based on tuple length

4. **Rewrote 4 parser functions** as thin wrappers calling `_apply_keyword_table`

5. **Converted `parse_solver_keyword`** from 4-way if/elif to
   `_SOLVER_KEYWORD_TABLES.get(solver)` dict lookup

### Design rationale

- **Dict dispatch**: Replaces O(n) linear if/elif scanning with O(1) dict
  lookup.  Each parser function is now a 1-line delegation.
- **Shared base tables**: `_BOXSUB_KEYWORDS` and `_UHF_BASE_KEYWORDS` use
  `**` spread to eliminate the triplication of boxsub entries and the
  duplication of UHF/HWAVE entries.
- **Two-format tuples**: The 2-tuple (scalar) vs 4-tuple (array element)
  distinction keeps the table entries compact while supporting both assignment
  patterns without sentinel markers.
- **Preserved public API**: All four `parse_*_keyword` functions retain their
  original signatures and are still importable by name.

### Test results

- **Baseline**: 775 passed
- **After**: 794 passed (775 + 19 new)
- **New tests**:
  - `TestKeywordTableStructure` (5 tests) — table shapes, types, lowercase keys
  - `TestUHFHWAVEKeywordSharing` (6 tests) — base sharing, extra keys, boxsub
  - `TestApplyKeywordTable` (8 tests) — scalar int/float/string, array elements,
    unknown keyword, empty table
- **Integration**: 5/5 pass
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`

**Files modified**: `python/keyword_parser.py`, `test/unit/test_keyword_parser.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Apply the same dict-dispatch pattern to
`parse_common_keyword()` (231-branch if/elif chain, ~498 lines), or continue
Phase 3 with `_check_conserved_quantities()` refactoring.

---

## Step 25 — `parse_common_keyword()` dict-dispatch with generator functions

**Date**: 2026-01-27
**Phase**: 3 — Pythonic control-flow consolidation
**Scope**: `python/keyword_parser.py`, `test/unit/test_keyword_parser.py`

### What changed

Replaced the 231-branch if/elif chain in `parse_common_keyword()` (~498 lines)
with a single `_COMMON_KEYWORDS` dict dispatch table (~110 lines).  Introduced
three generator functions to programmatically build repetitive table entries:

1. **`_j_matrix_keywords(prefix, scalar_field, matrix_field)`** — generates 10
   entries (1 isotropic scalar + 9 anisotropic matrix components) for each of
   12 J-exchange families: J, J0, J0', J0'', J1, J1', J1'', J2, J2', J2'',
   J', J''.  12 families × 10 entries = 120 entries from 12 one-line calls.

2. **`_cutoff_vec_keywords(prefix, vec_field)`** — generates 9 entries for each
   of 3 cutoff-vector families (cutoff_t, cutoff_u, cutoff_j): 3×3 grid of
   `a{0,1,2}{w,l,h}` suffixes.  3 families × 9 entries = 27 entries from
   3 one-line calls.

3. **`_box_keywords()`** — generates 9 entries for the supercell box array
   (`a0w`, `a0l`, … `a2h`).

4. **`_COMMON_KEYWORDS` dict** — built using `**` spread from the three
   generators plus manual entries for remaining keywords (model, lattice,
   hopping t, Coulomb V, phase, etc.).

5. **`parse_common_keyword()`** — reduced to a single-line delegation:
   `return _apply_keyword_table(_COMMON_KEYWORDS, keyword, value, StdI)`.

### Lines changed

- **Removed**: ~498 lines (231-branch if/elif chain in `parse_common_keyword`)
- **Added**: ~110 lines (`_COMMON_KEYWORDS` dict + 3 generator functions)
- **Net**: ~388 lines removed from keyword_parser.py

### Design rationale

- **Generator functions**: The 12 J-families each have identical structure
  (1 scalar + 9 matrix components), making them ideal for generation.  This
  eliminates 120 near-identical if/elif branches with 12 one-line calls.
- **Reuses `_apply_keyword_table`**: The same generic dispatch function from
  Step 24 handles both common and solver-specific keywords.
- **Consistent tuple format**: All entries use the same 2-tuple (scalar) or
  4-tuple (array element) format established in Step 24.
- **`nelec` alias**: The C code mapped `nelec` to the `ncond` field; this is
  preserved as `"nelec": (store_with_check_dup_i, "ncond")`.

### Test results

- **Baseline**: 794 passed
- **After**: 820 passed (794 + 26 new)
- **New tests**:
  - `TestCommonKeywordsTable` (11 tests) — table type, lowercase keys, tuple
    shapes, presence of all J families, all hopping/Coulomb keywords, phases,
    box entries, cutoffs, nelec alias
  - `TestJMatrixKeywords` (6 tests) — entry count, scalar entry, 9 components,
    index mapping, 4-tuple format, primed prefixes
  - `TestCutoffVecKeywords` (4 tests) — entry count, all keys, index mapping,
    4-tuple format
  - `TestBoxKeywords` (5 tests) — entry count, all keys, index mapping, int
    cast, target field
- **Integration**: 5/5 pass
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`

**Files modified**: `python/keyword_parser.py`, `test/unit/test_keyword_parser.py`,
`python/history/refactoring_log.md`

**Suggested next step**: Refactor `_check_conserved_quantities()` or
`_check_lattice_params()` in `stdface_main.py`, or begin Phase 4 with
the lattice/writer subpackage reorganization.

---

## Step 26 — `_check_conserved_quantities()` data-driven rules table

**Date**: 2026-01-27
**Phase**: 3 — Pythonic control-flow consolidation
**Scope**: `python/writer/common_writer.py`, `test/unit/test_common_writer.py`

### What changed

Replaced the 3-level nested if/elif chain in `_check_conserved_quantities()`
(branching on model × solver × lGC, ~37 lines of logic) with a declarative
rules table `_CONSERVED_QTY_RULES`.

1. **`_CONSERVED_QTY_RULES` dict** — 12 entries keyed by
   `(ModelType, is_hphi: bool, lGC: int)`, each mapping to a
   `(ncond_label, ncond_action, sz2_action)` 3-tuple.  Actions are string
   codes: `"required"`, `"not_used"`, `"default_0"`, or `None` (skip).

2. **`_check_conserved_quantities()`** — rewrote as a flat lookup + dispatch
   of the two actions.  The Spin+mVMC special case (`ncond = 0`) is handled
   as a post-validation assignment, preserving the original C ordering where
   `not_used_i` is called before the assignment.

### Design rationale

- **Declarative table**: All 12 combinations of (model, solver-class, GC-flag)
  are visible in one place, making it easy to audit and modify rules.
- **Flat dispatch**: Eliminates 3 levels of nesting, reducing cyclomatic
  complexity from 11 to 3.
- **Action codes**: String-based actions keep the table compact while
  supporting the three distinct validation operations.  `None` means no
  action, avoiding empty branches.
- **Preserved semantics**: The `ncond_label` field captures the name
  variation ("nelec" vs "ncond") that the original code used in different
  branches.

### Test results

- **Baseline**: 820 passed
- **After**: 830 passed (820 + 10 new)
- **New tests**:
  - `TestConservedQtyRulesTable` (10 tests) — entry count, all models present,
    both is_hphi flags per model, both lGC values per model, 3-tuple structure,
    label is string, valid action codes, Hubbard+HPhi uses "nelec",
    Hubbard+non-HPhi defaults Sz2, Spin rules identical for is_hphi flag
- **Integration**: 5/5 pass
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`

**Files modified**: `python/writer/common_writer.py`,
`test/unit/test_common_writer.py`, `python/history/refactoring_log.md`

**Suggested next step**: Refactor `_check_lattice_params()` in
`stdface_main.py`, vectorize `input_spin_nn()` conflict detection in
`lattice/input_params.py`, or begin Phase 4 subpackage reorganization.

---

## Step 27 — Vectorize conflict detection in `input_spin_nn()` / `input_spin()`

**Date**: 2026-01-27
**Phase**: 3 — Pythonic control-flow consolidation
**Scope**: `python/lattice/input_params.py`, `test/unit/test_input_params.py`

### What changed

Replaced the O(81) quadruple-nested loop and O(9) double-nested loops in
`input_spin_nn()` with vectorised NumPy conflict detection, and extracted
5 reusable helper functions:

1. **`_has_set_elements(mat)`** — `bool(np.any(~np.isnan(mat)))`.  Single
   vectorised operation replaces iterating over all 9 elements.

2. **`_first_set_index(mat)`** — finds `(row, col)` of the first non-NaN
   element using `np.argmax(~np.isnan(mat))`.  Used by error messages to
   report which specific element triggered the conflict.

3. **`_check_scalar_vs_matrix(scalar, mat, scalar_name, mat_name)`** —
   replaces 4 separate scalar-vs-matrix conflict checks in the original
   element loop (lines 81–92).  Each check is now a single function call.

4. **`_check_matrix_vs_matrix(mat_a, mat_b, name_a, name_b)`** — replaces
   the O(81) quadruple-nested loop (lines 94–101) with two `_first_set_index`
   calls.  If both matrices have any set element, it's a conflict.

5. **`_resolve_spin_matrix(J0, J0name, *, J, J0All, JAll)`** — extracted
   the value-resolution cascade (lines 104–118) into a standalone function
   with keyword-only arguments.  Reused by both `input_spin_nn()` and
   `input_spin()`, eliminating the duplicated resolution logic.

### Lines changed

- **Removed**: ~52 lines (nested loops in `input_spin_nn` + `input_spin`)
- **Added**: ~110 lines (5 helper functions with docstrings)
- **Net**: +58 lines (helpers are more verbose due to docstrings, but
  eliminate all nested loops and duplicated resolution logic)

### Performance improvement

- `_check_matrix_vs_matrix`: O(1) NumPy vectorised vs O(81) Python loop
- `_check_scalar_vs_matrix`: O(1) vs O(9) per call
- `_has_set_elements`: O(1) vectorised vs O(9) Python loop

### Design rationale

- **Vectorised checks**: NumPy's `np.any(~np.isnan(...))` replaces Python
  loops, leveraging C-level array operations.
- **Reusable helpers**: `_check_scalar_vs_matrix` is called 4 times in
  `input_spin_nn` and once in `input_spin`, replacing 5 separate loop blocks.
- **`_resolve_spin_matrix`**: Keyword-only `J`, `J0All`, `JAll` arguments make
  the cascade explicit.  `input_spin` reuses it with only `J0All`, avoiding
  the duplicated resolution code from the original.
- **Error messages**: Slightly more informative than the original (both
  conflicting names are always shown), but identical abort behavior.

### Test results

- **Baseline**: 830 passed
- **After**: 850 passed (830 + 20 new)
- **New tests**:
  - `TestHasSetElements` (4 tests) — all-NaN, one set, all set, zero counts
  - `TestFirstSetIndex` (4 tests) — all-NaN, single, first-of-multiple, corner
  - `TestCheckScalarVsMatrix` (3 tests) — NaN scalar, NaN matrix, conflict
  - `TestCheckMatrixVsMatrix` (3 tests) — both NaN, one NaN, both set
  - `TestResolveSpinMatrix` (6 tests) — zeros default, J0All diagonal, JAll
    diagonal, priority, J fallback, explicit preserved
- **Integration**: 5/5 pass
  - `fulldiag_hubbard_chain_longrange`
  - `fulldiag_hubbard_square_longrange`
  - `fulldiag_kondo_chain_longrange`
  - `fulldiag_spingcboost_chain_success`
  - `fulldiag_wannier90_alphaBETS`

**Files modified**: `python/lattice/input_params.py`,
`test/unit/test_input_params.py`, `python/history/refactoring_log.md`

**Suggested next step**: Refactor `_check_lattice_params()` in
`stdface_main.py`, or begin Phase 4 subpackage reorganization.

---

## Step 28 — Data-driven `_reset_vals()` common field resets

**Date**: 2026-01-27
**Phase**: 3 — Pythonic control-flow consolidation
**Scope**: `python/stdface_main.py`, `test/unit/test_reset_vals_dispatch.py`

### What changed

Converted `_reset_vals()` from ~100 individual field assignments into two
data-driven tables, consistent with the existing `_SOLVER_RESET_SCALARS` /
`_SOLVER_RESET_ARRAYS` pattern already used for solver-specific fields.

#### New data tables

1. **`_COMMON_RESET_SCALARS`** — 70-entry list of `(field_name, sentinel)`
   tuples covering all solver-independent scalar fields:
   - 41 float fields → `NaN_d`
   - 6 integer fields → `NaN_i`  (Height, L, W, S2, ncond, Sz2)
   - 12 complex fields → `NaN_c` (t, tp, tpp, t0, t0p, t0pp, ...)
   - 11 Wannier90 cutoff scalars → `NaN_d`

2. **`_COMMON_RESET_ARRAYS`** — 23-entry list of `(field_name, fill_value)`
   tuples covering all solver-independent array fields:
   - 3 lattice vectors (length, box, direct)
   - 12 spin coupling 3×3 matrices (J, Jp, Jpp, J0, ...)
   - 1 phase vector
   - 6 Wannier90 cutoff arrays (cutoff_tR, cutoff_tVec, ...)

#### Refactored `_reset_vals()`

Before (131 lines, ~100 individual `StdI.x = val` assignments):
```python
StdI.a = NaN_d
StdI.Gamma = NaN_d
StdI.Gamma_y = NaN_d
# ... 97 more assignments ...
```

After (37 lines, 2 loops + 3 special cases):
```python
for name, value in _COMMON_RESET_SCALARS:
    setattr(StdI, name, value)
for name, value in _COMMON_RESET_ARRAYS:
    getattr(StdI, name)[...] = value
# D matrix special case, solver-specific, lBoost
```

### Lines changed

- **Removed**: ~95 lines (individual field assignments)
- **Added**: ~100 lines (data table definitions) + ~37 lines (new function body)
- **Net**: function body reduced from 131 to 37 lines; tables are declarative
  data, not logic

### Design rationale

- **Consistency**: Extends the `_SOLVER_RESET_SCALARS`/`_SOLVER_RESET_ARRAYS`
  pattern (introduced in Step 19) to the solver-independent fields, making the
  entire `_reset_vals` function data-driven.
- **Single-point-of-truth**: Adding a new field requires only a table entry, not
  a new line of imperative code.
- **Grouping preserved**: Fields are grouped by category in the table (lattice,
  magnetic field, spin couplings, hopping, Coulomb, Wannier90 cutoffs) via
  comments, matching the original organisation.
- **Special cases inline**: The D matrix (zeros except D[2,2]=NaN) and `lBoost=0`
  remain as explicit statements since they don't fit the uniform table pattern.

### Test results

- **Baseline**: 850 passed
- **After**: 875 passed (850 + 25 new)
- **New tests**:
  - `TestCommonResetTableStructure` (10 tests) — lists non-empty, entries are
    tuples, fields exist on StdIntList, no duplicates, values are sentinels
  - `TestCommonResetTableContent` (6 tests) — all 12 hopping params present,
    all 13 Coulomb params, all 12 J*All scalars, all 12 J matrices,
    Wannier90 cutoff scalars, Wannier90 cutoff arrays
  - `TestResetVals` (9 tests) — pi/pi180 set, float/int/complex scalars NaN,
    arrays filled, D matrix special case, lBoost=0, solver fields set
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/stdface_main.py`,
`test/unit/test_reset_vals_dispatch.py`, `python/history/refactoring_log.md`

**Suggested next step**: Refactor `general_j()` in
`lattice/interaction_builder.py` to extract the exchange-term dispatch logic,
or tackle `_check_lattice_params()` in `stdface_main.py`.

---

## Step 29 — Extract `_resolve_string_param()` helper in `hphi_writer.py`

**Date**: 2026-01-27
**Phase**: 3 — Pythonic control-flow consolidation
**Scope**: `python/writer/hphi_writer.py`, `test/unit/test_hphi_writer.py`

### What changed

Extracted a reusable `_resolve_string_param()` helper function from 6
identical parameter-resolution blocks in `print_calc_mod()`.

#### Before (lines 213–301, 89 lines):

Six blocks following the same pattern:
```python
if StdI.FIELD == UNSET_STRING:
    StdI.FIELD = "default"
    print("... DEFAULT VALUE ...")
    iLocal = default_int
else:
    print(f"... = {StdI.FIELD}")
    iLocal = DISPATCH.get(StdI.FIELD)
    if iLocal is None:
        print(f"\n ERROR ! ...")
        exit_program(-1)
```

Repeated for: Restart, InitialVecType, EigenVecIO, HamIO, CalcSpec,
OutputExVec — each block ~12–15 lines, total ~89 lines.

#### After (20 lines):

```python
iRestart = _resolve_string_param(
    StdI, "Restart", "Restart", "none", 0, RESTART_TO_INT)
ivt_default = -1 if StdI.method in (TPQ, CTPQ) else 0
iInitialVecType = _resolve_string_param(
    StdI, "InitialVecType", "InitialVecType", "c",
    ivt_default, INITIAL_VEC_TYPE_TO_INT)
InputEigenVec, OutputEigenVec = _resolve_string_param(
    StdI, "EigenVecIO", "EigenVecIO", "none", (0, 0), EIGENVEC_IO_TO_FLAGS)
if StdI.method == MethodType.TIME_EVOLUTION:
    InputEigenVec = 1
iOutputHam, iInputHam = _resolve_string_param(
    StdI, "HamIO", "HamIO", "none", (0, 0), HAM_IO_TO_FLAGS)
iCalcSpec = _resolve_string_param(
    StdI, "CalcSpec", "CalcSpec", "none", 0, CALC_SPEC_TO_INT)
iOutputExVec = _resolve_string_param(
    StdI, "OutputExVec", "OutputExcitedVec", "none",
    0, OUTPUT_EX_VEC_TO_INT)
```

#### New helper: `_resolve_string_param(StdI, field, label, default, default_value, dispatch)`

Handles the uniform pattern:
1. If field == UNSET_STRING → set to default, print default message, return default_value
2. Else → look up in dispatch table, print value, abort on unknown

Supports both scalar (`int`) and tuple (`(int, int)`) return types, used by
EigenVecIO and HamIO which return flag pairs.

### Lines changed

- **Removed**: 89 lines (6 repetitive blocks)
- **Added**: 47 lines (helper function with docstring) + 20 lines (6 call sites)
- **Net**: -22 lines in `print_calc_mod()`, function body reduced from 187 to 148 lines

### Design rationale

- **DRY**: All 6 blocks shared identical logic; the only differences were
  field name, label, default, and dispatch table — exactly what parameters are for.
- **Polymorphic return**: `_resolve_string_param` returns `int | tuple[int, int]`,
  supporting both scalar codes (Restart, CalcSpec, etc.) and flag pairs (EigenVecIO,
  HamIO).
- **Special cases preserved**: InitialVecType's method-dependent default is computed
  before the call.  Time-evolution's `InputEigenVec = 1` override remains inline.
- **Stdout format**: Print format is slightly normalised (uniform field-width
  alignment) but this only affects console log output, not generated `.def` files.

### Test results

- **Baseline**: 875 passed
- **After**: 884 passed (875 + 9 new)
- **New tests**:
  - `TestResolveStringParam` (9 tests):
    - `test_unset_returns_default_value` — UNSET_STRING → default int
    - `test_unset_sets_field_to_default` — field overwritten
    - `test_set_value_returns_dispatched_int` — dispatch lookup
    - `test_set_value_field_unchanged` — no overwrite
    - `test_unknown_value_raises_exit` — SystemExit on bad input
    - `test_tuple_dispatch` — (int, int) return for EigenVecIO
    - `test_tuple_default` — (int, int) default for unset field
    - `test_prints_default_message` — capsys check for DEFAULT
    - `test_prints_set_value` — capsys check for set value
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/writer/hphi_writer.py`,
`test/unit/test_hphi_writer.py`, `python/history/refactoring_log.md`

**Suggested next step**: Apply the same `_resolve_string_param` pattern to
`print_excitation()` in `hphi_writer.py` (196 lines with repeated spectrum-type
dispatch), or refactor `general_j()` in `lattice/interaction_builder.py`.

---

## Step 30 — Extract `_configure_spectrum_ops()` from `print_excitation()`

**Date**: 2026-01-27
**Target**: `python/writer/hphi_writer.py` — `print_excitation()` function

### Problem

`print_excitation()` contained a 79-line if/elif chain (lines 394–473) that
dispatched on `SpectrumType` to configure operator counts, coefficients, and
spin indices.  Two branches were fully duplicated:

- The `UNSET_STRING` default branch (lines 394–412) was identical to the
  explicit `"szsz"` branch (lines 415–431) except for the print message.
- Five spectrum types (`"szsz"`, `"s+s-"`, `"density"`, `"up"`, `"down"`)
  each set `NumOp`, `coef[]`, `spin[][]`, and `SpectrumBody` inline.

### Solution

Extracted a `_configure_spectrum_ops(spectrum_type, model, S2, coef, spin)`
helper that:

1. Takes the resolved spectrum type string (after UNSET_STRING → "szsz"
   defaulting is handled by the caller).
2. Populates `coef` and `spin` arrays in place.
3. Returns `(NumOp, SpectrumBody)`.
4. Raises `SystemExit` for unknown types (matching original behavior).

The default-handling logic (UNSET_STRING → "szsz" with DEFAULT message) stays
in `print_excitation()` — just 5 lines instead of duplicating the full szsz
block.  The 79-line if/elif chain is replaced by a single call:

```python
NumOp, StdI.SpectrumBody = _configure_spectrum_ops(
    StdI.SpectrumType, StdI.model, StdI.S2, coef, spin)
```

Net reduction: 79 lines → 1 line in the caller + 72-line standalone helper
(including docstring).  Eliminated the 17-line szsz duplication entirely.

### Tests

- **Unit**: 903 passed (884 + 19 new in `TestConfigureSpectrumOps`)
  - `test_szsz_hubbard_numop`, `test_szsz_hubbard_body`,
    `test_szsz_hubbard_coef`, `test_szsz_hubbard_spin`
  - `test_szsz_spin_s2_1_numop`, `test_szsz_spin_s2_1_coef`
  - `test_szsz_spin_s2_2_numop`, `test_szsz_spin_s2_2_coef`
  - `test_spsm_hubbard_numop`, `test_spsm_hubbard_coef_and_spin`
  - `test_spsm_spin_s2_2_numop`, `test_spsm_spin_s2_2_coef`
  - `test_density_returns`, `test_density_coef`
  - `test_up_returns`, `test_up_spin`
  - `test_down_returns`, `test_down_spin`
  - `test_unknown_spectrum_type_exits`
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/writer/hphi_writer.py`,
`test/unit/test_hphi_writer.py`, `python/history/refactoring_log.md`

**Suggested next step**: Refactor `general_j()` in
`lattice/interaction_builder.py` (124 lines with nested conditionals for
J/J'/J'' coupling dispatch), or apply a similar extraction to `vector_potential()`
in `hphi_writer.py` (95 lines with 4-branch PumpType dispatch).

---

## Step 31 — Extract `_spin_ladder_factor()` and `_add_spin_half_terms()` from `general_j()`

**Date**: 2026-01-27
**Target**: `python/lattice/interaction_builder.py` — `general_j()` function

### Problem

`general_j()` was 125 lines with two code smells:

1. **Repeated expression**: `math.sqrt(S * (S + 1.0) - Sz * (Sz + 1.0))`
   appeared 6 times across branches (2)–(5) of the spin interaction loop.
   This is the standard spin-ladder matrix element and was computed inline
   each time.

2. **Inline spin-1/2 shortcut**: Lines 286–320 (35 lines) handled the special
   case where at least one site is spin-1/2, adding Hund, Cinter, Ex, and
   PairLift terms with nested solver/model conditionals.  This logic was
   mixed into the main function body, making the overall flow hard to follow.

### Solution

Extracted two helpers:

1. **`_spin_ladder_factor(S, Sz)`** — computes `sqrt(S(S+1) - Sz(Sz+1))`.
   Replaces 6 inline `math.sqrt(...)` calls with a named function that
   clarifies the physics (spin raising/lowering matrix element).

2. **`_add_spin_half_terms(StdI, J, isite, jsite)`** — handles the spin-1/2
   shortcut path.  Returns `(use_z, use_ex)` boolean flags indicating which
   general-path branches are still needed.  Replaced the `ZGeneral`/`ExGeneral`
   integer flags with descriptive boolean names.

The main `general_j()` function now:
- Calls `_add_spin_half_terms()` for the spin-1/2 case (1 line vs 35)
- Uses `_spin_ladder_factor()` for ladder coefficients (readable vs inline sqrt)
- Computes `fi` and `fj` once per loop iteration instead of repeating the sqrt

### Tests

- **Unit**: 920 passed (903 + 17 new)
  - `TestSpinLadderFactor` (7 tests): S=1/2, S=1, S=3/2 cases, boundary
    values (Sz=+S gives 0), match against inline formula
  - `TestAddSpinHalfTerms` (10 tests): Hund/Cinter added, site indices,
    return flags for diagonal/off-diagonal J, Ex/PairLift for diagonal J,
    kondo/mVMC sign, mVMC Jxx≠Jyy fallback, no Ex for off-diagonal J
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/lattice/interaction_builder.py`,
`test/unit/test_interaction_builder.py`, `python/history/refactoring_log.md`

**Suggested next step**: Refactor `vector_potential()` in `hphi_writer.py`
(107 lines with 4-branch PumpType dispatch and repeated time-stepping loops),
or tackle `_print_interactions()` in `writer/common_writer.py`.

---

## Step 32 — Extract pump-type handler functions and dispatch table from `vector_potential()`

**Date**: 2026-01-27
**Target**: `python/writer/hphi_writer.py` — `vector_potential()` function

### Problem

`vector_potential()` was 107 lines with a 4-branch `PumpType` if/elif chain
(lines 610–653).  Three branches (`"pulselaser"`, `"aclaser"`, `"dclaser"`)
each contained the same nested loop structure:

```python
for it in range(StdI.Lanczos_max):
    time = StdI.dt * float(it)
    for ii in range(3):
        StdI.At[it][ii] = <formula>
        Et[it][ii] = <formula>
```

Only the A(t) and E(t) formulas differed.  The `"pulselaser"` branch alone
was 13 lines of dense arithmetic with two multi-line expressions that were
hard to audit.

### Solution

1. Extracted three per-pump-type functions with uniform signature
   `(time, V, freq, tshift, tdump) -> (At, Et)`:
   - **`_pulselaser_At_Et()`** — Gaussian-enveloped cosine pulse.
     Uses intermediate variables (`dt_shift`, `gauss`, `cos_val`, `sin_val`)
     to make the physics readable.
   - **`_aclaser_At_Et()`** — sinusoidal AC laser.
   - **`_dclaser_At_Et()`** — linearly ramped DC laser.

2. Created **`_PUMP_TYPE_HANDLERS`** dispatch table mapping each PumpType
   string to `(PumpBody, handler_fn)`:
   - `"quench"` → `(2, None)` — no field computation needed.
   - `"pulselaser"` / `"aclaser"` / `"dclaser"` → `(1, <function>)`.

3. Rewrote `vector_potential()` to:
   - Resolve UNSET_STRING default (same as before).
   - Look up `(PumpBody, handler_fn)` from `_PUMP_TYPE_HANDLERS`.
   - Error-exit for unknown PumpType (single check vs. `else` branch).
   - Run a single generic time-stepping loop calling `handler_fn()`.

Net reduction: 44-line if/elif chain → 6-line dispatch + generic loop.
The three extracted functions total 63 lines (including docstrings) but
each is independently testable and the physics is clearly documented.

### Tests

- **Unit**: 938 passed (920 + 18 new)
  - `TestPulseLaserAtEt` (5 tests): peak at center, decay far away,
    E=0 at symmetric center, sign flip with negative V, match against
    original inline formula
  - `TestAcLaserAtEt` (3 tests): sin(0)=0 at tshift, E=V*freq at tshift,
    quarter-period gives A=1
  - `TestDcLaserAtEt` (3 tests): A=V*t linear, E=-V constant, A(0)=0
  - `TestPumpTypeHandlers` (7 tests): all 4 entries correct, unknown
    returns None, all one-body have handlers, unknown PumpType exits
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/writer/hphi_writer.py`,
`test/unit/test_hphi_writer.py`, `python/history/refactoring_log.md`

**Suggested next step**: Refactor `_print_interactions()` in
`writer/common_writer.py` (large function with repeated file-writing
patterns), or extract the keyword-parsing if/elif chain in
`keyword_parser.py` into a dispatch table.

---

## Step 33 — Extract `_find_cell_index()` to eliminate duplicated cell-search loops

**Date**: 2026-01-27
**Target**: `python/lattice/site_util.py`, `python/writer/mvmc_variational.py`

### Problem

A 5-line cell-search pattern was duplicated across the codebase:

```python
result = 0
for kCell in range(StdI.NCell):
    if (cellV[0] == StdI.Cell[kCell, 0]
            and cellV[1] == StdI.Cell[kCell, 1]
            and cellV[2] == StdI.Cell[kCell, 2]):
        result = kCell
```

This exact pattern appeared in:
- `site_util.py:find_site()` — 2 searches in a combined loop (10 lines)
- `mvmc_variational.py:generate_orb()` — 2 occurrences (12 lines)
- `mvmc_variational.py:proj()` — 1 occurrence (4 lines + gated inner block)

Total: **5 occurrences**, ~26 lines of duplicated search logic.

### Solution

Extracted `_find_cell_index(StdI, cellV) -> int` into `site_util.py`:
- Linear search returning first matching cell index
- Returns 0 on no match (matching original C behavior)
- Early return on first match (vs. original which always scanned all cells)

Applied to all 5 call sites:
1. `site_util.py:find_site()` — replaced 10-line dual-search loop with 2 calls
2. `mvmc_variational.py:generate_orb()` — replaced 2 × 6-line loops with 2 calls
3. `mvmc_variational.py:proj()` — replaced 4-line search + nested block with
   call + reduced nesting depth by 1 level

Net reduction: ~26 lines of duplicated loops → 5 one-line calls + 22-line
helper (including docstring).

### Tests

- **Unit**: 946 passed (938 + 8 new in `TestFindCellIndex`)
  - `test_finds_first_cell`, `test_finds_middle_cell`, `test_finds_last_cell`
  - `test_all_cells_found` (loop over full chain)
  - `test_no_match_returns_zero`
  - `test_2d_square`, `test_3d_cell` (non-chain lattices)
  - `test_early_return_on_match` (duplicate cell, verifies first-match)
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/lattice/site_util.py` (added helper + applied),
`python/writer/mvmc_variational.py` (applied in `generate_orb()` and `proj()`),
`test/unit/test_site_util.py` (8 new tests),
`python/history/refactoring_log.md`

**Suggested next step**: Refactor `generate_orb()` further — extract the
Kondo orbital assignment block (lines 266–302, 37 lines of repetitive
index arithmetic) into a helper, or refactor `print_jastrow()` in
`mvmc_variational.py` (109 lines with mixed algorithm/IO).

---

## Step 34 — Extract `_assign_orb_sector()` from `generate_orb()`

**Date**: 2026-01-27
**Baseline**: 946 passed (unit), 83/83 integration

### Problem

In `generate_orb()` (`writer/mvmc_variational.py`), the orbital assignment
block had ~50 lines of repetitive index arithmetic for 4 sectors:

1. **Base sector** `(i_off=0, j_off=0)` — non-Kondo code
2. **Kondo sectors** `(half,0)`, `(0,half)`, `(half,half)` — 3 blocks with
   identical structure differing only by site offsets

Each sector performs the same 2 operations:
- If first encounter (`CellDone == 0`), assign a fresh orbital index
- Copy the reference cell's orbital to the current cell

The 4 blocks were written inline with duplicated index arithmetic.

### Solution

Extracted `_assign_orb_sector(StdI, iOrb, anti_val, iCell, jCell, iCell2,
jCell2, isite, jsite, i_off, j_off, is_new) -> int`:
- Computes row/col indices from offsets and cell indices
- Handles both "new orbital" and "copy existing" paths
- Returns updated `iOrb` counter

Replaced the 50-line inline block with a data-driven sector loop:
```python
sectors = [(0, 0)]
if StdI.model == ModelType.KONDO:
    half = StdI.nsite // 2
    sectors += [(half, 0), (0, half), (half, half)]
for isite in range(StdI.NsiteUC):
    for jsite in range(StdI.NsiteUC):
        for i_off, j_off in sectors:
            iOrb = _assign_orb_sector(...)
```

Net reduction: ~50 lines of repetitive inline code → 8-line loop + 20-line
helper (including docstring).

### Tests

- **Unit**: 954 passed (946 + 8 new in `TestAssignOrbSector`)
  - `test_new_orbital_assignment` — verifies fresh iOrb is assigned and counter increments
  - `test_existing_orbital_copy` — verifies copy from reference cell when not new
  - `test_anti_orb_set` — verifies AntiOrb is set correctly
  - `test_kondo_offset` — verifies half-offset indexing for Kondo sectors
  - `test_multiple_sites` — verifies correct indices with NsiteUC > 1
  - `test_returns_incremented_iorb` — verifies return value semantics
  - `test_no_increment_when_not_new` — verifies iOrb unchanged for existing
  - `test_reference_cell_values_used` — verifies ref cell Orb values propagate
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/writer/mvmc_variational.py` (extracted helper +
applied in `generate_orb()`), `test/unit/test_mvmc_variational.py` (8 new
tests), `python/history/refactoring_log.md`

**Suggested next step**: Refactor `print_jastrow()` in `mvmc_variational.py`
(109 lines with mixed algorithm/IO and deep nesting), or refactor
`print_orb_para()` in `mvmc_writer.py` (152 lines, largest writer function).

---

## Step 35 — Extract `_merge_duplicate_terms()` from `print_trans()` and `print_pump()`

**Date**: 2026-01-27
**Baseline**: 954 passed (unit), 83/83 integration

### Problem

The "merge duplicate index quadruples + count non-negligible entries" pattern
appeared in two places:

1. `common_writer.py:print_trans()` — lines 134-148 (15 lines):
   merge transfer terms with identical `(i, s_i, j, s_j)` indices, then
   count entries with `abs > 1e-6`.

2. `hphi_writer.py:print_pump()` — lines 787-800 (14 lines):
   same merge-then-count for pump terms at each timestep.

Both use identical logic: O(n²) pairwise comparison of 4-element index
quadruples, summing amplitudes of duplicates and zeroing the second copy,
then counting survivors above the threshold.

### Solution

Extracted `_merge_duplicate_terms(indx, vals, n) -> int` into
`common_writer.py`:
- Takes generic index array `(n, 4)` and value array `(n,)`
- Merges duplicates in-place, zeros absorbed entries
- Returns count of entries with `abs(val) > 1e-6`

Applied at both call sites:
1. `print_trans()`: replaced 15-line merge+count block with single call
2. `print_pump()`: replaced 14-line merge+count block with single call
   (imported from `common_writer`)

Net reduction: ~29 lines of duplicated inline code → 2 one-line calls +
19-line helper (including docstring).

### Tests

- **Unit**: 964 passed (954 + 10 new in `TestMergeDuplicateTerms`)
  - `test_no_duplicates` — non-duplicate terms unchanged
  - `test_two_duplicates_merged` — duplicate pair merged correctly
  - `test_three_duplicates_merged` — triple merged into one
  - `test_cancellation_yields_zero_count` — cancelling terms give count 0
  - `test_negligible_entries_excluded` — below-threshold entries excluded
  - `test_empty_input` — zero-length input returns zero
  - `test_single_entry` — single entry processed correctly
  - `test_complex_amplitudes` — complex values summed correctly
  - `test_partial_index_match_not_merged` — partial matches not merged
  - `test_n_less_than_array_length` — only first n entries processed
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/writer/common_writer.py` (added helper, applied
in `print_trans()`), `python/writer/hphi_writer.py` (imported helper,
applied in `print_pump()`), `test/unit/test_common_writer.py` (10 new
tests), `python/history/refactoring_log.md`

**Suggested next step**: Refactor `print_jastrow()` in `mvmc_variational.py`
(109 lines with mixed algorithm/IO and deep nesting), or refactor
`print_orb_para()` in `mvmc_writer.py` (152 lines, largest writer function).
Alternatively, eliminate the dead `if has_anti / else` branch in the parallel
section of `print_orb_para()` (lines 213-226) where both branches are
identical.

---

## Step 36 — Extract `_has_anti_period()` and remove dead branch in `print_orb_para()`

**Date**: 2026-01-27
**Baseline**: 964 passed (unit), 83/83 integration

### Problem

Two issues in `mvmc_writer.py`:

1. **Duplicated anti-period check**: The expression
   `(StdI.AntiPeriod[0]==1 or StdI.AntiPeriod[1]==1 or StdI.AntiPeriod[2]==1)`
   was computed inline at 2 call sites: `print_orb()` (line 63-65) and
   `print_orb_para()` (line 186-188).

2. **Dead branch in parallel section**: In `print_orb_para()`'s
   `orbitalidxgen.def` parallel output section (lines 213-226), the
   `if has_anti` and `else` branches produced identical output — the
   `has_anti` flag does not affect the parallel section's format (both
   branches write the `reverse` column unconditionally). This dead branch
   added 8 lines of redundant code.

### Solution

1. Extracted `_has_anti_period(StdI) -> bool` helper (14 lines with
   docstring) that checks whether any direction has anti-periodic boundary.
   Applied at both call sites.

2. Removed the dead `if has_anti / else` branch in the parallel section,
   keeping the single unconditional code path.

Net reduction: 2 × 3-line inline computations + 8-line dead branch →
1 helper + 2 one-line calls. Total: ~14 lines removed, 14-line helper
added — net code is cleaner and the dead branch is eliminated.

### Tests

- **Unit**: 970 passed (964 + 6 new in `TestHasAntiPeriod`)
  - `test_no_anti_period` — all-zero returns falsy
  - `test_first_direction_anti` — AntiPeriod[0]==1 returns truthy
  - `test_second_direction_anti` — AntiPeriod[1]==1 returns truthy
  - `test_third_direction_anti` — AntiPeriod[2]==1 returns truthy
  - `test_all_directions_anti` — all-one returns truthy
  - `test_non_one_value_returns_false` — value 2 (not 1) returns falsy
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/writer/mvmc_writer.py` (added helper, applied
in `print_orb()` and `print_orb_para()`, removed dead branch),
`test/unit/test_mvmc_writer.py` (6 new tests + import),
`python/history/refactoring_log.md`

**Suggested next step**: Refactor `print_jastrow()` in `mvmc_variational.py`
(109 lines with mixed algorithm/IO and deep nesting), or split the 5 phases
of `print_orb_para()` into separate helpers (copy, symmetrise, renumber,
write-para, write-gen).

---

## Step 37 — Split `print_jastrow()` into computation helpers and I/O

**Date**: 2026-01-27
**Baseline**: 970 passed (unit), 83/83 integration

### Problem

`print_jastrow()` in `mvmc_variational.py` was 109 lines mixing two
distinct concerns:

1. **Computation** (80 lines): Two branches that build a `Jastrow` index
   matrix and compute `NJastrow`:
   - **Momentum-projected** (`abs(NMPTrans)==1` or unset): copy orbital
     indices → symmetrise → exclude local-spin sites → renumber
   - **Global optimization** (else): model-specific cell-based assignment
     (spin → single index; hubbard/kondo → per-cell-displacement)

2. **I/O** (20 lines): Write `jastrowidx.def` from the computed matrix.

The two computation branches were deeply nested inline, making them
impossible to test independently from file I/O.

### Solution

Extracted two computation helpers:

1. `_jastrow_momentum_projected(StdI, Jastrow) -> (NJastrow, Jastrow)`:
   30 lines — copy, symmetrise, exclude local spins, renumber. Returns
   both the count and the (possibly new) array since `Jastrow = -1 - Jastrow`
   creates a new object.

2. `_jastrow_global_optimization(StdI, Jastrow) -> NJastrow`:
   40 lines — spin/hubbard/kondo cell-displacement logic. Modifies Jastrow
   in-place, returns count.

`print_jastrow()` reduced from 109 lines to 25 lines: allocate matrix,
dispatch to one of the two helpers, write file.

### Tests

- **Unit**: 980 passed (970 + 10 new)
  - `TestJastrowMomentumProjected` (6 tests):
    - `test_returns_positive_njastrow` — NJastrow > 0
    - `test_off_diagonal_non_negative` — off-diagonal entries ≥ 0
    - `test_hubbard_starts_at_zero` — index 0 present in Hubbard
    - `test_spin_starts_at_minus_one` — Spin local-spin handling
    - `test_symmetrized_output` — output matrix is symmetric
    - `test_local_spin_sites_excluded` — local-spin rows/cols uniform
  - `TestJastrowGlobalOptimization` (4 tests):
    - `test_spin_model_single_index` — NJastrow=1, all zeros
    - `test_hubbard_positive_njastrow` — NJastrow > 0
    - `test_symmetric_jastrow` — symmetric output
    - `test_diagonal_excluded` — on-site terms remain 0
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/writer/mvmc_variational.py` (extracted 2
helpers, simplified `print_jastrow()`),
`test/unit/test_mvmc_variational.py` (10 new tests + imports),
`python/history/refactoring_log.md`

**Suggested next step**: Split the 5 phases of `print_orb_para()` in
`mvmc_writer.py` into separate helpers (copy, symmetrise, renumber,
write-para, write-gen), or refactor `wannier90()` in
`lattice/wannier90.py` (449 lines, 9-level nesting — the largest remaining
function).

---

## Step 38 — Extract `_compute_parallel_orbitals()` from `print_orb_para()`

**Date**: 2026-01-27
**Baseline**: 980 passed (unit), 83/83 integration

### Problem

`print_orb_para()` in `mvmc_writer.py` was 142 lines mixing computation
(3 phases: copy, symmetrise, renumber) with file I/O (2 files). The
computation phases (40 lines) could not be tested independently from
the I/O that writes `orbitalidxpara.def` and `orbitalidxgen.def`.

### Solution

Extracted `_compute_parallel_orbitals(nsite, NOrb, Orb, AntiOrb)` that
returns `(OrbGC, reverse, NOrbGC)`:

1. **Copy**: replicate `Orb` → `OrbGC` and `AntiOrb` → `reverse` as
   plain Python lists-of-lists.
2. **Symmetrise**: for each orbital index, set `OrbGC[j][i] = OrbGC[i][j]`
   and `reverse[j][i] = -reverse[i][j]`.
3. **Renumber**: walk lower triangle, assign negative temporaries, then
   invert to non-negative indices.

`print_orb_para()` now:
- Calls `_compute_parallel_orbitals()` to get `(OrbGC, reverse, NOrbGC)`
- Writes `orbitalidxpara.def` and `orbitalidxgen.def`

Reduced from 142 lines to 80 (I/O) + 62 (helper with docstring).

The helper takes plain arguments (`nsite`, `NOrb`, `Orb`, `AntiOrb`)
rather than `StdI`, making it easier to test without constructing a full
`StdIntList`.

### Tests

- **Unit**: 989 passed (980 + 9 new in `TestComputeParallelOrbitals`)
  - `test_returns_three_values` — correct return types
  - `test_orbgc_dimensions` — correct matrix sizes
  - `test_norbgc_non_negative` — NOrbGC ≥ 0
  - `test_off_diagonal_non_negative` — off-diagonal OrbGC ≥ 0
  - `test_off_diagonal_bounded_by_norbgc` — off-diagonal OrbGC < NOrbGC
  - `test_symmetrised` — OrbGC[i][j] == OrbGC[j][i]
  - `test_reverse_antisymmetric` — reverse[i][j] == -reverse[j][i]
  - `test_uniform_input_single_orbital` — all-zero Orb → NOrbGC=1
  - `test_does_not_modify_input` — input arrays unchanged
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/writer/mvmc_writer.py` (extracted helper,
simplified `print_orb_para()`), `test/unit/test_mvmc_writer.py`
(9 new tests + imports), `python/history/refactoring_log.md`

**Suggested next step**: Refactor `wannier90()` in `lattice/wannier90.py`
(449 lines, 9-level nesting — the largest remaining function), or
extract the `_init_site_sub()` / `_fold_site_sub()` sub-lattice
computation from `mvmc_variational.py` into a shared lattice utility.

---

### 2026-01-27 — Step 39: Extract `_parse_double_counting_mode()` from `wannier90()`

**Target**: `wannier90()` in `lattice/wannier90.py`

**Baseline**: 989 passed (unit), 83/83 integration

### Problem

`wannier90()` in `lattice/wannier90.py` is 449 lines — the largest
remaining function.  Lines 599–616 contained an inline `if/elif/else`
chain (18 lines) converting the `double_counting_mode` string to a
`_DCMode` enum value, with an error exit on invalid input.  This
validation logic was embedded in the middle of the function, making it
hard to test the string→enum conversion independently.

### Solution

Extracted `_parse_double_counting_mode(mode_str: str) -> _DCMode`:

- Accepts the mode string from the input file.
- Returns the corresponding `_DCMode` enum member for recognised values
  (`"none"`, `UNSET_STRING`, `"hartree"`, `"hartree_u"`, `"full"`).
- Calls `exit_program(-1)` for any unrecognised string.
- Full NumPy-style docstring with Parameters/Returns/Raises sections.

`wannier90()` now has a single line:
```python
idcmode = _parse_double_counting_mode(StdI.double_counting_mode)
```

Reduced `wannier90()` by 17 lines (replaced 18 inline lines with
1 call + 35 lines of well-documented, testable helper outside).

### Tests

- **Unit**: 996 passed (989 + 7 new in `TestParseDoubleCountingMode`)
  - `test_none_returns_notcorrect` — `"none"` → `NOTCORRECT`
  - `test_unset_string_returns_notcorrect` — sentinel → `NOTCORRECT`
  - `test_hartree_returns_hartree` — `"hartree"` → `HARTREE`
  - `test_hartree_u_returns_hartree_u` — `"hartree_u"` → `HARTREE_U`
  - `test_full_returns_full` — `"full"` → `FULL`
  - `test_invalid_mode_exits` — unrecognised string → `SystemExit`
  - `test_return_type_is_dcmode` — return is `_DCMode` instance
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/lattice/wannier90.py` (new
`_parse_double_counting_mode()` helper, simplified `wannier90()`),
`test/unit/test_wannier90.py` (7 new tests),
`python/history/refactoring_log.md`

**Suggested next step**: Continue refactoring `wannier90()` — extract
the three repetitive "read cutoff + `_read_w90`" blocks (hopping,
Coulomb, Hund) into a helper, or extract the 210-line transfer/interaction
loop body into sub-functions for each interaction type (hopping, Coulomb U,
Hund J).

---

### 2026-01-27 — Step 40: Extract per-cell interaction helpers from `wannier90()`

**Target**: `wannier90()` in `lattice/wannier90.py`

**Baseline**: 996 passed (unit), 83/83 integration

### Problem

`wannier90()` contained a ~210-line loop body (lines 783–995) that
processed three interaction types (hopping, Coulomb U, Hund J) inline,
with deep nesting (up to 7 levels).  Each interaction type had
distinct logic for local/non-local terms, model-dependent branching,
and double-counting corrections, all interleaved in a single loop
iteration.

### Solution

Extracted three helper functions, each handling one interaction type
for a single unit cell:

- **`_apply_hopping_terms()`** (~40 lines): Processes hopping (t)
  terms — local on-site energies (Hubbard) and non-local
  super-exchange (spin) or hopping integrals (Hubbard).

- **`_apply_coulomb_terms()`** (~70 lines): Processes Coulomb U
  terms — local intra-site Coulomb, non-local inter-site Coulomb,
  and Hartree/Hartree-Fock double-counting corrections.

- **`_apply_hund_terms()`** (~80 lines): Processes Hund J coupling
  terms — exchange, pair-hopping (Hubbard), spin exchange,
  and Hartree/Hartree-Fock double-counting corrections.

The main loop in `wannier90()` now reads:
```python
_apply_hopping_terms(StdI, kCell, iW, iL, iH, NtUJ, tUJ, tUJindx, Uspin)
_apply_coulomb_terms(StdI, kCell, iW, iL, iH, NtUJ, tUJ, tUJindx, idcmode, DenMat)
_apply_hund_terms(StdI, kCell, iW, iL, iH, NtUJ, tUJ, tUJindx, idcmode, DenMat)
```

Replaced ~200 inline lines with 3 function calls plus ~190 lines of
well-documented, independently testable helpers.  Each helper has a
full NumPy-style docstring with Parameters section.

### Tests

- **Unit**: 1007 passed (996 + 11 new)
  - `TestApplyHoppingTerms` (4 tests):
    - `test_none_indices_is_noop` — no-op when tUJindx[0] is None
    - `test_local_hopping_hubbard` — local term adds on-site transfer
    - `test_local_hopping_spin_is_noop` — local term skipped for spin
    - `test_zero_terms_is_noop` — no-op when NtUJ[0] is 0
  - `TestApplyCoulombTerms` (4 tests):
    - `test_none_indices_is_noop` — no-op when tUJindx[1] is None
    - `test_local_coulomb_adds_cintra` — local U adds Cintra
    - `test_local_coulomb_dc_adds_transfer` — DC correction adds transfer
    - `test_zero_terms_is_noop` — no-op when NtUJ[1] is 0
  - `TestApplyHundTerms` (3 tests):
    - `test_none_indices_is_noop` — no-op when tUJindx[2] is None
    - `test_local_hund_term_skipped` — on-site Hund skipped
    - `test_zero_terms_is_noop` — no-op when NtUJ[2] is 0
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/lattice/wannier90.py` (3 new helper
functions, simplified `wannier90()` loop body),
`test/unit/test_wannier90.py` (11 new tests + helpers),
`python/history/refactoring_log.md`

**Suggested next step**: Extract the three repetitive "set cutoff
parameters + call `_read_w90`" blocks (hopping, Coulomb, Hund — lines
650–718) into a helper, or continue reducing `wannier90()` by
extracting the parameter validation section (lines 727–770) into a
dedicated function.

---

### 2026-01-27 — Step 41: Extract `_read_w90_with_cutoff()` from `wannier90()`

**Target**: `wannier90()` in `lattice/wannier90.py`

**Baseline**: 1007 passed (unit), 83/83 integration

### Problem

`wannier90()` contained three repetitive blocks (~25 lines each, ~70
lines total) for reading hopping, Coulomb, and Hund interaction files.
Each block followed the same pattern:

1. Set cutoff threshold and cutoff length with `print_val_d`
2. Set cutoff R-vector defaults with `print_val_i`
3. Set cutoff Vec matrix from `box * 0.5` with `print_val_d`
4. Build filename and call `_read_w90`

The only differences were the parameter names (label case: `t`/`U`/`J`),
default values, file suffix, interaction index, and lambda scaling.

### Solution

Extracted `_read_w90_with_cutoff()` (~30 lines of logic, ~80 lines with
docstring) that parametrises all the differences:

- `label` / `label_suffix`: Controls the parameter name generation
  (e.g. `cutoff_t` vs `cutoff_UR`) to match the C code's inconsistent
  casing (hopping uses lowercase `t` throughout; Coulomb/Hund use
  lowercase threshold `u`/`j` but uppercase `U`/`J` for length/R/Vec).
- `cutoff_R_defaults`: Tuple of `(int | None)` — `None` skips a
  dimension (used by hopping when W/L/Height are unset).
- `cutoff_length_default`: Differs between hopping (`-1.0`) and
  Coulomb/Hund (`0.3`).

The three blocks in `wannier90()` are now three calls:
```python
StdI.cutoff_t, StdI.cutoff_length_t = _read_w90_with_cutoff(
    StdI, "t", "t", "_hr.dat", ...)
StdI.cutoff_u, StdI.cutoff_length_U = _read_w90_with_cutoff(
    StdI, "u", "U", "_ur.dat", ...)
StdI.cutoff_j, StdI.cutoff_length_J = _read_w90_with_cutoff(
    StdI, "j", "J", "_jr.dat", ...)
```

Reduced `wannier90()` by ~35 lines (replaced ~70 inline lines with
~35 lines including the hopping_R_defaults computation + 3 calls).

### Tests

- **Unit**: 1013 passed (1007 + 6 new in `TestReadW90WithCutoff`)
  - `test_returns_updated_cutoff_values` — returns resolved defaults
  - `test_sets_cutoff_R_defaults` — cutoff_R set from defaults
  - `test_skips_cutoff_R_for_none_default` — None skips dimension
  - `test_sets_cutoff_vec_from_box` — cutoff_Vec = box * 0.5
  - `test_reads_terms` — populates NtUJ/tUJ after reading
  - `test_missing_file_is_noop` — gracefully skips missing file
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/lattice/wannier90.py` (new
`_read_w90_with_cutoff()` helper, simplified `wannier90()`),
`test/unit/test_wannier90.py` (6 new tests),
`python/history/refactoring_log.md`

**Suggested next step**: Continue reducing `wannier90()` by extracting
the parameter validation section (lambda/alpha/double-counting checks,
~30 lines) or the Hamiltonian parameter setup (model check + locspinflag,
~30 lines) into dedicated helpers.  Alternatively, move to a different
target: `stdface_main.py` still has ~3,200 lines and could benefit from
further module extraction.

---

## Step 42 — Extract `set_local_spin_flags()` shared utility
**Date**: 2026-01-27
**Target**: `python/lattice/site_util.py` + all 10 lattice modules

### Problem

Every lattice module contained an identical ~13-line block that:
1. Optionally doubles `nsite` for KONDO models
2. Allocates `locspinflag` array
3. Fills the array based on model type (SPIN → S2, HUBBARD → 0, KONDO → first half S2 / second half 0)

This block was copy-pasted across 10 files:
- **Vectorized** (8 files): square, ladder, triangular, honeycomb, kagome, orthorhombic, fc_ortho, pyrochlore
- **Loop-based** (1 file): chain_lattice
- **Simplified** (1 file): wannier90 (no KONDO branch)

Total: ~130 lines of duplicated logic.

### Solution

Added `set_local_spin_flags(StdI, nsite_base)` to `lattice/site_util.py`.
The function takes `nsite_base` (site count before KONDO doubling) as a
parameter since each lattice computes it differently:
- Chain: `StdI.L`
- Ladder: `StdI.L * StdI.NsiteUC`
- Others: `StdI.NsiteUC * StdI.NCell`

Each lattice now calls a single line:
```python
set_local_spin_flags(StdI, StdI.NsiteUC * StdI.NCell)
```

The chain_lattice loop-based code was replaced with the vectorized
version (semantically identical). The wannier90 simplified version now
also goes through the shared function (KONDO branch is unreachable since
wannier90 rejects KONDO models earlier).

### Tests

- **Unit**: 1019 passed (1013 + 6 new in `TestSetLocalSpinFlags`)
  - `test_spin_model` — all flags set to S2
  - `test_hubbard_model` — all flags set to 0
  - `test_kondo_model` — nsite doubled, first half S2, second half 0
  - `test_kondo_with_s2_3` — KONDO with S2=3 (spin-3/2)
  - `test_single_site` — edge case with 1 site
  - `test_kondo_single_site` — KONDO doubles single site to 2
- **Integration**: 83/83 pass (HPhi 59 + mVMC 14 + UHF 2 + H-wave 8)

**Files modified**: `python/lattice/site_util.py` (new function),
`python/lattice/chain_lattice.py`, `python/lattice/square_lattice.py`,
`python/lattice/ladder.py`, `python/lattice/triangular_lattice.py`,
`python/lattice/honeycomb_lattice.py`, `python/lattice/kagome.py`,
`python/lattice/orthorhombic.py`, `python/lattice/fc_ortho.py`,
`python/lattice/pyrochlore.py`, `python/lattice/wannier90.py`
(all 10 updated to call shared function),
`test/unit/test_site_util.py` (6 new tests),
`python/history/refactoring_log.md`

**Suggested next step**: Look for other cross-lattice duplicated patterns
(e.g. memory allocation / ntransMax+nintrMax computation, or the
interaction loop patterns). Alternatively, continue with `wannier90.py`
validation extraction or `stdface_main.py` module splitting.

---

## Step 43 — Eliminate redundant `not_used_c` function

**Date**: 2026-01-27

**Target files**: `python/param_check.py`, all 9 lattice files, `python/stdface_model_util.py`,
`test/unit/test_param_check.py`, `test/unit/test_stdface_model_util.py`

**What changed and why**:
`not_used_c(valname, val)` was a redundant specialisation of `not_used_d(valname, val)`.
`not_used_d` already handles complex values (extracts `val.real` via `isinstance(val, complex)`)
and produces the identical error message and exit behaviour. Every call site that used
`not_used_c` has been replaced with `not_used_d`, and the `not_used_c` function itself has
been removed.

Changes made:
- **`param_check.py`**: Removed `not_used_c` function definition (15 lines) and its entry
  in the module docstring.
- **9 lattice files** (`chain_lattice.py`, `square_lattice.py`, `ladder.py`,
  `triangular_lattice.py`, `honeycomb_lattice.py`, `kagome.py`, `orthorhombic.py`,
  `fc_ortho.py`, `pyrochlore.py`): Replaced all `not_used_c(` calls with `not_used_d(`
  and removed `not_used_c` from imports.
- **`stdface_model_util.py`**: Removed `not_used_c` from re-exports.
- **`test/unit/test_param_check.py`**: Removed `TestNotUsedC` class (2 tests) and
  `not_used_c` from imports and backward-compatibility test.
- **`test/unit/test_stdface_model_util.py`**: Replaced `not_used_c` calls in tests
  with `not_used_d` calls.

**Test results**:
- Unit tests: 1017 passed (2 fewer — removed redundant `TestNotUsedC` tests, whose
  coverage is already provided by `TestNotUsedD.test_handles_complex_nan` and
  `TestNotUsedD.test_exits_for_complex_set`)
- Integration tests: 83/83 passed

**Suggested next step**: Phase 1 cross-lattice extraction is now essentially complete.
Consider beginning Phase 2 (introducing lightweight classes) or Phase 3 (Python idiom
improvements such as replacing sentinel `NaN_i` with `Optional[int]`).

---

## Step 44 — Replace magic number thresholds with named constants

**Date**: 2026-01-27

**Target files**: `python/stdface_vals.py` (new constants),
`python/writer/interaction_writer.py`, `python/writer/common_writer.py`,
`python/writer/hphi_writer.py`, `python/lattice/interaction_builder.py`,
`python/lattice/site_util.py`, `python/lattice/wannier90.py`,
`test/unit/test_stdface_vals.py`

**What changed and why**:
Two magic number thresholds (`0.000001` / `1e-6` / `1.0e-6` and `1.0e-12`) were
used throughout the codebase for amplitude filtering and zero-checking, respectively.
These have been replaced with named constants defined in `stdface_vals.py`:

- **`AMPLITUDE_EPS = 1e-6`** — threshold for treating an amplitude as non-zero in
  output files (transfer, interaction, pump, initial-guess writers and J-matrix
  off-diagonal checks). Replaces 11 occurrences of `0.000001`, `1e-6`, and `1.0e-6`
  across 6 files.
- **`ZERO_BODY_EPS = 1e-12`** — threshold for skipping zero-amplitude terms when
  accumulating one-body (`trans`) and two-body (`intr`) Hamiltonian entries.
  Replaces 2 occurrences of `1.0e-12` in `interaction_builder.py`.

Docstrings that referenced the literal values were also updated to reference the
constant names.

**Files modified**:
- `python/stdface_vals.py` — added `AMPLITUDE_EPS` and `ZERO_BODY_EPS` constants
- `python/writer/interaction_writer.py` — import + 4 replacements + docstring
- `python/writer/common_writer.py` — import + 2 replacements + 3 docstring updates
- `python/writer/hphi_writer.py` — import + 1 replacement + 1 docstring update
- `python/lattice/interaction_builder.py` — import + 2 `ZERO_BODY_EPS` + 2 `AMPLITUDE_EPS`
- `python/lattice/site_util.py` — import + 1 replacement
- `python/lattice/wannier90.py` — import + 2 replacements
- `test/unit/test_stdface_vals.py` — added 5 tests in `TestNumericalTolerances`

**Test results**:
- Unit tests: 1022 passed (+5 new tolerance tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 3 with another Python idiom improvement:
extract an action-dispatch dict in `common_writer.py::_check_conserved_quantities()`
or split a long function like `hphi_writer.py::print_calc_mod()` (129 lines) into
parameter-resolution helpers.

---

### Step 45 — Unify `_write_1idx_file` / `_write_2idx_file` into `_write_interaction_file`

**Date**: 2026-01-27
**Phase**: 3 — Python idiom improvements
**Target**: `python/writer/interaction_writer.py`

**What changed**:

Replaced two nearly identical functions `_write_1idx_file` (37 lines) and
`_write_2idx_file` (37 lines) with a single generic `_write_interaction_file`
(35 lines) that accepts an `n_indices` parameter (1 or 2).

The two old functions differed only in the number of site indices formatted per
output line — `_write_1idx_file` wrote one `{:5d}` field while `_write_2idx_file`
wrote two.  The unified version uses a generator expression:

```python
idx_str = " ".join(f"{indx[k][i]:5d}" for i in range(n_indices))
```

The `_process_interaction` driver was updated to call the unified function,
removing the `if n_indices == 1: … else: …` branching that previously selected
between the two writers.

**Files modified**:
- `python/writer/interaction_writer.py` — removed `_write_1idx_file` and
  `_write_2idx_file`, added `_write_interaction_file`, updated
  `_process_interaction` call site
- `test/unit/test_interaction_writer.py` — added `TestWriteInteractionFile` class
  with 4 tests (1-idx basic, 2-idx basic, all-zero, header format) and updated
  import list

**Test results**:
- Unit tests: 1026 passed (+4 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 3 — consider extracting a helper for the
repeated `tmp_path / monkeypatch.chdir` test pattern, or begin Phase 2 by
introducing a data class for interaction-type metadata (replacing the dict-of-str
`_INTERACTION_TYPES` table).

---

### Step 46 — Extract boost.def output helpers (`write_boost_mag_field`, `write_boost_j_full`, `write_boost_j_symmetric`)

**Date**: 2026-01-27
**Phase**: 3 — Python idiom improvements (DRY)
**Target**: `python/lattice/site_util.py`, 4 lattice files

**What changed**:

Extracted three shared helpers into `site_util.py` that eliminate repeated
`boost.def` formatting code across 4 lattice files:

- **`write_boost_mag_field(fp, StdI)`** — writes the magnetic-field comment and
  data line (`-0.5 * Gamma`, `-0.5 * Gamma_y`, `-0.5 * h`).  Replaces 4
  identical 2–3-line blocks (one per lattice).

- **`write_boost_j_full(fp, J, scale=0.25)`** — writes a full 3×3 J-coupling
  matrix row by row.  Used by chain lattice (2 calls, replacing 2 loop blocks).

- **`write_boost_j_symmetric(fp, J, scale=0.25)`** — writes the upper-triangle
  symmetric form of a 3×3 J-coupling matrix (rows use `J[0,0] J[0,1] J[0,2]` /
  `J[0,1] J[1,1] J[1,2]` / `J[0,2] J[1,2] J[2,2]`).  Used by honeycomb (3
  calls), ladder (5 calls), and kagome (3 calls), replacing 11 blocks of 3
  explicit `fp.write` lines each = 33 inline lines.

Net reduction: ~49 lines of duplicated formatting replaced by 13 one-line helper
calls plus ~50 lines of documented helper definitions — a net size wash but a
significant DRY improvement.

**Files modified**:
- `python/lattice/site_util.py` — added 3 helper functions and updated module
  docstring
- `python/lattice/chain_lattice.py` — import + replaced 2 mag-field lines and
  2 J-matrix loop blocks with helper calls
- `python/lattice/honeycomb_lattice.py` — import + replaced 1 mag-field line and
  3 × 3-line J-matrix blocks with helper calls
- `python/lattice/ladder.py` — import + replaced 1 mag-field line and 5 × 3-line
  J-matrix blocks with helper calls
- `python/lattice/kagome.py` — import + replaced 1 mag-field line and 3 × 3-line
  J-matrix blocks with helper calls
- `test/unit/test_site_util.py` — added 4 test classes (`TestWriteBoostMagField`,
  `TestWriteBoostJFull`, `TestWriteBoostJSymmetric`) with 10 tests

**Test results**:
- Unit tests: 1036 passed (+10 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 3 — convert the
`_parse_double_counting_mode` if/elif chain in `wannier90.py` to a dict dispatch,
or begin Phase 2 by introducing typed dataclasses for interaction-type metadata.

---

### Step 47 — Convert `_parse_double_counting_mode` if/elif to dict dispatch (2026-01-27)

**Phase**: 3 — Pythonic idioms  
**Target**: `python/lattice/wannier90.py`  
**Category**: dict dispatch (if/elif → table lookup)

**What changed**:

Replaced the 5-branch if/elif chain in `_parse_double_counting_mode()` with a
module-level `_DC_MODE_MAP` dictionary that maps mode strings (including the
`UNSET_STRING` sentinel) to `_DCMode` enum members.  The function now performs a
single `_DC_MODE_MAP.get(mode_str)` lookup instead of sequential string
comparisons, with the error path falling through when the key is absent.

**Why**: Dict dispatch is the standard Python idiom for fixed string-to-value
mapping.  It is more concise, easier to extend (add one line), and avoids
repeated `elif` boilerplate.

**Files modified**:
- `python/lattice/wannier90.py` — added `_DC_MODE_MAP` dict; rewrote
  `_parse_double_counting_mode()` to use it
- `test/unit/test_wannier90.py` — added 2 tests (`test_dc_mode_map_keys`,
  `test_dc_mode_map_values`) to validate the dispatch table structure

**Test results**:
- Unit tests: 1038 passed (+2 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: The codebase is now extensively refactored in Phase 3.
Most remaining if/elif chains are physics-specific branching (spectrum type,
Coulomb/Hund term application) not suited to simple dict dispatch.  Consider
beginning Phase 2 — introducing typed dataclasses for interaction-type metadata
or consolidating `StdIntList` fields into logical groups.

---

### Step 48 — Split `wannier90()` into smaller helper functions (2026-01-27)

**Phase**: 2 — Introduce structure (function splitting)  
**Target**: `python/lattice/wannier90.py`  
**Category**: Split long function into focused helpers

**What changed**:

The 202-line `wannier90()` main entry point was split into three extracted
helper functions plus the remaining orchestrator:

1. **`_validate_wannier_params(StdI)`** — checks and stores Hamiltonian
   parameters (`h`, `Gamma`, `Gamma_y`, `S2`/`mu`), reports unused parameters,
   and exits on unsupported Kondo model.  Extracted from 15 lines of inline code.

2. **`_build_wannier_interactions(StdI, NtUJ, tUJ, tUJindx, idcmode, DenMat)`**
   — computes upper-bound array sizes, calls `malloc_interactions`, sets up
   spin super-exchange `Uspin`, and runs the main cell loop applying local,
   hopping, Coulomb, and Hund terms.  Extracted from 48 lines of inline code.

3. **`_write_wan2site(StdI)`** — writes `wan2site.dat` with site-to-cell mapping.
   Extracted from 15 lines of inline code.

The remaining `wannier90()` function is now ~110 lines and reads as a clear
sequential pipeline: geometry → scaling → double-counting → file reads →
validation → interactions → output.

**Why**: The original 202-line function mixed geometry setup, parameter
validation, interaction computation, and file I/O in a single scope.  Splitting
it improves readability, testability, and sets the stage for further Phase 2
class introduction.

**Files modified**:
- `python/lattice/wannier90.py` — added 3 helper functions (`_validate_wannier_params`,
  `_build_wannier_interactions`, `_write_wan2site`); refactored `wannier90()` to call them
- `test/unit/test_wannier90.py` — added 3 test classes (`TestValidateWannierParams`,
  `TestBuildWannierInteractions`, `TestWriteWan2site`) with 7 tests

**Test results**:
- Unit tests: 1045 passed (+7 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by splitting other long functions
(e.g. `print_excitation` in `hphi_writer.py` at 126 lines, or `_read_w90` at
190 lines), or begin introducing a class for the Green-function index generation
in `common_writer.py`.

---

### Step 49 — Split `_read_w90()` post-processing into helpers (2026-01-27)

**Phase**: 2 — Introduce structure (function splitting)  
**Target**: `python/lattice/wannier90.py`  
**Category**: Split long function into focused helpers

**What changed**:

The 190-line `_read_w90()` function was split by extracting two post-processing
helpers, reducing it to ~95 lines:

1. **`_apply_boundary_weights(indx_tot, Weight_tot, nWSC, StdI)`** — computes
   the band-lattice extent per dimension and halves weights at the model
   super-cell boundary for even-sized periodic lattice dimensions.  Returns the
   `Band_lattice` array.  Extracted from 15 lines of inline code.

2. **`_count_and_store_terms(Mat_tot, indx_tot, Weight_tot, nWSC, NsiteUC,
   cutoff, itUJ, NtUJ, tUJindx, tUJ)`** — applies per-WSC weights, counts
   matrix elements above the cutoff, prints the effective-term summary, and
   packs surviving terms into the `tUJ` / `tUJindx` output arrays.  Extracted
   from 42 lines of inline code.

The remaining `_read_w90()` now focuses on file I/O, cutoff application, and
inversion-symmetry deduplication, delegating post-processing to the two helpers.

**Why**: The original 190-line function mixed file parsing, cutoff filtering,
weight adjustment, term counting, and term storage.  Splitting the self-contained
post-processing into named helpers improves readability, testability, and reuse.

**Files modified**:
- `python/lattice/wannier90.py` — added 2 helper functions
  (`_apply_boundary_weights`, `_count_and_store_terms`); refactored `_read_w90()`
  to call them; removed now-unused `Band_lattice` / `Model_lattice` local
  allocations
- `test/unit/test_wannier90.py` — added 2 test classes
  (`TestApplyBoundaryWeights` with 4 tests, `TestCountAndStoreTerms` with 5
  tests)

**Test results**:
- Unit tests: 1054 passed (+9 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by splitting `print_excitation` in
`hphi_writer.py` (126 lines), or begin introducing a `GreenFunctionIndices`
class for the cohesive group of Green-function index generators in
`common_writer.py`.

---

### Step 50 — Split `print_excitation()` into smaller helpers (2026-01-27)

**Phase**: 2 — Introduce structure (function splitting)  
**Target**: `python/writer/hphi_writer.py`  
**Category**: Split long function into focused helpers

**What changed**:

The 126-line `print_excitation()` function was split by extracting two helpers,
reducing the orchestrator to ~32 lines:

1. **`_compute_fourier_coefficients(StdI)`** — computes per-site Fourier
   coefficients :math:`\exp(2\pi i \mathbf{q} \cdot \mathbf{r})` for spectrum
   excitations, with Kondo-model duplication.  Returns ``(fourier_r, fourier_i)``
   lists.  Extracted from 18 lines of inline code.

2. **`_write_excitation_file(StdI, NumOp, coef, spin, fourier_r, fourier_i)`**
   — writes ``single.def`` (SpectrumBody=1) or ``pair.def`` (SpectrumBody=2)
   with the appropriate header and formatted data lines.  Handles Kondo half-site
   count for single excitations.  Extracted from 35 lines of inline code.

The remaining `print_excitation()` now reads as a clean pipeline: allocate arrays
→ resolve SpectrumType → configure operators → compute Fourier coefficients →
write file.

**Why**: The original 126-line function mixed parameter setup, numerical
computation (Fourier phases), and file I/O.  Splitting these concerns improves
readability and makes each piece independently testable.

**Files modified**:
- `python/writer/hphi_writer.py` — added 2 helper functions
  (`_compute_fourier_coefficients`, `_write_excitation_file`); refactored
  `print_excitation()` to call them; removed unused local Fourier arrays
- `test/unit/test_hphi_writer.py` — updated imports; added 2 test classes
  (`TestComputeFourierCoefficients` with 4 tests, `TestWriteExcitationFile` with
  4 tests)

**Test results**:
- Unit tests: 1062 passed (+8 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by introducing a `GreenFunctionIndices`
class in `common_writer.py` to encapsulate the cohesive group of Green-function
index generators (`_green1_indices_corr`, `_green1_indices_raw`,
`_green2_indices_corr`, `_green2_indices_raw`, and shared helpers `_spin_max`,
`_skip_local_spin_pair`, `_kondo_site`).

---

### Step 51 — Introduce `GreenFunctionIndices` class (2026-01-27)

**Phase**: 2 — Introduce structure (class introduction)
**Target**: `python/writer/common_writer.py`
**Category**: Encapsulate cohesive function group into a class

**What changed**:

Introduced `GreenFunctionIndices`, the first Phase 2 class, encapsulating all
Green-function index generation logic that previously lived as 7 module-level
functions:

1. **Class `GreenFunctionIndices`** — constructor takes 5 parameters (`nsite`,
   `NsiteUC`, `locspinflag`, `is_kondo`, `is_mvmc`).  Provides 7 methods:
   - `spin_max(site)` — max spin index for a site
   - `skip_local_spin_pair(site_a, site_b)` — local-spin pair filter
   - `kondo_site(isite)` — Kondo UC-to-physical site mapping
   - `green1_corr()` — one-body indices, correlation mode
   - `green1_raw()` — one-body indices, raw mode
   - `green2_corr()` — two-body indices, correlation mode
   - `green2_raw()` — two-body indices, raw mode

2. **Updated `print_1_green` and `print_2_green`** to instantiate
   `GreenFunctionIndices` and call methods instead of module-level functions.

3. **Kept module-level thin wrappers** (`_spin_max`, `_skip_local_spin_pair`,
   `_kondo_site`, `_green1_indices_corr`, `_green1_indices_raw`,
   `_green2_indices_corr`, `_green2_indices_raw`) for backward compatibility
   with existing tests and callers.  The wrappers now delegate to the class.

**Why**: The 7 functions all share the same set of parameters (`nsite`,
`NsiteUC`, `locspinflag`, Kondo/mVMC flags) and form a cohesive group — the
classic signal for introducing a class.  The class eliminates repeated parameter
threading and makes the relationship between helpers and generators explicit.

**Files modified**:
- `python/writer/common_writer.py` — added `GreenFunctionIndices` class; updated
  module docstring; converted `print_1_green`/`print_2_green` to use class;
  kept thin wrapper functions for backward compatibility
- `test/unit/test_common_writer.py` — added import for `GreenFunctionIndices`;
  added 4 test classes:
  - `TestGreenFunctionIndicesConstruction` (3 tests)
  - `TestGreenFunctionIndicesHelpers` (7 tests)
  - `TestGreenFunctionIndicesGreen1` (4 tests)
  - `TestGreenFunctionIndicesGreen2` (4 tests)

**Test results**:
- Unit tests: 1080 passed (+18 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by extracting a `KeywordParser` class
from `keyword_parser.py`, or split `print_jastrow()` in `mvmc_variational.py`
(109 lines with mixed algorithm/IO).

---

### Step 52 — Split `print_calc_mod()` into smaller helpers (2026-01-27)

**Phase**: 2 — Introduce structure (function splitting)
**Target**: `python/writer/hphi_writer.py`
**Category**: Split long function into focused helpers

**What changed**:

The 122-line `print_calc_mod()` function was split by extracting two helpers:

1. **`_validate_ngpu_scalapack(StdI)`** — validates the optional `NGPU` and
   `Scalapack` parameters, printing values and calling `exit_program` if
   out of range.  Extracted from 16 lines of inline validation code.

2. **`_write_calcmod_file(StdI, iCalcType, iCalcModel, ...)`** — writes
   `calcmod.def` with the 11 resolved integer parameters plus conditional
   NGPU/Scalapack lines.  Extracted from 25 lines of file-I/O code.

The remaining `print_calc_mod()` now reads as a clean pipeline: resolve method
→ resolve model → resolve 5 string parameters → validate → write file.

**Why**: The original 122-line function mixed three concerns: parameter
resolution (string→integer mappings), validation (NGPU/Scalapack bounds), and
file I/O.  Splitting these improves readability and makes each piece
independently testable.

**Files modified**:
- `python/writer/hphi_writer.py` — added 2 helper functions
  (`_validate_ngpu_scalapack`, `_write_calcmod_file`); refactored
  `print_calc_mod()` to call them
- `test/unit/test_hphi_writer.py` — updated imports; added 2 test classes
  (`TestValidateNGPUScalapack` with 6 tests, `TestWriteCalcmodFile` with
  5 tests)

**Test results**:
- Unit tests: 1091 passed (+11 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by splitting `print_jastrow()` in
`mvmc_variational.py` (109 lines with mixed algorithm/IO), or extract a
`KeywordParser` class from `keyword_parser.py`.

---

### Step 53 — Split `print_gutzwiller()` into smaller helpers (2026-01-27)

**Phase**: 2 — Introduce structure (function splitting)
**Target**: `python/writer/mvmc_writer.py`
**Category**: Split long function into focused helpers

**What changed**:

The 110-line `print_gutzwiller()` function was split by extracting three helpers:

1. **`_gutzwiller_momentum_projected(StdI, Gutz)`** — computes Gutzwiller
   indices in momentum-projected mode (``abs(NMPTrans) == 1`` or unset).
   Uses diagonal orbital indices, excludes local-spin sites, and renumbers
   with negative temporaries.  Returns ``NGutzwiller``.

2. **`_gutzwiller_global_optimization(StdI, Gutz)`** — computes Gutzwiller
   indices in global-optimisation mode (all other ``NMPTrans``).  Model-
   dependent assignment: Hubbard→site mod NsiteUC, Spin→all 0,
   Kondo→conduction 0 / localized isite+1.  Returns ``NGutzwiller``.

3. **`_write_gutzwiller_file(StdI, NGutzwiller, Gutz)`** — writes
   ``gutzwilleridx.def`` with header, per-site index data, and per-index
   optimisation flags.

The remaining `print_gutzwiller()` is now a 10-line orchestrator: allocate
array → dispatch to computation helper → write file.

**Why**: The original 110-line function mixed two distinct computation
algorithms (momentum-projected vs. global-optimisation) with file I/O.
Splitting separates concerns and makes each algorithm independently testable.

**Files modified**:
- `python/writer/mvmc_writer.py` — added 3 helper functions; refactored
  `print_gutzwiller()` to call them
- `test/unit/test_mvmc_writer.py` — updated imports; added 3 test classes:
  - `TestGutzwillerMomentumProjected` (4 tests)
  - `TestGutzwillerGlobalOptimization` (3 tests)
  - `TestWriteGutzwillerFile` (3 tests)

**Test results**:
- Unit tests: 1101 passed (+10 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by splitting `print_orb_para()` in
`mvmc_writer.py` (92 lines writing two files), or extract the InterAll 3-pass
merge logic in `interaction_writer.py` into a more structured form.

---

### Step 54 — Split `print_orb_para()` into file-writing helpers (2026-01-27)

**Phase**: 2 — Introduce structure (function splitting)
**Target**: `python/writer/mvmc_writer.py`
**Category**: Split long function into focused helpers

**What changed**:

The 92-line `print_orb_para()` function wrote two separate files inline.
Extracted two file-writing helpers:

1. **`_write_orbitalidxpara(nsite, ComplexType, OrbGC, reverse, NOrbGC)`** —
   writes ``orbitalidxpara.def`` with header, upper-triangle site pairs
   with orbital indices and sign-reversal flags, and per-index optimisation
   lines.

2. **`_write_orbitalidxgen(StdI, OrbGC, reverse, NOrbGC)`** — writes
   ``orbitalidxgen.def`` with header, anti-parallel section (using
   ``AntiOrb`` if anti-periodic boundaries exist), parallel section
   (upper triangle with spin-0 and spin-1 blocks), and optimisation lines.

The remaining `print_orb_para()` is now a 12-line orchestrator: compute
parallel orbitals → write para file → write gen file.

**Why**: The original function interleaved two complete file-write operations.
Splitting makes each file's format independently testable and reduces the
cognitive load of understanding what `print_orb_para` does.

**Files modified**:
- `python/writer/mvmc_writer.py` — added 2 helper functions
  (`_write_orbitalidxpara`, `_write_orbitalidxgen`); refactored
  `print_orb_para()` to call them
- `test/unit/test_mvmc_writer.py` — updated imports; added 2 test classes:
  - `TestWriteOrbitalidxpara` (4 tests)
  - `TestWriteOrbitalidxgen` (4 tests)

**Test results**:
- Unit tests: 1109 passed (+8 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by extracting the InterAll 3-pass
merge logic in `interaction_writer.py` into a more structured form, or
introduce a `KeywordParser` class from `keyword_parser.py`.

### Step 55 — Extract `_write_gnuplot_header()` from `init_site()` (2026-01-27)

**Phase**: 2 — Introduce structure (function splitting)
**Target**: `python/lattice/site_util.py`
**Category**: Extract focused helper from long function

**What changed**:

The `init_site()` function (~150 lines) contained an inline block (~38 lines)
that writes the gnuplot header for `lattice.gp` in 2-D lattices.  This block
computes the four corner positions of the super-cell, sets axis ranges,
styles, and draws boundary arrows.

1. **Extracted `_write_gnuplot_header(fp, StdI)`** — writes the full gnuplot
   header including commented PDF terminal setup, axis ranges, layout commands
   (`set size square`, `unset key/tics/border`), three style lines, and four
   boundary arrows connecting the corners 0→1→3→2→0.

2. **Simplified the arrow-writing code** — replaced a 4-case `if/elif` chain
   (each case writing one `set arrow` line with hard-coded corner indices) with
   a concise loop over `corners = [(0, 1), (1, 3), (3, 2), (2, 0)]`.

The remaining `init_site()` is reduced from ~150 to ~115 lines and delegates
to `_write_gnuplot_header(fp, StdI)` when `dim == 2`.

**Why**: The gnuplot header block was a self-contained piece of output logic
embedded in a function that also handles cell allocation, reciprocal-box
computation, and phase setup.  Extracting it makes both the header and the
remaining init_site easier to understand and test independently.

**Files modified**:
- `python/lattice/site_util.py` — added `_write_gnuplot_header(fp, StdI)`;
  replaced inline block in `init_site()` with a call to it; simplified
  4-case if/elif to loop over corner pairs
- `test/unit/test_site_util.py` — updated imports to include
  `_write_gnuplot_header`; added `TestWriteGnuplotHeader` class (6 tests):
  writes xrange, writes yrange, writes four arrows, writes style lines,
  corner positions for identity direct lattice, layout commands

**Test results**:
- Unit tests: 1115 passed (+6 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by targeting remaining long functions
or introducing classes for repeated patterns (e.g., `KeywordParser` class from
`keyword_parser.py`, or structuring the InterAll 3-pass merge logic in
`interaction_writer.py`).

### Step 56 — Split `_write_wannier90()` into matrix-building and body-writing helpers (2026-01-27)

**Phase**: 2 — Introduce structure (function splitting)
**Target**: `python/writer/export_wannier90.py`
**Category**: Split long function into focused helpers

**What changed**:

The 107-line `_write_wannier90()` function interleaved matrix construction
with file I/O.  Split into two focused helpers:

1. **`_build_wannier_matrix(nintr_table, intr_table, nsiteuc, nspin)`** —
   scans interaction items to determine coordinate half-ranges, allocates
   a flat complex matrix, populates it from the interaction table, and fills
   Hermitian-conjugate entries.  Returns ``(rr, nvol, matrix)``.

2. **`_write_wannier_body(fp, rr, nvol, nsiteuc, nspin, matrix)`** —
   iterates over all real-space cells and orbital/spin pairs, writing one
   line per entry.  Uses extended format (with ``s``, ``t`` columns) when
   ``nspin > 1``, compact format otherwise.

The remaining `_write_wannier90()` is now a ~25-line orchestrator: build
matrix → open file → write header → write body → close file.

**Why**: The original function mixed data transformation (range computation,
matrix allocation, Hermitian-conjugate filling) with file I/O (header format,
body iteration, spin-dependent column layout).  Separating these makes each
piece independently testable and easier to understand.

**Files modified**:
- `python/writer/export_wannier90.py` — added `_build_wannier_matrix()` and
  `_write_wannier_body()`; refactored `_write_wannier90()` to call them
- `test/unit/test_export_wannier90.py` — added `io` import; added 2 test
  classes:
  - `TestBuildWannierMatrix` (5 tests): rr computation, range expansion,
    value placement, Hermitian conjugate, spin matrix size
  - `TestWriteWannierBody` (4 tests): compact format, extended format,
    correct values, multi-volume iterations

**Test results**:
- Unit tests: 1124 passed (+9 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by splitting `_configure_spectrum_ops()`
(91 lines) in `hphi_writer.py` into per-spectrum-type handlers, or split
`_check_mod_para_mvmc()` (63 lines) in `common_writer.py` into parameter-group
helpers.

### Step 57 — Dict-dispatch `_configure_spectrum_ops()` via per-type handlers (2026-01-28)

**Phase**: 2/3 — Introduce structure + Python idioms (dict dispatch)
**Target**: `python/writer/hphi_writer.py`
**Category**: Replace if/elif chain with dict dispatch and per-type handler functions

**What changed**:

The 91-line `_configure_spectrum_ops()` function contained a five-branch
`if/elif` chain dispatching by spectrum-type string (`"szsz"`, `"s+s-"`,
`"density"`, `"up"`, `"down"`).  Refactored into:

1. **`_spectrum_szsz(model, S2, coef, spin)`** — handles SzSz spectrum;
   for SPIN models produces `S2+1` operators, otherwise 2 operators.

2. **`_spectrum_spsm(model, S2, coef, spin)`** — handles S+S- spectrum;
   for high-spin SPIN models produces ladder coefficients, otherwise 1
   operator.

3. **`_spectrum_density(model, S2, coef, spin)`** — handles density
   spectrum; always 2 operators with coefficient 1.

4. **`_spectrum_up(model, S2, coef, spin)`** — handles spin-up single-
   particle spectrum; 1 operator.

5. **`_spectrum_down(model, S2, coef, spin)`** — handles spin-down
   single-particle spectrum; 1 operator.

6. **`_SPECTRUM_HANDLERS`** — dict mapping spectrum-type strings to the
   handler functions above.

The remaining `_configure_spectrum_ops()` is now a 5-line dispatcher:
look up handler in `_SPECTRUM_HANDLERS`, call it if found, otherwise
print error and exit.

**Why**: The if/elif chain was a classic "switch on string" pattern that
maps naturally to a dict dispatch.  Each handler is now independently
testable, has its own docstring, and is easier to understand in isolation.
This also follows the same pattern used for pump-type handlers
(`_PUMP_TYPE_HANDLERS`) established in Step 32.

**Files modified**:
- `python/writer/hphi_writer.py` — added 5 handler functions
  (`_spectrum_szsz`, `_spectrum_spsm`, `_spectrum_density`, `_spectrum_up`,
  `_spectrum_down`) and `_SPECTRUM_HANDLERS` dispatch dict; simplified
  `_configure_spectrum_ops()` to use the dispatch table
- `test/unit/test_hphi_writer.py` — updated imports; added 5 test classes:
  - `TestSpectrumHandlersSzSz` (1 test): handler matches dispatcher
  - `TestSpectrumHandlersSpsm` (1 test): handler matches dispatcher
  - `TestSpectrumHandlersDensity` (1 test): handler matches dispatcher
  - `TestSpectrumHandlersUpDown` (2 tests): up/down match dispatcher
  - `TestSpectrumHandlersDispatch` (4 tests): all 5 types registered,
    correct function references

**Test results**:
- Unit tests: 1133 passed (+9 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by splitting
`_check_mod_para_mvmc()` (63 lines) in `common_writer.py` into
parameter-group helpers, or split `general_j()` (88 lines) in
`interaction_builder.py` into per-term helpers.

### Step 58 — Extract `_compute_reciprocal_box()` and `_enumerate_cells()` from `init_site()` (2026-01-28)

**Phase**: 2 — Introduce structure (function splitting)
**Target**: `python/lattice/site_util.py`
**Category**: Extract focused helpers from long function

**What changed**:

The `init_site()` function (~115 lines after Step 55) contained two
substantial inline blocks for computing the reciprocal lattice and
enumerating super-cell contents.  Extracted into:

1. **`_compute_reciprocal_box(StdI)`** (~25 lines) — computes ``NCell``
   as the determinant of ``box`` and ``rbox`` as the cofactor matrix.
   Flips sign if determinant is negative.  Exits on degenerate (zero-det)
   super-cell.

2. **`_enumerate_cells(StdI)`** (~25 lines) — finds the bounding box of
   all corner vertices, then iterates over the integer grid to find which
   coordinates fold back to ``nBox == [0,0,0]``, filling ``StdI.Cell``.

The remaining `init_site()` is now ~75 lines: parse L/W/H or box
parameters → fix 2-D direct matrix → compute phases → allocate tau →
call `_compute_reciprocal_box` → call `_enumerate_cells` → write gnuplot.

**Why**: The reciprocal-box computation and cell enumeration are
independent, self-contained algorithms.  Extracting them makes each
piece independently testable (especially the degenerate-cell error path
and tilted-box enumeration), and reduces `init_site()` to a cleaner
sequential orchestrator.

**Files modified**:
- `python/lattice/site_util.py` — added `_compute_reciprocal_box(StdI)`
  and `_enumerate_cells(StdI)`; replaced inline blocks (4) and (5) in
  `init_site()` with calls to them
- `test/unit/test_site_util.py` — updated imports; added 2 test classes:
  - `TestComputeReciprocalBox` (5 tests): identity box, 1D chain, 2D
    square, negative determinant flip, zero determinant exit
  - `TestEnumerateCells` (4 tests): 1D chain cells, 2D square cells,
    single cell at origin, tilted box distinctness

**Test results**:
- Unit tests: 1142 passed (+9 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by extracting helper functions
from the large lattice builders (chain, square, honeycomb, kagome, etc.)
which share a common pattern of geometry setup + model-type interaction
loops, or split `_export_transfer()` (90 lines) in `export_wannier90.py`.

---

### Step 59 — Extract `open_lattice_gp` / `close_lattice_gp` helpers
**Date**: 2026-01-28
**Phase**: 2 (function splitting) + 3 (Python idioms — DRY)

**Target files**:
- `python/lattice/site_util.py` (added helpers)
- `python/lattice/chain_lattice.py` (consumer)
- `python/lattice/square_lattice.py` (consumer)
- `python/lattice/triangular_lattice.py` (consumer)
- `python/lattice/honeycomb_lattice.py` (consumer)
- `python/lattice/kagome.py` (consumer)
- `python/lattice/ladder.py` (consumer)

**What changed and why**:
All six 2D lattice builder modules duplicated a 3-line open pattern
(`fp = None; if solver != HWAVE or lattice_gp == 1: fp = open(...)`)
and a 3–4 line close pattern (`if ...: fp.write(footer); fp.close();
print_geometry(StdI)`).  This was consolidated into two shared helpers:

- `open_lattice_gp(StdI)` — conditionally opens `lattice.gp`, returns
  file handle or `None`.
- `close_lattice_gp(fp, StdI)` — writes the gnuplot footer, closes the
  file (if open), and calls `print_geometry(StdI)`.
- `_LATTICE_GP_FOOTER` — module-level constant for the footer string.

Each lattice file was updated:
- 3-line open block → `fp = open_lattice_gp(StdI)` (1 line)
- 3–4 line close + geometry block → `close_lattice_gp(fp, StdI)` (1 line)
- Removed `from .geometry_output import print_geometry` (now in site_util)
- Removed `SolverType` import where no longer needed (5 of 6 files;
  chain keeps it for `chain_boost`)

Net line change: ~−30 across the 6 consumers, +40 in site_util (helpers
+ docstrings).

**New tests** (9 tests in `test/unit/test_site_util.py`):
- `TestOpenLatticeGp` (5 tests): HPhi returns file, mVMC returns file,
  HWAVE returns None, HWAVE+lattice_gp=1 returns file, file created on disk
- `TestCloseLatticeGp` (4 tests): footer written & closed, None fp safe,
  geometry.dat created (print_geometry called), footer constant content

**Test results**:
- Unit tests: 1151 passed (+9 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Continue Phase 2 by splitting
`_export_transfer()` (90 lines) in `export_wannier90.py`, or examine
remaining long functions like `_check_mod_para_mvmc()` (63 lines) or
the boost functions.

---

### Step 60 — Extract `_build_transfer_table()` from `_export_transfer()`
**Date**: 2026-01-28
**Phase**: 2 (function splitting — separate computation from I/O)

**Target file**: `python/writer/export_wannier90.py`

**What changed and why**:
`_export_transfer()` was 90 lines mixing computation (converting
accumulated absolute-index entries into relative-coordinate
`_IntrItem` entries with deduplication and sign-flip) with I/O
(writing the Wannier90 file).

Extracted the computation core into a new function:

- `_build_transfer_table(StdI, nintr, intr_index, intr_value, spin_dep)`
  — converts accumulated transfer entries from absolute site indices
  to deduplicated `_IntrItem` entries in relative coordinates.
  Handles sign-flip (`*= -1`), spin-dependent filtering, and
  consistency warnings for duplicate keys with different values.
  Returns `list[_IntrItem]`.

`_export_transfer()` is now a thin orchestrator (~20 lines):
accumulate → build table → write file.

**New tests** (9 tests in `test/unit/test_export_wannier90.py`):
- `TestBuildTransferTable` (9 tests): empty input, single entry,
  sign flip, deduplication, spin_dep=0 filters, spin_dep=0 keeps (0,0),
  spin_dep=1 keeps all, relative coordinates, inconsistent values warning

**Test results**:
- Unit tests: 1160 passed (+9 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Apply the same computation/I/O split to
`_export_inter()` (52 lines) and `_export_coulomb_intra()` (42 lines),
or examine the longer lattice builder functions (honeycomb 252 lines,
kagome 247 lines) for structural simplification.

---

### Step 61 — Extract `_build_inter_table()` and `_build_coulomb_intra_table()`
**Date**: 2026-01-28
**Phase**: 2 (function splitting — separate computation from I/O)

**Target file**: `python/writer/export_wannier90.py`

**What changed and why**:
Following Step 60's suggestion, applied the same computation/I/O
separation to the remaining two export functions:

1. `_export_inter()` (52 lines) mixed computation (converting
   accumulated absolute-index entries into relative-coordinate
   `_IntrItem` entries with deduplication) with file I/O.
   Extracted `_build_inter_table(StdI, nintr, intr_index, intr_value)`
   — like `_build_transfer_table` but without sign-flip, dedup by
   `(rr, a, b)` (no spin dimension), spin always `(0, 0)`.

2. `_export_coulomb_intra()` (42 lines) similarly mixed computation
   with I/O. Extracted `_build_coulomb_intra_table(StdI, nintr,
   intr_index, intr_value)` — on-site only (`rr=[0,0,0]`, `a==b`),
   dedup by `isite`.

Both `_export_inter()` and `_export_coulomb_intra()` are now thin
orchestrators (~15-20 lines): accumulate → build table → write file.

The three builder functions share a common structure but differ in
dedup keys (transfer: `(rr,a,b,s,t)`, inter: `(rr,a,b)`, coulomb:
`(a)`) making a single shared extraction non-trivial. Each keeps its
own focused builder for clarity.

**New tests** (15 tests in `test/unit/test_export_wannier90.py`):
- `TestBuildInterTable` (8 tests): empty input, single entry,
  no sign flip, deduplication, different sites kept, spin indices
  always zero, relative coordinate length, inconsistent values warning
- `TestBuildCoulombIntraTable` (7 tests): empty input, single entry,
  dedup by site, different UC sites kept, spin indices zero, a equals b,
  inconsistent values warning

**Test results**:
- Unit tests: 1175 passed (+15 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: The three Wannier90 builder functions now
share visible structural similarity (loop over entries, unfold, dedup,
warn). A possible Phase 3 step could unify them via a shared helper
with a key-extraction callback. Alternatively, examine the longer
lattice builder functions (honeycomb 252 lines, kagome 247 lines) for
structural simplification, or tackle `_check_mod_para_mvmc()` (63 lines).

---

### Step 62 — Extract `_write_gnuplot_bond()` from `set_label()`
**Date**: 2026-01-28
**Phase**: 2 (function splitting — eliminate duplication)

**Target file**: `python/lattice/site_util.py`

**What changed and why**:
`set_label()` contained two nearly identical I/O blocks (lines 473–484
and 499–510 before this change): each wrote two `set label` commands
(with single/double-digit formatting) and one `set arrow` command.
The only difference between the two call sites was the computed
coordinates and site indices — the I/O pattern was identical.

Extracted the repeated I/O into a new helper:

- `_write_gnuplot_bond(fp, isite, jsite, xi, yi, xj, yj, connect)`
  — writes two gnuplot label commands (using 1-digit or 2-digit
  format depending on site index) and, when `connect < 3`, an arrow
  command connecting the two positions.

`set_label()` now calls `_write_gnuplot_bond()` twice (once for the
reversed bond, once for the normal bond), replacing 24 lines of
duplicated conditional I/O with 2 one-line calls.

Net line change: +35 (new helper with docstring), −22 (removed
duplicated blocks) = +13 overall (docstring accounts for most of it).

**New tests** (8 tests in `test/unit/test_site_util.py`):
- `TestWriteGnuplotBond` (8 tests): single-digit sites, double-digit
  sites, arrow written for connect < 3, no arrow for connect >= 3,
  positions in output, exactly 3 lines with arrow, exactly 2 lines
  without arrow, mixed digit sites

**Test results**:
- Unit tests: 1183 passed (+8 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: The lattice builder model-dispatch pattern
(repeated `if model == SPIN: general_j() else: hopping(); coulomb()`
blocks across 6 lattice files) is a good Phase 2 target. Also consider
the longer lattice builders (honeycomb 252 lines, kagome 247 lines)
for structural simplification, or unify the three Wannier90 builder
functions via a shared helper with a key-extraction callback.

---

### Step 63 — Extract `add_neighbor_interaction()` for lattice builders
**Date**: 2026-01-28
**Phase**: 2 (function splitting — eliminate cross-file duplication)

**Target files**:
- `python/lattice/interaction_builder.py` (new function)
- `python/lattice/chain_lattice.py`
- `python/lattice/square_lattice.py`
- `python/lattice/triangular_lattice.py`
- `python/lattice/honeycomb_lattice.py`
- `python/lattice/kagome.py`
- `python/lattice/ladder.py`

**What changed and why**:
All 6 lattice builder functions repeated the same 5-line pattern for
each neighbor interaction bond:

```python
isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, dW, dL, ucI, ucJ, connect)
if StdI.model == ModelType.SPIN:
    general_j(StdI, J, StdI.S2, StdI.S2, isite, jsite)
else:
    hopping(StdI, Cphase * t, isite, jsite, dR)
    coulomb(StdI, V, isite, jsite)
```

This appeared 3 times (chain), 6 times (square), 9 times (triangular),
12 times (honeycomb), 10 times (kagome), and 5 times (ladder) — 45
instances total.

Added `add_neighbor_interaction(StdI, fp, iW, iL, diW, diL, isiteUC,
jsiteUC, connect, J, t, V)` to `interaction_builder.py` which combines
`set_label` with the model-dependent dispatch. Each lattice builder now
replaces the 5-line block with a single 2-line call.

Each of the 6 lattice files:
- Updated import from `interaction_builder` to include
  `add_neighbor_interaction` (removed `hopping`, `coulomb` imports)
- Replaced all neighbor blocks with `add_neighbor_interaction()` calls

Net line reduction: ~135 lines removed across the 6 lattice files,
+55 lines added (new function with full docstring + local import).

**New tests** (8 tests in `test/unit/test_interaction_builder.py`):
- `TestAddNeighborInteraction` (8 tests): returns site indices,
  spin model calls general_j, hubbard calls hopping, hubbard calls
  coulomb, spin model no hopping, hubbard model no general_j,
  fp=None no error, gnuplot output written

**Test results**:
- Unit tests: 1191 passed (+8 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: The local-term block (mag_field + general_j
for spin, hubbard_local + Kondo terms for electron) is also repeated
identically across all 6 lattice builders. Extracting an
`add_local_terms()` helper would further reduce each lattice builder.
Alternatively, unify the three Wannier90 builder functions, or examine
the remaining long functions (`_check_mod_para_mvmc` 63 lines,
`vector_potential` 86 lines).

---

### Step 64 — Extract `add_local_terms()` helper from lattice builders
**Date**: 2026-01-28
**Files changed**:
- `python/lattice/interaction_builder.py` (added function, updated docstring)
- `python/lattice/chain_lattice.py` (updated imports, replaced local block)
- `python/lattice/square_lattice.py` (updated imports, replaced local block)
- `python/lattice/triangular_lattice.py` (updated imports, replaced local block)
- `python/lattice/honeycomb_lattice.py` (updated imports, replaced local block)
- `python/lattice/kagome.py` (updated imports, replaced local block, also
  converted one remaining unconverted neighbor block from Step 63)
- `python/lattice/ladder.py` (updated imports, replaced local block)
- `test/unit/test_interaction_builder.py` (added import, 8 new tests)

**What changed**:
The on-site interaction block was repeated in all 6 lattice builders.
Each instance contained: for SPIN models, `mag_field()` + `general_j(D)`
for anisotropy; for electron models, `hubbard_local()` plus optional
Kondo coupling (`general_j(J)` + `mag_field()` on the localized spin).

Added `add_local_terms(StdI, isite, jsite_kondo)` to
`interaction_builder.py` which encapsulates this model-dependent dispatch.
Each lattice builder now replaces the 5–8 line block with a single call.

For honeycomb (NsiteUC=2) and kagome (NsiteUC=3), the previously
unrolled per-site blocks were converted to a loop over `NsiteUC`.

Also fixed kagome.py: one neighbor block ("Nearest neighbor intra cell
0 → 1") was not converted in Step 63 — converted it to use
`add_neighbor_interaction()` and removed the now-unused `set_label`
import.

**New tests** (8 tests in `test/unit/test_interaction_builder.py`):
- `TestAddLocalTerms`: spin adds mag_field, spin adds anisotropy,
  hubbard adds local terms, hubbard no Kondo coupling, kondo adds
  J-coupling, kondo adds localized mag_field, spin no Cintra terms,
  spin transverse field

**Test results**:
- Unit tests: 1199 passed (+8 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: The lattice builders are now significantly
slimmer. Consider unifying the Wannier90 builder functions
(`wannier90_read_lattice_data`, `wannier90_read_hr`, `wannier90`)
which share similar patterns, or examine the remaining long functions
(`_check_mod_para_mvmc` 63 lines, `vector_potential` 86 lines).
Alternatively, refactor the boost functions (`chain_boost`,
`honeycomb_boost`, `kagome_boost`, `ladder_boost`) which share
similar structure.

### Step 65 — Add `add_neighbor_interaction_3d()` for 3D lattice builders (2026-01-28)

**Target files**: `interaction_builder.py`, `orthorhombic.py`,
`fc_ortho.py`, `pyrochlore.py`

**What changed and why**: Steps 63–64 extracted `add_neighbor_interaction()`
(2D, using `set_label`) and `add_local_terms()` for the six 2D lattice
builders. The three 3D lattice builders (`orthorhombic`, `fc_ortho`,
`pyrochlore`) had the same repetitive pattern — calling `find_site` then
dispatching on model type to call `general_j` or `hopping`+`coulomb` —
but using `find_site` instead of `set_label`.

Added `add_neighbor_interaction_3d()` in `interaction_builder.py` as the
3D counterpart of `add_neighbor_interaction()`. It combines `find_site` +
model-dependent dispatch (spin → `general_j`; electron → `hopping` +
`coulomb`).

Applied it to all three 3D lattice builders:
- **orthorhombic.py**: Replaced 13 neighbor blocks + local-term block
  with `add_neighbor_interaction_3d()` + `add_local_terms()` calls.
  Removed direct imports of `mag_field`, `general_j`, `hubbard_local`,
  `hopping`, `coulomb` and `find_site`.
- **fc_ortho.py**: Replaced 9 neighbor blocks with
  `add_neighbor_interaction_3d()` calls. Kept local-term block manual
  because the C original intentionally omits `mag_field` on the Kondo
  `jsite` (unlike orthorhombic.py), so `add_local_terms()` would change
  behavior. Removed `hopping`, `coulomb`, `find_site` imports; kept
  `mag_field`, `general_j`, `hubbard_local`.
- **pyrochlore.py**: Replaced 12 neighbor blocks with
  `add_neighbor_interaction_3d()` calls. Kept local-term block manual
  because it loops over NsiteUC=4 sub-sites with a hard-coded `isite + 3`
  offset for the Kondo J-coupling. Removed `hopping`, `coulomb`,
  `find_site` imports; kept `mag_field`, `general_j`, `hubbard_local`.

**Line count changes**:
- orthorhombic.py: 217 → ~100 lines (removed ~117 lines)
- fc_ortho.py: 252 → ~198 lines (removed ~54 lines)
- pyrochlore.py: 286 → ~220 lines (removed ~66 lines)

**New tests** (7 tests in `test/unit/test_interaction_builder.py`):
- `TestAddNeighborInteraction3D`: returns site indices, spin model calls
  general_j, hubbard model calls hopping, hubbard model calls coulomb,
  spin model no hopping, hubbard model no general_j, neighbor along height

**Test results**:
- Unit tests: 1206 passed (+7 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: The lattice builders are now compact. Consider
refactoring the boost functions (`chain_boost`, `honeycomb_boost`,
`kagome_boost`, `ladder_boost`) which share similar patterns, or examine
remaining long functions (`_check_mod_para_mvmc`, `vector_potential`).
The Wannier90 builder functions also share common patterns that could
benefit from unification.

### Step 66 — Extract `write_boost_6spin_star()` / `write_boost_6spin_pair()` and compact boost data (2026-01-28)

**Target files**: `site_util.py`, `chain_lattice.py`, `ladder.py`,
`honeycomb_lattice.py`, `kagome.py`

**What changed and why**: All four boost functions (`chain_boost`,
`ladder_boost`, `honeycomb_boost`, `kagome_boost`) contained identical
6-spin write loops — one for `list_6spin_star` (~6 lines) and one for
`list_6spin_pair` (~7 lines). These were duplicated across all four files.

Additionally, `chain_boost` and `ladder_boost` used element-by-element
assignment for their `list_6spin_star` arrays (8+ lines per function),
while `honeycomb_boost` and `kagome_boost` already used compact slice
notation.

Changes made:

1. **Extracted** `write_boost_6spin_star()` and `write_boost_6spin_pair()`
   into `site_util.py` alongside existing `write_boost_mag_field()`,
   `write_boost_j_full()`, `write_boost_j_symmetric()`.

2. **Replaced** the duplicated write loops in all four boost functions
   with calls to the new helpers.

3. **Converted** `chain_boost`'s element-by-element `list_6spin_star`
   assignment (8 lines per pivot) to `[ipivot, :] = [...]` slice notation
   (1 line per pivot).

4. **Converted** `chain_boost`'s element-by-element `list_6spin_pair`
   assignment (66 lines, 7 assignments × 8 interactions) to slice notation
   (8 lines per pivot).

5. **Converted** `ladder_boost`'s element-by-element `list_6spin_star`
   assignment to slice notation.

**Line count changes**:
- chain_lattice.py: 357 → 280 lines (−77 lines)
- ladder.py: 316 → 294 lines (−22 lines)
- honeycomb_lattice.py: 320 → 308 lines (−12 lines)
- kagome.py: 332 → 320 lines (−12 lines)
- site_util.py: 658 → 700 lines (+42 lines for helpers + docstrings)
- Net reduction: ~81 lines

**New tests** (8 tests in `test/unit/test_site_util.py`):
- `TestWriteBoost6spinStar`: writes header, writes pivot labels, writes
  correct values, num_pivot_lines
- `TestWriteBoost6spinPair`: writes header, writes correct interaction
  count, writes correct values, multiple pivots

**Test results**:
- Unit tests: 1214 passed (+8 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Consider examining the remaining long functions
(`_check_mod_para_mvmc` in `mvmc_writer.py`, `vector_potential` in
`hphi_writer.py`), or unifying the Wannier90 builder functions. The
boost functions are now compact and use shared helpers for all I/O.

---

### Step 67 — Extract `compute_max_interactions()` for section-4 upper-limit computation
**Date**: 2026-01-28
**Target**: `lattice/interaction_builder.py`, 7 lattice builder files

**What changed and why**:
Each lattice builder has a "section 4" block that computes `ntransMax` and
`nintrMax` upper limits for array allocation.  Seven of the nine builders
(square, triangular, honeycomb, kagome, orthorhombic, fc_ortho, pyrochlore)
use an identical formula parameterised only by the total number of neighbor
bond types (`n_bonds`).  Extracted a shared helper
`compute_max_interactions(StdI, n_bonds)` into `interaction_builder.py` and
replaced the 8–12 line section-4 blocks in all seven files with a 2–3 line
call.

Chain and ladder were excluded because they use different base variables
(`StdI.L` instead of `StdI.nsite`/`StdI.NCell`) and have different Kondo
formulas.  Wannier90 was excluded (completely different structure).

**Changes**:

1. **Added** `compute_max_interactions(StdI, n_bonds)` to
   `lattice/interaction_builder.py` (lines 468–504).  Formula:
   - SPIN: `ntransMax = nsite*(S2+1+2*S2)`,
     `nintrMax = NCell*(NsiteUC+n_bonds)*(3*S2+1)^2`
   - ELEC: `ntransMax = NCell*2*(2*NsiteUC+2*n_bonds)`,
     `nintrMax = NCell*(NsiteUC+4*n_bonds)`
   - KONDO: adds `nsite//2 * (S2+1+2*S2)` and `nsite//2 * (3*S2+1)^2`

2. **Updated** 7 lattice files — added `compute_max_interactions` to imports
   and replaced section-4 blocks:
   - square_lattice.py: `n_bonds=2+2+2` (6)
   - triangular_lattice.py: `n_bonds=3+3+3` (9)
   - honeycomb_lattice.py: `n_bonds=3+6+3` (12)
   - kagome.py: `n_bonds=6+6` (12)
   - orthorhombic.py: `n_bonds=3+6+4` (13)
   - fc_ortho.py: `n_bonds=6+3` (9)
   - pyrochlore.py: `n_bonds=12`

3. **Updated** module docstring in `interaction_builder.py` to list the new
   function.

**New tests** (8 tests in `test/unit/test_interaction_builder.py`):
- `TestComputeMaxInteractions`:
  - `test_spin_ntrans`: verifies SPIN ntransMax formula
  - `test_spin_nintr`: verifies SPIN nintrMax formula
  - `test_spin_s1`: verifies S=1 (S2=2) formulas
  - `test_hubbard_ntrans`: verifies HUBBARD ntransMax formula
  - `test_hubbard_nintr`: verifies HUBBARD nintrMax formula
  - `test_kondo_adds_spin_terms`: verifies KONDO extra terms
  - `test_multi_site_uc`: verifies NsiteUC > 1 (pyrochlore-like)
  - `test_returns_tuple`: verifies return type

**Test results**:
- Unit tests: 1222 passed (+8 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Consider extracting the remaining duplicated
section-4 blocks in chain_lattice.py and ladder.py into their own helpers,
or target other cross-cutting patterns such as the Kondo local-term
variations across fc_ortho.py and pyrochlore.py.

### 2026-01-28 — Step 68: Remove dead file-handle pattern in 3D lattice files; add `close_lattice_xsf`

**Target**: `lattice/site_util.py`, `lattice/orthorhombic.py`,
`lattice/fc_ortho.py`, `lattice/pyrochlore.py`

**Problem**: The three 3D lattice files (orthorhombic, fc_ortho, pyrochlore)
all opened `lattice.xsf` for writing, passed the file handle `fp` to
`init_site(StdI, fp, 3)`, which **ignores `fp` when `dim=3`**, then closed
the file.  Immediately afterward, `print_xsf(StdI)` re-opened
`lattice.xsf` from scratch and wrote the actual content.  The initial
`open`/`close` was dead code — the file was created empty and instantly
overwritten.

The 2D lattice files already had a clean pattern: `open_lattice_gp()` /
`close_lattice_gp()`.  No 3D equivalent existed.

**Changes**:

1. **Added** `close_lattice_xsf(StdI)` to `lattice/site_util.py` — the 3D
   counterpart of `close_lattice_gp()`.  Calls `print_xsf(StdI)` and
   `print_geometry(StdI)`.

2. **Updated** 3 lattice files:
   - Removed `fp = open("lattice.xsf", "w")` (dead code)
   - Changed `init_site(StdI, fp, 3)` → `init_site(StdI, None, 3)`
   - Replaced `fp.close(); print_xsf(StdI); print_geometry(StdI)` →
     `close_lattice_xsf(StdI)`
   - Updated imports: added `close_lattice_xsf` from `.site_util`, removed
     now-unused `print_geometry, print_xsf` from `.geometry_output`

3. **Updated** module docstring in `site_util.py` to list `close_lattice_xsf`.

**New tests** (4 tests in `test/unit/test_site_util.py`):
- `TestCloseLatticeXsf`:
  - `test_creates_lattice_xsf`: verifies `lattice.xsf` is created
  - `test_creates_geometry_dat`: verifies `geometry.dat` is created
  - `test_lattice_xsf_contains_crystal`: verifies XCrySDen header
  - `test_geometry_dat_has_content`: verifies non-empty output

**Test results**:
- Unit tests: 1226 passed (+4 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Consider further cleanup of shared patterns across
lattice files, or target large domain-specific modules (wannier90.py,
export_wannier90.py, common_writer.py) for internal decomposition.

---

### 2026-01-28 — Step 69: Extract boost-output helpers into `boost_output.py`

**Target**: `python/lattice/site_util.py` (719 lines → 587 lines)

**What changed**: Extracted 5 `write_boost_*` functions (~132 lines) from
`site_util.py` into a new `lattice/boost_output.py` module.  These functions
write sections of HPhi's `boost.def` file and are logically separate from the
site initialisation and folding utilities in `site_util.py`.

**Rationale**: `site_util.py` mixed two unrelated concerns — super-cell site
management (init_site, find_site, set_label, etc.) and boost-output file
formatting.  Separating them improves cohesion and makes both modules easier
to navigate.

1. **Created** `python/lattice/boost_output.py` — new module with:
   - `write_boost_mag_field` — magnetic-field header line
   - `write_boost_j_full` — full 3×3 J-coupling matrix
   - `write_boost_j_symmetric` — upper-triangle symmetric J-coupling matrix
   - `write_boost_6spin_star` — list_6spin_star array
   - `write_boost_6spin_pair` — list_6spin_pair array

2. **Removed** 5 functions (~132 lines) from `site_util.py` and updated its
   module docstring to no longer list the boost functions.

3. **Updated imports** in 4 lattice files:
   - `chain_lattice.py` — split boost imports to `from .boost_output import ...`
   - `honeycomb_lattice.py` — same
   - `kagome.py` — same
   - `ladder.py` — same

4. **Created** `test/unit/test_boost_output.py` with 18 tests in 5 classes:
   - `TestWriteBoostMagField` (2 tests)
   - `TestWriteBoostJFull` (4 tests)
   - `TestWriteBoostJSymmetric` (4 tests)
   - `TestWriteBoost6spinStar` (4 tests)
   - `TestWriteBoost6spinPair` (4 tests)

**Test results**:
- Unit tests: 1226 passed (18 moved from test_site_util.py to test_boost_output.py)
- Integration tests: 83/83 passed

**Suggested next step**: Consider internal decomposition of large domain-specific
modules (wannier90.py at 1356 lines, common_writer.py at 1314 lines) or further
simplification of the lattice dispatch/builder pattern.

---

### 2026-01-28 — Step 70: Remove backward-compatibility thin wrappers from `common_writer.py`

**Target**: `python/writer/common_writer.py` (1314 lines → 1126 lines)

**What changed**: Removed 7 module-level thin wrapper functions that duplicated
`GreenFunctionIndices` class methods.  These wrappers (`_spin_max`,
`_skip_local_spin_pair`, `_kondo_site`, `_green1_indices_corr`,
`_green1_indices_raw`, `_green2_indices_corr`, `_green2_indices_raw`) were only
called by tests, never by production code — the actual callers (`print_1_green`,
`print_2_green`) already use the class directly.

**Rationale**: The wrappers were backward-compatibility shims added when the
`GreenFunctionIndices` class was introduced.  With the class now stable and
fully tested, the wrappers are dead code that bloats the module and creates
a maintenance burden (every class method change requires updating a parallel
wrapper function).

1. **Removed** 7 wrapper functions (~186 lines) from `common_writer.py`.

2. **Updated** `GreenFunctionIndices` class docstring to remove mention of
   the now-deleted wrappers.

3. **Updated** `test/unit/test_common_writer.py`:
   - Removed imports of the 7 deleted wrappers.
   - Removed 7 wrapper-based test classes (24 tests total):
     `TestSpinMax`, `TestSkipLocalSpinPair`, `TestKondoSite`,
     `TestGreen1IndicesCorr`, `TestGreen1IndicesRaw`,
     `TestGreen2IndicesCorr`, `TestGreen2IndicesRaw`.
   - Removed 4 "matches wrapper" tests from the class-based test classes.
   - Merged unique coverage from the deleted wrapper tests into the
     class-based test classes, adding 8 new tests:
     - `test_spin_max_local_spin_half`, `test_spin_max_local_spin_one`
     - `test_skip_local_spin_pair_one_itinerant`
     - `test_kondo_site_single_site_uc`
     - `test_green1_corr_2site_hubbard_basic`, `test_green1_corr_local_spin_pair_skipped`
     - `test_green1_raw_2site_itinerant`, `test_green1_raw_1site_itinerant`,
       `test_green1_raw_local_spin_pair_skipped`
     - `test_green2_corr_kondo_doubles_range`
     - `test_green2_raw_1site_itinerant`, `test_green2_raw_local_spin_filtering`

**Test results**:
- Unit tests: 1210 passed (net -16: removed 24 wrapper tests + 4 "matches wrapper"
  tests, added 12 new class-based tests)
- Integration tests: 83/83 passed

**Suggested next step**: Consider internal decomposition of the largest remaining
modules (wannier90.py at 1356 lines, common_writer.py at 1126 lines) or further
cleanup of the lattice dispatch/builder pattern.

---

## Step 71 — Convert bare `open()/close()` to `with` statements (2026-01-28)

**Target files**: `chain_lattice.py`, `ladder.py`, `honeycomb_lattice.py`, `kagome.py`,
`wannier90.py`, `stdface_main.py`, `export_wannier90.py`

**What**: Converted all bare `fp = open(...); ...; fp.close()` patterns to
idiomatic Python `with open(...) as fp:` context-manager blocks across 7 files
(10 instances total). This is a Phase 3 refactoring (Python idioms).

**Why**: The C-style `open()/close()` pattern leaks file handles if an exception
occurs between the open and close. The `with` statement guarantees cleanup
regardless of how the block exits. This is a standard Python best practice.

**Changes**:

1. **`python/lattice/chain_lattice.py`** — `chain_boost()`: wrapped `boost.def`
   write in `with` block.
2. **`python/lattice/ladder.py`** — `ladder_boost()`: wrapped `boost.def` write
   in `with` block.
3. **`python/lattice/honeycomb_lattice.py`** — `honeycomb_boost()`: wrapped
   `boost.def` write in `with` block.
4. **`python/lattice/kagome.py`** — `kagome_boost()`: wrapped `boost.def` write
   in `with` block.
5. **`python/lattice/wannier90.py`** — 4 instances:
   - `_geometry_w90()`: wrapped geometry file read in `with` block (renamed `fp`
     → `fp_geom`).
   - `_read_w90()`: wrapped `_hr.dat`/`_ur.dat`/`_jr.dat` read in `with` block
     (renamed `fp` → `fp_hr`).
   - `_read_density_matrix()`: wrapped `_dr.dat` read in `with` block (renamed
     `fp` → `fp_dr`).
   - `wannier90()`: wrapped `lattice.xsf` write in narrow `with` block around
     `init_site()` call; removed orphaned `fp_xsf.close()`.
6. **`python/stdface_main.py`** — `stdface_main()`: wrapped input file parse loop
   in `with` block (renamed `fp` → `fp_in`).
7. **`python/writer/export_wannier90.py`** — 2 instances:
   - `_write_geometry()`: wrapped geometry write in `with` block (renamed `fp` →
     `fp_out`).
   - `_write_wannier90()`: wrapped Wannier90 format write in `with` block
     (renamed `fp` → `fp_out`).

**Not changed**: `site_util.py`'s `close_lattice_gp()` — the file handle is
opened in `open_lattice_gp()` and passed through callers before being closed
in `close_lattice_gp()`. Converting this requires changing the caller pattern
across all 2D lattice files, which is a separate refactoring step.

**Test results**:
- Unit tests: 1210 passed (no change)
- Integration tests: 83/83 passed

**Suggested next step**: Convert the `open_lattice_gp()`/`close_lattice_gp()`
cross-function file handle pattern to use `contextlib.contextmanager` or
restructure the 2D lattice builders to use `with` statements directly.

---

## Step 72 — Replace `open_lattice_gp`/`close_lattice_gp` with `lattice_gp` context manager (2026-01-28)

**Target files**: `site_util.py`, `chain_lattice.py`, `ladder.py`,
`square_lattice.py`, `triangular_lattice.py`, `honeycomb_lattice.py`,
`kagome.py`, `test_site_util.py`

**What**: Replaced the two-function `open_lattice_gp()`/`close_lattice_gp()`
pattern with a single `lattice_gp()` context manager using
`contextlib.contextmanager`. Updated all 6 caller lattice files and the
corresponding unit tests.

**Why**: The previous pattern required callers to manually pair
`open_lattice_gp()`/`close_lattice_gp()` calls, which is error-prone (the
footer write and `print_geometry()` call could be skipped if an exception
occurs between them). A context manager ensures cleanup runs unconditionally,
and the `with` statement makes the intent clearer. This completes the
Phase 3 `open()/close()` → `with` conversion started in Step 71.

**Changes**:

1. **`python/lattice/site_util.py`**:
   - Added `from contextlib import contextmanager` and `Iterator` to typing
     imports.
   - Replaced `open_lattice_gp()` + `close_lattice_gp()` (2 functions, ~40
     lines) with `lattice_gp()` context manager (~30 lines) using
     `@contextmanager` decorator.
   - Updated module docstring to list `lattice_gp` instead of the old
     functions.
   - Updated `close_lattice_xsf` docstring cross-reference.

2. **6 lattice caller files** (chain, ladder, square, triangular, honeycomb,
   kagome):
   - Changed imports: `open_lattice_gp, close_lattice_gp` → `lattice_gp`.
   - Changed function bodies: `fp = open_lattice_gp(StdI) ... close_lattice_gp(fp, StdI)` → `with lattice_gp(StdI) as fp:` wrapping the
     entire lattice setup body.

3. **`test/unit/test_site_util.py`**:
   - Replaced `TestOpenLatticeGp` (5 tests) + `TestCloseLatticeGp` (4 tests)
     with `TestLatticeGp` (10 tests) covering: yields file for HPhi/mVMC,
     yields None for HWAVE, yields file for HWAVE+lattice_gp, creates file,
     writes footer on exit, None fp no error, calls print_geometry, footer
     constant content.
   - Updated import from `open_lattice_gp, close_lattice_gp` → `lattice_gp`.

**Test results**:
- Unit tests: 1210 passed (net +1 from 9→10 tests in this class)
- Integration tests: 83/83 passed

**Suggested next step**: The codebase now has zero bare `fp.close()` calls.
Potential next targets: replace remaining `TextIO | None` file-handle
pass-through patterns (e.g., in `init_site` and `set_label`), or move to
other Phase 3 idiom improvements such as replacing sentinel `NaN_i` checks
with `Optional` types in selected hot paths.

---

## Step 73 — Extract `_det_and_cofactor()` from duplicated reciprocal-box computation (2026-01-28)

**Target files**: `lattice/site_util.py`, `writer/mvmc_variational.py`,
`test/unit/test_site_util.py`

**What**: Extracted the 3x3 integer determinant + cofactor matrix computation
into a pure function `_det_and_cofactor(box)` in `site_util.py`. Both
`_compute_reciprocal_box()` (site_util.py) and `_init_site_sub()`
(mvmc_variational.py) now call this shared helper instead of duplicating
the 15-line Sarrus-rule + cofactor algorithm inline.

**Why**: The identical algorithm (determinant via Sarrus rule, cofactor via
minors, sign-flip for negative determinant) was duplicated between
`_compute_reciprocal_box()` for the main lattice and `_init_site_sub()` for
the mVMC sub-lattice. Extracting it into a pure function that takes a box
matrix and returns `(det, cofactor)`:

- Eliminates code duplication (~15 lines removed from each call site)
- Makes the math independently testable (no side effects)
- Improves readability — callers now express *intent* ("compute det & cofactor")
  rather than *mechanism* (nested loops with index arithmetic)

**Changes**:

1. **`python/lattice/site_util.py`**:
   - Added `_det_and_cofactor(box: np.ndarray) -> tuple[int, np.ndarray]`:
     pure function computing determinant and cofactor matrix of a 3x3 integer
     matrix, with sign-flip for negative determinant.
   - Simplified `_compute_reciprocal_box()`: replaced ~20 lines of inline
     determinant + cofactor computation with a single call to
     `_det_and_cofactor(StdI.box)`.

2. **`python/writer/mvmc_variational.py`**:
   - Added `_det_and_cofactor` to imports from `lattice.site_util`.
   - Simplified `_init_site_sub()`: replaced ~18 lines of inline determinant +
     cofactor computation with a single call to
     `_det_and_cofactor(StdI.boxsub)`.

3. **`test/unit/test_site_util.py`**:
   - Added `_det_and_cofactor` to imports.
   - Added `TestDetAndCofactor` class with 7 tests: identity matrix, diagonal
     matrix, negative determinant flip, zero determinant, 1D chain, return
     type verification, off-diagonal box.

**Test results**:
- Unit tests: 1217 passed (+7 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Extract the LWH-vs-box conflict validation pattern
that is duplicated between `init_site()` (L/W/Height vs box) and
`_init_site_sub()` (Lsub/Wsub/Hsub vs boxsub), or look at other
deduplication opportunities in `mvmc_variational.py`.

---

## Step 74 — Extract `_validate_box_params()` from duplicated LWH-vs-box validation

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`, `python/writer/mvmc_variational.py`

**Motivation**: Both `init_site()` and `_init_site_sub()` contain a
~30-line block that: (1) checks whether L/W/Height and box entries conflict,
(2) if L/W/Height specified, sets the box to `diag(W, L, Height)`, or
(3) if box entries specified (or neither), fills box entries from defaults
using `print_val_i`.  The two blocks differ only in field names (suffix
`""` vs `"sub"`) and default values (identity vs parent box).  Extracting
a parameterised helper eliminates the duplication.

**Changes**:

1. **`python/lattice/site_util.py`**:
   - Added `_validate_box_params(L, W, Height, box, suffix="", defaults=None)`
     helper before `_det_and_cofactor`.  Handles all three branches
     (conflict → exit, LWH → diagonal box, box/defaults → fill entries).
     The `suffix` parameter controls label names (`"Height"` vs `"Hsub"`)
     and the `defaults` parameter supplies fallback values (identity by
     default, parent box for sub-lattice).
   - Simplified `init_site()`: replaced ~25 lines of inline LWH-vs-box
     validation with a single call:
     `StdI.L, StdI.W, StdI.Height = _validate_box_params(StdI.L, StdI.W, StdI.Height, StdI.box)`

2. **`python/writer/mvmc_variational.py`**:
   - Added `_validate_box_params` to imports from `lattice.site_util`.
   - Simplified `_init_site_sub()`: replaced ~25 lines of inline validation
     with a single call:
     `StdI.Lsub, StdI.Wsub, StdI.Hsub = _validate_box_params(StdI.Lsub, StdI.Wsub, StdI.Hsub, StdI.boxsub, suffix="sub", defaults=StdI.box)`

3. **`test/unit/test_site_util.py`**:
   - Added `_validate_box_params` and `NaN_i` to imports.
   - Updated module docstring.
   - Added `TestValidateBoxParams` class with 10 tests: LWH fills diagonal,
     partial LWH defaults to 1, box with identity defaults, box with custom
     defaults, neither specified fills from defaults, conflict raises
     SystemExit, conflict with suffix, suffix affects labels, no-suffix
     height label, box modified in-place.

**Test results**:
- Unit tests: 1227 passed (+10 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Look for further deduplication opportunities in
`mvmc_variational.py` (e.g. the remaining `_init_site_sub` logic that
mirrors `init_site`'s reciprocal-box and cell-enumeration steps), or
continue with other Phase 3 idiom improvements.

---

## Step 75 — Extract `_fold_to_cell()` from duplicated site-folding logic

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`, `python/writer/mvmc_variational.py`

**Motivation**: `_fold_site()` in `site_util.py` and `_fold_site_sub()` in
`mvmc_variational.py` contained structurally identical ~20-line algorithms
for folding a site coordinate into a cell.  They differed only in which
fields they used: `_fold_site` used `(rbox, NCell, box)` while
`_fold_site_sub` used `(rboxsub, NCellsub, boxsub)`.  Extracting a
parameterised pure function eliminates the duplication.

**Changes**:

1. **`python/lattice/site_util.py`**:
   - Added `_fold_to_cell(rbox, ncell, box, iCellV)` pure function
     before `_fold_site`.  Takes the reciprocal matrix, cell count, box
     matrix, and coordinate as explicit parameters.  Contains the
     complete three-step folding algorithm (fractional transform, periodic
     image search, back-transform).
   - Simplified `_fold_site()`: replaced ~20-line inline algorithm with a
     single delegation: `return _fold_to_cell(StdI.rbox, StdI.NCell, StdI.box, iCellV)`.

2. **`python/writer/mvmc_variational.py`**:
   - Added `_fold_to_cell` to imports from `lattice.site_util`.
   - Simplified `_fold_site_sub()`: replaced ~20-line inline algorithm
     with a single delegation:
     `return _fold_to_cell(StdI.rboxsub, StdI.NCellsub, StdI.boxsub, iCellV)`.

3. **`test/unit/test_site_util.py`**:
   - Added `_fold_to_cell` to imports.
   - Updated module docstring.
   - Added `TestFoldToCell` class with 7 tests: identity cell inside,
     wraps at boundary, negative wraps, 2D square cell, matches
     `_fold_site`, sub-cell parameters, origin coordinate.

**Test results**:
- Unit tests: 1234 passed (+7 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Extract anti-periodic dot-product helper from
duplicated `sum(AntiPeriod[k] * nBox[k])` pattern in `proj()` and
`generate_orb()`, or look at other deduplication opportunities across
the writer modules.

---

## Step 76 — Extract `_anti_period_dot()` and `_parity_sign()` helpers

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/mvmc_variational.py`

**Motivation**: The anti-periodic boundary dot product
`AntiPeriod[0]*nBox[0] + AntiPeriod[1]*nBox[1] + AntiPeriod[2]*nBox[2]`
appeared in three places across `proj()` (lines 91-93 and 100-102) and
`generate_orb()` (line 246).  The parity-to-sign conversion
`1 if x % 2 == 0 else -1` appeared at two sites (lines 116 and 247).
Extracting named helpers improves readability and eliminates duplication.

**Changes**:

1. **`python/writer/mvmc_variational.py`**:
   - Added `_anti_period_dot(AntiPeriod, nBox) -> int`: computes the dot
     product of the anti-period flag vector and the super-cell image index.
   - Added `_parity_sign(value) -> int`: converts an integer's parity
     to +1 (even) or -1 (odd).
   - Simplified `proj()`: replaced the two 3-line anti-period dot products
     with a single `ap_dot = _anti_period_dot(StdI.AntiPeriod, nBox)`
     computed once per jCell.  Replaced the inline parity conversion with
     `a = _parity_sign(Anti[iSym][jsite])`.
   - Simplified `generate_orb()`: replaced inline dot-product + parity
     conversion with `_parity_sign(_anti_period_dot(...))`.

2. **`test/unit/test_mvmc_variational.py`**:
   - Added `_anti_period_dot` and `_parity_sign` to imports.
   - Updated module docstring.
   - Added `TestAntiPeriodDot` class with 5 tests: all-zero, single
     direction, all-ones, mixed signs, large nBox.
   - Added `TestParitySign` class with 6 tests: zero, one, even positive,
     odd positive, negative even, negative odd.

**Test results**:
- Unit tests: 1245 passed (+11 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Pythonize `_has_anti_period()` in `mvmc_writer.py`
to use `any()`, or look for other Phase 3 opportunities such as replacing
C-style 3-element loops with numpy operations or extracting the
commensurate-check logic from `_init_site_sub()`.

---

## Step 77 — Extract `_cell_vector()` for Cell-row extraction

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`, `python/writer/mvmc_variational.py`

**Motivation**: The pattern `[int(StdI.Cell[idx, k]) for k in range(3)]`
(extracting a cell-coordinate row as an integer list) appeared 15+ times
across `mvmc_variational.py`, with each occurrence manually converting
from float to int.  Extracting a named helper centralises the conversion,
gives the operation a semantic name, and enables cleaner cell-arithmetic
using `zip()` list comprehensions instead of `for k in range(3)` index
loops.

**Changes**:

1. **`python/lattice/site_util.py`**:
   - Added `_cell_vector(Cell, idx) -> list[int]`: extracts row `idx`
     from a `(NCell, 3)` float array as `[int, int, int]`.
   - Simplified `_find_cell_index()`: replaced 3-line element-wise
     equality check with `_cell_vector(StdI.Cell, k) == cellV`.

2. **`python/writer/mvmc_variational.py`**:
   - Added `_cell_vector` to imports from `lattice.site_util`.
   - Simplified `proj()`:
     - Extracted `iCV = _cell_vector(StdI.Cell, iCell)` once per outer loop.
     - Replaced inline equality check (3 lines) with `iCellV == iCV`.
     - Used `zip()` for jCell vector construction.
   - Simplified `generate_orb()`:
     - Extracted `iCV` and `jCV` per loop iteration.
     - Replaced inline list comprehensions with `zip()`-based arithmetic.
   - Simplified `_jastrow_global_optimization()`:
     - Extracted `dCV = _cell_vector(StdI.Cell, dCell)` once per outer loop.
     - Replaced 3-line `Cell[dCell, 0] == 0 and ...` check with
       `dCV == [0, 0, 0]`.
     - Replaced individual `int(StdI.Cell[...])` in `find_site()` calls
       with `dCV[i]` and `iCV[i]` indexing.

3. **`test/unit/test_site_util.py`**:
   - Added `_cell_vector` to imports.
   - Updated module docstring.
   - Added `TestCellVector` class with 7 tests: basic extraction, returns
     int, negative values, zero vector, returns list, equality comparison,
     used in find_cell_index.

**Test results**:
- Unit tests: 1252 passed (+7 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Apply `_cell_vector` to `export_wannier90.py`
(2 remaining call sites), Pythonize `_has_anti_period()` to use `any()`,
or look for other Phase 3 opportunities such as extracting the
commensurate-check from `_init_site_sub()`.

---

## Step 78 — Pythonize `_has_anti_period()` using `any()`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/mvmc_writer.py`

**Motivation**: The `_has_anti_period()` function used a C-style chained
`or` comparison checking `StdI.AntiPeriod[0] == 1 or StdI.AntiPeriod[1] == 1
or StdI.AntiPeriod[2] == 1`.  This pattern is more idiomatically expressed
in Python using `any()` with a generator expression.

**Changes**:

1. **`python/writer/mvmc_writer.py`**:
   - Simplified `_has_anti_period()`: replaced 3-way chained `or` comparison
     with `any(ap == 1 for ap in StdI.AntiPeriod)`.
   - The explicit `== 1` check is preserved (not just truthiness) because
     the function specifically detects the anti-periodic flag value of 1,
     as validated by the existing test `test_non_one_value_returns_false`.

**Test results**:
- Unit tests: 1252 passed (no new tests needed — 6 existing tests fully cover the function)
- Integration tests: 83/83 passed

**Suggested next step**: Apply `_cell_vector` to `export_wannier90.py`
(2 remaining call sites), extract the commensurate-check from
`_init_site_sub()`, or continue with other Phase 3 idiom improvements.

---

## Step 79 — Apply `_cell_vector` to `export_wannier90.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/export_wannier90.py`

**Motivation**: The file contained two occurrences of the C-style pattern
`[int(StdI.Cell[jcell, i] - StdI.Cell[icell, i]) for i in range(3)]`
for computing cell-coordinate differences.  This pattern duplicates
the float-to-int conversion that `_cell_vector()` centralises, and
uses index iteration instead of `zip()`-based subtraction.

**Changes**:

1. **`python/writer/export_wannier90.py`**:
   - Added import: `from lattice.site_util import _cell_vector`
   - Replaced cell-difference pattern at line ~590 (in `_export_trans`):
     ```python
     jCV = _cell_vector(StdI.Cell, jcell)
     iCV = _cell_vector(StdI.Cell, icell)
     rr = [j - i for j, i in zip(jCV, iCV)]
     ```
   - Replaced identical pattern at line ~748 (in `_export_coulomb`):
     same refactoring.

**Test results**:
- Unit tests: 1252 passed (no new tests needed — existing integration coverage is sufficient)
- Integration tests: 83/83 passed

**Suggested next step**: Apply similar `_cell_vector` usage to
`geometry_output.py` (4 remaining sites with Cell differences), extract
the commensurate-check from `_init_site_sub()`, or continue with other
Phase 3 idiom improvements.

---

## Step 80 — Extract `_cell_diff()` helper in `geometry_output.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/geometry_output.py`

**Motivation**: The file contained repetitive cell-coordinate difference
calculations written as 3-line f-string interpolations accessing
`StdI.Cell[iCell, k] - StdI.Cell[0, k]` for k=0,1,2.  Extracting a local
`_cell_diff()` helper centralises this pattern and improves readability.

**Note**: Initially attempted to import `_cell_vector` from `site_util`,
but this created a circular import (`site_util → geometry_output →
site_util`).  The solution was to create a local `_cell_diff()` helper
that encapsulates the cell difference calculation directly.

**Changes**:

1. **`python/lattice/geometry_output.py`**:
   - Added `_cell_diff(Cell, iCell, jCell) -> list[int]` helper function
     that computes the integer difference between two cell coordinate rows.
   - Simplified cell-output loops in `print_geometry()`: replaced 4 inline
     3-element f-string interpolations with calls to `_cell_diff()`.

2. **`test/unit/test_geometry_output.py`**:
   - Added import of `_cell_diff`.
   - Added `TestCellDiff` class with 5 tests: basic difference, returns
     integers, negative difference, same cell (zero), float-to-int conversion.

**Test results**:
- Unit tests: 1257 passed (+5 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Extract the commensurate-check from
`_init_site_sub()` in `mvmc_variational.py`, or look for other Phase 3
idiom improvements such as replacing remaining C-style index loops with
`enumerate()` or `zip()`.

---

## Step 81 — Extract `_check_commensurate()` from `_init_site_sub()`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/mvmc_variational.py`

**Motivation**: The `_init_site_sub()` function contained a 9-line nested
triple loop that checked whether the sublattice is commensurate with the
main lattice.  Extracting this into a pure helper function:
1. Improves readability by giving the algorithm a semantic name
2. Makes the logic unit-testable independently
3. Replaces a C-style accumulator pattern with a generator expression

**Changes**:

1. **`python/writer/mvmc_variational.py`**:
   - Added `_check_commensurate(rbox_sub, box, ncell_sub) -> bool` helper
     that returns `True` if the sublattice is commensurate.  Uses a
     generator expression `sum(... for kk in range(3))` instead of the
     C-style `prod = 0; for kk: prod +=` pattern.
   - Simplified `_init_site_sub()`: replaced 9-line triple nested loop
     with a single call to `_check_commensurate()`.

2. **`test/unit/test_mvmc_variational.py`**:
   - Added `_check_commensurate` to imports.
   - Added `TestCheckCommensurate` class with 5 tests: identity is
     commensurate, 2x2 in 4x4 commensurate, 3x3 in 4x4 incommensurate,
     same box commensurate, mismatched dimensions incommensurate.

**Test results**:
- Unit tests: 1262 passed (+5 new tests)
- Integration tests: 83/83 passed

**Suggested next step**: Look for other Phase 3 idiom improvements such as
replacing remaining C-style index loops with `enumerate()` or `zip()`,
or replacing `for i in range(3): ... [i]` patterns with direct iteration.

---

## Step 82 — Replace C-style index loops with direct iteration in `geometry_output.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/geometry_output.py`

**Motivation**: The file contained several C-style `for ii in range(3):`
loops that indexed arrays with `array[ii, ...]`.  Python allows iterating
directly over rows of numpy arrays and using `enumerate()` when both the
index and value are needed.

**Changes**:

1. **`print_geometry()` — direct vector lines**:
   - Before: `for ii in range(3): fp.write(f"{StdI.direct[ii, 0]:25.15e} ...")`
   - After: `for row in StdI.direct: fp.write(f"{row[0]:25.15e} ...")`

2. **`print_geometry()` — box lines**:
   - Before: `for ii in range(3): fp.write(f"{int(StdI.box[ii, 0])} ...")`
   - After: `for row in StdI.box: fp.write(f"{int(row[0])} ...")`

3. **`print_xsf()` — CONVVEC section**:
   - Before: `for ii in range(3): row = [0,0,0]; for jj: if ii==jj: row[jj] = length[ii]`
   - After: `for ii, length_val in enumerate(StdI.length): row = [0,0,0]; row[ii] = length_val`
   - This eliminates the inner conditional loop by directly setting the diagonal element.

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing tests cover output format)
- Integration tests: 83/83 passed

**Suggested next step**: Apply similar direct iteration patterns to
`print_xsf()` PRIMVEC section (matrix multiplication with nested loops),
or look at other files like `export_wannier90.py` or `wannier90.py` for
`range(3)` loop simplifications.

---

## Step 83 — Replace nested loops with numpy matrix multiplication in `print_xsf()`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/geometry_output.py`

**Motivation**: The `print_xsf()` function contained two triple-nested
loops that computed matrix-vector products.  These C-style patterns:
```python
for jj in range(3):
    for kk in range(3):
        vec[jj] += A[kk] * B[kk, jj]
```
can be replaced with numpy's `@` operator for clearer, more efficient code.

**Changes**:

1. **Added `import numpy as np`** at module level.

2. **PRIMVEC section** (was 6 lines, now 4 lines):
   - Before: Triple-nested loop computing `vec[jj] += box[ii, kk] * direct[kk, jj]`
   - After: `primvec = StdI.box @ StdI.direct; for vec in primvec: ...`

3. **PRIMCOORD section** (was 7 lines, now 4 lines):
   - Before: Double-nested loop computing `vec[jj] += (Cell + tau)[kk] * direct[kk, jj]`
   - After: `frac_coord = Cell[iCell, :] + tau[isite, :]; vec = frac_coord @ direct`

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing tests cover output format)
- Integration tests: 83/83 passed

**Suggested next step**: Apply similar numpy matrix operations to other
files like `export_wannier90.py` or `wannier90.py` that have nested
loops computing matrix products, or look for other `for i in range(3)`
patterns to simplify.

---

## Step 84 — Replace C-style loops with numpy/Pythonic patterns in `export_wannier90.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/export_wannier90.py`

**Motivation**: The file contained several C-style nested loops that could
be replaced with numpy matrix operations and Pythonic iteration patterns.

**Changes**:

1. **`_unfold_site()` — matrix-vector products** (was 13 lines, now 10 lines):
   - Before: `for i in range(3): for j in range(3): v[i] += rbox[i,j] * v_in[j]`
   - After: `v = (StdI.rbox.astype(int) @ np.array(v_in)).tolist()`
   - Similarly for the second matrix-vector product: `w = (v @ box) // NCell`

2. **`_write_geometry()` — direct row iteration** (was 8 lines, now 6 lines):
   - Before: `for ii in range(3): fp.write(f"{direct[ii, 0]} ...")`
   - After: `for row in StdI.direct: fp.write(f"{row[0]} ...")`
   - Similarly for tau array: `for tau_row in StdI.tau[:NsiteUC]:`

3. **`_write_interaction()` — min/max finding** (was 10 lines, now 9 lines):
   - Before: `rmin = [table[0].r[i] for i in range(3)]; for k: for i: ...`
   - After: `rmin = list(table[0].r); for item in table[1:]: for i, r in enumerate(item.r):`
   - Used `zip()` for final range computation: `[max(abs(lo), abs(hi)) for lo, hi in zip(rmin, rmax)]`

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing integration coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Apply similar numpy/Pythonic patterns to
`wannier90.py` (4+ occurrences of `for i in range(3)` loops), or continue
with other files.

---

## Step 85 — Replace C-style loops with numpy operations in `wannier90.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/wannier90.py`

**Motivation**: The file contained several C-style nested loops that could
be replaced with numpy vectorized operations and Pythonic patterns.

**Changes**:

1. **`_check_in_box()` — matrix-vector product** (was 9 lines, now 4 lines):
   - Before: `for i in range(3): for j in range(3): judge_vec[i] += rvec[j] * inv[j, i]`
   - After: `judge_vec = rvec @ inverse_matrix`
   - Before: `abs(judge_vec[0]) <= 1 and abs(judge_vec[1]) <= 1 and ...`
   - After: `np.all(np.abs(judge_vec) <= 1)`

2. **`_geometry_w90()` — direct row iteration for printing**:
   - Before: `for ii in range(3): print(f"... {StdI.direct[ii, 0]} ...")`
   - After: `for row in StdI.direct: print(f"... {row[0]} ...")`
   - Similarly for tau array with slicing: `for tau_row in StdI.tau[:StdI.NsiteUC]:`

3. **`_apply_boundary_weights()` — vectorized max and masking** (was 17 lines, now 11 lines):
   - Before: nested loop to find max: `for iWSC: for ii: if abs(...) > Band_lattice[ii]:`
   - After: `Band_lattice = np.max(np.abs(indx_tot[:nWSC]), axis=0).astype(int)`
   - Before: nested loop for weight update: `for iWSC: if abs(indx_tot[iWSC, ii]) == ...`
   - After: `mask = np.abs(indx_tot[:nWSC, ii]) == Model_lattice[ii]; Weight_tot[:nWSC][mask] *= 0.5`

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Continue with other `range(3)` loop patterns in
`wannier90.py` (lines 377-378, 490, 744, 750-751, 834), or move to other
files like `honeycomb_lattice.py` or `param_check.py`.

---

## Step 86 — Replace element-wise comparisons with numpy vectorized ops in `wannier90.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/wannier90.py`

**Motivation**: The file contained several patterns comparing array elements
one-by-one with `and`/`or` chains that could be replaced with `np.all()`
and `np.any()`, plus nested loops for array assignment that could use slicing.

**Changes**:

1. **Cutoff check** (was 4 lines, now 2 lines):
   - Before: `if abs(indx[0]) > cutoff_R[0] or abs(indx[1]) > ... or abs(indx[2]) > ...:`
   - After: `if np.any(np.abs(indx_tot[iWSC]) > cutoff_R):`

2. **Inversion symmetry check** (was 6 lines, now 2 lines):
   - Before: `if indx[iWSC, 0] == -indx[jWSC, 0] and ... and ...:`
   - After: `if np.all(indx_tot[iWSC] == -indx_tot[jWSC]):`
   - Before: Nested loop `for iWan: for jWan: Mat[iWSC, iWan, jWan] = 0`
   - After: `Mat_tot[iWSC, :, :] = 0.0`

3. **Origin check** (was 4 lines, now 2 lines):
   - Before: `if indx[0] == 0 and indx[1] == 0 and indx[2] == 0:`
   - After: `if np.all(indx_tot[iWSC] == 0):`
   - Before: `for iWan: for jWan in range(iWan): Mat[iWan, jWan] = 0`
   - After: `Mat_tot[iWSC, iWan, :iWan] = 0.0`

4. **DenMat population** (was 4 lines, now 2 lines):
   - Before: `key = (int(indx[0]), int(indx[1]), int(indx[2]))`
   - After: `key = tuple(indx_tot[iWSC].astype(int))`
   - Before: Nested loop `for iWan: for jWan: DenMat[key][iWan,jWan] = Mat[...]`
   - After: `DenMat[key][:, :] = Mat_tot[iWSC, :nWan, :nWan]`

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Continue with remaining `range(3)` patterns in
`wannier90.py` (lines 732-736, 738-744), or move to other files like
`honeycomb_lattice.py` or `param_check.py`.

---

## Step 87 — Replace nested matrix threshold check with `np.any()` in boost validation

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/honeycomb_lattice.py`, `python/lattice/kagome.py`

**Motivation**: Both files contained identical nested loops checking if any
element of the 3x3 `Jp` matrix exceeds a threshold. This is a classic use
case for `np.any()` with `np.abs()`.

**Changes**:

1. **`honeycomb_lattice.py` — `_write_boost_honey()`** (was 4 lines, now 2 lines):
   - Before: `for i1 in range(3): for i2 in range(3): if abs(StdI.Jp[i1, i2]) > 1.0e-8:`
   - After: `if np.any(np.abs(StdI.Jp) > 1.0e-8):`

2. **`kagome.py` — `_write_boost_kagome()`** (was 4 lines, now 2 lines):
   - Same pattern replaced with same numpy vectorized check.

**Note**: Similar patterns in `input_params.py` and `param_check.py` were
not changed because they require element-specific handling (printing
different suffixes for each element, applying different defaults based on
diagonal/off-diagonal position).

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for other simple loop patterns that can be
vectorized, or consider Phase 3 opportunities like replacing string
constants with enums, or reviewing the remaining large functions for
splitting opportunities.

---

## Step 88 — Replace `range(len(...))` with direct iteration in `export_wannier90.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/export_wannier90.py`

**Motivation**: The file contained 3 occurrences of the anti-pattern
`for j in range(len(intr_table))` followed by `intr_table[j].field`.
Python's direct iteration over lists is cleaner and more idiomatic.

**Changes**:

1. **`_export_trans()` search loop** (lines ~591):
   - Before: `for j in range(len(intr_table)): if intr_table[j].r == rr ...`
   - After: `for item in intr_table: if item.r == rr ...`

2. **`_export_coulomb()` search loop** (lines ~748):
   - Before: `for j in range(len(intr_table)): if intr_table[j].r == rr and ... and intr_table[j].s == ispin ...`
   - After: `for item in intr_table: if item.r == rr and ... and item.s == ispin ...`

3. **`_export_coulomb_intra_unique()` search loop** (lines ~858):
   - Before: `for j in range(len(intr_table)): if intr_table[j].a == isite ...`
   - After: `for item in intr_table: if item.a == isite ...`

**Benefits**:
- More readable: `item.v` is clearer than `intr_table[j].v`
- More Pythonic: direct iteration over sequences is preferred
- Less error-prone: no off-by-one risks with explicit indices

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for other `range(len(...))` patterns or
continue with Phase 3 idiom improvements in other files.

---

## Step 89 — Vectorize periodic folding loop in `_unfold_site()`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/export_wannier90.py`

**Motivation**: The `_unfold_site()` function had a loop folding each
coordinate independently to the [-N/2, N/2] range.  This can be done
in a single vectorized operation using `np.where()`.

**Changes**:

**`_unfold_site()` — periodic folding** (was 6 lines, now 4 lines):
- Before:
  ```python
  v = (StdI.rbox.astype(int) @ np.array(v_in)).tolist()
  for i in range(3):
      vv = v[i] / StdI.NCell
      if vv > 0.5:
          v[i] -= StdI.NCell
      elif vv <= -0.5:
          v[i] += StdI.NCell
  w = (np.array(v) @ StdI.box.astype(int)) // StdI.NCell
  ```
- After:
  ```python
  v = StdI.rbox.astype(int) @ np.array(v_in)
  vv = v / StdI.NCell
  v = np.where(vv > 0.5, v - StdI.NCell, v)
  v = np.where(vv <= -0.5, v + StdI.NCell, v)
  w = (v @ StdI.box.astype(int)) // StdI.NCell
  ```

**Benefits**:
- No explicit loop needed
- `v` stays as numpy array throughout (no `.tolist()` conversion)
- Clearer separation between computation and folding steps

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for similar element-wise conditional patterns
that can use `np.where()`, or continue with other Phase 3 improvements.

---

## Step 90 — Vectorize Model_lattice computation in `wannier90.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/wannier90.py`

**Motivation**: The Model_lattice computation in `wannier90_geometry()`
used a list comprehension with conditional logic. This can be expressed
more cleanly using `np.where()` for element-wise conditionals.

**Changes**:

**`wannier90_geometry()` — Model_lattice initialization**:
- Before:
  ```python
  dims = [StdI.W, StdI.L, StdI.Height]
  Model_lattice = np.array([d // 2 if d % 2 == 0 else 0 for d in dims], dtype=int)
  ```
- After:
  ```python
  dims = np.array([StdI.W, StdI.L, StdI.Height], dtype=int)
  Model_lattice = np.where(dims % 2 == 0, dims // 2, 0)
  ```

**Benefits**:
- Uses numpy's vectorized conditional instead of Python list comprehension
- `dims` is a numpy array from the start, enabling direct element-wise ops
- Cleaner expression of the conditional logic

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Continue looking for element-wise conditional
patterns or list comprehensions that can be vectorized with `np.where()`.
Other Phase 3 opportunities include replacing string constants with enums
or reviewing remaining large functions for splitting.

---

## Step 91 — Vectorize `_fold_to_cell()` in `site_util.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The `_fold_to_cell()` function had four separate C-style
loop blocks performing matrix-vector multiplication and element-wise
floor division. These can all be replaced with numpy operations.

**Changes**:

**`_fold_to_cell()` — coordinate folding** (was 18 lines, now 10 lines):
- Before:
  ```python
  iCellV_frac = [0, 0, 0]
  for ii in range(3):
      for jj in range(3):
          iCellV_frac[ii] += rbox[ii, jj] * iCellV[jj]

  nBox = [0, 0, 0]
  for ii in range(3):
      nBox[ii] = (iCellV_frac[ii] + ncell * 1000) // ncell - 1000

  for ii in range(3):
      iCellV_frac[ii] -= ncell * nBox[ii]

  iCellV_fold = [0, 0, 0]
  for ii in range(3):
      for jj in range(3):
          iCellV_fold[ii] += box[jj, ii] * iCellV_frac[jj]
      iCellV_fold[ii] = (iCellV_fold[ii] + ncell * 1000) // ncell - 1000

  return nBox, iCellV_fold
  ```
- After:
  ```python
  iCellV_arr = np.asarray(iCellV)
  iCellV_frac = rbox @ iCellV_arr

  nBox = (iCellV_frac + ncell * 1000) // ncell - 1000

  iCellV_frac = iCellV_frac - ncell * nBox

  iCellV_fold = (box.T @ iCellV_frac + ncell * 1000) // ncell - 1000

  return nBox.astype(int).tolist(), iCellV_fold.astype(int).tolist()
  ```

**Benefits**:
- Four nested `for` loops replaced with numpy matrix operations
- Clear separation of each transformation step
- Uses `@` operator for matrix-vector multiplication
- Uses `box.T` for transpose instead of explicit indexing
- Preserves the return type (list of int) for backward compatibility

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for similar nested loop patterns in
`site_util.py` (there are several more at lines 332, 361, 370, 424, 477, 539)
that can be vectorized with numpy matrix operations.

---

## Step 92 — Vectorize `_det_and_cofactor()` in `site_util.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The `_det_and_cofactor()` function had C-style loops for
computing the 3x3 determinant (Sarrus rule) and cofactor matrix. These
can be replaced with numpy's `np.linalg.det` and vectorized indexing.

**Changes**:

**`_det_and_cofactor()` — determinant and cofactor** (was 15 lines, now 12 lines):
- Before:
  ```python
  det = 0
  for ii in range(3):
      det += (int(box[0, ii])
              * int(box[1, (ii + 1) % 3])
              * int(box[2, (ii + 2) % 3])
              - int(box[0, ii])
              * int(box[1, (ii + 2) % 3])
              * int(box[2, (ii + 1) % 3]))

  cofactor = np.zeros((3, 3), dtype=float)
  for ii in range(3):
      for jj in range(3):
          cofactor[ii, jj] = (int(box[(ii + 1) % 3, (jj + 1) % 3])
                              * int(box[(ii + 2) % 3, (jj + 2) % 3])
                              - int(box[(ii + 1) % 3, (jj + 2) % 3])
                              * int(box[(ii + 2) % 3, (jj + 1) % 3]))
  ```
- After:
  ```python
  det = int(round(np.linalg.det(box.astype(float))))

  # Vectorized cofactor using index arrays
  idx = np.array([1, 2, 0])   # (i+1) % 3
  idx2 = np.array([2, 0, 1])  # (i+2) % 3
  cofactor = (box[idx][:, idx] * box[idx2][:, idx2]
              - box[idx][:, idx2] * box[idx2][:, idx])
  ```

**Benefits**:
- Uses numpy's optimized `np.linalg.det` instead of manual Sarrus rule
- Cofactor computation uses vectorized fancy indexing instead of nested loops
- Cleaner and more readable
- `idx` and `idx2` arrays encode the cyclic permutation pattern

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Continue vectorizing loops in `site_util.py`.
Check remaining patterns at lines 424, 477, 539 for similar opportunities.

---

## Step 93 — Vectorize ExpPhase/AntiPeriod computation in `init_site()`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The `init_site()` function had a C-style loop computing
`ExpPhase[ii] = exp(i * pi180 * phase[ii])` and conditionally setting
`AntiPeriod[ii]` based on whether `ExpPhase[ii] ≈ -1`. This can be
vectorized using numpy's `np.exp` and `np.where`.

**Changes**:

**`init_site()` — phase factor computation** (was 7 lines, now 2 lines):
- Before:
  ```python
  for ii in range(3):
      StdI.ExpPhase[ii] = (math.cos(StdI.pi180 * StdI.phase[ii])
                           + 1j * math.sin(StdI.pi180 * StdI.phase[ii]))
      if abs(StdI.ExpPhase[ii] + 1.0) < AMPLITUDE_EPS:
          StdI.AntiPeriod[ii] = 1
      else:
          StdI.AntiPeriod[ii] = 0
  ```
- After:
  ```python
  StdI.ExpPhase = np.exp(1j * StdI.pi180 * StdI.phase)
  StdI.AntiPeriod = np.where(np.abs(StdI.ExpPhase + 1.0) < AMPLITUDE_EPS, 1, 0)
  ```

**Benefits**:
- Uses `np.exp(1j * ...)` instead of explicit `cos + i*sin` (Euler's formula)
- Uses `np.where` for vectorized conditional assignment
- Reduces 7 lines to 2 lines
- More readable and idiomatic numpy code

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Continue looking for C-style loops in `site_util.py`
or other files that can be vectorized with numpy operations.

---

## Step 94 — Vectorize `find_site()` in `site_util.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The `find_site()` function had C-style patterns for
computing the distance vector `dR` (element-by-element assignment) and
the boundary phase `Cphase` (loop with multiplication). Both can be
vectorized with numpy operations.

**Changes**:

**`find_site()` — dR and Cphase computation** (was 7 lines, now 4 lines):
- Before:
  ```python
  dR = np.zeros(3)
  dR[0] = -float(diW) + StdI.tau[isiteUC, 0] - StdI.tau[jsiteUC, 0]
  dR[1] = -float(diL) + StdI.tau[isiteUC, 1] - StdI.tau[jsiteUC, 1]
  dR[2] = -float(diH) + StdI.tau[isiteUC, 2] - StdI.tau[jsiteUC, 2]

  jCellV = [iW + diW, iL + diL, iH + diH]
  nBox, jCellV = _fold_site(StdI, jCellV)
  Cphase = 1.0 + 0j
  for ii in range(3):
      Cphase *= StdI.ExpPhase[ii] ** nBox[ii]
  ```
- After:
  ```python
  di = np.array([diW, diL, diH], dtype=float)
  dR = -di + StdI.tau[isiteUC, :] - StdI.tau[jsiteUC, :]

  jCellV = [iW + diW, iL + diL, iH + diH]
  nBox, jCellV = _fold_site(StdI, jCellV)
  Cphase = np.prod(StdI.ExpPhase ** np.array(nBox))
  ```

**Benefits**:
- `dR` computed with single vectorized expression instead of 4 lines
- `Cphase` uses `np.prod(a ** b)` instead of explicit loop
- More concise and idiomatic numpy code

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Continue looking for C-style loops in other files
(e.g., `boost_output.py`, `wannier90.py`, `hphi_writer.py`) that can be
vectorized with numpy operations.

---

## Step 95 — Replace `for i in range(3)` with direct row iteration in `boost_output.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/boost_output.py`

**Motivation**: The `write_boost_j_full()` function had a C-style
`for i in range(3)` loop indexing into a numpy array. This can be
replaced with direct iteration over rows, which is more Pythonic.

**Changes**:

**`write_boost_j_full()` — matrix row output** (was 5 lines, now 3 lines):
- Before:
  ```python
  for i in range(3):
      fp.write(
          f"{scale * J[i, 0]:25.15e} "
          f"{scale * J[i, 1]:25.15e} "
          f"{scale * J[i, 2]:25.15e}\n"
      )
  ```
- After:
  ```python
  for row in J:
      scaled = scale * row
      fp.write(f"{scaled[0]:25.15e} {scaled[1]:25.15e} {scaled[2]:25.15e}\n")
  ```

**Benefits**:
- Direct iteration over numpy array rows instead of index-based access
- Pre-computes `scaled` array for cleaner formatting
- More Pythonic and readable

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for similar `for i in range(n)` patterns in
`wannier90.py` or `hphi_writer.py` that iterate over array rows and can
use direct iteration instead.

---

## Step 96 — Replace dot product loop and trig with numpy in `interaction_builder.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/interaction_builder.py`

**Motivation**: The `_general_hopping()` function had a C-style loop
computing a dot product and explicit `cos/sin` for a complex exponential.
Both can be replaced with numpy operations.

**Changes**:

**`_general_hopping()` — phase computation** (was 4 lines, now 2 lines):
- Before:
  ```python
  for it in range(StdI.Lanczos_max):
      Cphase = 0.0
      for ii in range(3):
          Cphase += StdI.At[it][ii] * dR[ii]
      coef = math.cos(Cphase) + 1j * math.sin(-Cphase)
  ```
- After:
  ```python
  for it in range(StdI.Lanczos_max):
      Cphase = np.dot(StdI.At[it], dR)
      coef = np.exp(-1j * Cphase)
  ```

**Benefits**:
- `np.dot()` replaces explicit element-wise loop for dot product
- `np.exp(-1j * Cphase)` replaces `cos + i*sin(-...)` (Euler's formula)
- Cleaner, more mathematical expression of the phase factor
- Reduces 4 lines to 2 lines

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for other `math.cos/sin` pairs that compute
complex exponentials and can use `np.exp(1j * ...)` instead. Also check
for remaining explicit dot product loops.

---

## Step 97 — Replace sum-comprehension dot product and simplify list comparison in `_enumerate_cells()`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The `_enumerate_cells()` function had:
1. A generator expression `sum(nBox[jj] * int(StdI.box[jj, ii]) for jj in range(3))`
   that computes a dot product - can use numpy `@` operator
2. A verbose check `nBox[0] == 0 and nBox[1] == 0 and nBox[2] == 0`
   that can be simplified to `nBox == [0, 0, 0]`

**Changes**:

**`_enumerate_cells()` — bounding box and cell enumeration**:
- Before:
  ```python
  for ii in range(3):
      for n2 in range(2):
          for n1 in range(2):
              for n0 in range(2):
                  nBox = [n0, n1, n2]
                  edge = sum(nBox[jj] * int(StdI.box[jj, ii]) for jj in range(3))
                  ...
  ...
  if nBox[0] == 0 and nBox[1] == 0 and nBox[2] == 0:
  ```
- After:
  ```python
  box_int = StdI.box.astype(int)
  for ii in range(3):
      for n2 in range(2):
          for n1 in range(2):
              for n0 in range(2):
                  nBox = np.array([n0, n1, n2])
                  edge = nBox @ box_int[:, ii]
                  ...
  ...
  if nBox == [0, 0, 0]:
  ```

**Benefits**:
- `nBox @ box_int[:, ii]` is cleaner than explicit sum-comprehension
- `box_int` conversion hoisted outside nested loops for efficiency
- `nBox == [0, 0, 0]` is more readable than three separate comparisons

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for remaining sum-comprehension patterns that
compute dot products, or consider Phase 3 opportunities in other areas like
replacing string constants with enums.

---

## Step 98 — Vectorize `_check_commensurate()` in `mvmc_variational.py`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/mvmc_variational.py`

**Motivation**: The `_check_commensurate()` function had nested loops with a
sum-comprehension computing dot products. Since it checks all 9 elements of a
3x3 product matrix, it can be fully vectorized using numpy matrix multiplication.

**Changes**:

**`_check_commensurate()` — commensurate check** (was 6 lines, now 3 lines):
- Before:
  ```python
  for ii in range(3):
      for jj in range(3):
          prod = sum(int(rbox_sub[ii, kk]) * int(box[jj, kk]) for kk in range(3))
          if prod % ncell_sub != 0:
              return False
  return True
  ```
- After:
  ```python
  # Compute all dot products: prod[i,j] = rbox_sub[i,:] · box[j,:]
  prod = rbox_sub.astype(int) @ box.astype(int).T
  return bool(np.all(prod % ncell_sub == 0))
  ```

**Benefits**:
- Nested loops replaced with single matrix multiplication
- `np.all()` checks divisibility across entire matrix at once
- `bool()` wrapper ensures Python `True/False` (not numpy `np.True_`)
- Much more concise and mathematically clear

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Most C-style loop patterns have been addressed.
Consider reviewing Phase 3 opportunities like adding enums for remaining
string constants, or look at Phase 2 opportunities for introducing classes.

---

## Step 99 — Replace if/elif chains with dict dispatch in `_check_conserved_quantities()`

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/common_writer.py`

**Motivation**: The `_check_conserved_quantities()` function had two if/elif
chains dispatching on action strings ("required", "not_used", "default_0").
This is a classic case for dict dispatch, making the code more extensible
and reducing conditional branches.

**Changes**:

Added two dispatch dicts:
```python
_NCOND_ACTION_DISPATCH: dict[str, callable] = {
    "required": required_val_i,
    "not_used": not_used_i,
}

_SZ2_ACTION_DISPATCH: dict[str, callable] = {
    "required": lambda StdI: required_val_i("2Sz", StdI.Sz2),
    "not_used": lambda StdI: not_used_i("2Sz", StdI.Sz2),
    "default_0": lambda StdI: setattr(StdI, 'Sz2', print_val_i("2Sz", StdI.Sz2, 0)),
}
```

**`_check_conserved_quantities()` — action dispatch**:
- Before:
  ```python
  if ncond_action == "required":
      required_val_i(ncond_label, StdI.ncond)
  elif ncond_action == "not_used":
      not_used_i(ncond_label, StdI.ncond)
  ...
  if sz2_action == "required":
      required_val_i("2Sz", StdI.Sz2)
  elif sz2_action == "not_used":
      not_used_i("2Sz", StdI.Sz2)
  elif sz2_action == "default_0":
      StdI.Sz2 = print_val_i("2Sz", StdI.Sz2, 0)
  ```
- After:
  ```python
  if ncond_action is not None:
      _NCOND_ACTION_DISPATCH[ncond_action](ncond_label, StdI.ncond)
  ...
  if sz2_action is not None:
      _SZ2_ACTION_DISPATCH[sz2_action](StdI)
  ```

**Benefits**:
- If/elif chains replaced with single dict lookups
- Adding new actions only requires adding dict entries
- Separates action definitions from dispatch logic
- Consistent with other dispatch patterns in the codebase

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for remaining if/elif chains that dispatch on
string values, or consider Phase 2 class introduction opportunities.

---

## Step 100 — Consolidate `SPIN_SUFFIXES` constant to eliminate duplication

**Date**: 2026-01-28
**Phase**: 3 — Leverage Python idioms
**Target**: `python/param_check.py`, `python/lattice/input_params.py`

**Motivation**: The 3x3 spin-interaction suffix matrix (["x","xy","xz"], etc.)
was duplicated: defined locally in `not_used_j()` in `param_check.py` and as
`_SUFFIXES` in `input_params.py`. Consolidating to a single source of truth
improves maintainability.

**Changes**:

1. **`param_check.py`**: Added module-level constant `SPIN_SUFFIXES`:
   ```python
   SPIN_SUFFIXES: list[list[str]] = [
       ["x", "xy", "xz"],
       ["yx", "y", "yz"],
       ["zx", "zy", "z"],
   ]
   ```
   Updated `not_used_j()` to use `SPIN_SUFFIXES` instead of local variable.

2. **`input_params.py`**: Removed local `_SUFFIXES` definition, now imports:
   ```python
   from param_check import exit_program, SPIN_SUFFIXES
   _SUFFIXES = SPIN_SUFFIXES  # Alias for backward compatibility
   ```

**Benefits**:
- Single source of truth for spin-interaction suffixes
- `param_check.py` is the canonical location (low-level utility module)
- Backward compatible: `_SUFFIXES` alias preserved in `input_params.py`
- Reduces maintenance burden when suffix naming conventions change

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Reached Step 100! Consider reviewing the overall
refactoring progress and identifying remaining high-impact opportunities,
or shift focus to Phase 2 class introduction (e.g., SolverWriter hierarchy).

---

## Step 101 — Vectorize Fourier coefficient computation in `hphi_writer.py`

**Date**: 2026-01-29
**Phase**: 3 — Leverage Python idioms
**Target**: `python/writer/hphi_writer.py`

**Motivation**: The `_compute_fourier_coefficients()` function had a
nested C-style loop computing `exp(2πi q·r)` for each site individually
using `math.cos/sin` and an explicit dot product sum. This is a natural
candidate for numpy vectorization using broadcasting.

**Changes**:

1. **Added `import numpy as np`** to `hphi_writer.py` (was not imported).

2. **`_compute_fourier_coefficients()` — vectorized** (was 12 lines, now 11 lines):
   - Before:
     ```python
     isite = 0
     for icell in range(StdI.NCell):
         for itau in range(StdI.NsiteUC):
             Cphase = ((StdI.Cell[icell][0] + StdI.tau[itau][0]) * StdI.SpectrumQ[0]
                       + (StdI.Cell[icell][1] + StdI.tau[itau][1]) * StdI.SpectrumQ[1]
                       + (StdI.Cell[icell][2] + StdI.tau[itau][2]) * StdI.SpectrumQ[2])
             fourier_r[isite] = math.cos(2.0 * StdI.pi * Cphase)
             fourier_i[isite] = math.sin(2.0 * StdI.pi * Cphase)
             isite += 1
     ```
   - After:
     ```python
     if n_computed > 0:
         # Broadcasting: (NCell,1,3) + (1,NsiteUC,3) -> (NCell,NsiteUC,3)
         positions = StdI.Cell[:NCell].astype(float)[:, np.newaxis, :] + \
                     StdI.tau[:NsiteUC][np.newaxis, :, :]
         Cphase_flat = (2.0 * StdI.pi * (positions @ StdI.SpectrumQ)).ravel()
         fourier_r[:n_computed] = np.cos(Cphase_flat)
         fourier_i[:n_computed] = np.sin(Cphase_flat)
     ```

3. **Kondo duplication** also simplified from explicit loop to array slicing:
   `fourier_r[half:] = fourier_r[:half]`

**Benefits**:
- Nested loops replaced with numpy broadcasting and matrix multiplication
- `np.cos/np.sin` applied to entire array at once
- Handles `NCell == 0` gracefully (guard for `Cell is None`)
- ~10x faster for large lattices due to vectorized operations

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for other nested loop patterns in
`hphi_writer.py` that can be vectorized, or consider Phase 2 class
introduction opportunities.

---

## Step 102 — Vectorize bounding box computation in `_enumerate_cells()`

**Date**: 2026-01-29
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The bounding box computation in `_enumerate_cells()` used
four nested loops (`for ii, n2, n1, n0`) to check all 8 cube corners against
each dimension. This can be fully vectorized using numpy matrix multiplication
and `min/max` along axes.

**Changes**:

1. **Added `import itertools`** to `site_util.py`.

2. **`_enumerate_cells()` — bounding box** (was 10 lines, now 6 lines):
   - Before:
     ```python
     bound = [[0, 0], [0, 0], [0, 0]]
     box_int = StdI.box.astype(int)
     for ii in range(3):
         for n2 in range(2):
             for n1 in range(2):
                 for n0 in range(2):
                     nBox = np.array([n0, n1, n2])
                     edge = nBox @ box_int[:, ii]
                     if edge < bound[ii][0]:
                         bound[ii][0] = edge
                     if edge > bound[ii][1]:
                         bound[ii][1] = edge
     ```
   - After:
     ```python
     box_int = StdI.box.astype(int)
     corners = np.array(list(itertools.product(range(2), repeat=3)))
     edges = corners @ box_int  # shape (8, 3)
     bound = list(zip(edges.min(axis=0).tolist(), edges.max(axis=0).tolist()))
     ```

**Benefits**:
- Four nested loops replaced with single matrix multiplication
- `min/max` along axis replaces manual tracking of min/max
- `itertools.product` generates cube corners cleanly
- Much more concise and mathematically clear

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Consider replacing the cell enumeration triple loop
below (lines 421-428) with `itertools.product`, or look for other nested
loop patterns to simplify.

---

## Step 103 — Replace cell enumeration triple loop with `itertools.product`

**Date**: 2026-01-29
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The cell enumeration in `_enumerate_cells()` had three nested
`for` loops iterating over the bounding box ranges. This is a natural fit
for `itertools.product`, reducing nesting depth from 3 to 1.

**Changes**:

**`_enumerate_cells()` — cell enumeration** (reduced nesting):
- Before:
  ```python
  for ic2 in range(bound[2][0], bound[2][1] + 1):
      for ic1 in range(bound[1][0], bound[1][1] + 1):
          for ic0 in range(bound[0][0], bound[0][1] + 1):
              iCellV = [ic0, ic1, ic2]
              ...
  ```
- After:
  ```python
  for ic2, ic1, ic0 in itertools.product(
      range(bound[2][0], bound[2][1] + 1),
      range(bound[1][0], bound[1][1] + 1),
      range(bound[0][0], bound[0][1] + 1),
  ):
      iCellV = [ic0, ic1, ic2]
      ...
  ```

**Important**: The iteration order `ic2, ic1, ic0` (outermost to innermost)
must match the original C code to preserve the cell enumeration order, which
affects output file ordering. An initial attempt with `ic0, ic1, ic2` caused
33 integration test failures.

**Benefits**:
- Reduces nesting depth from 3 levels to 1
- More Pythonic — `itertools.product` is the standard tool for Cartesian products
- Preserves exact iteration order

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Look for other triple-nested loops that can use
`itertools.product`, or consider other Phase 3 opportunities.

---

## Step 104 — Replace element-wise matrix product with numpy `@` in `_write_gnuplot_header()`

**Date**: 2026-01-29
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The `_write_gnuplot_header()` function computed supercell
corner positions using 4 lines of element-by-element matrix multiplication.
This is a standard 2x2 matrix product that can be expressed as `box @ direct`.

**Changes**:

**`_write_gnuplot_header()` — corner positions** (was 5 lines, now 3 lines):
- Before:
  ```python
  pos = np.zeros((4, 2))
  pos[1, 0] = StdI.direct[0, 0] * StdI.box[0, 0] + StdI.direct[1, 0] * StdI.box[0, 1]
  pos[1, 1] = StdI.direct[0, 1] * StdI.box[0, 0] + StdI.direct[1, 1] * StdI.box[0, 1]
  pos[2, 0] = StdI.direct[0, 0] * StdI.box[1, 0] + StdI.direct[1, 0] * StdI.box[1, 1]
  pos[2, 1] = StdI.direct[0, 1] * StdI.box[1, 0] + StdI.direct[1, 1] * StdI.box[1, 1]
  pos[3, :] = pos[1, :] + pos[2, :]
  ```
- After:
  ```python
  pos = np.zeros((4, 2))
  pos[1:3, :] = StdI.box[:2, :2] @ StdI.direct[:2, :2]
  pos[3, :] = pos[1, :] + pos[2, :]
  ```

**Benefits**:
- 4 lines of scalar arithmetic replaced with 1 matrix multiplication
- Mathematically clear: corner positions = box vectors × direct lattice
- More concise and less error-prone

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize the 2D position computations in
`set_label()` (lines 612-635 of site_util.py), which compute xi, yi, xj, yj
using element-wise dot products with `direct`.

---

## Step 105 — Vectorize 2D position computations in `set_label()`

**Date**: 2026-01-29
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The `set_label()` function computed 2D site positions
(`xi, yi, xj, yj`) using 8 lines of element-by-element dot products with
the `direct` lattice matrix. Each pair `(x, y) = frac @ direct[:2, :2]`
is a standard matrix-vector product.

**Changes**:

**`set_label()` — 2D position computation** (was 16 lines, now 8 lines):
- Before (each block):
  ```python
  xi = (StdI.direct[0, 0] * (iW + StdI.tau[jsiteUC, 0])
        + StdI.direct[1, 0] * (iL + StdI.tau[jsiteUC, 1]))
  yi = (StdI.direct[0, 1] * (iW + StdI.tau[jsiteUC, 0])
        + StdI.direct[1, 1] * (iL + StdI.tau[jsiteUC, 1]))
  xj = ...  # similar 2 lines
  yj = ...  # similar 2 lines
  ```
- After (each block):
  ```python
  D = StdI.direct[:2, :2]
  frac_i = np.array([iW + StdI.tau[jsiteUC, 0], iL + StdI.tau[jsiteUC, 1]])
  frac_j = np.array([iW - diW + StdI.tau[isiteUC, 0], iL - diL + StdI.tau[isiteUC, 1]])
  xi, yi = frac_i @ D
  xj, yj = frac_j @ D
  ```

**Benefits**:
- 8 lines of scalar arithmetic per block replaced with 3 lines
- Uses numpy `@` operator for matrix-vector product
- `D` matrix hoisted and reused across both blocks
- More mathematically clear: position = fractional_coords × direct_lattice

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Most element-by-element matrix operations in
`site_util.py` have been vectorized. Consider looking at other files
for similar patterns, or review Phase 2 class introduction opportunities.

---

## Step 106 — Replace element-wise direct-matrix assignments with numpy slicing in `init_site()`

**Date**: 2026-01-29
**Phase**: 3 — Leverage Python idioms
**Target**: `python/lattice/site_util.py`

**Motivation**: The `init_site()` function had 5 individual scalar
assignments to zero out z-components and set the third lattice vector
for 2D lattices. These can be expressed more concisely with numpy slicing.

**Changes**:

**`init_site()` — 2D direct matrix setup** (was 5 lines, now 2 lines):
- Before:
  ```python
  StdI.direct[0, 2] = 0.0
  StdI.direct[1, 2] = 0.0
  StdI.direct[2, 0] = 0.0
  StdI.direct[2, 1] = 0.0
  StdI.direct[2, 2] = 1.0
  ```
- After:
  ```python
  StdI.direct[:2, 2] = 0.0              # zero z-component of first two vectors
  StdI.direct[2, :] = [0.0, 0.0, 1.0]   # third vector = unit z
  ```

**Benefits**:
- 5 scalar assignments reduced to 2 slice assignments
- Self-documenting with inline comments
- More idiomatic numpy

**Test results**:
- Unit tests: 1262 passed (no new tests needed — existing coverage)
- Integration tests: 83/83 passed

**Suggested next step**: Review remaining files for similar element-wise
numpy assignments that can use slicing, or consider Phase 2 class
introduction opportunities (e.g., SolverWriter hierarchy).

---

## Step 107 — Remove unused `import math` from `site_util.py`

**Date**: 2026-01-29
**Phase**: 3 — Leverage Python idioms (cleanup)
**Target**: `python/lattice/site_util.py`

**Motivation**: After Steps 91-106 replaced `math.cos/sin` with `np.exp`
and element-wise computations with numpy operations, the `math` module
is no longer used in `site_util.py`. Removing the unused import keeps
the module clean.

**Changes**:

- Removed `import math` from `site_util.py` (no remaining `math.*` calls)
- Verified all other files with `import math` still use it

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: The Phase 3 "Leverage Python idioms" work on
`site_util.py` is now complete (Steps 91-107). Consider targeting other
files for remaining Phase 3 improvements, or begin Phase 2 class
introduction (SolverWriter hierarchy, LatticeBuilder, etc.).

---

## Step 108 — Vectorize `large_value()` in `hphi_writer.py`

**Date**: 2026-01-29
**File**: `python/writer/hphi_writer.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: The `large_value()` function used 7 separate `for` loops to
sum absolute values of interaction arrays. These are numpy arrays, so the
loops can be replaced with vectorized `np.sum(np.abs(...))` calls.

**Changes**:

- Replaced 7 `for ... in range(N): large_value0 += abs(...)` loops with a
  single expression using `np.sum(np.abs(array[:count]))` for each term
- Reduces ~20 lines to ~8 lines while improving performance

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize `_read_w90()` in `wannier90.py` — replace
manual double-loop matrix multiply (lines 366-371) with `@` operator and
`np.linalg.norm()`.

---

## Step 109 — Vectorize matrix multiply in `_read_w90()` in `wannier90.py`

**Date**: 2026-01-29
**File**: `python/lattice/wannier90.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: The `_read_w90()` function used a manual double `for` loop
(lines 366-371) to compute `dR[ii] += direct[jj, ii] * tau_diff[jj]`, which
is a matrix-vector product `direct.T @ tau_diff`. The subsequent length
computation used `math.sqrt(dR[0]**2 + dR[1]**2 + dR[2]**2)` instead of
`np.linalg.norm()`.

**Changes**:

- Replaced 4-line double loop with `tau_diff = ...; dR = StdI.direct.T @ tau_diff`
- Replaced `math.sqrt(dR[0]**2 + dR[1]**2 + dR[2]**2)` with `np.linalg.norm(dR)`

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize `_geometry_w90()` in `wannier90.py` —
replace element-by-element coordinate reading loop with numpy array assignment.

---

## Step 110 — Simplify coordinate reading in `_geometry_w90()` in `wannier90.py`

**Date**: 2026-01-29
**File**: `python/lattice/wannier90.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: Two loops in `_geometry_w90()` assigned array elements one at
a time (`arr[i, 0] = ...; arr[i, 1] = ...; arr[i, 2] = ...`). Using row
slice assignment `arr[i, :] = [...]` is more concise and idiomatic.

**Changes**:

- `StdI.direct[ii, :]` row assignment replaces 3 separate element assignments (lines 125-129)
- `StdI.tau[isite, :]` row assignment replaces 3 separate element assignments (lines 137-141)

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize remaining C-style patterns in `wannier90.py`,
such as the Wigner-Seitz cell enumeration loops or the inversion symmetry loop.

---

## Step 111 — Vectorize `_count_and_store_terms()` in `wannier90.py`

**Date**: 2026-01-29
**File**: `python/lattice/wannier90.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: `_count_and_store_terms()` had two triple-nested loops. The
first multiplied each matrix element by its weight one-by-one — this is a
broadcast operation. The second manually stored surviving terms element by
element into pre-allocated arrays — this can use numpy masking.

**Changes**:

- Replaced per-element weight multiplication with broadcasting:
  `Mat_tot[:nWSC] *= Weight_tot[:nWSC, np.newaxis, np.newaxis]`
- Replaced second triple loop (store terms) with `np.nonzero(mask)` and
  `np.column_stack()` to extract indices and values in one shot
- Kept the printing loop unchanged (output order matters for compatibility)

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize the inversion symmetry loop in `_read_w90()`
in `wannier90.py`, or target `_merge_duplicate_terms()` in `common_writer.py`.

---

## Step 112 — Simplify index/Rmin/Rmax assignments in `wannier90.py`

**Date**: 2026-01-29
**File**: `python/lattice/wannier90.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: Both `_read_w90()` and `_read_density_matrix()` assigned
`indx_tot[iWSC, 0/1/2]` element-by-element. The density matrix reader also
used a `for ii in range(3)` loop with if-chains to track `Rmin`/`Rmax`.

**Changes**:

- `_read_w90()`: replaced 3 element-wise `indx_tot` assignments with row
  slice `indx_tot[iWSC, :] = [...]`
- `_read_density_matrix()`: same row-slice pattern for `indx_tot`, plus
  replaced `for ii in range(3)` min/max loop with `np.minimum()`/`np.maximum()`

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Replace the triple-nested `DenMat` dictionary
construction loop in `_read_density_matrix()` with `itertools.product`, or
target `_merge_duplicate_terms()` in `common_writer.py`.

---

## Step 113 — Replace triple loop with `itertools.product` in `_read_density_matrix()`

**Date**: 2026-01-29
**File**: `python/lattice/wannier90.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: The `DenMat` dictionary construction used three nested
`for` loops over `range(Rmin[i], Rmax[i]+1)`. Using `itertools.product`
flattens this to a single loop, improving readability.

**Changes**:

- Added `import itertools` to `wannier90.py`
- Replaced triple-nested loop with `itertools.product(range(...), range(...), range(...))`

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Target `_merge_duplicate_terms()` in `common_writer.py`
(O(n²) nested comparison loop), or vectorize the `_print_uhf_initial()` loops
in `wannier90.py`.

---

## Step 115 — Vectorize `_cell_diff()` in `geometry_output.py`

**Date**: 2026-01-29
**File**: `python/lattice/geometry_output.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: `_cell_diff()` used a list comprehension with `range(3)` to
compute element-wise differences. Since `Cell` is a numpy array, row
subtraction with `.astype(int).tolist()` is more idiomatic.

**Changes**:

- Replaced `[int(Cell[iCell, k] - Cell[jCell, k]) for k in range(3)]`
  with `(Cell[iCell] - Cell[jCell]).astype(int).tolist()`

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize the pump potential computation loop
(`for ii in range(3)`) in `vector_potential()` in `hphi_writer.py`, or
look for remaining C-style patterns in `wannier90.py`.

---

## Step 116 — Replace diagonal loop with `np.diag()` in `_apply_hopping_terms()`

**Date**: 2026-01-29
**File**: `python/lattice/wannier90.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: A `for ii in range(3)` loop set diagonal elements of a 3×3
matrix to the same computed value. Using `np.diag([val, val, val])` is more
concise and clearly expresses the intent of creating a diagonal matrix.

**Changes**:

- Replaced `np.zeros((3,3))` + loop setting `Jtmp[ii, ii]` with
  `np.diag([diag_val, diag_val, diag_val])`
- Extracted the repeated scalar computation into `diag_val`

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize the pump potential computation loop
(`for ii in range(3)`) in `vector_potential()` in `hphi_writer.py`, or
target remaining `for i in range(3): for j in range(3):` patterns in
`wannier90.py` (line 710) or `site_util.py` (line 322).

---

## Step 117 — Replace nested `range(3)` double loops with `itertools.product`

**Date**: 2026-01-29
**Files**: `python/param_check.py`, `python/lattice/input_params.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: Both `not_used_j()` in `param_check.py` and
`_resolve_spin_matrix()` in `input_params.py` used `for i1 in range(3):
for i2 in range(3):` nested loops to iterate over the 3×3 spin interaction
matrix. `itertools.product(range(3), repeat=2)` flattens this to a single
loop.

**Changes**:

- `param_check.py`: added `import itertools`; replaced double loop in
  `not_used_j()` with `itertools.product`
- `input_params.py`: added `import itertools`; replaced double loop in
  `_resolve_spin_matrix()` with `itertools.product`; fixed indentation

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Phase 3 improvements on small loops are nearing
exhaustion. Consider beginning Phase 2 class introduction (e.g., SolverWriter
hierarchy or LatticeBuilder), or audit remaining files for further idiom
opportunities.

---

## Step 118 — Replace list-of-lists with numpy arrays for `At` and `Et`

**Date**: 2026-01-29
**File**: `python/writer/hphi_writer.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: `StdI.At` and `Et` were allocated as list-of-lists
(`[[0.0] * 3 for _ in range(N)]`), a C-style pattern. Since they are
used as 2D numeric arrays (indexed as `[it][ii]`, passed to `np.dot()`),
`np.zeros((N, 3))` is more appropriate and enables future vectorization
of the pump computation loops.

**Changes**:

- Replaced `[[0.0] * 3 for _ in range(StdI.Lanczos_max)]` with
  `np.zeros((StdI.Lanczos_max, 3))` for both `StdI.At` and `Et`

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize the `for ii in range(3)` inner loop
in `vector_potential()` now that `At`/`Et` are numpy arrays, or look
for other list-of-lists patterns to convert.

---

## Step 119 — Replace inner `range(3)` loop with row assignment in `vector_potential()`

**Date**: 2026-01-29
**File**: `python/writer/hphi_writer.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: The inner `for ii in range(3)` loop in `vector_potential()`
called `handler_fn` per component and assigned to `At[it][ii]` and
`Et[it][ii]` individually. Now that `At`/`Et` are numpy arrays (Step 118),
we can collect the 3 results in a list comprehension and assign entire rows
with `At[it, :] = ...` and `Et[it, :] = ...`.

**Changes**:

- Replaced `for ii in range(3): At[it][ii], Et[it][ii] = handler_fn(...)`
  with list comprehension + row slice assignment

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Look for other list-of-lists or C-style allocation
patterns, or begin more substantial Phase 2/3 improvements such as converting
`StdI.npump` from list to numpy array.

---

## Step 120 — Remove `_SUFFIXES` alias in `input_params.py`

**Date**: 2026-01-29
**File**: `python/lattice/input_params.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: Step 100 consolidated `SPIN_SUFFIXES` into `param_check.py`
and added a backward-compatibility alias `_SUFFIXES = SPIN_SUFFIXES` in
`input_params.py`. Since this alias is purely internal, it can be removed
in favour of using `SPIN_SUFFIXES` directly everywhere in the module.

**Changes**:

- Removed `_SUFFIXES = SPIN_SUFFIXES` alias declaration
- Replaced all 7 uses of `_SUFFIXES` with `SPIN_SUFFIXES` throughout the file

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Look for other backward-compatibility shims or
unused aliases to remove, or begin Phase 2 work on introducing a
`KeywordParser` class.

---

## Step 114 — Replace O(n²) merge loop with dict-based O(n) in `_merge_duplicate_terms()`

**Date**: 2026-01-29
**File**: `python/writer/common_writer.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: `_merge_duplicate_terms()` used an O(n²) nested loop to find
duplicate index quadruples and merge their amplitudes. A dict keyed by the
4-tuple tracks first-occurrence indices, reducing complexity to O(n). The
counting loop was also replaced with a generator expression.

**Changes**:

- Replaced O(n²) double loop with single-pass dict-based merge
- Replaced counting loop with `sum(1 for k in range(n) if abs(vals[k]) > AMPLITUDE_EPS)`

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize `_print_uhf_initial()` loops in
`wannier90.py`, or target the `cos`/`sin` patterns remaining in
`interaction_builder.py`.

---

## Step 115 — Vectorize and deduplicate `_print_uhf_initial()` in `wannier90.py`

**Date**: 2026-01-29
**File**: `python/lattice/wannier90.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: `_print_uhf_initial()` had three issues: (1) an O(n²) Python
loop for counting nonzero entries, replaced with `np.count_nonzero`; (2) a
second O(n²) loop for writing, replaced with `np.nonzero` iteration; (3)
near-identical Coulomb (U) and Exchange (J) loop blocks, merged into a single
loop over `idx in (1, 2)`.

**Changes**:

- Replaced double counting loop with `np.abs(IniGuess) > AMPLITUDE_EPS` mask + `np.count_nonzero`
- Replaced double output loop with `np.nonzero(mask)` iteration
- Merged duplicate U/J loops into single parameterized loop
- Precomputed `0.5 * IniGuess[isite, jsite]` to avoid repeated multiplication

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Vectorize remaining `cos`/`sin` patterns in
`hphi_writer.py` or `interaction_builder.py`.

---

## Step 116 — Replace O(n²) accumulation loop with dict in `_accumulate_list()`

**Date**: 2026-01-29
**File**: `python/writer/export_wannier90.py`
**Phase**: 3 — Leverage Python idioms

**Motivation**: `_accumulate_list()` used an O(n²) linear scan (`_is_equal_key`
on every existing entry) to find duplicate keys and merge values. Replaced with
a dict keyed by tuple for O(1) lookup. The reverse-iteration zero-elimination
loop was also replaced with a single forward-pass filter.

**Changes**:

- Replaced O(n²) linear-scan deduplication with `dict[tuple[int, ...], complex]`
- Replaced reverse-iteration zero-elimination with forward-pass filter
- Preserved insertion order via `key_order` list

**Test results**:
- Unit tests: 1262 passed
- Integration tests: 83/83 passed

**Suggested next step**: Replace `_is_equal_key` (now unused for accumulation)
or refactor `_build_inter_table` / `_build_transfer_table` which also have O(n²)
linear scans for deduplication.
