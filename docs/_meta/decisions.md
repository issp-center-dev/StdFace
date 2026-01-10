# Documentation Decisions (StdFace)

This file is append-only. Add dated entries.

## Global rules
- Headers define the public API.
- No inference from implementation; if unclear, write `unspecified in the current code`.
- Developer docs first; user docs derived from dev docs.

## Tooling policy
- Serena is mandatory for symbol verification.
- context7 is allowed only for Sphinx/reST syntax and conventions.

---

## 2026-01-10: Phase 0 Bootstrap Validation

**Decision**: Accept skeleton structure as complete despite minor RST underline issue.

**Rationale**: The `docs/index.rst` title underline is 1 character short (20 vs 21), but this is a cosmetic warning, not an error. Sphinx will still build successfully. The issue is documented in `docs/_meta/issues/2026-01-10_issues.md` for future correction.

**Alternatives considered**: Fix immediately vs. defer. Chose to defer since it's non-blocking and Phase 0 goal was validation, not correction.

---

## 2026-01-10: Phase 1 - Headers in src/, not include/

**Decision**: Document headers in their actual location (`src/`) rather than expecting `include/`.

**Rationale**: The project does not follow the typical C convention of placing public headers in `include/`. All 5 headers are in `src/`. Per claude.md ground rules, we document what exists without inference. The header_map.yaml explicitly notes this deviation.

**Alternatives considered**:
1. Treat all src/ headers as internal (rejected - they contain the public API)
2. Request project restructuring (rejected - out of scope for documentation)
3. Document as-is with explicit note (chosen)

---

## 2026-01-10: Phase 1 - Conditional API documentation

**Decision**: Include conditionally-compiled functions in header_map.yaml with explicit `conditional` field.

**Rationale**: `export_wannier90.h` contains functions only available when `_HWAVE` is defined. These are part of the public API for HWAVE mode and should be documented, but with clear indication of the condition.

**Alternatives considered**:
1. Omit conditional functions (rejected - they are public API for one mode)
2. Document without noting condition (rejected - misleading)
3. Document with conditional field (chosen)
