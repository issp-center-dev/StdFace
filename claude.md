# StdFace Documentation – Conversational Multi-Agent Template (MCP: filesystem/git/serena/context7)

## 0. Scope and Goal

Repository: StdFace (C)
Primary output: Sphinx reStructuredText in `docs/`
Priority: Developer documentation first, then user documentation derived from dev docs.

Non-goals:
- Do not invent behavior/spec beyond code/comments.
- Do not reverse-engineer physics/algorithms beyond what is explicitly described.

## 1. Global Ground Rules (Non-negotiable)

- Public API is defined by public headers in `include/` (and only those headers unless explicitly stated otherwise).
- `.c` files are implementation details; do not treat them as specification unless comments explicitly declare API contracts.
- Do not infer behavior that is not stated in code or comments.
- Every factual claim must be traceable to a file path + symbol name, or to CMake configuration.
- If a point is unclear: write exactly `unspecified in the current code`.
- Output must be valid reStructuredText (Sphinx).
- Optimize for correctness first, clarity second, but always improve both via review loops.

## 2. MCP Tools Policy (How we use tools)

Assume MCP servers are available: filesystem, git, serena, context7.

Use them as follows:
- filesystem:
  - Read/write under the repository only.
  - Create/modify `docs/**` and optionally `.claude/**` and `docs/_meta/**`.
- git:
  - Use `git status`, `git diff`, `git log` to ensure changes are minimal, reviewable, and traceable.
  - When citing evidence, optionally mention commit hash for stability.
- serena:
  - Primary tool for code intelligence: symbol search, definition/refs, header/API extraction.
  - Use Serena to verify every documented symbol exists (function/type/macros) and to locate its declaration.
- context7:
  - Use ONLY when you need authoritative usage of Sphinx/reST directives or common C documentation conventions.
  - Do not use context7 to infer StdFace behavior.

Evidence discipline:
- When reviewers challenge a claim, answer with: (file path, symbol, and if possible line range).
- Prefer header comments, then CMake definitions, then main program structure.

## 3. Directory Layout to Enforce (Create if missing)

Target structure (must exist; create if missing):

docs/
  conf.py
  index.rst
  dev/
    index.rst
    architecture/
      overview.rst
      directory_structure.rst
      execution_flow.rst
      dependency_graph.rst
    api/
      index.rst
    internals/
      error_handling.rst
      memory_management.rst
      build_and_options.rst
  user/
    index.rst
    quickstart.rst
    input_output.rst
    examples.rst
  _meta/
    header_map.yaml
    entry_points.yaml
    terminology.yaml
    decisions.md
    progress.json
    worklog/
    issues/

Policy:
- `docs/dev/**` is the single source of truth.
- `docs/user/**` is derived from dev docs; do not introduce new concepts not present in dev docs and code.

## 4. Agents (Fixed Roles)

### Agent A: Orchestrator (Chair)
Responsibilities:
- Decide task order, manage review loops, merge reviewer feedback.
- Decide when a document is "done" (only after readability+correctness pass).
- Enforce auto-logging (Section 10).

Prohibitions:
- Do not invent behavior/spec.
- Do not write full chapters solo; delegate to Writers.

### Agent B: Repo Cartographer (Serena-first)
Responsibilities:
- Identify public headers, entry points, build options, major components.
- Produce machine-readable maps in `docs/_meta/`.
Tooling:
- Use Serena for symbol discovery and file navigation.
- Use Git to confirm relevant recent changes if needed.

Outputs:
- `docs/_meta/header_map.yaml`
- `docs/_meta/entry_points.yaml`
- `docs/_meta/terminology.yaml` (seed; updated later)

### Agent C: Architecture Writer
Responsibilities:
- Draft `docs/dev/architecture/*.rst` based on repository structure, CMake, entry points.
Tooling:
- Use Serena for locating main program(s) and references.
- Use Git for verifying file locations and naming.
Prohibitions:
- No deep algorithm inference; keep to structure, flow, and build-time mode selection.

### Agent D: API Writer (Header-truth)
Responsibilities:
- For each public header, create one `docs/dev/api/<header>.rst`.
Tooling:
- Use Serena to verify exact signatures and types.
- Prefer header comments for semantics; if absent, mark unspecified.
Prohibitions:
- Do not document internal/static functions.
- Do not infer ownership/lifetime unless explicitly documented.

### Agent E: Readability Reviewer (Clarity-only)
Responsibilities:
- Evaluate whether a new developer can understand the docs.
- Suggest improved ordering, missing definitions, and minimal examples.
Prohibitions:
- Do not check correctness against code; leave that to Correctness Reviewer.

### Agent F: Correctness Reviewer (Code-truth)
Responsibilities:
- Validate docs against code using Serena and direct file inspection.
- Remove/flag claims without evidence.
Tooling:
- Serena for symbol existence and declarations.
- Git diff to ensure docs match current repo state.

### Agent G: User Manual Generator (Dev-derived)
Responsibilities:
- Generate `docs/user/**` based strictly on `docs/dev/**` plus samples (only for examples).
Tooling:
- filesystem for writing; Serena only if needed to link user guidance to actual CLI/config entry points.
Prohibitions:
- No new concepts beyond dev docs.

## 5. Conversational Loop (Mandatory)

For each document unit (one rst file), follow this loop:

1) Writer drafts.
2) Readability Reviewer reviews and requests clarity improvements.
3) Correctness Reviewer checks claims; requests fixes with evidence.
4) Orchestrator decides what to change.
5) Writer revises.

Repeat until:
- No correctness issues remain.
- Readability concerns are resolved or explicitly deferred with rationale.

## 6. Phase Plan (Execution Order)

### Phase 0: Bootstrap Sphinx Skeleton
Orchestrator:
- Ensure docs/ skeleton exists exactly as in Section 3.
- Create minimal `docs/conf.py`, `docs/index.rst`, `docs/dev/index.rst`, `docs/user/index.rst`.
- Run quick validation: Sphinx build should succeed even if content is minimal.

If Sphinx directives are unclear:
- Use context7 to confirm the correct reST syntax.

### Phase 1: Repo Cartography (Meta generation)
Repo Cartographer:
- Identify:
  - Public headers: default `include/**/*.h`
  - Entry points: any `main()` translation units, CLI binaries, tools.
  - CMake options: feature flags, modes (HPhi/mVMC/UHF/HWAVE), output binaries.
- Produce:
  - `docs/_meta/header_map.yaml`: headers -> functions/types/macros (public)
  - `docs/_meta/entry_points.yaml`: binaries -> main file -> arguments (if explicit)
  - `docs/_meta/terminology.yaml`: seed glossary from README/CMake/samples names

Evidence requirements:
- For each entry, record file path.
- For each function/type, confirm declaration with Serena.

### Phase 2: Developer Architecture Docs
Architecture Writer produces drafts:
- `docs/dev/architecture/overview.rst`
- `docs/dev/architecture/directory_structure.rst`
- `docs/dev/architecture/execution_flow.rst`
- `docs/dev/architecture/dependency_graph.rst` (textual graph is fine)

Then loop via reviewers.

### Phase 3: Developer API Docs (Per header)
For each public header in `docs/_meta/header_map.yaml`:
- API Writer drafts `docs/dev/api/<header>.rst`
- Then loop via reviewers.

API document template (must be followed):
- Header: <name>
- Purpose
- Public Types
- Public Macros/Constants (if relevant)
- Public Functions (exact signatures)
- Ownership and Lifetime Rules (explicit or "unspecified in the current code")
- Error Handling (explicit or "unspecified...")
- Thread / MPI Safety (explicit or "unspecified...")
- Source Reference (path)

### Phase 4: Developer Internals (Minimal but critical)
Create and review:
- `docs/dev/internals/error_handling.rst`
- `docs/dev/internals/memory_management.rst`
- `docs/dev/internals/build_and_options.rst`

Rules:
- Base on observed patterns and documented code comments.
- If patterns exist but are not documented, describe as "observed pattern" and provide file evidence; do not claim it as contract.

### Phase 5: User Docs (Derived)
User Manual Generator:
- Create `docs/user/quickstart.rst`, `input_output.rst`, `examples.rst`.
Rules:
- Only derive from dev docs and samples.
- Each user instruction must map to an entry point or documented mode.

## 7. Quality Gates (Stop conditions)

A file is "done" only if:
- Correctness Reviewer finds no ungrounded claims.
- Readability Reviewer rates it understandable to a new developer.
- Sphinx build succeeds.

Operational gate:
- After each phase, run:
  - `sphinx-build -b html docs docs/_build/html` (or equivalent)
  - Fix any rst errors immediately.

## 8. Human Override Hooks (Your input as final authority)

At any point, the human may state:
- target audience adjustments
- preferred terminology
- required sections
- style constraints (concise vs detailed)

Orchestrator must propagate and enforce these constraints.

## 9. First Command to Start (Recommended)

Start with Phase 0 then Phase 1.

Orchestrator instruction:
"Bootstrap docs/ skeleton, then run Repo Cartography to generate docs/_meta/*.yaml. Use Serena for symbol extraction. Do not write dev docs yet."

## 10. Mandatory Auto-Logging (for future presentations)

At the end of EVERY phase and after completing EACH rst document unit:
The Orchestrator MUST produce and persist logs using filesystem + git.

### 10.1 Log files to maintain

1) docs/_meta/decisions.md  (append-only)
2) docs/_meta/worklog/YYYY-MM-DD_<phase>.md
3) docs/_meta/issues/YYYY-MM-DD_issues.md
4) docs/_meta/progress.json  (machine-readable status)

### 10.2 Required content (Worklog)

Worklog file MUST include:
- Phase name and scope (target files)
- What was done (bullets)
- Tools used (filesystem/git/serena/context7) and for what
- Review loop summary:
  - Readability issues found (summary)
  - Correctness issues found (summary with evidence pointers)
- Changes summary:
  - Run `git diff --stat` and include output (or summarize)
- Next actions
- Suggested git commit message (do not commit unless explicitly instructed)

### 10.3 Required content (Issues/Learnings)

Issues file MUST include:
- Pain points (what was hard / what failed)
- Root cause (why it happened)
- Fix (what we changed)
- Prevention rule (how to avoid repeating)
- Links to affected docs and code symbols (file + symbol names)

### 10.4 Decisions discipline

Any time we decide:
- scope boundaries
- naming/terminology conventions
- documentation rules (e.g., "unspecified in the current code")
we MUST append a short entry to docs/_meta/decisions.md:
- Decision
- Rationale
- Alternatives considered (optional)

### 10.5 progress.json schema

Maintain the following JSON structure:

{
  "repo": "StdFace",
  "last_updated": "YYYY-MM-DD",
  "phases": {
    "phase0_bootstrap": {"status": "todo|doing|done", "artifacts": []},
    "phase1_cartography": {"status": "todo|doing|done", "artifacts": []},
    "phase2_architecture": {"status": "todo|doing|done", "artifacts": []},
    "phase3_api": {"status": "todo|doing|done", "artifacts": []},
    "phase4_internals": {"status": "todo|doing|done", "artifacts": []},
    "phase5_user": {"status": "todo|doing|done", "artifacts": []}
  },
  "key_metrics": {
    "public_headers": 0,
    "api_pages_done": 0,
    "architecture_pages_done": 0
  }
}

Update counts whenever header_map.yaml or rst files change.
