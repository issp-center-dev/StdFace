# StdFace Lattice Catalog

`lattice_catalog/` is a set of YAML files that transcribe every combination of the lattices
(chain, ladder, square, triangular, honeycomb, kagome, orthorhombic, fc_ortho,
pyrochlore) and models (Spin, Hubbard, Kondo) that StdFace supports out of the box,
plus one example of a wannier90 conversion, into the geometry / system / model
three-layer structure defined by the reference specification "Lattice Definition
Specification (draft)", 2026/07/27 (hereafter the "draft spec"). Each YAML is a
**source-cross-checked transcription**, read and copied directly from StdFace's
C/Python implementation (the `_BONDS` table, the operator-generation logic, the
parameter-resolution rules), and it passes the linter checks in `tools/
lint_catalog.py` (schema, bond-duplicate detection, param-reference consistency,
operator/type consistency, and cross-checking against `manifest.yaml`).
Matters the draft spec leaves unspecified (fermion site_dof, the single-site
operator vocabulary, `{param, scale, default}` references, the `catalog:` header,
etc.) are filled in with an experimental extension dialect explicitly marked as
`dialect: experimental`. Note, however, that what this catalog guarantees is
"static agreement with the reference implementation's source"; it does not go so
far as to verify runtime equivalence (an oracle comparison) -- that is, whether
actually running the generated Hamiltonian numerically produces results matching
StdFace's output (see `manual.md` Section 1.1 / chapter 7 for details).

## Directory layout

```
lattice_catalog/
  README.md                 # this file
  CONVENTIONS.md             # authoring conventions (normative document)
  manual.md                  # guide (how to read it, examples, full keyword mapping table)
  manifest.yaml               # consistency-check ledger (bonds_per_uc/coordination, etc., machine-readable)
  tools/
    lint_catalog.py           # semantic-check linter (C1-C12)
    keyword_inventory.py      # generates an inventory of all StdFace keywords
    test_tools.py              # automated tests for the tools themselves (development use only)
  chain/        chain_{spin,hubbard,kondo}.yaml
  ladder/       ladder_w2_{spin,hubbard,kondo}.yaml, ladder_w3_{spin,hubbard,kondo}.yaml
  square/       square_{spin,hubbard,kondo}.yaml
  triangular/   triangular_{spin,hubbard,kondo}.yaml
  honeycomb/    honeycomb_{spin,hubbard,kondo}.yaml
  kagome/       kagome_{spin,hubbard,kondo}.yaml
  orthorhombic/ orthorhombic_{spin,hubbard,kondo}.yaml
  fc_ortho/     fc_ortho_{spin,hubbard,kondo}.yaml
  pyrochlore/   pyrochlore_{spin,hubbard,kondo}.yaml
  wannier90/    example_hubbard.yaml
```

With 9 lattices × 3 models, that is 30 files, plus 1 wannier90 Hubbard example
file, for a total of **31 files** that are checked by `lint_catalog.py`
(`manifest.yaml` is a data file and is not counted among the checked files).

## Using the linter

Development-time dependency: this linter requires PyYAML (`pyyaml`) (a development-only
tool, not included in the runtime dependencies of `python/pyproject.toml`.
If it is not installed, the linter exits with a clear message prompting
`pip install pyyaml`).

Run it from the repository root:

```bash
python3 lattice_catalog/tools/lint_catalog.py
```

For all 31 files, it performs schema checks (including deep type checks of
geometry/sites/bonds/site_dof), label consistency, duplicate detection under bond
reversal equivalence, param-reference consistency, J-tensor component
completeness, operator/site_dof type consistency, the sign convention
(`CONVENTIONS.md` §6.2), and cross-checking against `manifest.yaml` (a "counting
expansion" check that actually expands the bonds onto a torus of the specified
size and compares the measured `bonds_per_uc`/`coordination` against the
expected values, plus a strictness check of `min_size_for_check` itself). If
there are no errors, it prints

```
31 files, 0 errors
```

For details on the individual checks (C1-C12), see `CONVENTIONS.md` and
`manual.md` Section 1.1.

`tools/` also contains `keyword_inventory.py`, which scans and lists the
keywords of every StdFace solver (HPhi/HWAVE/UHF/mVMC), and `test_tools.py`
(development use only, a self-contained test runner that does not use pytest),
which checks the behavior of the linter and inventory tools themselves.

## Guide to manual.md

How to read this catalog, examples, and the full mapping table against every
keyword are collected in `manual.md`. Its structure is as follows:

1. **Overview and how to read it** -- a summary of the three-layer specification,
   the positioning of `dialect: experimental`, this catalog's verification level
   (source cross-checking, the linter, and manifest consistency checks -- not a
   numeric oracle comparison), the directory layout, and supplementary
   explanation of the main conventions in CONVENTIONS.md.
2. **Keyword mapping table** -- a table that assigns **every one** of the
   StdFace keywords reported by `keyword_inventory.py` (351 items: 313
   keywords + 26 lattice_alias + 12 model_alias) to one of 7 categories:
   geometry / system / bonds+couplings / onsite / site_dof / wannier90 / out of
   scope.
3. **Per-lattice notes** -- for each of the 9 lattices, its geometry, bond
   definition table, manifest consistency-check basis, and provenance.
4. **Per-model operator mapping** -- the derivation of the sign table (the
   verification chain), Spin's J-tensor resolution rules, Hubbard/Kondo's
   operators and signs, and the GC variants.
5. **wannier90 conversion specification** -- the conversion rules from
   RESPACK/Wannier90 output (the sign, cutoff, and Hermitian canonical-pair
   selection for the H/U channels, automatic generation of Spin's
   superexchange, and a prose description of the J channel).
6. **Specification extension proposals** -- an organization, as a proposed
   addendum to the draft spec, of the 6 extensions introduced by
   `dialect: experimental` (fermion site_dof, the single-site operator
   vocabulary, named two-body operators, `{param, scale, default}`
   references, the `catalog:` header, and the sketch of general 4-fermion
   terms).
7. **Known limitations and unhandled items** -- the oracle comparison not yet
   performed, the W=2/3 limitation of ladder, the scope of wannier90 support,
   known inconsistencies in the upstream implementation, and so on.

The normative document (which takes precedence in case of conflict with this
one) is `CONVENTIONS.md`; `manual.md` is positioned as its explanation, examples,
and index.

## Reference specification

The geometry / system / model three-layer structure of this catalog conforms to
the reference specification "Lattice Definition Specification (draft)",
2026/07/27 (the basis for the design decisions organized in
`docs/superpowers/specs/2026-08-15-lattice-catalog-design.md`). Matters the
draft does not explicitly specify are documented and adopted as an
"experimental extension dialect" in this README, `CONVENTIONS.md` §7, and
`manual.md` chapter 6.
