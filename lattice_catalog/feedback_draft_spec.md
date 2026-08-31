# Feedback on the Lattice Definition Specification (draft 2026/07/27) — from cataloging and hands-on trials

- Date: 2026-08-15
- Basis: (1) cataloging every StdFace lattice and interaction (the 31 YAML
  files in this directory, cross-checked against source and lint-verified),
  (2) hands-on trials of 5 lattices/models outside the catalog
  (see [`examples_experimental/`](examples_experimental/))
- Purpose: proposals for the next revision of the draft specification "Lattice
  Definition Specification (draft)" (2026/07/27). Each item is given a
  priority (High/Medium/Low) and cites the sample it is based on.

## 0. Overall assessment

**The bond description via (from, to, R, type) and the three-layer
separation (geometry / system / model) hold up strongly in practice** — across
all 9 StdFace lattices x 3 models, plus Shastry-Sutherland, the alternating-bond
(dimerized) chain, breathing-type lattices, DM interactions, and cylindrical
boundaries, there was not a single case where the geometry/connectivity
description hit a wall. The improvement proposals concentrate mainly in two
areas: **(A) a parameter-declaration mechanism** and **(B) the operator
vocabulary**.

## 1. Proposals for the draft specification (in priority order)

### 1.1 [High] Introduce an in-file `parameters:` declaration block

- **Problem**: the draft writes coefficients as numeric literals (§7.3). The
  catalog dialect introduced `{param: <name>}` references, but bound the set
  of referenceable names to the StdFace keyword inventory. As a result, a
  model with parameters that don't exist in StdFace (e.g. the alternating-bond
  chain's `J_intra`/`J_inter`, or DM's `Dz`) can only be written by "bending
  the meaning of an existing name and borrowing it."
  (Evidence: `ssh_dimer_chain_spin.yaml`, `kagome_dm_spin.yaml` — wiped out by
  C7)
- **Proposal**: introduce a declaration block at the top of the file.

  ```yaml
  parameters:
    J_intra: {default: 1.0, type: real}
    J_inter: {default: 0.5, type: real}
    Dz:      {default: 0.1, type: real}
  ```

  The couplings side would reference it as `{param: J_intra}`. The validator
  would only need to check "is the referenced name declared?" — and for a
  StdFace-compatible catalog, the declaration block could be auto-generated
  from the inventory (retaining backward compatibility with the current
  catalog).

### 1.2 [High] Give the isotropic spin operator a place in the vocabulary (consistency with the draft itself)

- **Problem**: draft §4.3/§7.3 uses `operator: "S_i . S_j"`, but the catalog
  dialect canonicalized only the 9-component tensor_terms expansion, to
  support anisotropy. As a result, **the draft's own flagship example
  (Shastry-Sutherland) does not pass the catalog dialect's validator**.
  (Evidence: `shastry_sutherland_spin.yaml` — C9: unknown operator)
- **Proposal**: make `"S_i . S_j"` (isotropic spin-spin) an official
  vocabulary entry, and specify a mechanical expansion rule into
  tensor_terms (the same coefficient on the 3 diagonal components) in the
  spec. Isotropic as shorthand, anisotropic as tensor_terms — defining these
  roles at the spec level lets both coexist. At the same time, we recommend
  also adding `"s_i . s_j"` (fermion-fermion spin product), needed for
  multi-orbital Hund coupling, to the vocabulary (see 1.4).

### 1.3 [High] Bond-orientation convention: prefer source-order (definition-order) preservation over canonicalization

- **Problem**: draft §4.2 proposes an "orientation-canonicalization
  convention (e.g. lexicographic order)", but for DM interactions the bond's
  orientation *is* the sign of the coefficient itself (J_xy = +D_z,
  J_yx = -D_z requires a transpose under bond reversal). In a design where a
  single coupling definition is shared by multiple bonds, rewriting the
  orientation via canonicalization **becomes inexpressible the moment bonds
  that need reversal and bonds that don't are mixed together**.
- **Proposal**: adopt into the spec: "each bond is kept in its defined
  orientation and written exactly once. Uniqueness is checked via reversal
  equivalence (type, i, j, R) = (type, j, i, -R). The coefficient
  transformation under reversal (complex conjugate for hopping, transpose
  for the exchange tensor) is the consumer's responsibility." (We have
  confirmed the catalog can describe DM interactions without breakdown under
  this scheme — `kagome_dm_spin.yaml`)

### 1.4 [High] Multi-orbital support: orbital-indexed operators, or a general 4-fermion term

- **Problem**: `site_dof` can declare an orbital count, but there is no way
  for the operator side to refer to an orbital. Orbital-resolved hopping
  t_ab, Hund coupling, and pair hopping cannot be written. Even the
  "orbital = separate label" workaround only reaches inter-orbital
  density-density (U').
  (Evidence: `two_orbital_hubbard_chain.yaml`)
- **Proposal**: either of the following (can be combined).
  1. Specify a **general 4-fermion term** — an ordered sequence of
     creation/annihilation operators with site/orbital/spin bindings (a
     sketch already exists in this catalog's manual §6.6).
  2. As a near-term, practical scope: allow orbital indices on named
     operators (e.g. `hop[a,b]`, `density-density[a,b]`, `s_i . s_j`).
  The same mechanism is also needed to describe wannier90's J channels
  (Hund, exchange, pair-hopping; manual §5.5), so this is high priority.

### 1.5 [Medium] Allow `{param, scale}` in tensor_terms' coeff

- **Problem**: for antisymmetric exchange (DM), two components share the
  same parameter with opposite sign (`Jxy = +Dz`, `Jyx = -Dz`). Forcing an
  independent parameter name per component invites users to break the
  convention by "entering the same value in two places."
- **Proposal**: allow the `{param: Dz, scale: -1.0}` form in coeff.
  (The catalog dialect's C8/C12 enforce StdFace-compatible 9-component
  naming; this would be demoted to "an additional StdFace-catalog-only
  convention," not a general rule.)

### 1.6 [Medium] Make onsite (one-body terms) optional

- **Problem**: requiring an empty `onsite:` even for models with no onsite
  term is descriptive noise. (Evidence: `shastry_sutherland_spin.yaml` and
  others — C2)
- **Proposal**: make `onsite` optional, and define "omitted" as "no
  one-body term."

### 1.7 [Medium] Consistency-check semantics for open-boundary systems

- **Problem**: direction-dependent boundary conditions (e.g. cylindrical)
  are straightforward to write, but the coordination-number consistency
  check (this catalog's counting-expansion check) assumes a torus, and
  cannot express or verify the coordination deficit at edge sites along an
  open direction. (Evidence: `square_cylinder_spin.yaml`)
- **Proposal**: define a verification spec with expected values held on two
  layers — "bulk coordination number + edge deficit count" (this could also
  double as an acceptance test for the expansion engine).

### 1.8 [Low] Clarify that aperiodic modulation is out of scope

- Aperiodic modulation — a single impurity, or a coupling that differs only
  at the edge — is, in principle, out of scope for the three-layer design
  (unit cell + translation); a supercell approximation is the only recourse.
  This is a reasonable trade-off, but stating "out of scope" explicitly in
  the spec would help manage users' expectations.

## 2. Catalog-side tooling improvements (a TODO for this repository; no spec revision needed)

| Item | Description | Status |
|---|---|---|
| standalone mode | add `--standalone` to the linter (skip C10/C11's manifest cross-check, allow paths outside `lattice_catalog/`) so it can be used for syntax/consistency checks on user-authored lattices | not started (proposed) |
| onsite optional | once 1.6 is adopted into the draft, update C2 to match | pending spec |
| parameters declaration | once 1.1 is adopted, change C7 to a "declared-reference" check, and for StdFace-compatible catalogs keep the inventory cross-check as an additional check | pending spec |

## 3. Summary (one page for the next meeting)

- **Strengths (no change needed)**: the (from,to,R,type) bond description,
  the three-layer separation, direction-dependent boundaries. No breakdown
  in geometry description across all StdFace models plus the 5 non-trivial
  lattice trials.
- **To decide now (High)**: the parameter-declaration block / giving the
  isotropic operator a vocabulary entry / switching the orientation
  convention to source-order preservation / a path for multi-orbital
  (4-fermion terms).
- **Good to settle (Medium)**: coeff's {param, scale} / making onsite
  optional / the consistency-check format for open boundaries.
- **To document explicitly (Low)**: aperiodic modulation is out of scope.
