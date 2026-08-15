# lattice_catalog Authoring Conventions (stdface-catalog/0.1)

This document defines the normative conventions that the YAML catalogs under `lattice_catalog/` must follow.
The reference specification is "Lattice Definition Specification (draft)", 2026/07/27 (hereafter the "draft spec").
Where this document and the draft spec conflict, the draft spec takes precedence for the parts it covers,
and this document is normative for matters the draft does not address (experimental extensions).
For the design rationale, see §4/§6 of `docs/superpowers/specs/2026-08-15-lattice-catalog-design.md`.

---

## 1. Positioning and the dialect declaration

Every YAML file must begin with the following header.

```yaml
catalog:
  schema: stdface-catalog/0.1
  dialect: experimental
  lattice: <lattice name>   # e.g. chain, square, kagome, ...
  model: <model name>     # e.g. spin, hubbard, kondo
```

- `schema` is the version identifier for this convention document. Fixed value `stdface-catalog/0.1`.
- `dialect` is a fixed value, `experimental`, indicating an extension dialect relative to the draft spec
  (see §7).
- `lattice` / `model` are the source used for the manifest cross-check (linter check C10).
  They must match the `lattice` / `model` of the corresponding entry in
  `lattice_catalog/manifest.yaml`.
- The parts that conform to the draft spec (the geometry / system / model three-layer structure itself)
  and the extension dialect items enumerated in §7 of this document (fermion site_dof, the operator
  vocabulary, `{param, scale, default}` references, the `catalog:` header, etc.)
  are treated as distinct. For matters the draft does not specify, the conventions in this document
  are the sole normative reference.

## 2. File structure

Every YAML file is self-contained and has the following four elements as top-level keys.

```yaml
catalog:  { ... }   # §1
geometry: { ... }   # §3
system:   { ... }   # §4
model:    { ... }   # §5, §6
```

References to other files (via `include` etc.) are not used. 1 file = 1 lattice × 1 model.

## 3. geometry conventions

- `geometry.sites` is an array, and **the array order matches the "site number" of draft spec
  §4.2**. That is, `sites[0]` corresponds to site 0 in the spec, `sites[1]` to site 1, and this
  correspondence is maintained mechanically.
- Site labels (`sites[*].label`) use the conventional names from the spec document or the source
  implementation (`A`, `B`, `A3`, etc.) and must not be duplicated within a file (C2).
- `geometry` holds the unit cell's lattice vectors `geometry.lattice_vectors`
  (a dictionary of named vectors `a1`, `a2`, `a3`, ... for as many dimensions as there are;
  the matrix `A` in §5 is these stacked as row vectors), and each site's
  fractional coordinates `frac`.

## 4. system conventions

- Repeat counts such as `W, L, Height` are mapped to `system.size` (an integer array).
  `system.size` only indicates a representative value for stan.in and is independent of the
  linter's check size (the manifest's `min_size_for_check`). `min_size_for_check`
  is a check size derived independently from the maximum `|R|` component of each bond,
  and for each direction it is the smallest odd number exceeding `2 * max|R|`.
- `phase0`-`phase2` are mapped to `system.boundary` (an array whose elements are each
  `{twist: {param: phaseN}}`), and the phase factor for a path crossing the boundary
  n times is interpreted by the consumer as `exp(i · n · π · θ / 180)`
  (**in degrees**). `{param: phaseN}` itself is expressed in the general form of §6.3.
- The chain lattice is logically one-dimensional, but the current implementation internally
  uses a two-dimensional representation with `W=1`. The transcription from `phase0` to the
  internal `phase[1]` is an implementation detail, explained in manual chapter 3 (it does not
  affect the YAML's semantics).
- `box` (the supercell transformation matrix) corresponds to `system.supercell` (row-vector
  convention `A_super = S · A`, `det(S) ≠ 0`). `supercell` and `W/L/Height`
  (= `size`) are used exclusively of each other. The details of this correspondence are
  documented only in the manual; this catalog's YAML uses the `size` representation as the base.

## 5. bonds conventions

- **Definition of the cell difference R**: `R = cell(to) − cell(from)` (an integer vector,
  whose dimension matches `dimension`).
- **Definition of the displacement δ**: `δ = (frac_to − frac_from) + R · A`
  (`A` is the matrix formed by stacking the row vectors `a1`, `a2`, ... of
  `geometry.lattice_vectors`).
- **type names**: Use StdFace's keyword names as-is
  (`J0`, `J0'`, `t0`, etc.). Names containing a prime are written in YAML with
  quotes (`"J0'"`). Splitting primes into a mechanical identifier is a future task
  (manual chapter 7).
- **Direction (source-order preservation)**: The `from` → `to` direction of a bond
  preserves the source order of the reference implementation as-is.
  - Ordinary bonds: the `site_i → site_j` order in the `_BONDS` table.
  - Kondo coupling: the argument order of `general_j(..., isite, jsite)`
    (**the itinerant site is the first argument**).
  Reordering into a canonical form is not performed. Reordering would require the complex
  conjugate of complex hopping and the transpose of the exchange tensor, which cannot be
  expressed with a single `couplings` key
  (design §4.2, Round 2 decision).
- **Uniqueness check**: Each bond is described exactly once. The linter detects duplicates
  under **reversal equivalence** `(type, i, j, R) ≡ (type, j, i, −R)` (C6).
  It is valid for different `type`s (e.g. `t0` and `V0`) to coexist on the same
  geometric bond; this is distinguished by including `type` in the equivalence
  key.
- **Coefficient transformation under reversal (normative for the consumer)**: When reading
  `(i,j,R)` as `(j,i,−R)`,
  - hopping (complex number): take the complex conjugate (`t_ij = conj(t_ji)`).
  - exchange tensor (J family, 3×3): take the transpose (`J_ij = J_ji^T`).
  Since the catalog itself is described uniquely in source order, this transformation is a
  norm needed only when the **consumer of the catalog** (the resolver/expansion engine)
  performs a reversal reinterpretation; catalog authors only need to write things in source
  order.

## 6. couplings / onsite conventions

### 6.1 value semantics

`model.couplings[<type>].value` (for scalar/vector-family couplings: `hop` /
`density-density` / `s_i . S_j`) and
`model.onsite[<site>][<term>].value` are

```
H = Σ_bonds value · operator + Σ_onsite value · operator
```

that is, they are the **coefficient of the physical Hamiltonian**. Note that this is
**distinct** from the coefficient convention of HPhi's solver output files such as
`trans.def` (where the sign is implicitly reversed, as in
`H_trans = −Σ t c†c`). The correspondence between the two is cross-checked
via the following verification chain (documented with derivations in manual chapter 4):

```
StdFace parameters → builder call (interaction_builder.py)
→ trans/intr coefficients → solver output convention → sign of the physical Hamiltonian
```

**J-family / onsite exception**: The J-family exchange interaction
(`model.couplings[<type>].operator.tensor_terms`) has no shared `value`;
each component term individually carries `coeff: {param: ...}` (§6.4). Each onsite
term (`model.onsite[<site>][<term>]`) likewise takes the form
`operator: {tensor_terms: [{ops: [...], coeff: <number>}]}`, but here
`tensor_terms` is a single-element array whose `coeff` is a **literal number**
representing sign only (`+1.0` / `-1.0`, etc.). The actual external parameter
reference is carried by that same term's `value` (`{param, scale, default}`, §6.3).

**Values sourced from external data (wannier90 example)**: For concrete numeric values
read from external data such as `_ur.dat` (e.g. measured values of `H_mn(R)`), it is
acceptable to write a numeric literal directly in `value` (the general form
`{param, scale, default}` is not required). Since C7 scans only dict-typed `value`s
that carry a `param` key via `_find_param_refs`, a non-dict `value` (a numeric
literal) is outside the scope of that check
(see `lattice_catalog/wannier90/example_hubbard.yaml`).

### 6.2 Sign table (normative)

| StdFace | Contribution to physical Hamiltonian | YAML representation |
|---|---|---|
| t family (hopping) | −t Σ_σ (c†c + h.c.) | `couplings[t0]: {operator: hop, value: {param: t0, scale: -1.0}}` |
| mu (chemical potential) | −mu N | onsite `operator.tensor_terms: [{ops: [N], coeff: -1.0}]`, `value: {param: mu}` |
| U (onsite Coulomb) | +U n↑n↓ | onsite `operator.tensor_terms: [{ops: [NupNdn], coeff: 1.0}]`, `value: {param: U}` |
| V family (intersite Coulomb) | +V n_i n_j | `couplings[V0]: {operator: density-density, value: {param: V0}}` (scale may be omitted, default +1.0) |
| J family (exchange interaction) | +Σ_ab J_ab S^a S^b | `couplings[J0]: {operator: {tensor_terms: [...]}}` (§6.4) |
| h / Gamma / Gamma_y (magnetic field) | −h Sz − Γ Sx − Γy Sy | onsite `operator.tensor_terms: [{ops: [Sz/Sx/Sy], coeff: -1.0}]`, `value: {param: h/Gamma/Gamma_y}` |
| D (single-ion anisotropy) | +D (Sz)² | onsite `operator.tensor_terms: [{ops: [Szz], coeff: 1.0}]`, `value: {param: D}` |
| Kondo J (s·S coupling) | +J s·S | `couplings[J]: {operator: "s_i . S_j", value: {param: J}}` (+1.0 if scale is omitted) |

**The sign must always be carried in the data (`scale` / `coeff`), never written as a
comment.** This is a required convention to guarantee verifiability (mechanical cross-checking
by the linter and downstream processing).

### 6.3 General form of param references

`{param, scale, default}` is an extension dialect (§7) with the following general form.

```yaml
value: {param: <name>, scale: <real number, default 1.0>, default: <number, default 0>}
```

- **If param is specified (by the consumer's input)**: value = `scale × param`.
- **If param is not specified**: value = `default` (**treated as the final value, with
  scale not applied**).

That is, `scale` is "a multiplier that applies only when param is given," and
`default` is "the final value when param is absent," not the product with `scale`.
This branching is itself normative and is fixed by the C7 check and the tool tests
(`test_tools.py`).

Example: `{param: 2S, scale: 0.5, default: 0.5}`
- If `2S=1` is given → S = 0.5 × 1 = 0.5.
- If `2S` is not given → S = 0.5 (`default` is taken as-is,
  `scale` is not applied).

Type constraints: `scale` is a real number, `default` is a number (integer or real, either
is acceptable). Type constraints for individual parameters follow the semantics of each
keyword (e.g. `2S` is a positive integer).

The wannier90 hop channel follows a different rule from the above: `wannier90.py::
_apply_hopping_terms` passes the nonlocal term as `hopping(StdI, -Cphase*tUJ[0][it],
jsite, isite, dR)` (the builder call itself carries a sign flip), and this
cancels with the solver-side sign flip of trans.def (§6.1), with the result that
**the physical hopping coefficient is `+H_mn(R)` itself (no sign flip)**
(the flip `H_mn → −H_mn` does not occur). The verified derivation is documented
separately in Section 5.2 of the manual (this is outside the scope of this catalog's
`{param, scale, default}` general form).

### 6.4 Canonical form of the J-family 9-component tensor_terms

The J-family exchange interaction is expressed as the 9 components of
`model.couplings[<type>].operator.tensor_terms` (it is not reduced to a single parameter
such as an isotropic `J0` -- all 9 terms are always written out in full). Each term
consists of a product of two spin operators `ops: [S?, S?]` and its
coefficient `coeff: {param: <name>}`.

Correspondence table between ops pairs and param suffixes:

| ops pair | suffix | example (prefix=J0) |
|---|---|---|
| `[Sx, Sx]` | `x` | `J0x` |
| `[Sy, Sy]` | `y` | `J0y` |
| `[Sz, Sz]` | `z` | `J0z` |
| `[Sx, Sy]` | `xy` | `J0xy` |
| `[Sx, Sz]` | `xz` | `J0xz` |
| `[Sy, Sx]` | `yx` | `J0yx` |
| `[Sy, Sz]` | `yz` | `J0yz` |
| `[Sz, Sx]` | `zx` | `J0zx` |
| `[Sz, Sy]` | `zy` | `J0zy` |

Note that for the diagonal components (`x`, `y`, `z`), both ends of `ops` are the same
operator, while for off-diagonal components (`xy`, etc.), they are an ordered product of
different operators (`xy` and `yx` are independent parameters).

### 6.5 Parameter resolution order (implementation-conformant)

The value of each component is resolved with the same priority order as
`input_params.py::_resolve_spin_matrix` (the catalog provides this resolution rule as a
norm; the actual resolution is performed by the **catalog consumer** (the resolver). The
catalog's output contract is "resolved numeric coefficients per bond"):

1. Component-local (e.g. `J0xy`)
2. Component-global (e.g. `Jxy`)
3. Scalar-local, diagonal only (e.g. `J0`, applies only to the `x=y=z` components)
4. Scalar-global, diagonal only (e.g. `J`, applies only to the `x=y=z` components)
5. 0 (default)

**Conflict rules** (identical to `input_spin_nn`. All of the following combinations are
detected as errors):

- Scalar vs. scalar: `J` (global scalar) and `J0` (local scalar) specified at the same time.
- Scalar vs. matrix (all 4 combinations): `J` vs. `J` (matrix),
  `J` vs. `J0` (matrix), `J0` (scalar) vs. `J` (matrix),
  `J0` (scalar) vs. `J0` (matrix).
- Matrix vs. matrix: components of `J0` and components of `J` specified at the same time.

**Exception for the primed series** (identical to `input_spin`): primed series such as
`J0'`, `J0''`, `J1'`, etc. have **no global fallback**. That is, only the conflict
between the scalar version `J0'` and the component version `J0'xy` is checked for
`J0'xy`, and the corresponding global variable (`J'`, etc.) never appears in the
fallback chain.

Resolution table by prefix (summary):

| prefix series | local component | local scalar | global component fallback | global scalar fallback |
|---|---|---|---|---|
| `J0`, `J1`, `J2` (nearest-neighbor series) | yes | yes | `J` components | `J` (scalar) |
| `J0'`, `J0''`, `J1'`, ... (primed series) | yes | yes | none | none |
| `J` (the global itself) | -- | -- | -- | -- |

### 6.6 Operator semantics and endpoint order

- `hop` = Σ_σ (c†_iσ c_jσ + h.c.) (unsigned. The sign is carried in
  `value`'s `scale`/`coeff` -- §6.2).
- `density-density` = n_i n_j.
- `s_i . S_j` (Kondo coupling): **the first endpoint (i) is the itinerant electron spin,
  the second endpoint (j) is the localized spin**. Endpoint order is meaningful, and the
  order is guaranteed by the source-order-preservation convention of §5.
- onsite operator vocabulary: `N`, `Nup`, `Ndn`, `NupNdn`, `Sx`, `Sy`, `Sz`, `Szz`.
  These appear as elements of
  `model.onsite[<site>][<term>].operator.tensor_terms[0].ops` (a single operator;
  see the exception in §6.1).

Type-compatibility table between operators and `site_dof` (checked by linter C9):

| operator | required site_dof types (endpoint order) |
|---|---|
| `hop` (couplings, `operator` string) | fermion - fermion |
| `density-density` (couplings, `operator` string) | fermion - fermion |
| `s_i . S_j` (couplings, `operator` string) | fermion - spin (in this order, the first endpoint is fermion) |
| J tensor (couplings, `operator.tensor_terms`, ops are `S?`) | spin - spin |
| onsite spin operators (`Sx`/`Sy`/`Sz`/`Szz`) | spin, or fermion (as the electron spin) |
| onsite `N`/`Nup`/`Ndn`/`NupNdn` | fermion |

Applicable onsite scope by model (same content as the Global Constraints table):

| model | terms that may appear onsite |
|---|---|
| Spin | magnetic field (`h`, `Gamma`, `Gamma_y`), single-ion anisotropy `D` |
| Hubbard | `mu` (chemical potential), `U`, magnetic field (`h`/`Gamma`/`Gamma_y` for the electron spin) |
| Kondo | the `_c` label carries the full Hubbard set (`mu`, `U`, magnetic field), the `_s` label carries the magnetic field, and both labels carry the magnetic field (§6.7) |

### 6.7 Kondo's two-label convention

The Kondo model represents 1 physical site as **two geometry labels**
(`<X>_c` = the fermion degree of freedom of the itinerant electron, `<X>_s` = the spin
degree of freedom of the localized spin). Both labels share **the same fractional
coordinate**, but physically represent 2 degrees of freedom on 1 site.

- The full Hubbard set (`t`, `U`, `mu`) applies only to the `_c` label.
- The magnetic field (`h`, `Gamma`, `Gamma_y`) applies to **both** the `_c` and `_s`
  labels (per the sign table in §6.2).
- Kondo `J` (`s_i . S_j`) is in the order `_c → _s` (the itinerant one is the first
  endpoint, §6.6).
- **Pyrochlore exception**: the current implementation (common to C/Python) generates J
  in an asymmetric form, "the itinerant site of sublattice 3 ↔ the localized spins of
  all sublattices" (`GeneralJ` in `src/Pyrochlore.c` /
  `python/stdface/lattice/pyrochlore.py::_local`). This catalog faithfully reproduces the
  behavior of the current implementation as-is (`from` is fixed to the itinerant site
  `A3_c`). This may be a bug in the upstream implementation, but since this catalog's
  verification level is "cross-checking against the current implementation's source,"
  the upstream behavior takes precedence. A note will be added to the manual that if the
  upstream is fixed in the future, the catalog will need versioning (a schema version
  bump).

## 7. List of extension dialects (experimental extensions)

Under `dialect: experimental`, this catalog uses the following extensions that the draft
spec does not explicitly define. These are also documented as a proposal to the draft in
manual chapter 6.

1. **fermion site_dof**: expresses a fermionic degree of freedom in the form
   `{fermion: {orbitals: n}}` (distinguished from a spin degree of freedom; see the
   type-compatibility table in §6.6).
2. **Single-site operator vocabulary**: `N`, `Nup`, `Ndn`, `NupNdn`, `Sx`, `Sy`, `Sz`,
   `Szz` and their expanded forms (§6.6).
3. **Named two-body operators**: `hop`, `density-density`, `s_i . S_j`, together with the
   convention that endpoint order is meaningful (§6.6).
4. **param references** `{param, scale, default}` (§6.3). The same general form is used
   for all of `value` / `coeff` / `spin` (`2S`) / `twist` (`phaseN`).
5. **The `catalog:` header** (`schema` / `dialect` / `lattice` / `model`, §1).
6. **(sketch only -- future task) general 4-fermion terms**: a term representation with an
   ordered sequence of creation/annihilation operators and site/orbital bindings. Named
   operators (`hop`, etc.) would then be positioned as an abbreviated notation for this.
   This will be needed in the future for describing the J channel of wannier90 (Hund,
   exchange, pair-hopping), but is out of scope for this catalog
   (design §2 out of scope, manual Section 5.5).

## 8. Format for recording consistency checks and provenance

- The consistency-check values for each YAML entry (`n_sites_uc`, `bonds_per_uc`,
  `coordination`, `min_size_for_check`, etc.) are recorded in
  `lattice_catalog/manifest.yaml`, and the linter (C10, C11) cross-checks them
  automatically.
- Provenance is recorded as a set of "file path + function name + referenced commit hash"
  (a line number may additionally be noted as supplementary information, but the commit
  hash is the unique reference point).

```yaml
source:
  file: python/stdface/lattice/chain_lattice.py
  func: chain
  commit: <hash>
```

- For the manifest-side format definition and the list of keys, see the comment at the top
  of `manifest.yaml`.
