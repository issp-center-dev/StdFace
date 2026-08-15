# examples_experimental — hands-on trial samples

Unlike the main catalog (31 files), the YAML in this directory are
**experimental samples meant to demonstrate the limits of the new format's
(catalog dialect's) expressive power**. The design friction points observed
while writing lattices/models that don't exist in the catalog are
summarized in
[`../feedback_draft_spec.md`](../feedback_draft_spec.md).

## Important: relationship to the linter

Some of these files **intentionally fail to pass the linter** (that failure
is itself an observation). For this reason, this directory is excluded from
the linter's default scan
(`_DISCOVERY_EXCLUDED_DIRS` in `tools/lint_catalog.py`).
To reproduce the individual diagnostics, pass an explicit path:

```bash
python3 lattice_catalog/tools/lint_catalog.py \
    lattice_catalog/examples_experimental/kagome_dm_spin.yaml
```

## File list and key observations

| File | Subject | Key observation |
|---|---|---|
| `ssh_dimer_chain_spin.yaml` | Alternating-bond (dimerized) chain | Geometry is straightforward. Free-form parameter names (e.g. `J_intra`) are rejected by C7 — the closedness of the parameter vocabulary |
| `shastry_sutherland_spin.yaml` | Shastry-Sutherland (draft §7.3) | The bond structure can be written almost isomorphically to the draft. The draft's own isotropic operator `"S_i . S_j"` has no place in the dialect's vocabulary and is rejected by C9 |
| `kagome_dm_spin.yaml` | Kagome with DM interaction | DM can be expressed via antisymmetric tensor components. However, sharing a parameter across components (`Dz` with opposite sign in 2 components) is rejected by C8/C12. Confirms that the source-order-preservation convention is essential for preserving the DM sign |
| `two_orbital_hubbard_chain.yaml` | Two-orbital Hubbard chain | The orbital = separate-label scheme reaches as far as U'/t_ab, but Hund coupling and pair hopping (4-fermion terms) cannot be written with the current vocabulary |
| `square_cylinder_spin.yaml` | Cylindrical boundary (for DMRG) | `boundary: [periodic, open]` is straightforward to write. However, since the C11 consistency check assumes a torus, an open boundary cannot be verified |

Details of the observations for each sample are recorded in that file's
leading comment.
