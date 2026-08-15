#!/usr/bin/env python3
"""lattice_catalog/tools の自動テスト(開発時専用)。

pytest 非依存、assert ベース。`python3 test_tools.py` で完結し、
成功時は最後に `all tools tests passed` を出力する。

対象:
- keyword_inventory.py: canonical 集約・安定ソート・重複なし・
  ソルバーレジストリ件数下限・prime 付き成分の存在・
  model_alias がレジストリ(MODEL_ALIASES + MODEL_ALIASES_HPHI_BOOST)と
  完全一致すること(A1)
- lint_catalog.py: C1-C12 の正例・負例(インライン YAML / dict fixture)。
  C11 min_size_for_check の厳密性(A2)、C2 のスキーマ深部検査(A3)、
  C12 符号規約検査(A4)を含む
- param 参照の意味論(resolve_param): scale×param / default そのまま
- 不正 YAML・null 文書での非例外動作
- トップレベル失敗の非例外動作(A6): 破損 manifest.yaml
  (FatalLintError)、孤立 manifest エントリ検出

開発時依存: lint_catalog.py 経由で PyYAML (``pyyaml``) を必要とする
(開発時専用ツール一式であり、``python/pyproject.toml`` の実行時依存
には含めない — A5)。

Run
---
python3 lattice_catalog/tools/test_tools.py
"""
from __future__ import annotations

import copy
import sys
from pathlib import Path

TOOLS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS_DIR))

import keyword_inventory as ki  # noqa: E402
import lint_catalog as lc  # noqa: E402

_PASS_COUNT = 0


def check(cond: bool, msg: str) -> None:
    global _PASS_COUNT
    if not cond:
        raise AssertionError(f"FAILED: {msg}")
    _PASS_COUNT += 1


def has_prefix(errs: list[str], prefix: str) -> bool:
    return any(e.startswith(prefix) or f"{prefix}" in e for e in errs)


def only_ids(errs: list[str]) -> set[str]:
    """Extract the leading C<N> id from each diagnostic (best-effort)."""
    ids = set()
    for e in errs:
        for tok in e.replace(":", " ").split():
            if tok.startswith("C") and tok[1:].isdigit():
                ids.add(tok)
                break
    return ids


# ===========================================================================
#  keyword_inventory.py
# ===========================================================================

def test_inventory_2S_from_solver_table() -> None:
    entries = ki.build_inventory()
    e = [x for x in entries if x["keyword_canonical"] == "2S" and x["kind"] == "keyword"]
    check(len(e) == 1, "'2S' should appear exactly once as a keyword entry")
    check("HPhi" in e[0]["sources"], "'2S' should be sourced from the HPhi solver plugin table")


def test_inventory_prime_j_component() -> None:
    entries = ki.build_inventory()
    names = {x["keyword_canonical"] for x in entries if x["kind"] == "keyword"}
    check("J0'x" in names, "prime-bearing J component \"J0'x\" must be present as canonical \"J0'x\"")
    check("J0''xy" in names, "double-prime J component \"J0''xy\" must be present")


def test_inventory_scalar_j_canonicalised_uppercase() -> None:
    """The bare Kondo scalar exchange keyword ``j`` must canonicalise to
    ``J`` (matching its sibling ``v`` -> ``V``; both are single-letter
    scalar members of the ``J``/``V`` families and appear as ``"J"`` /
    ``"V"`` in the underlying keyword tables, e.g.
    ``_j_matrix_keywords("j", "JAll", "J")``). Regression for a Task 1
    finding: ``_CANON_EXCEPTIONS`` had ``"v": "V"`` but was missing the
    analogous ``"j": "J"`` entry, so catalog ``{param: J}`` references
    (chain/chain_kondo.yaml) failed C7 against the (wrongly lower-cased)
    inventory.
    """
    entries = ki.build_inventory()
    names = {x["keyword_canonical"] for x in entries if x["kind"] == "keyword"}
    check("J" in names, "scalar keyword 'j' must canonicalise to 'J'")
    check("j" not in names, "lower-case 'j' must not remain as a separate canonical keyword")


def test_inventory_stable_sort_no_dup() -> None:
    e1 = ki.build_inventory()
    e2 = ki.build_inventory()
    check(e1 == e2, "build_inventory() must be deterministic/stably sorted across runs")
    keys = [(e["keyword_canonical"], e["kind"]) for e in e1]
    check(len(keys) == len(set(keys)), "no duplicate (keyword_canonical, kind) pairs allowed")
    check(keys == sorted(keys), "entries must be sorted by (keyword_canonical, kind)")


def test_inventory_solver_plugin_count() -> None:
    names = [name for name, _table in ki.iter_solver_plugins()]
    check(len(names) >= 4, f"expected >=4 distinct solver plugins, got {len(names)}: {names}")
    check(len(names) == len(set(names)), "solver plugin names must be de-duplicated by identity")


def test_inventory_lattice_and_model_aliases() -> None:
    entries = ki.build_inventory()
    lat = {x["keyword_canonical"] for x in entries if x["kind"] == "lattice_alias"}
    mdl = {x["keyword_canonical"] for x in entries if x["kind"] == "model_alias"}
    check("chain" in lat, "lattice alias 'chain' must be present")
    check({"spin", "hubbard", "kondo"} <= mdl, "model aliases spin/hubbard/kondo must be present")


def test_inventory_model_alias_set_matches_registry() -> None:
    """A1: model_alias must be generated from the *actual* alias registry
    (``MODEL_ALIASES`` + ``MODEL_ALIASES_HPHI_BOOST`` in
    ``stdface.core.stdface_main``), not a hardcoded {spin, hubbard, kondo}
    set of the 3 canonical ``ModelType`` values.
    """
    entries = ki.build_inventory()
    mdl = {x["keyword_canonical"] for x in entries if x["kind"] == "model_alias"}
    expected = set(ki._MODEL_ALIASES) | set(ki._MODEL_ALIASES_HPHI_BOOST)
    check(mdl == expected,
          f"inventory model_alias set must equal the registry-derived set, "
          f"got {mdl}, expected {expected}")
    check(len(expected) > 3,
          "the real alias registry has more than the 3 canonical model names "
          f"(GC/Boost variants included); got {expected}")

    targets = {x["keyword_canonical"]: x["canonical_target"]
               for x in entries if x["kind"] == "model_alias"}
    for alias, cfg in ki._MODEL_ALIASES.items():
        check(targets[alias] == str(cfg.model.value),
              f"model_alias {alias!r} canonical_target must be {cfg.model.value!r}, "
              f"got {targets[alias]!r}")
    for alias, cfg in ki._MODEL_ALIASES_HPHI_BOOST.items():
        check(targets[alias] == str(cfg.model.value),
              f"model_alias {alias!r} (HPhi Boost) canonical_target must be "
              f"{cfg.model.value!r}, got {targets[alias]!r}")

    sources = {x["keyword_canonical"]: set(x["sources"])
               for x in entries if x["kind"] == "model_alias"}
    for alias in ki._MODEL_ALIASES_HPHI_BOOST:
        check("model_registry_hphi_boost" in sources[alias],
              f"HPhi-Boost-only alias {alias!r} must be sourced as "
              "'model_registry_hphi_boost'")


# ===========================================================================
#  resolve_param semantics (C7-pinned)
# ===========================================================================

def test_resolve_param_semantics() -> None:
    ref = {"param": "2S", "scale": 0.5, "default": 0.5}
    check(lc.resolve_param(ref, 1) == 0.5, "param given: value = scale * param")
    check(lc.resolve_param(ref, None) == 0.5, "param absent: value = default AS-IS (no scale applied)")

    ref2 = {"param": "t0", "scale": -1.0}  # default omitted -> 0
    check(lc.resolve_param(ref2, 2.0) == -2.0, "scale * param with explicit scale")
    check(lc.resolve_param(ref2, None) == 0, "default omitted -> 0, scale not applied")

    ref3 = {"param": "V0"}  # scale omitted -> 1.0, default omitted -> 0
    check(lc.resolve_param(ref3, 3.0) == 3.0, "scale defaults to 1.0")
    check(lc.resolve_param(ref3, None) == 0, "default defaults to 0")


# ===========================================================================
#  lint_catalog.py fixtures
# ===========================================================================

INVENTORY = ki.build_inventory()
KEYWORDS = {e["keyword_canonical"] for e in INVENTORY if e["kind"] == "keyword"}

_J0_TERMS = [
    {"ops": ["Sx", "Sx"], "coeff": {"param": "J0x"}},
    {"ops": ["Sx", "Sy"], "coeff": {"param": "J0xy"}},
    {"ops": ["Sx", "Sz"], "coeff": {"param": "J0xz"}},
    {"ops": ["Sy", "Sy"], "coeff": {"param": "J0y"}},
    {"ops": ["Sy", "Sx"], "coeff": {"param": "J0yx"}},
    {"ops": ["Sy", "Sz"], "coeff": {"param": "J0yz"}},
    {"ops": ["Sz", "Sz"], "coeff": {"param": "J0z"}},
    {"ops": ["Sz", "Sx"], "coeff": {"param": "J0zx"}},
    {"ops": ["Sz", "Sy"], "coeff": {"param": "J0zy"}},
]


def base_doc() -> dict:
    """A minimal, fully-valid chain/spin catalog document (1 bond type)."""
    return {
        "catalog": {
            "schema": "stdface-catalog/0.1",
            "dialect": "experimental",
            "lattice": "chain",
            "model": "spin",
        },
        "geometry": {
            "dimension": 1,
            "lattice_vectors": {"a1": [1.0]},
            "sites": [{"label": "A", "frac": [0.0]}],
        },
        "system": {
            "size": [7],
            "boundary": [{"twist": {"param": "phase0"}}],
        },
        "model": {
            "site_dof": {
                "A": {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}},
            },
            "bonds": [
                {"type": "J0", "from": "A", "to": "A", "R": [1]},
            ],
            "couplings": {
                "J0": {"operator": {"tensor_terms": copy.deepcopy(_J0_TERMS)}},
            },
            "onsite": {
                "A": {
                    "aniso_z": {
                        "operator": {"tensor_terms": [{"ops": ["Szz"], "coeff": 1.0}]},
                        "value": {"param": "D", "scale": 1.0, "default": 0},
                    },
                },
            },
        },
    }


def base_manifest() -> dict:
    return {
        "chain/chain_spin.yaml": {
            "lattice": "chain",
            "model": "spin",
            "dimension": 1,
            "n_sites_uc": 1,
            "bonds_per_uc": {"J0": 1},
            "coordination": {"A": {"J0": 2}},
            "min_size_for_check": [3],  # = 2*max|R|+1 = 2*1+1 (bond R=[1])
            "source": {
                "file": "python/stdface/lattice/chain_lattice.py",
                "func": "chain",
                "commit": "deadbeef",
            },
        }
    }


REL = "chain/chain_spin.yaml"


def lint(doc: dict, manifest: dict | None = None) -> list[str]:
    m = manifest if manifest is not None else base_manifest()
    return lc.lint_document(doc, REL, KEYWORDS, m)


# ---------------------------------------------------------------------------
#  baseline sanity
# ---------------------------------------------------------------------------

def test_base_doc_is_valid() -> None:
    errs = lint(base_doc())
    check(errs == [], f"base_doc() should be fully valid, got: {errs}")


# ---------------------------------------------------------------------------
#  C1
# ---------------------------------------------------------------------------

def test_c1_positive() -> None:
    errs = lint(base_doc())
    check("C1" not in only_ids(errs), "valid catalog header must not raise C1")


def test_c1_negative_schema() -> None:
    d = base_doc()
    d["catalog"]["schema"] = "wrong/0.0"
    errs = lint(d)
    check("C1" in only_ids(errs), f"bad schema must raise C1, got {errs}")


def test_c1_negative_missing_dialect() -> None:
    d = base_doc()
    del d["catalog"]["dialect"]
    errs = lint(d)
    check("C1" in only_ids(errs), f"missing dialect must raise C1, got {errs}")


def test_c1_negative_missing_lattice_model() -> None:
    d = base_doc()
    del d["catalog"]["lattice"]
    del d["catalog"]["model"]
    errs = lint(d)
    check("C1" in only_ids(errs), f"missing lattice/model must raise C1, got {errs}")


# ---------------------------------------------------------------------------
#  C2
# ---------------------------------------------------------------------------

def test_c2_positive() -> None:
    errs = lint(base_doc())
    check("C2" not in only_ids(errs), "valid schema must not raise C2")


def test_c2_negative_missing_section() -> None:
    d = base_doc()
    del d["system"]
    errs = lint(d)
    check("C2" in only_ids(errs), f"missing system section must raise C2, got {errs}")


def test_c2_negative_bad_dimension() -> None:
    d = base_doc()
    d["geometry"]["dimension"] = 0
    errs = lint(d)
    check("C2" in only_ids(errs), f"dimension=0 must raise C2, got {errs}")


def test_c2_negative_duplicate_labels() -> None:
    d = base_doc()
    d["geometry"]["sites"].append({"label": "A", "frac": [0.5]})
    d["model"]["site_dof"]["A"] = {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}}
    errs = lint(d)
    check("C2" in only_ids(errs), f"duplicate site labels must raise C2, got {errs}")


def test_c2_negative_non_int_size() -> None:
    d = base_doc()
    d["system"]["size"] = [7.5]
    errs = lint(d)
    check("C2" in only_ids(errs), f"non-integer size must raise C2, got {errs}")


def test_c2_negative_lattice_vectors_key_count_mismatch() -> None:
    d = base_doc()
    d["geometry"]["lattice_vectors"] = {"a1": [1.0], "a2": [0.0, 1.0]}
    errs = lint(d)
    check("C2" in only_ids(errs),
          f"lattice_vectors with wrong key count for dimension=1 must raise C2, got {errs}")


def test_c2_negative_lattice_vectors_dim_mismatch() -> None:
    d = base_doc()
    d["geometry"]["lattice_vectors"] = {"a1": [1.0, 0.0]}  # dimension=1 expects length-1 vectors
    errs = lint(d)
    check("C2" in only_ids(errs),
          f"lattice_vectors[a1] with wrong vector length must raise C2, got {errs}")


def test_c2_negative_frac_dim_mismatch() -> None:
    d = base_doc()
    d["geometry"]["sites"][0]["frac"] = [0.0, 0.0]  # dimension=1 expects length-1 frac
    errs = lint(d)
    check("C2" in only_ids(errs), f"frac with wrong dimension must raise C2, got {errs}")


def test_c2_negative_non_string_label() -> None:
    d = base_doc()
    d["geometry"]["sites"][0]["label"] = 42
    errs = lint(d)
    check("C2" in only_ids(errs), f"non-string site label must raise C2, got {errs}")


def test_c2_negative_boundary_length_mismatch() -> None:
    d = base_doc()
    d["system"]["boundary"] = [{"twist": {"param": "phase0"}}, "periodic"]  # dimension=1 expects length 1
    errs = lint(d)
    check("C2" in only_ids(errs), f"boundary length != dimension must raise C2, got {errs}")


def test_c2_negative_site_dof_both_spin_and_fermion() -> None:
    d = base_doc()
    d["model"]["site_dof"]["A"] = {
        "spin": {"param": "2S", "scale": 0.5, "default": 0.5},
        "fermion": {"orbitals": 1},
    }
    errs = lint(d)
    check("C2" in only_ids(errs),
          f"site_dof entry with both spin and fermion must raise C2, got {errs}")


def test_c2_negative_mixed_type_duplicate_labels_no_crash() -> None:
    """C2 robustness: a duplicated *and* heterogeneous-type labels list
    (two string 'A' duplicates plus two int `5` duplicates) must not crash
    `sorted()` on mixed types and fall through to the generic top-level
    'exception while checking' diagnostic -- the string duplicate must
    still be reported cleanly."""
    d = base_doc()
    d["geometry"]["sites"] = [
        {"label": "A", "frac": [0.0]},
        {"label": "A", "frac": [0.0]},
        {"label": 5, "frac": [0.0]},
        {"label": 5, "frac": [0.0]},
    ]
    errs = lint(d)
    check("C2" in only_ids(errs), f"mixed-type duplicate labels must raise C2, got {errs}")
    check(not any("exception while checking" in e for e in errs),
          f"mixed-type duplicate labels must not fall through to a generic exception, got {errs}")
    check(any("duplicate site labels" in e and "'A'" in e for e in errs),
          f"the string duplicate 'A' must still be reported, got {errs}")


def test_c2_negative_non_string_label_sanitized_downstream_no_crash() -> None:
    """Sanitization: a non-string site label (5) must be diagnosed by C2
    and then dropped from the context passed to C4-C12 -- it must not
    reach `sorted()` on a mixed str/int set (labels vs. site_dof keys) in
    check_c4, which previously raised TypeError and collapsed all
    diagnostics into a single generic 'exception while checking' entry.
    """
    d = base_doc()
    d["geometry"]["sites"] = [{"label": "B", "frac": [0.0]}, {"label": 5, "frac": [0.0]}]
    d["model"]["site_dof"] = {"A": {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}}}
    errs = lint(d)
    check("C2" in only_ids(errs), f"non-string label must raise C2, got {errs}")
    check("C4" in only_ids(errs), f"sites/site_dof mismatch must raise C4, got {errs}")
    check(not any("exception while checking" in e for e in errs),
          f"must not fall through to a generic exception diagnostic, got {errs}")


def test_c2_negative_malformed_bond_sanitized_downstream_no_crash() -> None:
    """Sanitization: a bond with a non-string `type` (a list) and a
    non-string `from` (a dict) must be diagnosed by C2 and then dropped
    from the bonds passed to C4-C12 -- it must not reach the unhashable
    set-building in check_c5 (`{b["type"] for b in bonds}`) or the
    unhashable tuple-keying in check_c6's reversal_dup(), either of which
    previously raised TypeError and collapsed all diagnostics into a
    single generic 'exception while checking' entry.
    """
    d = base_doc()
    d["model"]["bonds"][0]["type"] = ["J0"]
    d["model"]["bonds"][0]["from"] = {"x": 1}
    errs = lint(d)
    check("C2" in only_ids(errs), f"malformed bond fields must raise C2, got {errs}")
    check(not any("exception while checking" in e for e in errs),
          f"must not fall through to a generic exception diagnostic, got {errs}")


# ---------------------------------------------------------------------------
#  C3
# ---------------------------------------------------------------------------

def test_c3_positive() -> None:
    errs = lint(base_doc())
    check("C3" not in only_ids(errs), "matching R/dimension/size must not raise C3")


def test_c3_negative_bad_R_length() -> None:
    d = base_doc()
    d["model"]["bonds"][0]["R"] = [1, 0]
    errs = lint(d)
    check("C3" in only_ids(errs), f"R length mismatch must raise C3, got {errs}")


def test_c3_negative_bad_size_length() -> None:
    d = base_doc()
    d["system"]["size"] = [7, 7]
    errs = lint(d)
    check("C3" in only_ids(errs), f"size length != dimension must raise C3, got {errs}")


# ---------------------------------------------------------------------------
#  C4
# ---------------------------------------------------------------------------

def test_c4_positive() -> None:
    errs = lint(base_doc())
    check("C4" not in only_ids(errs), "consistent labels must not raise C4")


def test_c4_negative_bond_unknown_label() -> None:
    d = base_doc()
    d["model"]["bonds"][0]["to"] = "B"
    errs = lint(d)
    check("C4" in only_ids(errs), f"bond referencing unknown label must raise C4, got {errs}")


def test_c4_negative_onsite_unknown_label() -> None:
    d = base_doc()
    d["model"]["onsite"] = {"Z": d["model"]["onsite"].pop("A")}
    errs = lint(d)
    check("C4" in only_ids(errs), f"onsite referencing unknown label must raise C4, got {errs}")


def test_c4_negative_site_dof_mismatch() -> None:
    d = base_doc()
    d["model"]["site_dof"]["Extra"] = {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}}
    errs = lint(d)
    check("C4" in only_ids(errs), f"site_dof/sites mismatch must raise C4, got {errs}")


def test_c3_c4_negative_malformed_bond_preserves_original_index() -> None:
    """Regression: C2 sanitization drops malformed bonds from ctx["bonds"],
    which must not renumber the bonds that follow. A malformed bond at
    index 0 (type is a list) is dropped; the retained bond at index 1
    fails both C3 (bad R length) and C4 (undefined endpoint) and must
    still be cited as model.bonds[1], not model.bonds[0], by both."""
    d = base_doc()
    d["model"]["bonds"] = [
        {"type": ["J0"], "from": "A", "to": "A", "R": [1]},  # malformed: type is a list -> dropped by C2
        {"type": "J0", "from": "A", "to": "Z", "R": [1, 0]},  # retained: bad R length + undefined "to" label
    ]
    errs = lint(d)
    check("C2" in only_ids(errs), f"malformed bond at index 0 must raise C2, got {errs}")
    check(any(e.startswith("C3:") and "model.bonds[1]" in e for e in errs),
          f"C3 must cite model.bonds[1] for the retained bond, got {errs}")
    check(any(e.startswith("C4:") and "model.bonds[1]" in e for e in errs),
          f"C4 must cite model.bonds[1] for the retained bond, got {errs}")
    check(not any("model.bonds[0]" in e for e in errs if e.startswith("C3:") or e.startswith("C4:")),
          f"C3/C4 must not misreport the retained bond as model.bonds[0], got {errs}")
    check(not any("exception while checking" in e for e in errs),
          f"must not fall through to a generic exception diagnostic, got {errs}")


# ---------------------------------------------------------------------------
#  C5
# ---------------------------------------------------------------------------

def test_c5_positive() -> None:
    errs = lint(base_doc())
    check("C5" not in only_ids(errs), "matching bond type / couplings keys must not raise C5")


def test_c5_negative_undefined_type() -> None:
    d = base_doc()
    d["model"]["bonds"].append({"type": "J1", "from": "A", "to": "A", "R": [2]})
    errs = lint(d)
    check("C5" in only_ids(errs), f"bond type without couplings def must raise C5, got {errs}")


def test_c5_negative_unused_definition() -> None:
    d = base_doc()
    d["model"]["couplings"]["J1"] = {"tensor_terms": copy.deepcopy(_J0_TERMS)}
    errs = lint(d)
    check("C5" in only_ids(errs), f"unused couplings def must raise C5, got {errs}")


# ---------------------------------------------------------------------------
#  C6
# ---------------------------------------------------------------------------

def test_c6_positive_t_v_coexist() -> None:
    """Same geometric bond, different types (t0 & V0) -- legitimate, no C6."""
    d = base_doc()
    d["model"]["site_dof"]["A"] = {"fermion": {"orbitals": 1}}
    d["model"]["bonds"] = [
        {"type": "t0", "from": "A", "to": "A", "R": [1]},
        {"type": "V0", "from": "A", "to": "A", "R": [1]},
    ]
    d["model"]["couplings"] = {
        "t0": {"operator": "hop", "value": {"param": "t0", "scale": -1.0}},
        "V0": {"operator": "density-density", "value": {"param": "V0", "scale": 1.0}},
    }
    d["model"]["onsite"] = {}
    m = {REL: {**base_manifest()[REL], "bonds_per_uc": {"t0": 1, "V0": 1},
               "coordination": {"A": {"t0": 2, "V0": 2}}}}
    errs = lint(d, m)
    check("C6" not in only_ids(errs), f"t/V coexistence on same geometric bond must not raise C6, got {errs}")


def test_c6_negative_reversal_duplicate_same_type() -> None:
    d = base_doc()
    d["model"]["bonds"].append({"type": "J0", "from": "A", "to": "A", "R": [-1]})
    errs = lint(d)
    check("C6" in only_ids(errs), f"reversal-equivalent duplicate (same type) must raise C6, got {errs}")


# ---------------------------------------------------------------------------
#  C7
# ---------------------------------------------------------------------------

def test_c7_positive() -> None:
    errs = lint(base_doc())
    check("C7" not in only_ids(errs), "valid param references must not raise C7")


def test_c7_negative_unknown_param() -> None:
    d = base_doc()
    d["model"]["onsite"]["A"]["aniso_z"]["value"]["param"] = "NotAKeyword"
    errs = lint(d)
    check("C7" in only_ids(errs), f"unknown param name must raise C7, got {errs}")


def test_c7_negative_bad_scale_type() -> None:
    d = base_doc()
    d["model"]["onsite"]["A"]["aniso_z"]["value"]["scale"] = "big"
    errs = lint(d)
    check("C7" in only_ids(errs), f"non-numeric scale must raise C7, got {errs}")


def test_c7_negative_bad_default_type() -> None:
    d = base_doc()
    d["model"]["onsite"]["A"]["aniso_z"]["value"]["default"] = "zero"
    errs = lint(d)
    check("C7" in only_ids(errs), f"non-numeric default must raise C7, got {errs}")


# ---------------------------------------------------------------------------
#  C8
# ---------------------------------------------------------------------------

def test_c8_positive() -> None:
    errs = lint(base_doc())
    check("C8" not in only_ids(errs), "9 well-formed tensor terms must not raise C8")


def test_c8_negative_missing_term() -> None:
    d = base_doc()
    d["model"]["couplings"]["J0"]["operator"]["tensor_terms"].pop()
    errs = lint(d)
    check("C8" in only_ids(errs), f"8 terms (missing one) must raise C8, got {errs}")


def test_c8_negative_duplicate_ops() -> None:
    d = base_doc()
    d["model"]["couplings"]["J0"]["operator"]["tensor_terms"][-1] = {
        "ops": ["Sx", "Sx"], "coeff": {"param": "J0x"},
    }
    errs = lint(d)
    check("C8" in only_ids(errs), f"duplicate ops pair must raise C8, got {errs}")


def test_c8_negative_extra_term() -> None:
    d = base_doc()
    d["model"]["couplings"]["J0"]["operator"]["tensor_terms"].append(
        {"ops": ["Sx", "Sx"], "coeff": {"param": "J0x"}}
    )
    errs = lint(d)
    check("C8" in only_ids(errs), f"10 terms must raise C8, got {errs}")


def test_c8_negative_bad_param_name() -> None:
    d = base_doc()
    d["model"]["couplings"]["J0"]["operator"]["tensor_terms"][0]["coeff"]["param"] = "J0wrong"
    errs = lint(d)
    check("C8" in only_ids(errs), f"ops<->param suffix mismatch must raise C8, got {errs}")


# ---------------------------------------------------------------------------
#  C9
# ---------------------------------------------------------------------------

def _fermion_hop_doc() -> dict:
    """A fully-valid fermion-hop fixture (catalog.model consistent with its
    manifest, no double sign inversion — A6-iii)."""
    d = base_doc()
    d["catalog"]["model"] = "hubbard"
    d["model"]["site_dof"] = {"A": {"fermion": {"orbitals": 1}}}
    d["model"]["bonds"] = [{"type": "t0", "from": "A", "to": "A", "R": [1]}]
    d["model"]["couplings"] = {"t0": {"operator": "hop", "value": {"param": "t0", "scale": -1.0}}}
    d["model"]["onsite"] = {
        "A": {
            "chemical_potential": {
                "operator": {"tensor_terms": [{"ops": ["N"], "coeff": -1.0}]},
                "value": {"param": "mu"},
            },
        },
    }
    return d


def _fermion_hop_manifest() -> dict:
    return {REL: {"lattice": "chain", "model": "hubbard", "dimension": 1, "n_sites_uc": 1,
                  "bonds_per_uc": {"t0": 1}, "coordination": {"A": {"t0": 2}},
                  "min_size_for_check": [3],  # = 2*max|R|+1 = 2*1+1 (bond R=[1])
                  "source": {"file": "x", "func": "y", "commit": "z"}}}


def test_c9_positive_fermion_fermion_hop() -> None:
    errs = lint(_fermion_hop_doc(), _fermion_hop_manifest())
    check(errs == [], f"hop on fermion-fermion must be fully valid, got {errs}")


def test_c9_negative_hop_on_spin() -> None:
    d = _fermion_hop_doc()
    d["model"]["site_dof"]["A"] = {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}}
    errs = lint(d, _fermion_hop_manifest())
    check("C9" in only_ids(errs), f"hop on spin-spin must raise C9, got {errs}")


def test_c9_positive_kondo_fermion_spin() -> None:
    d = base_doc()
    d["catalog"]["model"] = "kondo"
    d["model"]["site_dof"] = {
        "A_c": {"fermion": {"orbitals": 1}},
        "A_s": {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}},
    }
    d["geometry"]["sites"] = [{"label": "A_c", "frac": [0.0]}, {"label": "A_s", "frac": [0.0]}]
    d["model"]["bonds"] = [{"type": "J", "from": "A_c", "to": "A_s", "R": [0]}]
    d["model"]["couplings"] = {"J": {"operator": "s_i . S_j", "value": {"param": "J", "scale": 1.0}}}
    d["model"]["onsite"] = {}
    m = {REL: {"lattice": "chain", "model": "kondo", "dimension": 1, "n_sites_uc": 2,
               "bonds_per_uc": {"J": 1},
               "coordination": {"A_c": {"J": 1}, "A_s": {"J": 1}},
               "min_size_for_check": [1],  # = 2*max|R|+1 = 2*0+1 (bond R=[0])
               "source": {"file": "x", "func": "y", "commit": "z"}}}
    errs = lint(d, m)
    check(errs == [], f"s_i.S_j fermion(1st)-spin(2nd) must be fully valid, got {errs}")


def test_c9_negative_kondo_wrong_order() -> None:
    d = base_doc()
    d["catalog"]["model"] = "kondo"
    d["model"]["site_dof"] = {
        "A_c": {"fermion": {"orbitals": 1}},
        "A_s": {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}},
    }
    d["geometry"]["sites"] = [{"label": "A_c", "frac": [0.0]}, {"label": "A_s", "frac": [0.0]}]
    d["model"]["bonds"] = [{"type": "J", "from": "A_s", "to": "A_c", "R": [0]}]
    d["model"]["couplings"] = {"J": {"operator": "s_i . S_j", "value": {"param": "J", "scale": 1.0}}}
    d["model"]["onsite"] = {}
    m = {REL: {"lattice": "chain", "model": "kondo", "dimension": 1, "n_sites_uc": 2,
               "bonds_per_uc": {"J": 1},
               "coordination": {"A_c": {"J": 1}, "A_s": {"J": 1}},
               "min_size_for_check": [1],  # = 2*max|R|+1 = 2*0+1 (bond R=[0])
               "source": {"file": "x", "func": "y", "commit": "z"}}}
    errs = lint(d, m)
    check("C9" in only_ids(errs), f"s_i.S_j with spin-first/fermion-second order must raise C9, got {errs}")


def test_c9_negative_j_tensor_on_fermion() -> None:
    d = base_doc()
    d["model"]["site_dof"]["A"] = {"fermion": {"orbitals": 1}}
    errs = lint(d)
    check("C9" in only_ids(errs), f"J tensor on fermion-fermion must raise C9, got {errs}")


def test_c9_negative_onsite_type_mismatch() -> None:
    d = base_doc()  # site A is spin
    d["model"]["onsite"]["A"]["aniso_z"] = {
        "operator": {"tensor_terms": [{"ops": ["N"], "coeff": -1.0}]},
        "value": {"param": "mu", "scale": -1.0},
    }
    errs = lint(d)
    check("C9" in only_ids(errs), f"N onsite op on spin site_dof must raise C9, got {errs}")


# ---------------------------------------------------------------------------
#  C10
# ---------------------------------------------------------------------------

def test_c10_positive() -> None:
    errs = lint(base_doc())
    check("C10" not in only_ids(errs), "matching manifest must not raise C10")


def test_c10_negative_missing_key() -> None:
    m = base_manifest()
    del m[REL]["source"]
    errs = lint(base_doc(), m)
    check("C10" in only_ids(errs), f"manifest entry missing a required key must raise C10, got {errs}")


def test_c10_negative_dimension_mismatch() -> None:
    m = base_manifest()
    m[REL]["dimension"] = 2
    errs = lint(base_doc(), m)
    check("C10" in only_ids(errs), f"dimension mismatch must raise C10, got {errs}")


def test_c10_negative_bonds_per_uc_mismatch() -> None:
    m = base_manifest()
    m[REL]["bonds_per_uc"] = {"J0": 2}
    errs = lint(base_doc(), m)
    check("C10" in only_ids(errs), f"bonds_per_uc mismatch must raise C10, got {errs}")


def test_c10_negative_no_entry() -> None:
    errs = lint(base_doc(), {})
    check("C10" in only_ids(errs), f"missing manifest entry entirely must raise C10, got {errs}")


# ---------------------------------------------------------------------------
#  C11
# ---------------------------------------------------------------------------

def test_c11_positive() -> None:
    errs = lint(base_doc())
    check("C11" not in only_ids(errs), "correct coordination must not raise C11")


def test_c11_positive_expand_and_count_direct() -> None:
    """C11: expand_and_count() genuinely expands the torus (not decorative).

    chain nn bond A-A, R=[1], on a size=[7] ring: each site instance has
    exactly one bond as 'from' and one as 'to' -> coordination 2, and no
    folding-degeneracy errors (7 cells * 1 bond = 7 distinct instance-pairs,
    which is exactly what a 1D ring of length 7 with a nn bond gives).
    """
    bonds = [{"type": "J0", "from": "A", "to": "A", "R": [1]}]
    measured, errs = lc.expand_and_count(1, ["A"], bonds, [7])
    check(measured == {"A": {"J0": 2}}, f"expected coordination A:{{J0:2}}, got {measured}")
    check(errs == [], f"valid min_size_for_check must not raise C11 expansion errors, got {errs}")


def test_c11_negative_wrong_coordination() -> None:
    m = base_manifest()
    m[REL]["coordination"] = {"A": {"J0": 4}}
    errs = lint(base_doc(), m)
    check("C11" in only_ids(errs), f"wrong coordination must raise C11, got {errs}")


def test_c11_negative_folding_degeneracy() -> None:
    """C11: an artificially small min_size_for_check must be caught.

    The same chain nn bond (A-A, R=[1]) on a size=[2] ring folds: the two
    directed bonds (cell 0 -> cell 1) and (cell 1 -> cell 0, i.e. cell 1 ->
    cell 0 via R=1 wrap) land on the *same* unordered instance pair
    {(cell 0, A), (cell 1, A)}. This is exactly the min_size_for_check-too-
    small case the field exists to guard against, so expand_and_count()
    must report it directly.
    """
    bonds = [{"type": "J0", "from": "A", "to": "A", "R": [1]}]
    measured, errs = lc.expand_and_count(1, ["A"], bonds, [2])
    check(any("C11" in e and "folding" in e for e in errs),
          f"min_size_for_check=[2] must raise a C11 folding-degeneracy error, got {errs}")

    # ...and the same must hold end-to-end through lint_document() when the
    # manifest declares an undersized min_size_for_check.
    m = base_manifest()
    m[REL]["min_size_for_check"] = [2]
    errs = lint(base_doc(), m)
    check("C11" in only_ids(errs),
          f"undersized manifest.min_size_for_check must raise C11 via lint(), got {errs}")


# ---------------------------------------------------------------------------
#  C11 (A2): min_size_for_check strict validation
# ---------------------------------------------------------------------------

def test_c11_negative_min_size_non_int() -> None:
    m = base_manifest()
    m[REL]["min_size_for_check"] = [3.5]
    errs = lint(base_doc(), m)
    check("C11" in only_ids(errs), f"non-int min_size_for_check must raise a clean C11 diagnostic, got {errs}")


def test_c11_negative_min_size_zero() -> None:
    m = base_manifest()
    m[REL]["min_size_for_check"] = [0]
    errs = lint(base_doc(), m)
    check("C11" in only_ids(errs), f"min_size_for_check=0 must raise a clean C11 diagnostic, got {errs}")


def test_c11_negative_min_size_negative() -> None:
    m = base_manifest()
    m[REL]["min_size_for_check"] = [-3]
    errs = lint(base_doc(), m)
    check("C11" in only_ids(errs), f"negative min_size_for_check must raise a clean C11 diagnostic, got {errs}")


def test_c11_negative_min_size_wrong_equality() -> None:
    """A positive odd value that is nonetheless not 2*max|R|+1 must be
    rejected (base_doc()'s single bond has R=[1] -> expected [3], not [5])."""
    m = base_manifest()
    m[REL]["min_size_for_check"] = [5]
    errs = lint(base_doc(), m)
    check("C11" in only_ids(errs),
          f"min_size_for_check not equal to 2*max|R|+1 must raise C11, got {errs}")


def test_c11_manifest_min_size_matches_all_catalog_files() -> None:
    """A2: every real manifest.yaml entry's min_size_for_check must already
    satisfy the strict equality (2*max|R component|+1 per direction,
    computed from that file's actual model.bonds) -- i.e. linting the real
    catalog must not raise any C11 min_size_for_check diagnostic."""
    manifest_files = lc.load_manifest()
    check(len(manifest_files) > 0, "real manifest.yaml must be loadable and non-empty")
    for path in lc._discover_default_files():
        rel_path = path.resolve().relative_to(lc.ROOT).as_posix()
        errs = lc.lint_file(path, KEYWORDS, manifest_files)
        min_size_errs = [e for e in errs if e.startswith("C11") and "min_size_for_check" in e]
        check(min_size_errs == [],
              f"{rel_path}: min_size_for_check must already satisfy the strict "
              f"A2 equality, got {min_size_errs}")


# NOTE: a per-label instance mismatch (translation non-invariance) is not
# reachable through the public expand_and_count() interface: for a fixed
# bond list, every cell contributes exactly one 'from'-incidence to
# (cell, b['from']) and exactly one 'to'-incidence to ((cell+R) mod size,
# b['to']) per bond, and cell -> (cell+R) mod size is a bijection of the
# torus for any fixed R. So the per-instance touch count for a given label
# is a sum of per-bond contributions that is, by construction, identical
# for every cell -- independent of any folding collisions. Folding degrades
# the *distinctness* of expanded bonds (covered above), not the coordination
# sum. There is therefore no way to reach the "mismatch" branch by
# constructing a bonds/size input; it is retained purely as a defensive
# invariant check. Per the task instructions we skip a dedicated negative
# test for it (unreachable by construction).


# ---------------------------------------------------------------------------
#  C12 (A4): sign-convention checks (CONVENTIONS.md §6.2)
# ---------------------------------------------------------------------------

def test_c12_positive() -> None:
    errs = lint(base_doc())
    check("C12" not in only_ids(errs), f"correctly-signed base_doc() must not raise C12, got {errs}")


def test_c12_positive_fermion_hop() -> None:
    errs = lint(_fermion_hop_doc(), _fermion_hop_manifest())
    check("C12" not in only_ids(errs), f"correctly-signed hop fixture must not raise C12, got {errs}")


def test_c12_negative_inverted_hop_scale() -> None:
    d = _fermion_hop_doc()
    d["model"]["couplings"]["t0"]["value"]["scale"] = 1.0  # should be -1.0
    errs = lint(d, _fermion_hop_manifest())
    check("C12" in only_ids(errs), f"hop with scale=+1.0 (inverted) must raise C12, got {errs}")


def test_c12_negative_doubled_sign_density_density() -> None:
    """density-density with an extra scale=-1.0 is a doubled sign inversion
    relative to the CONVENTIONS §6.2 table (scale absent or +1.0)."""
    d = base_doc()
    d["model"]["site_dof"]["A"] = {"fermion": {"orbitals": 1}}
    d["model"]["bonds"] = [{"type": "V0", "from": "A", "to": "A", "R": [1]}]
    d["model"]["couplings"] = {
        "V0": {"operator": "density-density", "value": {"param": "V0", "scale": -1.0}},
    }
    d["model"]["onsite"] = {}
    m = {REL: {**base_manifest()[REL], "bonds_per_uc": {"V0": 1},
               "coordination": {"A": {"V0": 2}}}}
    errs = lint(d, m)
    check("C12" in only_ids(errs),
          f"density-density with scale=-1.0 (doubled sign) must raise C12, got {errs}")


def test_c12_negative_j_coeff_with_scale() -> None:
    d = base_doc()
    d["model"]["couplings"]["J0"]["operator"]["tensor_terms"][0]["coeff"] = {
        "param": "J0x", "scale": 2.0,
    }
    errs = lint(d)
    check("C12" in only_ids(errs), f"J tensor_terms coeff with a scale key must raise C12, got {errs}")


def test_c12_negative_missing_value() -> None:
    d = _fermion_hop_doc()
    del d["model"]["couplings"]["t0"]["value"]
    errs = lint(d, _fermion_hop_manifest())
    check("C12" in only_ids(errs), f"hop coupling missing 'value' must raise C12, got {errs}")


def test_c12_wannier90_plain_number_exempt() -> None:
    """wannier90-dialect couplings may use a plain numeric value with no
    scale requirement enforced (CONVENTIONS §6.3 Rule W1 exception) --
    must not raise C12."""
    d = base_doc()
    d["catalog"]["lattice"] = "wannier90"
    d["model"]["site_dof"]["A"] = {"fermion": {"orbitals": 1}}
    d["model"]["bonds"] = [{"type": "t0", "from": "A", "to": "A", "R": [1]}]
    d["model"]["couplings"] = {"t0": {"operator": "hop", "value": -1.0}}
    d["model"]["onsite"] = {}
    m = {REL: {**base_manifest()[REL], "lattice": "wannier90",
               "bonds_per_uc": {"t0": 1}, "coordination": {"A": {"t0": 2}}}}
    errs = lint(d, m)
    check("C12" not in only_ids(errs), f"wannier90 plain-number hop value must not raise C12, got {errs}")


def _wannier90_hop_doc_and_manifest(value) -> tuple[dict, dict]:
    d = base_doc()
    d["catalog"]["lattice"] = "wannier90"
    d["model"]["site_dof"]["A"] = {"fermion": {"orbitals": 1}}
    d["model"]["bonds"] = [{"type": "t0", "from": "A", "to": "A", "R": [1]}]
    d["model"]["couplings"] = {"t0": {"operator": "hop", "value": value}}
    d["model"]["onsite"] = {}
    m = {REL: {**base_manifest()[REL], "lattice": "wannier90",
               "bonds_per_uc": {"t0": 1}, "coordination": {"A": {"t0": 2}}}}
    return d, m


def test_c12_wannier90_negative_dict_without_param() -> None:
    """The wannier90 plain-number exception (Rule W1) covers only *bare
    numeric literals*, not dicts. ``value: {"scale": -1.0}`` is neither a
    valid ``{param, ...}`` reference nor a plain number, so it must still
    raise C12 even under the wannier90 dialect -- pins the existing
    ``_is_number``/dict-shape branching in check_c12()."""
    d, m = _wannier90_hop_doc_and_manifest({"scale": -1.0})
    errs = lint(d, m)
    check("C12" in only_ids(errs),
          f"wannier90 hop value={{'scale': -1.0}} (dict without param) must raise C12, got {errs}")


def test_c12_wannier90_negative_boolean_value() -> None:
    """The wannier90 plain-number exception excludes booleans: ``_is_number``
    explicitly rejects ``bool`` even though ``True == 1`` in Python, so
    ``value: true`` must still raise C12 under the wannier90 dialect --
    pins the existing ``_is_number`` bool-exclusion in check_c12()."""
    d, m = _wannier90_hop_doc_and_manifest(True)
    errs = lint(d, m)
    check("C12" in only_ids(errs),
          f"wannier90 hop value=True (bool, not a plain number) must raise C12, got {errs}")


def test_c12_negative_plain_number_non_wannier90() -> None:
    d = base_doc()
    d["model"]["site_dof"]["A"] = {"fermion": {"orbitals": 1}}
    d["model"]["bonds"] = [{"type": "t0", "from": "A", "to": "A", "R": [1]}]
    d["model"]["couplings"] = {"t0": {"operator": "hop", "value": -1.0}}
    d["model"]["onsite"] = {}
    m = {REL: {**base_manifest()[REL], "bonds_per_uc": {"t0": 1}, "coordination": {"A": {"t0": 2}}}}
    errs = lint(d, m)
    check("C12" in only_ids(errs),
          f"plain-number coupling value outside wannier90 dialect must raise C12, got {errs}")


def test_c12_negative_onsite_coeff_wrong_sign() -> None:
    d = base_doc()
    d["model"]["onsite"]["A"]["aniso_z"]["operator"]["tensor_terms"][0]["coeff"] = -1.0  # should be +1.0
    errs = lint(d)
    check("C12" in only_ids(errs), f"aniso_z coeff=-1.0 (should be +1.0) must raise C12, got {errs}")


def test_c12_negative_hop_value_dict_without_param() -> None:
    """A non-J coupling `value` dict lacking a `param` key (e.g. only
    `scale`) must be a C12 error -- it is not a valid `{param, ...}`
    reference and is not the wannier90 plain-number exception either."""
    d = _fermion_hop_doc()
    d["model"]["couplings"]["t0"]["value"] = {"scale": -1.0}
    errs = lint(d, _fermion_hop_manifest())
    check("C12" in only_ids(errs),
          f"hop value dict without 'param' must raise C12, got {errs}")


def test_c12_negative_onsite_coeff_boolean() -> None:
    """`coeff: true` must be rejected even though `True == 1.0` in Python --
    the literal sign-convention coeff must be a real number, not a bool."""
    d = base_doc()
    d["model"]["onsite"]["A"]["aniso_z"]["operator"]["tensor_terms"][0]["coeff"] = True
    errs = lint(d)
    check("C12" in only_ids(errs), f"onsite coeff=True (bool) must raise C12, got {errs}")


def test_c12_all_31_catalog_files_clean() -> None:
    """A4: linting the real catalog must not raise any C12 diagnostic
    (kept as a targeted subset check; see
    test_all_31_catalog_files_fully_clean() for the full-errs version)."""
    manifest_files = lc.load_manifest()
    for path in lc._discover_default_files():
        rel_path = path.resolve().relative_to(lc.ROOT).as_posix()
        errs = lc.lint_file(path, KEYWORDS, manifest_files)
        c12_errs = [e for e in errs if "C12" in e]
        check(c12_errs == [], f"{rel_path}: must have no C12 diagnostics, got {c12_errs}")


def test_all_31_catalog_files_fully_clean() -> None:
    """Strengthened whole-catalog check: exactly 31 catalog YAML files are
    discovered, and each one lints with *zero* diagnostics of any kind
    (not just no-C12)."""
    manifest_files = lc.load_manifest()
    files = lc._discover_default_files()
    check(len(files) == 31, f"expected exactly 31 catalog YAML files, got {len(files)}: {files}")
    for path in files:
        rel_path = path.resolve().relative_to(lc.ROOT).as_posix()
        errs = lc.lint_file(path, KEYWORDS, manifest_files)
        check(errs == [], f"{rel_path}: must lint fully clean, got {errs}")


# ---------------------------------------------------------------------------
#  A6: robustness (top-level failure handling, orphan manifest entries)
# ---------------------------------------------------------------------------

def test_a6_load_manifest_bad_yaml_raises_fatal_lint_error(tmp_path=None) -> None:
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        bad = Path(d) / "manifest.yaml"
        bad.write_text("files: [unterminated\n  - broken", encoding="utf-8")
        try:
            lc.load_manifest(bad)
            raised = False
        except lc.FatalLintError:
            raised = True
        except Exception as e:  # noqa: BLE001
            raise AssertionError(f"expected FatalLintError, got {type(e).__name__}: {e}")
        check(raised, "load_manifest() with malformed YAML must raise FatalLintError, not crash")


def test_a6_load_manifest_missing_file_returns_empty() -> None:
    missing = Path("/nonexistent/path/manifest.yaml")
    check(lc.load_manifest(missing) == {}, "load_manifest() on a nonexistent path must return {}")


def test_a6_orphan_manifest_entry_detected() -> None:
    files = [lc.ROOT / REL]  # a real file that exists on disk
    manifest_files = {REL: {}, "nonexistent/ghost.yaml": {}}
    errs = lc.find_orphan_manifest_entries(files, manifest_files)
    check(any("ghost.yaml" in e for e in errs),
          f"a manifest entry with no corresponding file must be reported as an orphan, got {errs}")
    check(not any("chain_spin.yaml" in e for e in errs),
          f"a manifest entry that DOES have a corresponding file must not be reported, got {errs}")


def test_a6_no_orphans_in_real_manifest() -> None:
    """A6-ii end-to-end: the real manifest.yaml and the real discovered
    catalog files must be in exact set-equality (no orphans either way)."""
    manifest_files = lc.load_manifest()
    files = lc._discover_default_files()
    errs = lc.find_orphan_manifest_entries(files, manifest_files)
    check(errs == [], f"real manifest.yaml must have no orphan entries, got {errs}")


# ---------------------------------------------------------------------------
#  malformed YAML / null document
# ---------------------------------------------------------------------------

def test_malformed_yaml_does_not_raise() -> None:
    text = "catalog: [unterminated\n  - broken"
    errs = lc.lint_text(text, REL, KEYWORDS, base_manifest())
    check(isinstance(errs, list) and len(errs) >= 1, "malformed YAML must yield diagnostics, not raise")
    check(any("C2" in e for e in errs), f"malformed YAML diagnostic should be tagged C2, got {errs}")


def test_null_document_does_not_raise() -> None:
    errs = lc.lint_text("", REL, KEYWORDS, base_manifest())
    check(isinstance(errs, list) and len(errs) >= 1, "empty/null YAML document must yield diagnostics, not raise")
    check(any("C2" in e for e in errs), f"null document diagnostic should be tagged C2, got {errs}")


def test_missing_keys_do_not_raise() -> None:
    doc = {"catalog": {"schema": "stdface-catalog/0.1", "dialect": "experimental",
                        "lattice": "chain", "model": "spin"}}
    errs = lc.lint_document(doc, REL, KEYWORDS, base_manifest())
    check(isinstance(errs, list) and len(errs) >= 1, "docs with missing sections must yield diagnostics, not raise")


# ===========================================================================
#  Runner
# ===========================================================================

def main() -> None:
    tests = [obj for name, obj in sorted(globals().items())
             if name.startswith("test_") and callable(obj)]
    for t in tests:
        t()
    print(f"ran {len(tests)} test functions, {_PASS_COUNT} assertions")
    print("all tools tests passed")


if __name__ == "__main__":
    main()
