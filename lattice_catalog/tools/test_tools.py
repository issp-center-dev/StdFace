#!/usr/bin/env python3
"""lattice_catalog/tools の自動テスト(開発時専用)。

pytest 非依存、assert ベース。`python3 test_tools.py` で完結し、
成功時は最後に `all tools tests passed` を出力する。

対象:
- keyword_inventory.py: canonical 集約・安定ソート・重複なし・
  ソルバーレジストリ件数下限・prime 付き成分の存在
- lint_catalog.py: C1-C11 の正例・負例(インライン YAML / dict fixture)
- param 参照の意味論(resolve_param): scale×param / default そのまま
- 不正 YAML・null 文書での非例外動作

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
            "cell": [[1.0]],
            "sites": [{"label": "A", "frac": [0.0]}],
        },
        "system": {
            "size": [7],
        },
        "model": {
            "site_dof": {
                "A": {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}},
            },
            "bonds": [
                {"type": "J0", "from": "A", "to": "A", "R": [1]},
            ],
            "couplings": {
                "J0": {"tensor_terms": copy.deepcopy(_J0_TERMS)},
            },
            "onsite": [
                {"site": "A", "ops": ["Szz"], "coeff": {"param": "D", "scale": 1.0, "default": 0}},
            ],
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
            "min_size_for_check": [7],
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
    d["model"]["onsite"][0]["site"] = "Z"
    errs = lint(d)
    check("C4" in only_ids(errs), f"onsite referencing unknown label must raise C4, got {errs}")


def test_c4_negative_site_dof_mismatch() -> None:
    d = base_doc()
    d["model"]["site_dof"]["Extra"] = {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}}
    errs = lint(d)
    check("C4" in only_ids(errs), f"site_dof/sites mismatch must raise C4, got {errs}")


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
    d["model"]["onsite"] = []
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
    d["model"]["onsite"][0]["coeff"]["param"] = "NotAKeyword"
    errs = lint(d)
    check("C7" in only_ids(errs), f"unknown param name must raise C7, got {errs}")


def test_c7_negative_bad_scale_type() -> None:
    d = base_doc()
    d["model"]["onsite"][0]["coeff"]["scale"] = "big"
    errs = lint(d)
    check("C7" in only_ids(errs), f"non-numeric scale must raise C7, got {errs}")


def test_c7_negative_bad_default_type() -> None:
    d = base_doc()
    d["model"]["onsite"][0]["coeff"]["default"] = "zero"
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
    d["model"]["couplings"]["J0"]["tensor_terms"].pop()
    errs = lint(d)
    check("C8" in only_ids(errs), f"8 terms (missing one) must raise C8, got {errs}")


def test_c8_negative_duplicate_ops() -> None:
    d = base_doc()
    d["model"]["couplings"]["J0"]["tensor_terms"][-1] = {
        "ops": ["Sx", "Sx"], "coeff": {"param": "J0x"},
    }
    errs = lint(d)
    check("C8" in only_ids(errs), f"duplicate ops pair must raise C8, got {errs}")


def test_c8_negative_extra_term() -> None:
    d = base_doc()
    d["model"]["couplings"]["J0"]["tensor_terms"].append(
        {"ops": ["Sx", "Sx"], "coeff": {"param": "J0x"}}
    )
    errs = lint(d)
    check("C8" in only_ids(errs), f"10 terms must raise C8, got {errs}")


def test_c8_negative_bad_param_name() -> None:
    d = base_doc()
    d["model"]["couplings"]["J0"]["tensor_terms"][0]["coeff"]["param"] = "J0wrong"
    errs = lint(d)
    check("C8" in only_ids(errs), f"ops<->param suffix mismatch must raise C8, got {errs}")


# ---------------------------------------------------------------------------
#  C9
# ---------------------------------------------------------------------------

def _fermion_hop_doc() -> dict:
    d = base_doc()
    d["model"]["site_dof"] = {"A": {"fermion": {"orbitals": 1}}}
    d["model"]["bonds"] = [{"type": "t0", "from": "A", "to": "A", "R": [1]}]
    d["model"]["couplings"] = {"t0": {"operator": "hop", "value": {"param": "t0", "scale": -1.0}}}
    d["model"]["onsite"] = [{"site": "A", "ops": ["N"], "coeff": {"param": "mu", "scale": -1.0}}]
    return d


def _fermion_hop_manifest() -> dict:
    return {REL: {"lattice": "chain", "model": "hubbard", "dimension": 1, "n_sites_uc": 1,
                  "bonds_per_uc": {"t0": 1}, "coordination": {"A": {"t0": 2}},
                  "min_size_for_check": [7],
                  "source": {"file": "x", "func": "y", "commit": "z"}}}


def test_c9_positive_fermion_fermion_hop() -> None:
    errs = lint(_fermion_hop_doc(), _fermion_hop_manifest())
    check("C9" not in only_ids(errs), f"hop on fermion-fermion must not raise C9, got {errs}")


def test_c9_negative_hop_on_spin() -> None:
    d = _fermion_hop_doc()
    d["model"]["site_dof"]["A"] = {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}}
    errs = lint(d, _fermion_hop_manifest())
    check("C9" in only_ids(errs), f"hop on spin-spin must raise C9, got {errs}")


def test_c9_positive_kondo_fermion_spin() -> None:
    d = base_doc()
    d["model"]["site_dof"] = {
        "A_c": {"fermion": {"orbitals": 1}},
        "A_s": {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}},
    }
    d["geometry"]["sites"] = [{"label": "A_c", "frac": [0.0]}, {"label": "A_s", "frac": [0.0]}]
    d["model"]["bonds"] = [{"type": "J", "from": "A_c", "to": "A_s", "R": [0]}]
    d["model"]["couplings"] = {"J": {"operator": "s_i . S_j", "value": {"param": "J", "scale": 1.0}}}
    d["model"]["onsite"] = []
    m = {REL: {"lattice": "chain", "model": "kondo", "dimension": 1, "n_sites_uc": 2,
               "bonds_per_uc": {"J": 1},
               "coordination": {"A_c": {"J": 1}, "A_s": {"J": 1}},
               "min_size_for_check": [7],
               "source": {"file": "x", "func": "y", "commit": "z"}}}
    errs = lint(d, m)
    check("C9" not in only_ids(errs), f"s_i.S_j fermion(1st)-spin(2nd) must not raise C9, got {errs}")


def test_c9_negative_kondo_wrong_order() -> None:
    d = base_doc()
    d["model"]["site_dof"] = {
        "A_c": {"fermion": {"orbitals": 1}},
        "A_s": {"spin": {"param": "2S", "scale": 0.5, "default": 0.5}},
    }
    d["geometry"]["sites"] = [{"label": "A_c", "frac": [0.0]}, {"label": "A_s", "frac": [0.0]}]
    d["model"]["bonds"] = [{"type": "J", "from": "A_s", "to": "A_c", "R": [0]}]
    d["model"]["couplings"] = {"J": {"operator": "s_i . S_j", "value": {"param": "J", "scale": 1.0}}}
    d["model"]["onsite"] = []
    m = {REL: {"lattice": "chain", "model": "kondo", "dimension": 1, "n_sites_uc": 2,
               "bonds_per_uc": {"J": 1},
               "coordination": {"A_c": {"J": 1}, "A_s": {"J": 1}},
               "min_size_for_check": [7],
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
    d["model"]["onsite"][0] = {"site": "A", "ops": ["N"], "coeff": {"param": "mu", "scale": -1.0}}
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


def test_c11_negative_wrong_coordination() -> None:
    m = base_manifest()
    m[REL]["coordination"] = {"A": {"J0": 4}}
    errs = lint(base_doc(), m)
    check("C11" in only_ids(errs), f"wrong coordination must raise C11, got {errs}")


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
