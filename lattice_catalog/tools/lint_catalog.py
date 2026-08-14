#!/usr/bin/env python3
"""lattice_catalog の意味検査リンタ。

`lattice_catalog/CONVENTIONS.md` (= スペック §4/§6) が定める規約に対して、
`lattice_catalog/**/*.yaml`(manifest.yaml を除く)を検査する。

チェック項目(各診断メッセージは対応する ID を接頭辞に持つ):

- C1  catalog ヘッダ(schema/dialect の値、lattice/model の存在)
- C2  schema 検査(必須キー・型・サイトラベル重複・dimension/R/size の整数性)
- C3  R の長さ = dimension = size の長さ
- C4  ラベル整合(bonds の from/to、onsite のラベル、
      geometry.sites と site_dof の集合一致)
- C5  bonds の type ↔ couplings キー整合
- C6  反転同値での重複ボンド検出
- C7  param 参照の目録整合(scale/default の型検査を含む)
- C8  J coupling の 9 成分完全性・重複なし・ops↔接尾辞対応
- C9  演算子と site_dof の型整合、tensor_terms/onsite の ops 長
- C10 manifest 全項目突合
- C11 計数展開による配位数の実測比較

引数なしで実行すると `lattice_catalog/` 以下の全 `*.yaml`(manifest.yaml
を除く)を検査する。引数を渡すとそのファイル(群)のみを検査する
(ROOT = このスクリプトの一つ上のディレクトリ = lattice_catalog/。
ROOT 外のパスはエラーとして報告する)。

Run
---
python3 lattice_catalog/tools/lint_catalog.py
"""
from __future__ import annotations

import itertools
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = Path(__file__).resolve().parent

_REQUIRED_MANIFEST_KEYS = (
    "lattice", "model", "dimension", "n_sites_uc",
    "bonds_per_uc", "coordination", "min_size_for_check", "source",
)

_ONSITE_FERMION_OPS = {"N", "Nup", "Ndn", "NupNdn"}
_ONSITE_SPIN_OPS = {"Sx", "Sy", "Sz", "Szz"}

_J_COMPONENTS = {  # C8: ops 対 ↔ param 接尾辞
    ("Sx", "Sx"): "x",  ("Sy", "Sy"): "y",  ("Sz", "Sz"): "z",
    ("Sx", "Sy"): "xy", ("Sx", "Sz"): "xz", ("Sy", "Sx"): "yx",
    ("Sy", "Sz"): "yz", ("Sz", "Sx"): "zx", ("Sz", "Sy"): "zy",
}


# ---------------------------------------------------------------------------
#  Inventory / manifest loading
# ---------------------------------------------------------------------------

def load_inventory_keywords() -> set[str]:
    """Run ``keyword_inventory.py`` as a subprocess and return the set of
    canonical *keyword*-kind names (param references are checked against
    this set — lattice/model aliases are a separate namespace).

    Returns
    -------
    set[str]
        Canonical keyword names, e.g. ``{"t0", "J0x", "2S", "phase0", ...}``.
    """
    proc = subprocess.run(
        [sys.executable, str(TOOLS_DIR / "keyword_inventory.py")],
        capture_output=True, text=True, check=True,
    )
    entries = json.loads(proc.stdout)
    return {e["keyword_canonical"] for e in entries if e.get("kind") == "keyword"}


def load_manifest(manifest_path: Path | None = None) -> dict[str, Any]:
    """Load ``manifest.yaml`` and return its ``files`` mapping.

    Parameters
    ----------
    manifest_path : Path or None, optional
        Path to the manifest file. Defaults to ``ROOT / "manifest.yaml"``.

    Returns
    -------
    dict
        The ``files`` mapping (empty dict if absent/empty).
    """
    path = manifest_path if manifest_path is not None else ROOT / "manifest.yaml"
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as f:
        doc = yaml.safe_load(f)
    if not isinstance(doc, dict):
        return {}
    files = doc.get("files")
    return files if isinstance(files, dict) else {}


# ---------------------------------------------------------------------------
#  Core algorithms (C6, C8, C11) — kept verbatim to the reference spec
# ---------------------------------------------------------------------------

def reversal_dup(bonds: list[dict]) -> list[str]:
    """C6: 反転同値 (type,i,j,R)~(type,j,i,-R) での重複検出。

    キーに type を含める — 同一幾何ボンド上の t/V 共存は正当。
    """
    seen, errs = set(), []
    for b in bonds:
        k = (b["type"], b["from"], b["to"], tuple(b["R"]))
        rk = (b["type"], b["to"], b["from"], tuple(-x for x in b["R"]))
        if k in seen or rk in seen:
            errs.append(f"C6: (逆向き)重複ボンド {b}")
        seen.add(k)
    return errs


def check_j_coupling(type_name: str, tensor_terms: list[dict]) -> list[str]:
    """C8: 総項数 9・重複なし・ops↔接尾辞対応。

    重複 ops 対は交換係数の二重加算になるため明示的にエラー。
    """
    errs, seen = [], set()
    if len(tensor_terms) != 9:
        errs.append(f"C8: {type_name}: 項数 {len(tensor_terms)} != 9")
    for i, tt in enumerate(tensor_terms):
        pair = tuple(tt.get("ops", []))
        suffix = _J_COMPONENTS.get(pair)
        if suffix is None:
            errs.append(f"C8: {type_name}[{i}]: 不正な ops 対 {pair}")
            continue
        if pair in seen:
            errs.append(f"C8: {type_name}[{i}]: ops 対 {pair} の重複")
        seen.add(pair)
        coeff = tt.get("coeff")
        expected = f"{type_name}{suffix}"
        got = coeff.get("param") if isinstance(coeff, dict) else None
        if got != expected:
            errs.append(f"C8: {type_name}[{i}]: param {got} != {expected}")
    return errs


def expand_and_count(dim: int, labels: list[str], bonds: list[dict],
                      size: list[int]) -> dict[str, dict[str, int]]:
    """C11: min_size トーラス上でボンドを展開し、
    ラベル×type の配位数(そのラベルのサイト 1 個に接続する本数)を返す。

    from/to 双方の接続を数える。R=0 の同一ラベル自己ボンドは現行
    StdFace に存在しない前提とし、呼び出し側が事前検査する。
    """
    ncells = 1
    for s in size:
        ncells *= s
    touch: dict[str, dict[str, int]] = {lb: {} for lb in labels}
    for _cell in itertools.product(*[range(s) for s in size]):
        for b in bonds:
            t = b["type"]
            touch[b["from"]][t] = touch[b["from"]].get(t, 0) + 1
            touch[b["to"]][t] = touch[b["to"]].get(t, 0) + 1
    return {lb: {t: n // ncells for t, n in d.items()} for lb, d in touch.items()}


# ---------------------------------------------------------------------------
#  Small helpers
# ---------------------------------------------------------------------------

def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def resolve_param(ref: dict, param_value: float | None) -> float:
    """Reference resolver for the ``{param, scale, default}`` general form
    (CONVENTIONS.md §6.3). Not used by the structural C7 check itself
    (resolution is the *consumer's* responsibility — see §6.5), but kept
    here as the pinned, testable definition of the semantics.

    Parameters
    ----------
    ref : dict
        A ``{param, scale, default}``-shaped reference. ``scale`` defaults
        to ``1.0``, ``default`` defaults to ``0``.
    param_value : float or None
        The value supplied for ``ref["param"]`` by the (hypothetical)
        input, or ``None`` if not supplied.

    Returns
    -------
    float
        ``scale * param_value`` if *param_value* is given, else
        ``default`` as-is (``scale`` is **not** applied to ``default``).
    """
    if param_value is not None:
        return ref.get("scale", 1.0) * param_value
    return ref.get("default", 0)


def _find_param_refs(obj: Any, path: str = "$") -> list[tuple[str, dict]]:
    """Recursively find every dict containing a ``param`` key.

    Returns
    -------
    list[tuple[str, dict]]
        ``(path, dict)`` pairs for diagnostics.
    """
    found = []
    if isinstance(obj, dict):
        if "param" in obj:
            found.append((path, obj))
        for k, v in obj.items():
            found.extend(_find_param_refs(v, f"{path}.{k}"))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            found.extend(_find_param_refs(v, f"{path}[{i}]"))
    return found


def _site_dof_kind(entry: Any) -> str | None:
    """Return ``"fermion"``, ``"spin"``, or ``None`` (unrecognised) for a
    single ``site_dof`` entry (``{fermion: {...}}`` or ``{spin: {...}}``)."""
    if not isinstance(entry, dict):
        return None
    if "fermion" in entry:
        return "fermion"
    if "spin" in entry:
        return "spin"
    return None


_MISSING = object()  # sentinel: "not J-family shaped", distinct from an
                      # explicitly-present-but-malformed tensor_terms value.


def _coupling_operator_tensor_terms(c: Any) -> Any:
    """Extract the J-family ``tensor_terms`` list from a ``couplings`` entry.

    ``model.couplings[*]`` entries come in two shapes (CONVENTIONS.md §6,
    fixed by the chain catalog templates — Task 1):

    - scalar/vector-valued (``hop`` / ``density-density`` / ``s_i . S_j``):
      ``{operator: <str>, value: {param, scale, default}}``.
    - J-family (9-component exchange tensor): ``{operator: {tensor_terms:
      [...]}}`` — no shared top-level ``value``; each term carries its own
      ``coeff: {param: ...}``.

    Parameters
    ----------
    c : Any
        A single ``model.couplings`` value.

    Returns
    -------
    Any
        ``c["operator"]["tensor_terms"]`` (whatever its type — the caller
        validates it is a list) if *c* is J-family shaped, else the
        ``_MISSING`` sentinel.
    """
    if not isinstance(c, dict):
        return _MISSING
    op = c.get("operator")
    if isinstance(op, dict) and "tensor_terms" in op:
        return op["tensor_terms"]
    return _MISSING


# ---------------------------------------------------------------------------
#  Individual checks (C1-C11)
# ---------------------------------------------------------------------------

def check_c1(doc: dict) -> list[str]:
    errs = []
    catalog = doc.get("catalog")
    if not isinstance(catalog, dict):
        return ["C1: 'catalog' section missing or not a mapping"]
    if catalog.get("schema") != "stdface-catalog/0.1":
        errs.append(f"C1: catalog.schema != 'stdface-catalog/0.1' (got {catalog.get('schema')!r})")
    if catalog.get("dialect") != "experimental":
        errs.append(f"C1: catalog.dialect != 'experimental' (got {catalog.get('dialect')!r})")
    if not catalog.get("lattice"):
        errs.append("C1: catalog.lattice missing")
    if not catalog.get("model"):
        errs.append("C1: catalog.model missing")
    return errs


def check_c2(doc: dict) -> tuple[list[str], dict]:
    """Schema validation. Returns (errors, parsed-context).

    The returned context carries best-effort extracted values (dimension,
    site labels, ...) so later checks can proceed defensively even when
    parts of the document are malformed.
    """
    errs: list[str] = []
    ctx: dict[str, Any] = {}

    for key in ("catalog", "geometry", "system", "model"):
        if not isinstance(doc.get(key), dict):
            errs.append(f"C2: top-level key '{key}' missing or not a mapping")

    geometry = doc.get("geometry") if isinstance(doc.get("geometry"), dict) else {}
    system = doc.get("system") if isinstance(doc.get("system"), dict) else {}
    model = doc.get("model") if isinstance(doc.get("model"), dict) else {}

    dimension = geometry.get("dimension")
    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension <= 0:
        errs.append(f"C2: geometry.dimension must be a positive integer (got {dimension!r})")
        dimension = None
    ctx["dimension"] = dimension

    sites = geometry.get("sites")
    labels: list[str] = []
    if not isinstance(sites, list) or not sites:
        errs.append("C2: geometry.sites missing or empty")
    else:
        for i, s in enumerate(sites):
            if not isinstance(s, dict) or "label" not in s:
                errs.append(f"C2: geometry.sites[{i}] missing 'label'")
                continue
            labels.append(s["label"])
        if len(labels) != len(set(labels)):
            dups = sorted({lb for lb in labels if labels.count(lb) > 1})
            errs.append(f"C2: duplicate site labels in geometry.sites: {dups}")
    ctx["labels"] = labels

    size = system.get("size")
    if not isinstance(size, list) or not all(
        isinstance(x, int) and not isinstance(x, bool) for x in size
    ):
        errs.append(f"C2: system.size must be a list of integers (got {size!r})")
        size = None
    ctx["size"] = size

    site_dof = model.get("site_dof")
    if not isinstance(site_dof, dict) or not site_dof:
        errs.append("C2: model.site_dof missing or empty")
        site_dof = {}
    ctx["site_dof"] = site_dof

    bonds = model.get("bonds")
    if not isinstance(bonds, list):
        errs.append("C2: model.bonds missing or not a list")
        bonds = []
    for i, b in enumerate(bonds):
        if not isinstance(b, dict) or not {"type", "from", "to", "R"} <= b.keys():
            errs.append(f"C2: model.bonds[{i}] missing required keys (type/from/to/R)")
            continue
        if not isinstance(b["R"], list) or not all(
            isinstance(x, int) and not isinstance(x, bool) for x in b["R"]
        ):
            errs.append(f"C2: model.bonds[{i}].R must be a list of integers (got {b['R']!r})")
    ctx["bonds"] = bonds

    couplings = model.get("couplings")
    if not isinstance(couplings, dict):
        errs.append("C2: model.couplings missing or not a mapping")
        couplings = {}
    ctx["couplings"] = couplings

    onsite_raw = model.get("onsite")
    onsite: list[dict] = []
    if not isinstance(onsite_raw, dict):
        errs.append("C2: model.onsite missing or not a mapping of site -> {term_name: {...}}")
    else:
        for site, terms in onsite_raw.items():
            if not isinstance(terms, dict) or not terms:
                errs.append(f"C2: model.onsite[{site!r}] must be a non-empty mapping of term_name -> term")
                continue
            for term_name, term in terms.items():
                if not isinstance(term, dict):
                    errs.append(f"C2: model.onsite[{site!r}][{term_name!r}] must be a mapping")
                    continue
                op = term.get("operator")
                tt = op.get("tensor_terms") if isinstance(op, dict) else None
                ops_list: Any = None
                if isinstance(tt, list) and len(tt) == 1 and isinstance(tt[0], dict):
                    ops_list = tt[0].get("ops")
                else:
                    errs.append(
                        f"C2: model.onsite[{site!r}][{term_name!r}].operator.tensor_terms "
                        "must be a single-element list"
                    )
                onsite.append({
                    "site": site,
                    "term": term_name,
                    "ops": ops_list if isinstance(ops_list, list) else [],
                })
    ctx["onsite"] = onsite

    return errs, ctx


def check_c3(ctx: dict) -> list[str]:
    errs = []
    dimension, size, bonds = ctx["dimension"], ctx["size"], ctx["bonds"]
    if dimension is None or size is None:
        return errs  # already reported by C2
    if len(size) != dimension:
        errs.append(f"C3: len(system.size)={len(size)} != dimension={dimension}")
    for i, b in enumerate(bonds):
        R = b.get("R")
        if isinstance(R, list) and len(R) != dimension:
            errs.append(f"C3: model.bonds[{i}].R length={len(R)} != dimension={dimension}")
    return errs


def check_c4(ctx: dict) -> list[str]:
    errs = []
    labels = set(ctx["labels"])
    site_dof = ctx["site_dof"]
    bonds = ctx["bonds"]
    onsite = ctx["onsite"]

    if labels != set(site_dof.keys()):
        errs.append(
            "C4: geometry.sites labels != model.site_dof keys "
            f"(sites-only: {sorted(labels - set(site_dof))}, "
            f"site_dof-only: {sorted(set(site_dof) - labels)})"
        )

    for i, b in enumerate(bonds):
        for end in ("from", "to"):
            if b.get(end) not in labels:
                errs.append(f"C4: model.bonds[{i}].{end}={b.get(end)!r} is not a defined site label")

    for i, o in enumerate(onsite):
        site = o.get("site") if isinstance(o, dict) else None
        if site not in labels:
            errs.append(f"C4: model.onsite[{i}].site={site!r} is not a defined site label")

    return errs


def check_c5(ctx: dict) -> list[str]:
    errs = []
    used = {b["type"] for b in ctx["bonds"] if "type" in b}
    defined = set(ctx["couplings"].keys())
    for t in sorted(used - defined):
        errs.append(f"C5: bond type {t!r} has no couplings definition")
    for t in sorted(defined - used):
        errs.append(f"C5: couplings key {t!r} is not referenced by any bond")
    return errs


def check_c6(ctx: dict) -> list[str]:
    bonds = ctx["bonds"]
    if not bonds or any(not {"type", "from", "to", "R"} <= b.keys() for b in bonds):
        return []  # malformed bonds already reported by C2
    return reversal_dup(bonds)


def check_c7(doc: dict, inventory: set[str]) -> list[str]:
    errs = []
    for path, ref in _find_param_refs(doc):
        name = ref.get("param")
        if not isinstance(name, str) or name not in inventory:
            errs.append(f"C7: {path}: param {name!r} not found in keyword inventory")
        if "scale" in ref and not _is_number(ref["scale"]):
            errs.append(f"C7: {path}: scale must be a real number (got {ref['scale']!r})")
        if "default" in ref and not _is_number(ref["default"]):
            errs.append(f"C7: {path}: default must be numeric (got {ref['default']!r})")
    return errs


def check_c8(ctx: dict) -> list[str]:
    errs = []
    for type_name, c in ctx["couplings"].items():
        tt = _coupling_operator_tensor_terms(c)
        if tt is _MISSING:
            continue
        if not isinstance(tt, list):
            errs.append(f"C8: {type_name}: tensor_terms is not a list")
            continue
        errs.extend(check_j_coupling(type_name, tt))
    return errs


def check_c9(ctx: dict) -> list[str]:
    errs = []
    site_dof = ctx["site_dof"]
    kinds = {lb: _site_dof_kind(v) for lb, v in site_dof.items()}
    couplings = ctx["couplings"]

    for b in ctx["bonds"]:
        t = b.get("type")
        c = couplings.get(t)
        if not isinstance(c, dict):
            continue
        frm, to = b.get("from"), b.get("to")
        kf, kt = kinds.get(frm), kinds.get(to)
        tt = _coupling_operator_tensor_terms(c)
        if tt is not _MISSING:
            if kf != "spin" or kt != "spin":
                errs.append(f"C9: bond type {t!r} ({frm}->{to}): J tensor requires spin-spin site_dof (got {kf}-{kt})")
            for i, term in enumerate(tt if isinstance(tt, list) else []):
                ops = term.get("ops", []) if isinstance(term, dict) else []
                if len(ops) != 2:
                    errs.append(f"C9: {t}.tensor_terms[{i}]: ops length={len(ops)} != 2")
        else:
            op = c.get("operator")
            if op in ("hop", "density-density"):
                if kf != "fermion" or kt != "fermion":
                    errs.append(f"C9: bond type {t!r} ({frm}->{to}): operator {op!r} requires fermion-fermion site_dof (got {kf}-{kt})")
            elif op == "s_i . S_j":
                if kf != "fermion" or kt != "spin":
                    errs.append(f"C9: bond type {t!r} ({frm}->{to}): operator 's_i . S_j' requires fermion-spin site_dof in this order (got {kf}-{kt})")
            else:
                errs.append(f"C9: bond type {t!r}: unknown operator {op!r}")

    for i, o in enumerate(ctx["onsite"]):
        if not isinstance(o, dict):
            continue
        ops = o.get("ops", [])
        if len(ops) != 1:
            errs.append(f"C9: model.onsite[{i}]: ops length={len(ops)} != 1")
            continue
        op = ops[0]
        kind = kinds.get(o.get("site"))
        if op in _ONSITE_FERMION_OPS:
            if kind != "fermion":
                errs.append(f"C9: model.onsite[{i}]: op {op!r} requires fermion site_dof (got {kind})")
        elif op in _ONSITE_SPIN_OPS:
            if kind not in ("spin", "fermion"):
                errs.append(f"C9: model.onsite[{i}]: op {op!r} requires spin or fermion site_dof (got {kind})")
        else:
            errs.append(f"C9: model.onsite[{i}]: unknown onsite operator {op!r}")

    return errs


def check_c10_c11(doc: dict, ctx: dict, rel_path: str,
                   manifest_files: dict) -> list[str]:
    errs: list[str] = []
    entry = manifest_files.get(rel_path)
    if entry is None:
        return [f"C10: no manifest entry for {rel_path!r}"]
    if not isinstance(entry, dict):
        return [f"C10: manifest entry for {rel_path!r} is not a mapping"]

    missing = [k for k in _REQUIRED_MANIFEST_KEYS if k not in entry]
    if missing:
        errs.append(f"C10: manifest entry for {rel_path!r} missing keys: {missing}")
        return errs  # can't safely compare further

    catalog = doc.get("catalog", {}) if isinstance(doc.get("catalog"), dict) else {}
    if entry["lattice"] != catalog.get("lattice"):
        errs.append(f"C10: manifest.lattice={entry['lattice']!r} != catalog.lattice={catalog.get('lattice')!r}")
    if entry["model"] != catalog.get("model"):
        errs.append(f"C10: manifest.model={entry['model']!r} != catalog.model={catalog.get('model')!r}")
    if ctx["dimension"] is not None and entry["dimension"] != ctx["dimension"]:
        errs.append(f"C10: manifest.dimension={entry['dimension']!r} != geometry.dimension={ctx['dimension']!r}")
    if entry["n_sites_uc"] != len(ctx["labels"]):
        errs.append(f"C10: manifest.n_sites_uc={entry['n_sites_uc']!r} != len(geometry.sites)={len(ctx['labels'])}")

    bonds = ctx["bonds"]
    actual_bonds_per_uc: dict[str, int] = {}
    for b in bonds:
        t = b.get("type")
        actual_bonds_per_uc[t] = actual_bonds_per_uc.get(t, 0) + 1
    if entry["bonds_per_uc"] != actual_bonds_per_uc:
        errs.append(
            f"C10: manifest.bonds_per_uc={entry['bonds_per_uc']!r} != "
            f"actual per-unit-cell bond counts={actual_bonds_per_uc!r}"
        )

    if errs:
        return errs  # avoid noisy cascading C11 errors on top of C10 mismatches

    # --- C11: counting expansion -------------------------------------
    dimension = ctx["dimension"]
    labels = ctx["labels"]
    min_size = entry["min_size_for_check"]
    if not isinstance(min_size, list) or len(min_size) != dimension:
        errs.append(f"C11: manifest.min_size_for_check length must equal dimension={dimension} (got {min_size!r})")
        return errs

    self_loop = any(
        b.get("from") == b.get("to") and all(x == 0 for x in b.get("R", []))
        for b in bonds
    )
    if self_loop:
        errs.append("C2: self-loop bond (from == to, R == 0) is not supported")
        return errs

    if any(b.get("from") not in labels or b.get("to") not in labels for b in bonds):
        # Already reported by C4; skip the expansion to avoid a KeyError.
        return errs

    measured = expand_and_count(dimension, labels, bonds, min_size)
    expected = entry["coordination"]
    if measured != expected:
        errs.append(f"C11: coordination mismatch: manifest={expected!r} != measured={measured!r}")

    return errs


# ---------------------------------------------------------------------------
#  Per-file / per-document orchestration
# ---------------------------------------------------------------------------

def lint_document(doc: Any, rel_path: str, inventory: set[str],
                   manifest_files: dict) -> list[str]:
    """Run all checks (C1-C11) against a parsed YAML document.

    Never raises: any unexpected exception is converted to a C2 diagnostic.
    """
    if not isinstance(doc, dict):
        return ["C2: document is not a mapping (null or invalid top-level YAML)"]

    try:
        errs: list[str] = []
        errs.extend(check_c1(doc))
        c2_errs, ctx = check_c2(doc)
        errs.extend(c2_errs)
        errs.extend(check_c3(ctx))
        errs.extend(check_c4(ctx))
        errs.extend(check_c5(ctx))
        errs.extend(check_c6(ctx))
        errs.extend(check_c7(doc, inventory))
        errs.extend(check_c8(ctx))
        errs.extend(check_c9(ctx))
        errs.extend(check_c10_c11(doc, ctx, rel_path, manifest_files))
        return errs
    except Exception as e:  # noqa: BLE001 — must never crash the linter
        return [f"C2: exception while checking {rel_path!r}: {e!r}"]


def lint_text(text: str, rel_path: str, inventory: set[str],
              manifest_files: dict) -> list[str]:
    """Parse *text* as YAML and lint it. Never raises on bad YAML."""
    try:
        doc = yaml.safe_load(text)
    except yaml.YAMLError as e:
        return [f"C2: YAML parse error in {rel_path!r}: {e}"]
    return lint_document(doc, rel_path, inventory, manifest_files)


def lint_file(path: Path, inventory: set[str], manifest_files: dict) -> list[str]:
    try:
        rel_path = path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        rel_path = path.name
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as e:
        return [f"C2: cannot read {rel_path!r}: {e}"]
    return lint_text(text, rel_path, inventory, manifest_files)


# ---------------------------------------------------------------------------
#  CLI
# ---------------------------------------------------------------------------

def _discover_default_files() -> list[Path]:
    return sorted(
        p for p in ROOT.rglob("*.yaml")
        if p.name != "manifest.yaml"
    )


def _resolve_arg_paths(args: list[str]) -> tuple[list[Path], list[str]]:
    """Resolve CLI path arguments against ROOT.

    Returns
    -------
    tuple[list[Path], list[str]]
        Valid in-ROOT paths, and top-level error messages for paths that
        could not be resolved into ROOT.
    """
    files: list[Path] = []
    top_errors: list[str] = []
    for arg in args:
        resolved = Path(arg).resolve()
        try:
            resolved.relative_to(ROOT)
        except ValueError:
            top_errors.append(f"C2: path outside ROOT ({ROOT}): {resolved}")
            continue
        if resolved.name == "manifest.yaml":
            continue
        files.append(resolved)
    return files, top_errors


def main() -> None:
    args = sys.argv[1:]
    if args:
        files, top_errors = _resolve_arg_paths(args)
    else:
        files, top_errors = _discover_default_files(), []

    inventory = load_inventory_keywords()
    manifest_files = load_manifest()

    total_errors = list(top_errors)
    n_files = 0
    for path in files:
        n_files += 1
        errs = lint_file(path, inventory, manifest_files)
        for e in errs:
            print(f"{path}: {e}")
        total_errors.extend(errs)

    for e in top_errors:
        print(e)

    print(f"{n_files} files, {len(total_errors)} errors")
    sys.exit(1 if total_errors else 0)


if __name__ == "__main__":
    main()
