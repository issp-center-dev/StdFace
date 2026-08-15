#!/usr/bin/env python3
"""lattice_catalog の意味検査リンタ。

`lattice_catalog/CONVENTIONS.md` (= スペック §4/§6) が定める規約に対して、
`lattice_catalog/**/*.yaml`(manifest.yaml を除く)を検査する。

チェック項目(各診断メッセージは対応する ID を接頭辞に持つ):

- C1  catalog ヘッダ(schema/dialect の値、lattice/model の存在)
- C2  schema 検査(必須キー・型・サイトラベル重複・dimension/R/size の整数性、
      geometry.lattice_vectors/sites/system.boundary/site_dof の深い型検査)
- C3  R の長さ = dimension = size の長さ
- C4  ラベル整合(bonds の from/to、onsite のラベル、
      geometry.sites と site_dof の集合一致)
- C5  bonds の type ↔ couplings キー整合
- C6  反転同値での重複ボンド検出
- C7  param 参照の目録整合(scale/default の型検査を含む)
- C8  J coupling の 9 成分完全性・重複なし・ops↔接尾辞対応
- C9  演算子と site_dof の型整合、tensor_terms/onsite の ops 長
- C10 manifest 全項目突合
- C11 計数展開による配位数の実測比較、および min_size_for_check の
      厳密性検査(各成分が正の奇数であり、かつ
      `2 * max|R 成分| + 1`(その方向で bonds から実測)に厳密一致するか)
- C12 CONVENTIONS.md §6.2 の符号規約検査(dialect: experimental のみ対象)。
      couplings の hop/density-density/s_i.S_j の value.scale、onsite の
      hubbard_u/aniso_z/chemical_potential/field_* の tensor_terms coeff、
      J 族 tensor_terms coeff の形(`{param: ...}` のみで scale/default 不可)、
      非 J couplings の value の形(`{param, ...}` dict、または wannier90
      方言に限り生の数値リテラル)を検査する。

引数なしで実行すると `lattice_catalog/` 以下の全 `*.yaml`(manifest.yaml
を除く)を検査する。引数を渡すとそのファイル(群)のみを検査する
(ROOT = このスクリプトの一つ上のディレクトリ = lattice_catalog/。
ROOT 外のパスはエラーとして報告する)。引数なし実行時は、発見した YAML
パス集合と manifest.yaml のキー集合の完全一致も検査する(孤立した
manifest エントリは診断として報告する — A6-ii)。

開発時依存: 本ツールは PyYAML (``pyyaml``) を必要とする(開発時専用の
lint ツールであり、``python/pyproject.toml`` の実行時依存には含めない
— A5)。未インストールの場合は分かりやすいメッセージで終了する。

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

try:
    import yaml
except ImportError:
    raise SystemExit(
        "lint_catalog: PyYAML が必要です(開発時専用ツール)。"
        "pip install pyyaml を実行してください。"
    )

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

class FatalLintError(RuntimeError):
    """Raised for top-level (non-per-file) failures that must abort the run
    with a clean diagnostic instead of a traceback (A6-i)."""


def load_inventory_keywords() -> set[str]:
    """Run ``keyword_inventory.py`` as a subprocess and return the set of
    canonical *keyword*-kind names (param references are checked against
    this set — lattice/model aliases are a separate namespace).

    Returns
    -------
    set[str]
        Canonical keyword names, e.g. ``{"t0", "J0x", "2S", "phase0", ...}``.

    Raises
    ------
    FatalLintError
        If the subprocess fails (non-zero exit) or its stdout is not valid
        JSON — converted from a raw ``CalledProcessError``/``JSONDecodeError``
        traceback into a clean top-level diagnostic (A6-i).
    """
    try:
        proc = subprocess.run(
            [sys.executable, str(TOOLS_DIR / "keyword_inventory.py")],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError as e:
        raise FatalLintError(
            f"keyword_inventory.py failed (exit {e.returncode}): "
            f"{e.stderr.strip() or e.stdout.strip()}"
        ) from e
    try:
        entries = json.loads(proc.stdout)
    except json.JSONDecodeError as e:
        raise FatalLintError(f"keyword_inventory.py produced invalid JSON: {e}") from e
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

    Raises
    ------
    FatalLintError
        If the manifest exists but cannot be read or parsed
        (``OSError``/``yaml.YAMLError``) — converted into a clean top-level
        diagnostic instead of a traceback (A6-i).
    """
    path = manifest_path if manifest_path is not None else ROOT / "manifest.yaml"
    if not path.exists():
        return {}
    try:
        with path.open(encoding="utf-8") as f:
            doc = yaml.safe_load(f)
    except OSError as e:
        raise FatalLintError(f"cannot read manifest {path}: {e}") from e
    except yaml.YAMLError as e:
        raise FatalLintError(f"YAML parse error in manifest {path}: {e}") from e
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


def expand_and_count(
    dim: int, labels: list[str], bonds: list[dict], size: list[int],
) -> tuple[dict[str, dict[str, int]], list[str]]:
    """C11: min_size トーラス上に実際にボンドを展開し、ラベル×type の
    配位数(実測)を計数する。

    各セル ``c``(``0 <= c[k] < size[k]``)と各ボンド ``b`` について、
    始点インスタンス ``(c, b['from'])`` と終点インスタンス
    ``((c + b['R']) mod size, b['to'])`` の間に type=``b['type']`` の
    展開ボンドを 1 本張る(``mod`` はトーラスの周期境界条件)。

    - **配位数**: 各ラベルについて、そのラベルを持つ任意の 1 インスタンス
      に接続する展開ボンド数を type 別に数える。トーラスは並進対称なので
      理論上は全インスタンスが同じ値を持つはずであり、それを実際に
      全インスタンスについて検証する(食い違えば C11 エラー)。
    - **折り畳み縮退**: 展開ボンドは ``(type, {始点インスタンス, 終点
      インスタンス})``(順序なしペア)をキーとすると、本来
      ``len(bonds) * ncells`` 個の相異なるキーを持つはずである。2 本の
      展開ボンドが同じキーに縮退した場合、``min_size_for_check`` が
      小さすぎて異なる R が同一サイト対に折り畳まれていることを意味し、
      C11 エラーとして報告する。

    Parameters
    ----------
    dim : int
        次元(``len(size)`` と一致)。
    labels : list of str
        単位胞内のサイトラベル一覧(``geometry.sites`` の順序)。
    bonds : list of dict
        ``{type, from, to, R}`` を持つボンド定義のリスト。
    size : list of int
        トーラスの各方向のセル数(manifest の ``min_size_for_check``)。

    Returns
    -------
    tuple[dict, list[str]]
        ``(measured, errs)``。``measured`` はラベル→{type: 配位数}。
        ``errs`` は折り畳み縮退・並進非対称が見つかった場合の C11
        診断メッセージのリスト(通常は空)。
    """
    errs: list[str] = []
    cells = list(itertools.product(*[range(s) for s in size]))
    ncells = len(cells)

    touch: dict[tuple, dict[str, int]] = {(c, lb): {} for c in cells for lb in labels}
    pair_counts: dict[tuple, int] = {}

    for c in cells:
        for b in bonds:
            t = b["type"]
            r = b["R"]
            a_inst = (c, b["from"])
            b_cell = tuple((c[k] + r[k]) % size[k] for k in range(dim))
            b_inst = (b_cell, b["to"])

            touch[a_inst][t] = touch[a_inst].get(t, 0) + 1
            touch[b_inst][t] = touch[b_inst].get(t, 0) + 1

            endpoints = (a_inst, b_inst) if a_inst <= b_inst else (b_inst, a_inst)
            key = (t, endpoints)
            pair_counts[key] = pair_counts.get(key, 0) + 1

    n_expected = len(bonds) * ncells
    n_distinct = len(pair_counts)
    if n_distinct != n_expected:
        errs.append(
            "C11: folding degeneracy — expanded bonds collide on the "
            f"min_size_for_check torus (size={size!r}): {n_distinct} distinct "
            f"instance-pairs found, expected {n_expected} (= len(bonds) * ncells); "
            "min_size_for_check is too small"
        )

    measured: dict[str, dict[str, int]] = {}
    for lb in labels:
        rep_cell = cells[0]
        rep = touch[(rep_cell, lb)]
        for c in cells[1:]:
            got = touch[(c, lb)]
            if got != rep:
                errs.append(
                    f"C11: coordination is not translation-invariant for label "
                    f"{lb!r}: instance at cell {c} has {got!r}, expected {rep!r} "
                    f"(same as cell {rep_cell})"
                )
        measured[lb] = rep

    return measured, errs


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

    # --- geometry.lattice_vectors: mapping of exactly `dimension` named ---
    #     vectors a1..aN, each a numeric list of length `dimension` (A3).
    lattice_vectors = geometry.get("lattice_vectors")
    if not isinstance(lattice_vectors, dict):
        errs.append(f"C2: geometry.lattice_vectors missing or not a mapping (got {lattice_vectors!r})")
    elif dimension is not None:
        expected_names = [f"a{i}" for i in range(1, dimension + 1)]
        if set(lattice_vectors.keys()) != set(expected_names):
            errs.append(
                f"C2: geometry.lattice_vectors must have exactly the named vectors "
                f"{expected_names} for dimension={dimension} "
                f"(got {sorted(lattice_vectors.keys(), key=repr)})"
            )
        for name, vec in lattice_vectors.items():
            if not (isinstance(vec, list) and len(vec) == dimension and all(_is_number(x) for x in vec)):
                errs.append(
                    f"C2: geometry.lattice_vectors[{name!r}] must be a numeric list of "
                    f"length {dimension} (got {vec!r})"
                )

    sites = geometry.get("sites")
    labels: list[str] = []
    if not isinstance(sites, list) or not sites:
        errs.append("C2: geometry.sites missing or empty")
    else:
        for i, s in enumerate(sites):
            if not isinstance(s, dict) or "label" not in s:
                errs.append(f"C2: geometry.sites[{i}] missing 'label'")
                continue
            label = s["label"]
            if not isinstance(label, str):
                errs.append(f"C2: geometry.sites[{i}].label must be a string (got {label!r})")
            labels.append(label)
            frac = s.get("frac")
            if dimension is not None and not (
                isinstance(frac, list) and len(frac) == dimension and all(_is_number(x) for x in frac)
            ):
                errs.append(
                    f"C2: geometry.sites[{i}].frac must be a numeric list of length "
                    f"{dimension} (got {frac!r})"
                )
        # Only string labels are eligible for the duplicate check: a
        # heterogeneous-type labels list (e.g. a stray non-string label,
        # already reported above) must not reach `set()`/`sorted()` on
        # mixed types, which can raise TypeError and mask the real
        # diagnostics behind the generic top-level exception handler.
        str_labels = [lb for lb in labels if isinstance(lb, str)]
        if len(str_labels) != len(set(str_labels)):
            dups = sorted({lb for lb in str_labels if str_labels.count(lb) > 1})
            errs.append(f"C2: duplicate site labels in geometry.sites: {dups}")
        # Sanitize: a non-string label was already diagnosed above; it must
        # not reach C4-C12, which do further set()/sorted() operations on
        # site labels (heterogeneous types there raise TypeError and
        # collapse the diagnostics into a generic top-level exception).
        labels = str_labels
    ctx["labels"] = labels

    size = system.get("size")
    if not isinstance(size, list) or not all(
        isinstance(x, int) and not isinstance(x, bool) for x in size
    ):
        errs.append(f"C2: system.size must be a list of integers (got {size!r})")
        size = None
    ctx["size"] = size

    # --- system.boundary: list of length `dimension` (A3) -----------------
    boundary = system.get("boundary")
    if not isinstance(boundary, list):
        errs.append(f"C2: system.boundary missing or not a list (got {boundary!r})")
    elif dimension is not None and len(boundary) != dimension:
        errs.append(f"C2: len(system.boundary)={len(boundary)} != dimension={dimension}")

    site_dof = model.get("site_dof")
    if not isinstance(site_dof, dict) or not site_dof:
        errs.append("C2: model.site_dof missing or empty")
        site_dof = {}
    else:
        # --- site_dof values: mapping with exactly one of spin/fermion ----
        #     (exclusive — resolves the fermion+spin ambiguity note, A3) ---
        for label, entry in site_dof.items():
            if not isinstance(entry, dict):
                errs.append(f"C2: model.site_dof[{label!r}] must be a mapping (got {entry!r})")
                continue
            present = {k for k in ("spin", "fermion") if k in entry}
            if len(present) != 1:
                errs.append(
                    f"C2: model.site_dof[{label!r}] must have exactly one of "
                    f"'spin'/'fermion' (got {sorted(present)})"
                )
    ctx["site_dof"] = site_dof

    bonds_raw = model.get("bonds")
    if not isinstance(bonds_raw, list):
        errs.append("C2: model.bonds missing or not a list")
        bonds_raw = []
    # Sanitize: only bonds whose type/from/to are strings and whose R is a
    # list of ints are passed on to C3-C12. A malformed bond (e.g. type as
    # a list, from as a dict) is already diagnosed here; letting it reach
    # downstream checks risks unhashable-field TypeErrors (C6's tuple keys,
    # C5's set-building) that would collapse the diagnostics into a
    # generic top-level exception.
    bonds: list[dict] = []
    for i, b in enumerate(bonds_raw):
        if not isinstance(b, dict) or not {"type", "from", "to", "R"} <= b.keys():
            errs.append(f"C2: model.bonds[{i}] missing required keys (type/from/to/R)")
            continue
        valid = True
        for key in ("type", "from", "to"):
            if not isinstance(b[key], str):
                errs.append(f"C2: model.bonds[{i}].{key} must be a string (got {b[key]!r})")
                valid = False
        if not isinstance(b["R"], list) or not all(
            isinstance(x, int) and not isinstance(x, bool) for x in b["R"]
        ):
            errs.append(f"C2: model.bonds[{i}].R must be a list of integers (got {b['R']!r})")
            valid = False
        if valid:
            bonds.append(b)
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
            f"(sites-only: {sorted(labels - set(site_dof), key=repr)}, "
            f"site_dof-only: {sorted(set(site_dof) - labels, key=repr)})"
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
    for t in sorted(used - defined, key=repr):
        errs.append(f"C5: bond type {t!r} has no couplings definition")
    for t in sorted(defined - used, key=repr):
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


# onsite term_name -> required literal tensor_terms[0].coeff sign (C12,
# CONVENTIONS.md §6.2). Term names not listed here are outside this table
# (e.g. any future/unrecognised onsite term) and are left unchecked.
_ONSITE_COEFF_SIGN = {
    "hubbard_u": 1.0,
    "aniso_z": 1.0,
    "chemical_potential": -1.0,
    "field_z": -1.0,
    "field_x": -1.0,
    "field_y": -1.0,
}


def check_c12(doc: dict, ctx: dict) -> list[str]:
    """C12: CONVENTIONS.md §6.2 sign-convention check.

    Only applies to catalog files (``catalog.dialect == "experimental"`` —
    the only dialect this catalog currently defines; kept as an explicit
    guard rather than assumed). Checks, purely as diagnostics (never
    raises):

    - couplings ``hop`` requires ``value.scale == -1.0``.
    - couplings ``density-density`` / ``s_i . S_j`` require ``value.scale``
      absent or ``+1.0``.
    - onsite ``hubbard_u`` / ``aniso_z`` require literal
      ``tensor_terms[0].coeff == +1.0``.
    - onsite ``chemical_potential`` / ``field_z`` / ``field_x`` /
      ``field_y`` require literal ``tensor_terms[0].coeff == -1.0``.
    - J-family ``tensor_terms[*].coeff`` must be a ``{param: ...}`` dict
      *without* ``scale``/``default`` keys (those belong to the resolver,
      not the catalog — CONVENTIONS.md §6.4).
    - non-J coupling ``value`` must be either a ``{param, ...}`` dict, or
      (wannier90 exception only) a plain number.
    """
    errs: list[str] = []
    catalog = doc.get("catalog") if isinstance(doc.get("catalog"), dict) else {}
    if catalog.get("dialect") != "experimental":
        return errs
    lattice = catalog.get("lattice")

    for type_name, c in ctx["couplings"].items():
        if not isinstance(c, dict):
            continue
        tt = _coupling_operator_tensor_terms(c)
        if tt is not _MISSING:
            if isinstance(tt, list):
                for i, term in enumerate(tt):
                    if not isinstance(term, dict):
                        continue
                    coeff = term.get("coeff")
                    if not (
                        isinstance(coeff, dict) and "param" in coeff
                        and "scale" not in coeff and "default" not in coeff
                    ):
                        errs.append(
                            f"C12: {type_name}.tensor_terms[{i}].coeff must be a "
                            f"{{param: ...}} dict without scale/default (got {coeff!r})"
                        )
            continue

        op = c.get("operator")
        value = c.get("value")
        if isinstance(value, dict) and "param" in value:
            scale = value.get("scale", 1.0)
            if op == "hop" and scale != -1.0:
                errs.append(
                    f"C12: {type_name}: operator 'hop' requires value.scale == -1.0 (got {scale!r})"
                )
            elif op in ("density-density", "s_i . S_j") and scale != 1.0:
                errs.append(
                    f"C12: {type_name}: operator {op!r} requires value.scale absent or "
                    f"+1.0 (got {scale!r})"
                )
        elif _is_number(value):
            if lattice != "wannier90":
                errs.append(
                    f"C12: {type_name}.value is a plain number but catalog.lattice != "
                    "'wannier90' (only the wannier90 dialect may use a non-dict value)"
                )
        else:
            errs.append(
                f"C12: {type_name}.value must be a {{param, ...}} dict, or (wannier90 "
                f"exception) a plain number (got {value!r})"
            )

    onsite_raw = doc.get("model", {}).get("onsite") if isinstance(doc.get("model"), dict) else None
    if isinstance(onsite_raw, dict):
        for site, terms in onsite_raw.items():
            if not isinstance(terms, dict):
                continue
            for term_name, term in terms.items():
                expected = _ONSITE_COEFF_SIGN.get(term_name)
                if expected is None or not isinstance(term, dict):
                    continue
                op = term.get("operator")
                tt = op.get("tensor_terms") if isinstance(op, dict) else None
                if not (isinstance(tt, list) and len(tt) == 1 and isinstance(tt[0], dict)):
                    continue
                coeff = tt[0].get("coeff")
                if not (_is_number(coeff) and coeff == expected):
                    errs.append(
                        f"C12: model.onsite[{site!r}][{term_name!r}].operator."
                        f"tensor_terms[0].coeff must be {expected} (got {coeff!r})"
                    )

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

    # --- C11 (A2): each component must be a positive odd int, and must ----
    #     exactly equal 2*max|R component in that direction|+1 as computed
    #     from this file's bonds (the smallest odd integer > 2*max|R|).
    bad_component = False
    for i, comp in enumerate(min_size):
        if not isinstance(comp, int) or isinstance(comp, bool) or comp <= 0 or comp % 2 == 0:
            errs.append(
                f"C11: manifest.min_size_for_check[{i}]={comp!r} must be a positive odd integer"
            )
            bad_component = True
    if bad_component:
        return errs

    max_abs_r = [0] * dimension
    for b in bonds:
        R = b.get("R")
        if isinstance(R, list) and len(R) == dimension:
            for k, x in enumerate(R):
                if isinstance(x, int) and not isinstance(x, bool):
                    max_abs_r[k] = max(max_abs_r[k], abs(x))
    expected_min_size = [2 * m + 1 for m in max_abs_r]
    if min_size != expected_min_size:
        errs.append(
            f"C11: manifest.min_size_for_check={min_size!r} != expected {expected_min_size!r} "
            "(= 2*max|R component|+1 per direction, computed from model.bonds)"
        )
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

    measured, expand_errs = expand_and_count(dimension, labels, bonds, min_size)
    errs.extend(expand_errs)
    expected = entry["coordination"]
    if measured != expected:
        errs.append(f"C11: coordination mismatch: manifest={expected!r} != measured={measured!r}")

    return errs


# ---------------------------------------------------------------------------
#  Per-file / per-document orchestration
# ---------------------------------------------------------------------------

def lint_document(doc: Any, rel_path: str, inventory: set[str],
                   manifest_files: dict) -> list[str]:
    """Run all checks (C1-C12) against a parsed YAML document.

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
        errs.extend(check_c12(doc, ctx))
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


def find_orphan_manifest_entries(files: list[Path], manifest_files: dict) -> list[str]:
    """A6-ii: set-equality check between discovered YAML paths and
    ``manifest.yaml`` keys.

    Parameters
    ----------
    files : list[Path]
        Discovered catalog YAML files (as returned by
        :func:`_discover_default_files`).
    manifest_files : dict
        The ``files`` mapping from ``manifest.yaml``.

    Returns
    -------
    list[str]
        One ``C10`` diagnostic per manifest entry that has no corresponding
        discovered file (an "orphan" manifest entry), sorted by path.
    """
    discovered_rel: set[str] = set()
    for p in files:
        try:
            discovered_rel.add(p.resolve().relative_to(ROOT).as_posix())
        except ValueError:
            pass
    orphans = sorted(set(manifest_files) - discovered_rel)
    return [
        f"C10: manifest.yaml entry {rel!r} has no corresponding catalog YAML file"
        for rel in orphans
    ]


def main() -> None:
    args = sys.argv[1:]
    if args:
        files, top_errors = _resolve_arg_paths(args)
    else:
        files, top_errors = _discover_default_files(), []

    try:
        inventory = load_inventory_keywords()
        manifest_files = load_manifest()
    except FatalLintError as e:
        # A6-i: top-level failures get a clean diagnostic, not a traceback.
        print(f"C2: {e}")
        print("0 files, 1 errors")
        sys.exit(1)

    if not args:
        # A6-ii: default (whole-catalog) run also verifies set equality
        # between discovered YAML paths and manifest.yaml keys.
        top_errors.extend(find_orphan_manifest_entries(files, manifest_files))

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
