# 格子カタログ(新フォーマット対応一式)実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 現行 StdFace の全格子・全相互作用を新三層仕様(geometry/system/model + 拡張方言)の YAML カタログ + 日本語マニュアル + 意味検査リンタとして `lattice_catalog/` に作成する。

**Architecture:** Task 0 で規約・schema・リンタ・manifest を先に確定し、Task 1 (chain) で雛形を固めてから残り格子へ展開。YAML は 31 件((8格子 + ladder W=2/W=3)×3模型 + wannier90 Hubbard 例)。全ファイルがリンタ + manifest 検算 + ソース突合で検証される。

**Tech Stack:** YAML(手書き)、Markdown、Python 3(tools/ の開発時スクリプトのみ。PyYAML 使用可、リポジトリ依存には追加しない)。

**Spec:** `docs/superpowers/specs/2026-08-15-lattice-catalog-design.md`(Codex Round 1 反映版)

## Global Constraints

- 成果物はすべて `lattice_catalog/` 配下。リポジトリの既存コードは一切変更しない。
- マニュアル・YAML コメントは日本語。
- 各 YAML の先頭に `catalog: {schema: stdface-catalog/0.1, dialect: experimental}` を置く。
- ボンド type 名は StdFace キーワード名そのまま。プライムを含む名前(`J0'` 等)は YAML では必ず引用符付き。
- **R の定義**: `R = cell(to) − cell(from)`。実変位 δ = (frac_to − frac_from) + R·A。
- **向きの一意化**: from のサイト番号(`geometry.sites` の配列順)≦ to。同一サイトなら R 辞書順で正(最初の非零成分が正)。各ボンドは一度だけ書く。
- **反転変換**: `(i,j,R) → (j,i,−R)` の際、ホッピング係数は複素共役、スピン交換テンソルは転置(`J_ab → J_ba`)。
- 係数は `value: {param: <StdFaceキーワード>}` / `coeff: {param: ...}` で参照(wannier90 例のみ実数値可)。
- J 族は 9 成分 tensor_terms 正準形で書く(スペック §4.3)。パラメータ解決規則(成分 > スカラー対角 > 大域別名 > 0、同時指定エラー)は CONVENTIONS.md に規範として記載。
- `hop` に符号を埋め込まない(ハミルトニアン寄与は −t·hop、CONVENTIONS に係数規約を明記)。
- twist は度単位、境界 n 回横断で `exp(i·n·π·θ/180)`。`{param: phase0}` 等で参照。
- 出典記録は「ファイル名 + 関数名 + 参照コミットハッシュ」(行番号は補助)。
- 各タスクの完了条件に `python3 lattice_catalog/tools/lint_catalog.py` の成功を含む。
- コミットは Task ごと。コミットメッセージは中立的な表現(内部ツール名を出さない)。
- 参照実装: `python/stdface/lattice/<lattice>.py` の `_BONDS` テーブル。形式は
  `(dW, dL, site_i, site_j, nn_level, J系, t系, V系)`(3D は `(dW, dL, dH, site_i, site_j, J系, t系, V系)`)。
  `dW, dL(, dH)` → R、`site_i/site_j` → from/to ラベル、J系/t系/V系の変数名 → type 名。
- onsite 項の模型別適用範囲(`python/stdface/core/model_plugin.py` で確認済みの規範):
  - Spin: 磁場(−h Sz, −Gamma Sx, −Gamma_y Sy)+ D(Szz)
  - Hubbard: mu(−mu N)+ U(NupNdn)+ 磁場(電子スピンに −h/−Gamma/−Gamma_y)
  - Kondo: Hubbard 一式を `_c` に + **磁場 3 成分を `_s` にも**(両側適用)+ J 結合。D は非 Spin 模型では不使用。

---

### Task 0: 規約・schema・リンタ・manifest・キーワード目録

**Files:**
- Create: `lattice_catalog/CONVENTIONS.md`
- Create: `lattice_catalog/manifest.yaml`(空の骨格 + 形式定義コメント)
- Create: `lattice_catalog/tools/lint_catalog.py`
- Create: `lattice_catalog/tools/keyword_inventory.py`

**Interfaces:**
- Consumes: スペック §4(規約)、`python/stdface/core/keyword_parser.py` の
  `_COMMON_KEYWORDS`、`python/stdface/solvers/*/_plugin.py` のキーワードテーブル
- Produces: 全後続タスクが従う規約と検査基盤。
  `lint_catalog.py` は引数なしで `lattice_catalog/` 全体を検査し、
  違反があれば非零終了。`keyword_inventory.py` は全キーワードの
  JSON 目録を stdout に出力。

- [ ] **Step 1: CONVENTIONS.md を書く**

スペック §4.1–4.5 の規約を規範文書として転記・整理する。章立て:

```markdown
# lattice_catalog 記述規約 (stdface-catalog/0.1)

1. 位置づけと dialect 宣言
   - 各 YAML 先頭: catalog: {schema: stdface-catalog/0.1, dialect: experimental}
   - draft 仕様準拠部分と拡張方言(§7)の区別
2. ファイル構成(自己完結: catalog / geometry / system / model)
3. geometry 規約(サイト順序 = sites 配列順 = 仕様のサイト番号)
4. system 規約(size 代表値、twist: 度単位、exp(i·n·π·θ/180)、{param: phaseN})
5. bonds 規約
   - R = cell(to) − cell(from)、δ = (frac_to − frac_from) + R·A
   - type 名 = StdFace キーワード名(プライムは引用符付き)
   - 向きの一意化(from 番号 ≦ to、R 辞書順で正)
   - 反転変換(hopping: 複素共役、交換テンソル: 転置)
6. couplings / onsite 規約
   - J 族 9 成分 tensor_terms 正準形(成分キーワード名一覧)
   - パラメータ解決規則(成分 > スカラー対角 > 大域別名 > 0、同時指定エラー)
   - 演算子意味論: hop(符号なし、寄与は −t·hop)、density-density、
     s_i . S_j、onsite 語彙(N, Nup, Ndn, NupNdn, Sx, Sy, Sz, Szz)
   - 模型別 onsite 適用範囲(Global Constraints の表と同内容)
   - Kondo 2 ラベル(_c/_s、同一分率座標、物理的には 1 サイト 2 自由度)
7. 拡張方言一覧(fermion site_dof、演算子語彙、{param}参照、catalog ヘッダ)
8. 検算・出典記録の書式(manifest 参照、関数名+コミットハッシュ)
```

- [ ] **Step 2: manifest.yaml の形式を定義し骨格を書く**

```yaml
# lattice_catalog 検算台帳
# 各エントリはリンタが YAML 本体と突合する期待値。
# min_size_for_check: 異なる R が同一サイト対に折り畳まれない最小サイズ
#   (各方向 L > 2 * max|R成分| を満たす奇数を記載)
files: {}
# 記入例(Task 1 で実エントリを追加):
# chain/chain_spin.yaml:
#   lattice: chain
#   model: spin
#   dimension: 1
#   n_sites_uc: 1
#   bonds_per_uc: {J0: 1, "J0'": 1, "J0''": 1}
#   coordination: {J0: 2, "J0'": 2, "J0''": 2}
#   min_size_for_check: [7]
#   source: {file: python/stdface/lattice/chain_lattice.py, func: chain, commit: <hash>}
```

- [ ] **Step 3: lint_catalog.py を書く**

```python
#!/usr/bin/env python3
"""lattice_catalog 意味検査リンタ(開発時専用)。

チェック項目:
  C1  catalog ヘッダ(schema/dialect)の存在と値
  C2  必須トップレベルキー(geometry/system/model)と型
  C3  R の次元 = geometry.dimension、size の次元一致
  C4  bonds の from/to が site_dof・geometry.sites に存在
  C5  bonds の type ↔ couplings のキー整合(未定義参照・未使用定義)
  C6  向き規約(from番号 ≦ to、同一サイトは R 辞書順で正)
  C7  逆向き重複(反転変換 (i,j,R)~(j,i,-R) を同値として重複検出)
  C8  {param: ...} 参照名が keyword_inventory の目録に存在
  C9  tensor_terms の ops 長 = ボンド(2)/onsite(1) のサイト数
  C10 manifest.yaml との突合(bonds_per_uc・n_sites_uc・dimension)
使い方: python3 lattice_catalog/tools/lint_catalog.py [file.yaml ...]
        引数なしで lattice_catalog/**/*.yaml 全件。違反があれば exit 1。
"""
from __future__ import annotations
import glob, json, subprocess, sys
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parent.parent

def load_inventory() -> set[str]:
    out = subprocess.run(
        [sys.executable, str(ROOT / "tools" / "keyword_inventory.py")],
        capture_output=True, text=True, check=True)
    return {e["keyword_canonical"] for e in json.loads(out.stdout)}

def bond_key(b: dict) -> tuple:
    return (b["from"], b["to"], tuple(b["R"]))

def reversed_key(b: dict) -> tuple:
    return (b["to"], b["from"], tuple(-x for x in b["R"]))

def canonical_ok(b: dict, site_order: dict[str, int]) -> bool:
    i, j = site_order[b["from"]], site_order[b["to"]]
    if i < j:
        return True
    if i > j:
        return False
    nz = [x for x in b["R"] if x != 0]
    return bool(nz) and nz[0] > 0

def collect_params(node) -> list[str]:
    found = []
    if isinstance(node, dict):
        if set(node) >= {"param"}:
            found.append(node["param"])
        for v in node.values():
            found += collect_params(v)
    elif isinstance(node, list):
        for v in node:
            found += collect_params(v)
    return found

def lint_file(path: Path, inventory: set[str], manifest: dict) -> list[str]:
    errs = []
    doc = yaml.safe_load(path.read_text())
    rel = str(path.relative_to(ROOT))
    cat = doc.get("catalog", {})
    if cat.get("schema") != "stdface-catalog/0.1":                     # C1
        errs.append("C1: catalog.schema がない/不正")
    geom, model = doc.get("geometry", {}), doc.get("model", {})
    dim = geom.get("dimension")
    labels = [s["label"] for s in geom.get("sites", [])]
    site_order = {lb: i for i, lb in enumerate(labels)}
    if not (dim and labels and model):                                  # C2
        errs.append("C2: geometry/model の必須キー不足")
    if len(doc.get("system", {}).get("size", [])) != dim:               # C3
        errs.append("C3: system.size の次元が dimension と不一致")
    bonds = model.get("bonds", [])
    seen = set()
    for b in bonds:
        if len(b["R"]) != dim:                                          # C3
            errs.append(f"C3: R 次元不一致 {b}")
        if b["from"] not in site_order or b["to"] not in site_order:    # C4
            errs.append(f"C4: 未定義ラベル {b}")
            continue
        if not canonical_ok(b, site_order):                             # C6
            errs.append(f"C6: 向き規約違反 {b}")
        k, rk = bond_key(b), reversed_key(b)
        if k in seen or rk in seen:                                     # C7
            errs.append(f"C7: (逆向き)重複ボンド {b}")
        seen.add(k)
    types_used = {b["type"] for b in bonds}
    types_def = set(model.get("couplings", {}))
    for t in types_used - types_def:                                    # C5
        errs.append(f"C5: coupling 未定義の type {t}")
    for t in types_def - types_used:                                    # C5
        errs.append(f"C5: 未使用の coupling {t}")
    for p in collect_params(doc):                                       # C8
        if p not in inventory:
            errs.append(f"C8: 目録にないパラメータ参照 {p}")
    for scope, arity in (("couplings", 2), ("onsite", 1)):              # C9
        node = model.get(scope, {})
        entries = node.values() if scope == "couplings" else \
            (e for site in node.values() for e in site.values())
        for e in entries:
            op = e.get("operator")
            if isinstance(op, dict):
                for tt in op.get("tensor_terms", []):
                    if len(tt["ops"]) != arity:
                        errs.append(f"C9: ops 長 {len(tt['ops'])} != {arity}")
    m = manifest.get("files", {}).get(rel.removeprefix("lattice_catalog/"))
    if m:                                                               # C10
        if m["n_sites_uc"] != len(labels) or m["dimension"] != dim:
            errs.append("C10: manifest とサイト数/次元が不一致")
        per_uc: dict[str, int] = {}
        for b in bonds:
            per_uc[b["type"]] = per_uc.get(b["type"], 0) + 1
        if per_uc != m["bonds_per_uc"]:
            errs.append(f"C10: bonds_per_uc 不一致 {per_uc} != {m['bonds_per_uc']}")
    else:
        errs.append("C10: manifest 未登録")
    return [f"{rel}: {e}" for e in errs]

def main(argv: list[str]) -> int:
    files = [Path(a) for a in argv] or \
        [Path(p) for p in glob.glob(str(ROOT / "**" / "*.yaml"), recursive=True)
         if "manifest" not in p]
    inventory = load_inventory()
    manifest = yaml.safe_load((ROOT / "manifest.yaml").read_text())
    all_errs = []
    for f in files:
        all_errs += lint_file(f, inventory, manifest)
    for e in all_errs:
        print(e)
    print(f"{len(files)} files, {len(all_errs)} errors")
    return 1 if all_errs else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
```

- [ ] **Step 4: keyword_inventory.py を書く**

```python
#!/usr/bin/env python3
"""StdFace パーサーレジストリ全体からキーワード目録(JSON)を生成する。

母集合: core の _COMMON_KEYWORDS + 各ソルバープラグインのテーブル。
出力: [{"keyword_canonical": "J0x", "keyword_lower": "j0x",
        "source": "common|hphi|mvmc|uhf|hwave"}, ...]
keyword_canonical はカタログの {param: ...} 参照で使う表記
(StdFace ドキュメント慣例: J0x, t0, U, 2S など。lower 形との対応表を持つ)。
"""
from __future__ import annotations
import json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

from stdface.core.keyword_parser import _COMMON_KEYWORDS  # noqa: E402
from stdface.solvers.hphi._plugin import HPhiPlugin        # noqa: E402
from stdface.solvers.mvmc._plugin import MVMCPlugin        # noqa: E402
from stdface.solvers.uhf._plugin import UHFPlugin          # noqa: E402
from stdface.solvers.hwave._plugin import HWavePlugin      # noqa: E402

# 小文字キーワード → カタログ正準表記(機械変換 + 例外表)
_CANON_EXCEPTIONS = {"2s": "2S", "2sz": "2Sz", "gamma": "Gamma",
                     "gamma_y": "Gamma_y", "u": "U", "v": "V", "d": "D",
                     "k": "K", "l": "L", "w": "W", "h": "h", "mu": "mu"}

def canon(kw: str) -> str:
    if kw in _CANON_EXCEPTIONS:
        return _CANON_EXCEPTIONS[kw]
    if kw and kw[0] in "jv" and len(kw) > 1:      # J/V 族は先頭大文字
        return kw[0].upper() + kw[1:]
    return kw

def main() -> None:
    entries = []
    def add(table: dict, source: str) -> None:
        for kw in table:
            entries.append({"keyword_canonical": canon(kw),
                            "keyword_lower": kw, "source": source})
    add(_COMMON_KEYWORDS, "common")
    for plugin_cls, name in ((HPhiPlugin, "hphi"), (MVMCPlugin, "mvmc"),
                             (UHFPlugin, "uhf"), (HWavePlugin, "hwave")):
        table = getattr(plugin_cls, "keyword_table", None) or \
            getattr(plugin_cls(), "keyword_table", {})
        add(dict(table), name)
    json.dump(entries, sys.stdout, ensure_ascii=False, indent=1)

if __name__ == "__main__":
    main()
```

注意: ソルバープラグインのキーワードテーブル属性名は実装を確認して
合わせること(`solvers/hphi/_plugin.py` の `2s` を含むテーブル。
属性名が異なる場合はこのスクリプト側を修正する。既存コードは変更しない)。

- [ ] **Step 5: 動作確認(壊れ例で失敗、直して成功)**

Run: `python3 lattice_catalog/tools/keyword_inventory.py | python3 -c "import json,sys; d=json.load(sys.stdin); print(len(d), 'keywords'); assert any(e['keyword_canonical']=='2S' for e in d)"`
Expected: キーワード件数表示、assert 成功(`2S` がソルバー表から取れている)

一時ファイル `lattice_catalog/chain/chain_spin.yaml` に**故意に壊した**
最小 YAML(couplings 未定義 type、逆向き重複、未知 param)を置き:
Run: `python3 lattice_catalog/tools/lint_catalog.py`
Expected: C5/C7/C8/C10 の違反が報告され exit 1。確認後、一時ファイルを削除。

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/CONVENTIONS.md lattice_catalog/manifest.yaml lattice_catalog/tools/
git commit -m "Add lattice catalog conventions, manifest format, and semantic linter"
```

---

### Task 1: chain(3模型)— 雛形の確定

**Files:**
- Create: `lattice_catalog/chain/chain_spin.yaml`
- Create: `lattice_catalog/chain/chain_hubbard.yaml`
- Create: `lattice_catalog/chain/chain_kondo.yaml`
- Modify: `lattice_catalog/manifest.yaml`(chain 3 エントリ追加)

**Interfaces:**
- Consumes: Task 0 の CONVENTIONS.md・リンタ。
  `python/stdface/lattice/chain_lattice.py::chain`(幾何と `_BONDS`)
- Produces: 後続格子タスクが踏襲する 3 模型分の完全な雛形

chain の抽出済み事実:
- 1 サイト/単位胞、tau=(0,0,0)。論理 1 次元(内部は W=1 の 2D。
  `phase0` 入力は内部で `phase[1]` に転写される — manual 3 章に記載)
- ボンド: `J0/t0/V0` = R=[1]、`J0'/t0'/V0'` = R=[2]、`J0''/t0''/V0''` = R=[3]
- 別名: `J→J0, J'→J0', J''→J0'', t→t0, ..., V''→V0''`

- [ ] **Step 1: chain_spin.yaml を書く**

```yaml
# StdFace chain / Spin 模型
# stan.in 対応例:  model="Spin" / lattice="chain" / L=16 / J=1.0 (J は J0 の別名)
# 検算は manifest.yaml の chain/chain_spin.yaml エントリを参照。
# 出典: python/stdface/lattice/chain_lattice.py::chain (_BONDS), commit <hash>

catalog: {schema: stdface-catalog/0.1, dialect: experimental}

geometry:
  dimension: 1
  lattice_vectors:
    a1: [1.0]                    # StdFace: a (既定 1.0)
  sites:
    - {label: A, frac: [0.0]}

system:
  finite: true
  size: [16]                     # StdFace: L (代表値)
  boundary:
    - {twist: {param: phase0}}   # 度単位。phase0=0 で periodic と等価

model:
  site_dof:
    A: {spin: {param: 2S, default: 0.5}}   # S = 2S/2。既定 2S=1

  bonds:
    - {from: A, to: A, R: [1], type: J0}     # 最近接   (別名 J)
    - {from: A, to: A, R: [2], type: "J0'"}  # 次近接   (別名 J')
    - {from: A, to: A, R: [3], type: "J0''"} # 三次近接 (別名 J'')

  couplings:
    # J 族は 9 成分正準形。パラメータ解決: 成分 > スカラー対角 > 大域別名 > 0
    J0:
      operator:
        tensor_terms:
          - {ops: [Sx, Sx], coeff: {param: J0x}}
          - {ops: [Sy, Sy], coeff: {param: J0y}}
          - {ops: [Sz, Sz], coeff: {param: J0z}}
          - {ops: [Sx, Sy], coeff: {param: J0xy}}
          - {ops: [Sx, Sz], coeff: {param: J0xz}}
          - {ops: [Sy, Sx], coeff: {param: J0yx}}
          - {ops: [Sy, Sz], coeff: {param: J0yz}}
          - {ops: [Sz, Sx], coeff: {param: J0zx}}
          - {ops: [Sz, Sy], coeff: {param: J0zy}}
    "J0'":
      operator:
        tensor_terms:
          - {ops: [Sx, Sx], coeff: {param: "J0'x"}}
          - {ops: [Sy, Sy], coeff: {param: "J0'y"}}
          - {ops: [Sz, Sz], coeff: {param: "J0'z"}}
          - {ops: [Sx, Sy], coeff: {param: "J0'xy"}}
          - {ops: [Sx, Sz], coeff: {param: "J0'xz"}}
          - {ops: [Sy, Sx], coeff: {param: "J0'yx"}}
          - {ops: [Sy, Sz], coeff: {param: "J0'yz"}}
          - {ops: [Sz, Sx], coeff: {param: "J0'zx"}}
          - {ops: [Sz, Sy], coeff: {param: "J0'zy"}}
    "J0''":
      operator:
        tensor_terms:
          - {ops: [Sx, Sx], coeff: {param: "J0''x"}}
          - {ops: [Sy, Sy], coeff: {param: "J0''y"}}
          - {ops: [Sz, Sz], coeff: {param: "J0''z"}}
          - {ops: [Sx, Sy], coeff: {param: "J0''xy"}}
          - {ops: [Sx, Sz], coeff: {param: "J0''xz"}}
          - {ops: [Sy, Sx], coeff: {param: "J0''yx"}}
          - {ops: [Sy, Sz], coeff: {param: "J0''yz"}}
          - {ops: [Sz, Sx], coeff: {param: "J0''zx"}}
          - {ops: [Sz, Sy], coeff: {param: "J0''zy"}}

  onsite:
    A:
      field_z:                   # -h Sz
        operator: {tensor_terms: [{ops: [Sz], coeff: -1.0}]}
        value: {param: h}
      field_x:                   # -Gamma Sx
        operator: {tensor_terms: [{ops: [Sx], coeff: -1.0}]}
        value: {param: Gamma}
      field_y:                   # -Gamma_y Sy
        operator: {tensor_terms: [{ops: [Sy], coeff: -1.0}]}
        value: {param: Gamma_y}
      aniso_z:                   # D (Sz)^2
        operator: {tensor_terms: [{ops: [Szz], coeff: 1.0}]}
        value: {param: D}
```

- [ ] **Step 2: chain_hubbard.yaml を書く**

```yaml
# StdFace chain / Hubbard 模型
# stan.in 対応例:  model="Hubbard" / lattice="chain" / L=16 / t=1.0 / U=4.0
# 検算は manifest.yaml を参照。ハミルトニアン寄与は -t·hop (CONVENTIONS 6章)
# 出典: python/stdface/lattice/chain_lattice.py::chain (_BONDS), commit <hash>

catalog: {schema: stdface-catalog/0.1, dialect: experimental}

geometry:
  dimension: 1
  lattice_vectors:
    a1: [1.0]                    # StdFace: a
  sites:
    - {label: A, frac: [0.0]}

system:
  finite: true
  size: [16]                     # StdFace: L
  boundary:
    - {twist: {param: phase0}}   # 度単位

model:
  site_dof:
    A: {fermion: {orbitals: 1}}  # 拡張方言 (CONVENTIONS 7章)

  bonds:
    - {from: A, to: A, R: [1], type: t0}     # 別名 t
    - {from: A, to: A, R: [2], type: "t0'"}  # 別名 t'
    - {from: A, to: A, R: [3], type: "t0''"} # 別名 t''
    - {from: A, to: A, R: [1], type: V0}     # 別名 V
    - {from: A, to: A, R: [2], type: "V0'"}  # 別名 V'
    - {from: A, to: A, R: [3], type: "V0''"} # 別名 V''

  couplings:
    t0:    {operator: "hop", value: {param: t0}}     # 複素可。h.c. は反転変換で共役
    "t0'": {operator: "hop", value: {param: "t0'"}}
    "t0''": {operator: "hop", value: {param: "t0''"}}
    V0:    {operator: "density-density", value: {param: V0}}
    "V0'": {operator: "density-density", value: {param: "V0'"}}
    "V0''": {operator: "density-density", value: {param: "V0''"}}

  onsite:
    A:
      hubbard_u:                 # U n↑n↓
        operator: {tensor_terms: [{ops: [NupNdn], coeff: 1.0}]}
        value: {param: U}
      chemical_potential:        # -mu N
        operator: {tensor_terms: [{ops: [N], coeff: -1.0}]}
        value: {param: mu}
      field_z:                   # -h Sz (電子スピン)
        operator: {tensor_terms: [{ops: [Sz], coeff: -1.0}]}
        value: {param: h}
      field_x:
        operator: {tensor_terms: [{ops: [Sx], coeff: -1.0}]}
        value: {param: Gamma}
      field_y:
        operator: {tensor_terms: [{ops: [Sy], coeff: -1.0}]}
        value: {param: Gamma_y}
```

- [ ] **Step 3: chain_kondo.yaml を書く**

```yaml
# StdFace chain / Kondo 格子模型
# stan.in 対応例:  model="Kondo" / lattice="chain" / L=16 / t=1.0 / J=1.0
# 磁場 (h/Gamma/Gamma_y) は遍歴電子と局在スピンの両方に適用される
# (python/stdface/core/model_plugin.py::KondoModel)。
# 現行実装のサイト倍加 (前半=局在 S2, 後半=遍歴) との対応は manual 4章。
# 出典: python/stdface/lattice/chain_lattice.py::chain (_BONDS), commit <hash>

catalog: {schema: stdface-catalog/0.1, dialect: experimental}

geometry:
  dimension: 1
  lattice_vectors:
    a1: [1.0]                    # StdFace: a
  sites:
    - {label: A_c, frac: [0.0]}  # 遍歴電子 (1物理サイト上の自由度その1)
    - {label: A_s, frac: [0.0]}  # 局在スピン (同・その2)

system:
  finite: true
  size: [16]                     # StdFace: L
  boundary:
    - {twist: {param: phase0}}   # 度単位

model:
  site_dof:
    A_c: {fermion: {orbitals: 1}}
    A_s: {spin: {param: 2S, default: 0.5}}

  bonds:
    - {from: A_c, to: A_s, R: [0], type: J}   # Kondo 結合 (セル内)
    - {from: A_c, to: A_c, R: [1], type: t0}
    - {from: A_c, to: A_c, R: [2], type: "t0'"}
    - {from: A_c, to: A_c, R: [3], type: "t0''"}
    - {from: A_c, to: A_c, R: [1], type: V0}
    - {from: A_c, to: A_c, R: [2], type: "V0'"}
    - {from: A_c, to: A_c, R: [3], type: "V0''"}

  couplings:
    t0:    {operator: "hop", value: {param: t0}}
    "t0'": {operator: "hop", value: {param: "t0'"}}
    "t0''": {operator: "hop", value: {param: "t0''"}}
    V0:    {operator: "density-density", value: {param: V0}}
    "V0'": {operator: "density-density", value: {param: "V0'"}}
    "V0''": {operator: "density-density", value: {param: "V0''"}}
    J:     {operator: "s_i . S_j", value: {param: J}}

  onsite:
    A_c:
      hubbard_u:
        operator: {tensor_terms: [{ops: [NupNdn], coeff: 1.0}]}
        value: {param: U}
      chemical_potential:
        operator: {tensor_terms: [{ops: [N], coeff: -1.0}]}
        value: {param: mu}
      field_z:
        operator: {tensor_terms: [{ops: [Sz], coeff: -1.0}]}
        value: {param: h}
      field_x:
        operator: {tensor_terms: [{ops: [Sx], coeff: -1.0}]}
        value: {param: Gamma}
      field_y:
        operator: {tensor_terms: [{ops: [Sy], coeff: -1.0}]}
        value: {param: Gamma_y}
    A_s:
      field_z:                   # 局在側にも同じ 3 成分
        operator: {tensor_terms: [{ops: [Sz], coeff: -1.0}]}
        value: {param: h}
      field_x:
        operator: {tensor_terms: [{ops: [Sx], coeff: -1.0}]}
        value: {param: Gamma}
      field_y:
        operator: {tensor_terms: [{ops: [Sy], coeff: -1.0}]}
        value: {param: Gamma_y}
```

注意: chain_kondo の bonds で `A_c→A_s` R=[0] はサイト番号 0→1 なので
向き規約に適合。J は等方スカラーのみ(Kondo の `input_spin` はスカラー
入力)だが、成分キーワード(`Jx` 等)が共通テーブルに存在するため
`{param: J}` のスカラー参照とし、その旨コメントする。

- [ ] **Step 4: manifest.yaml に chain 3 エントリを追加**

commit ハッシュは `git rev-parse HEAD` の値を使用:

```yaml
files:
  chain/chain_spin.yaml:
    lattice: chain
    model: spin
    dimension: 1
    n_sites_uc: 1
    bonds_per_uc: {J0: 1, "J0'": 1, "J0''": 1}
    coordination: {J0: 2, "J0'": 2, "J0''": 2}
    min_size_for_check: [7]
    source: {file: python/stdface/lattice/chain_lattice.py, func: chain, commit: <hash>}
  chain/chain_hubbard.yaml:
    lattice: chain
    model: hubbard
    dimension: 1
    n_sites_uc: 1
    bonds_per_uc: {t0: 1, "t0'": 1, "t0''": 1, V0: 1, "V0'": 1, "V0''": 1}
    coordination: {t0: 2, V0: 2}
    min_size_for_check: [7]
    source: {file: python/stdface/lattice/chain_lattice.py, func: chain, commit: <hash>}
  chain/chain_kondo.yaml:
    lattice: chain
    model: kondo
    dimension: 1
    n_sites_uc: 2
    bonds_per_uc: {J: 1, t0: 1, "t0'": 1, "t0''": 1, V0: 1, "V0'": 1, "V0''": 1}
    coordination: {t0: 2, J: 1}
    min_size_for_check: [7]
    source: {file: python/stdface/lattice/chain_lattice.py, func: chain, commit: <hash>}
```

- [ ] **Step 5: リンタ実行**

Run: `python3 lattice_catalog/tools/lint_catalog.py`
Expected: `3 files, 0 errors` で exit 0

- [ ] **Step 6: ソース突合**

`chain_lattice.py::chain` の `_BONDS` 3 行と YAML bonds を突合
((0,1,...)→R=[1] 等)。onsite は `model_plugin.py` の
SpinModel/HubbardModel/KondoModel と突合(特に Kondo 両側磁場)。
不一致があれば YAML を修正して Step 5 からやり直す。

- [ ] **Step 7: Commit**

```bash
git add lattice_catalog/chain/ lattice_catalog/manifest.yaml
git commit -m "Add chain lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 2: square(3模型)

**Files:**
- Create: `lattice_catalog/square/square_{spin,hubbard,kondo}.yaml`
- Modify: `lattice_catalog/manifest.yaml`

**Interfaces:**
- Consumes: Task 0 規約、Task 1 の chain 雛形(全面踏襲)
- Produces: なし(独立成果物)

- [ ] **Step 1: `python/stdface/lattice/square_lattice.py` から抽出**

確定する項目: `NsiteUC`/`tau`/`direct` 既定(幾何)、`_BONDS` 全行
(R とサイト対と type 名)、模型別 `not_used_*`(YAML に含めない項)。
J0/J1 の方向割当(W方向/L方向)はソースコメントと
`input_spin_nn` の名前引数で確認。

- [ ] **Step 2: 3 つの YAML を書く**(chain 雛形と同一構造)

相違点: `dimension: 2`、`lattice_vectors` 2 本(ソース既定値)、
`system.size: [4, 4]  # StdFace: W, L`、
`boundary: [{twist: {param: phase0}}, {twist: {param: phase1}}]`、
bonds/couplings は抽出結果を全数記載(J 族は 9 成分正準形)。
Kondo は A_c/A_s の 2 ラベル、磁場両側適用。

- [ ] **Step 3: manifest エントリ追加**(3 件。coordination: 最近接 4
(J0+J1 各 2)、対角 `J'` 系 4。min_size_for_check は R 成分最大値から算出)

- [ ] **Step 4: リンタ実行**

Run: `python3 lattice_catalog/tools/lint_catalog.py`
Expected: `6 files, 0 errors`

- [ ] **Step 5: ソース突合**(`_BONDS` 全行、onsite)

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/square/ lattice_catalog/manifest.yaml
git commit -m "Add square lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 3: triangular(3模型)

**Files:**
- Create: `lattice_catalog/triangular/triangular_{spin,hubbard,kondo}.yaml`
- Modify: `lattice_catalog/manifest.yaml`

**Interfaces:** Task 2 と同じ(Consumes: Task 0/1)

- [ ] **Step 1: `python/stdface/lattice/triangular_lattice.py` から抽出**
(Task 2 Step 1 と同じ項目。格子ベクトル既定 a1=(1,0), a2=(1/2,√3/2) を
ソースで確認。t/J/V は 0,1,2 系 + p/pp)

- [ ] **Step 2: 3 つの YAML を書く**(検算: 最近接配位数 6 = J0/J1/J2 各 2)

- [ ] **Step 3: manifest エントリ追加**

- [ ] **Step 4: リンタ実行** — Expected: `9 files, 0 errors`

- [ ] **Step 5: ソース突合**

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/triangular/ lattice_catalog/manifest.yaml
git commit -m "Add triangular lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 4: honeycomb(3模型)

**Files:**
- Create: `lattice_catalog/honeycomb/honeycomb_{spin,hubbard,kondo}.yaml`
- Modify: `lattice_catalog/manifest.yaml`

**Interfaces:** Task 2 と同じ

- [ ] **Step 1: `python/stdface/lattice/honeycomb_lattice.py` から抽出**
(2 サイト A/B。J0/J1/J2 = 3 方向ボンド(Kitaev 慣例)の対応をコメントに
明記。Kondo は A_c/A_s/B_c/B_s の 4 ラベル)

- [ ] **Step 2: 3 つの YAML を書く**(検算: 最近接配位数 3)

- [ ] **Step 3: manifest エントリ追加**

- [ ] **Step 4: リンタ実行** — Expected: `12 files, 0 errors`

- [ ] **Step 5: ソース突合**

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/honeycomb/ lattice_catalog/manifest.yaml
git commit -m "Add honeycomb lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 5: kagome(3模型)

**Files:**
- Create: `lattice_catalog/kagome/kagome_{spin,hubbard,kondo}.yaml`
- Modify: `lattice_catalog/manifest.yaml`

**Interfaces:** Task 2 と同じ

- [ ] **Step 1: `python/stdface/lattice/kagome.py` から抽出**
(確認済み: 3 サイト、tau=(0,0),(1/2,0),(0,1/2)。`_BONDS` は
`kagome.py::kagome`。p 系まで、pp なし。検算: 最近接配位数 4。
Kondo は 6 ラベル)

- [ ] **Step 2: 3 つの YAML を書く**

- [ ] **Step 3: manifest エントリ追加**

- [ ] **Step 4: リンタ実行** — Expected: `15 files, 0 errors`

- [ ] **Step 5: ソース突合**

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/kagome/ lattice_catalog/manifest.yaml
git commit -m "Add kagome lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 6: ladder(W=2 / W=3 × 3模型 = 6件)

**Files:**
- Create: `lattice_catalog/ladder/ladder_w2_{spin,hubbard,kondo}.yaml`
- Create: `lattice_catalog/ladder/ladder_w3_{spin,hubbard,kondo}.yaml`
- Modify: `lattice_catalog/manifest.yaml`

**Interfaces:** Task 2 と同じ

- [ ] **Step 1: `python/stdface/lattice/ladder.py` から抽出**

`_BONDS` は W に依存して動的構築(`ladder.py::ladder`)。ループを読み、
脚内(J1, J1')・ラング(J0)・斜め(J2, J2')の割当を W の関数として
文書化した上で、W=2 と W=3 の具体的ボンド集合を書き下す。
単位胞 = 1 ラング(W サイト、ラベル A0..A{W-1})、`size: [16]  # StdFace: L`
の論理 1 次元表現。rung 方向は open(セル内結合のみ)である根拠を
ソースで確認して記載。

- [ ] **Step 2: 6 つの YAML を書く**

W 一般化規則(ラベルとボンドの生成則)を各ファイル冒頭コメントに記載。
検算(W=2): 脚方向配位数 2、ラング 1。

- [ ] **Step 3: manifest エントリ追加**(6 件)

- [ ] **Step 4: リンタ実行** — Expected: `21 files, 0 errors`

- [ ] **Step 5: ソース突合**(W=2/W=3 それぞれループを手で展開して照合)

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/ladder/ lattice_catalog/manifest.yaml
git commit -m "Add ladder (W=2, W=3) lattice catalog definitions"
```

---

### Task 7: orthorhombic(3模型)

**Files:**
- Create: `lattice_catalog/orthorhombic/orthorhombic_{spin,hubbard,kondo}.yaml`
- Modify: `lattice_catalog/manifest.yaml`

**Interfaces:**
- Consumes: Task 0/1
- Produces: 3D 格子の雛形(Task 8, 9 が参照)

- [ ] **Step 1: `python/stdface/lattice/orthorhombic.py` から抽出**
(3D: `dimension: 3`、`size: [4,4,4]  # W, L, Height`、boundary 3 方向
`{twist: {param: phase0/1/2}}`。`_BONDS` は `(dW,dL,dH,...)` 形式。
J0/J1/J2 = 3 軸 nn、pp の割当をソースで確認。検算: 最近接配位数 6)

- [ ] **Step 2: 3 つの YAML を書く**

- [ ] **Step 3: manifest エントリ追加**

- [ ] **Step 4: リンタ実行** — Expected: `24 files, 0 errors`

- [ ] **Step 5: ソース突合**

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/orthorhombic/ lattice_catalog/manifest.yaml
git commit -m "Add orthorhombic lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 8: fc_ortho(3模型)

**Files:**
- Create: `lattice_catalog/fc_ortho/fc_ortho_{spin,hubbard,kondo}.yaml`
- Modify: `lattice_catalog/manifest.yaml`

**Interfaces:** Consumes: Task 0/1/7(3D 雛形)

- [ ] **Step 1: `python/stdface/lattice/fc_ortho.py` から抽出**
(面心系の格子ベクトル既定を確認。J0/J1/J2 の方向割当と p/pp の範囲)

- [ ] **Step 2: 3 つの YAML を書く**

- [ ] **Step 3: manifest エントリ追加**

- [ ] **Step 4: リンタ実行** — Expected: `27 files, 0 errors`

- [ ] **Step 5: ソース突合**

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/fc_ortho/ lattice_catalog/manifest.yaml
git commit -m "Add face-centered orthorhombic lattice catalog definitions"
```

---

### Task 9: pyrochlore(3模型)

**Files:**
- Create: `lattice_catalog/pyrochlore/pyrochlore_{spin,hubbard,kondo}.yaml`
- Modify: `lattice_catalog/manifest.yaml`

**Interfaces:** Consumes: Task 0/1/7

- [ ] **Step 1: `python/stdface/lattice/pyrochlore.py` から抽出**
(4 副格子、FCC 格子ベクトル既定。四面体内(R=0)と四面体間の区別を
コメントに明記。検算: 最近接配位数 6)

- [ ] **Step 2: 3 つの YAML を書く**

**Kondo の特記事項(必須)**: 現行実装(C/Python 共通、
`src/Pyrochlore.c` の GeneralJ 呼び出しおよび
`python/stdface/lattice/pyrochlore.py::_local`)は Kondo 結合を
**副格子 3 の遍歴サイトと全副格子(0–3)の局在スピン**の間に生成する。
カタログはこの現行動作を忠実に再現する:

```yaml
  bonds:
    # Kondo 結合: 現行実装は A3_c と全局在スピンを結合する
    # (他格子の副格子ごと対応と異なる。上流実装のバグの可能性あり —
    #  manual 7章参照。ここでは現行動作を忠実に記載)
    - {from: A0_s, to: A3_c, R: [0, 0, 0], type: J}
    - {from: A1_s, to: A3_c, R: [0, 0, 0], type: J}
    - {from: A2_s, to: A3_c, R: [0, 0, 0], type: J}
    - {from: A3_c, to: A3_s, R: [0, 0, 0], type: J}
```

(from/to の向きはサイト番号順の規約に合わせて調整。`s_i . S_j` の
第 1 引数が遍歴側になるよう couplings 側の定義と整合させること。)

- [ ] **Step 3: manifest エントリ追加**

- [ ] **Step 4: リンタ実行** — Expected: `30 files, 0 errors`

- [ ] **Step 5: ソース突合**(Kondo 結合 4 本の対応を C/Python 両方で確認)

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/pyrochlore/ lattice_catalog/manifest.yaml
git commit -m "Add pyrochlore lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 10: wannier90 Hubbard 小例

**Files:**
- Create: `lattice_catalog/wannier90/example_hubbard.yaml`
- Modify: `lattice_catalog/manifest.yaml`

**Interfaces:**
- Consumes: Task 0 規約。`python/stdface/lattice/wannier90.py`
  (変換アルゴリズム)、`test/wannier90_data/` または `samples/`(入力データ)
- Produces: manual 5 章(Task 12)が参照する具体例

- [ ] **Step 1: 変換アルゴリズムを精読して変換規則を確定**

`wannier90.py` から以下を確定(manual 5 章の下書きを兼ねる):
- H_mn(R) → hop の符号規約(現行実装の符号反転の有無と位置)
- R=0 & m=n の対角要素 → onsite 一体項としての分離
- Hermite 対 `(m,n,R)/(n,m,−R)` の正準対選択(現行実装がどちらを
  採用しているか)と縮退重みの扱い
- cutoff(`cutoff_t/u/j`, `cutoff_*R`, `cutoff_length_*`)の適用順
- `lambda_U/lambda_J/alpha/doublecounting` の作用点

- [ ] **Step 2: 入力データを選ぶ**

`test/wannier90_data/` と `samples/` から軌道数最小の Hubbard 用
データセット(_hr.dat / _geom / _ur.dat)を選ぶ。

- [ ] **Step 3: example_hubbard.yaml を書く**

- `geometry`: _geom の格子ベクトルと Wannier 中心 → sites(W1, W2, …)
- `bonds` + `couplings`: 選んだ cutoff の**範囲内全要素**を Step 1 の
  規則で変換(上位 N 件抽出ではない)。正準対のみ列挙し、
  排除した Hermite 対の扱いをコメントで明記。R=0 対角は onsite へ。
- 値は実数値を直接記載(外部データ由来。param 参照不可の旨コメント)。
- 各要素に元ファイルの行番号をコメントで記録。

- [ ] **Step 4: manifest エントリ追加 + リンタ実行**

Run: `python3 lattice_catalog/tools/lint_catalog.py`
Expected: `31 files, 0 errors`
(リンタの C8 は数値直書きを許すこと — Task 0 の実装で
`{param: ...}` 形式のみ検査対象とし、数値はスキップされる)

- [ ] **Step 5: 検算** — 記載全要素を元データと突合、変換規則(符号・
正準対・onsite 分離)の適用を 2 要素以上で手計算検証しコメントに残す

- [ ] **Step 6: Commit**

```bash
git add lattice_catalog/wannier90/ lattice_catalog/manifest.yaml
git commit -m "Add wannier90 Hubbard mapping example for the lattice catalog"
```

---

### Task 11: manual.md 前半(1–2章: 概要・キーワード対応表)

**Files:**
- Create: `lattice_catalog/manual.md`(1–2 章)

**Interfaces:**
- Consumes: `tools/keyword_inventory.py` の出力、CONVENTIONS.md
- Produces: manual.md の章立て(Task 12, 13 が追記)

- [ ] **Step 1: 1 章(概要と読み方)を書く**

内容: 三層仕様の要約、本カタログの位置づけ(experimental dialect、
「実行検証済みではなくソース突合+リンタ検査」)、ファイル構成一覧、
CONVENTIONS.md の規約の詳説、`{param: ...}` とパラメータ解決規則の読み方。

- [ ] **Step 2: 2 章(キーワード対応表)を書く**

Run: `python3 lattice_catalog/tools/keyword_inventory.py > /tmp/inv.json`
の出力を母集合とし、**全エントリ**を分類した表を作る:

| StdFace キーワード | 出典テーブル | 新フォーマットでの対応 | 備考 |

分類: (a) geometry(`a`, `wx…`, `wlength…`)、(b) system(`W/L/Height`,
`phase0-2`, `box`→supercell(S) 写像)、(c) bonds+couplings(t/J/V 族全
変種、J 成分 9 種×全 prefix)、(d) onsite(`U`,`mu`,`D`,`h`,`Gamma`,
`Gamma_y`)、(e) site_dof(`2S`)、(f) wannier90(`cutoff_*`,`lambda*`,
`alpha`,`doublecounting` → 5 章参照)、(g) 対象外(`model`,`lattice`,
`outputmode`,`ncond/nelec`,`2Sz`,`K` ほかソルバー固有の計算条件 —
理由を各行に明記)

- [ ] **Step 3: 網羅性チェック**

表の行数と inventory のエントリ数を照合(スクリプト出力件数と一致、
欠落ゼロ)。カタログ YAML 中の全 `{param: ...}` が表の (a)–(f) の
いずれかに載っていることをリンタ C8 で再確認。

- [ ] **Step 4: Commit**

```bash
git add lattice_catalog/manual.md
git commit -m "Add lattice catalog manual: overview and keyword correspondence table"
```

---

### Task 12: manual.md 中盤(3–5章: 格子・模型・wannier90)

**Files:**
- Modify: `lattice_catalog/manual.md`(3–5 章を追記)

**Interfaces:**
- Consumes: Task 1–10 の全 YAML・manifest・各タスクの突合記録
- Produces: なし

- [ ] **Step 1: 3 章(格子ごとの解説)を書く**

9 格子それぞれ 1 節: 幾何(格子ベクトル・副格子表)、ボンド定義表
(type / from–to / R / 意味)、検算(manifest 値の根拠)、出典
(ファイル+関数+コミット)。chain の内部 2D 表現(phase0→phase[1])、
ladder の W 一般化規則もここに記載。

- [ ] **Step 2: 4 章(模型ごとの演算子対応)を書く**

- Spin: 9 成分正準形とパラメータ解決規則(成分キーワード 9 種の表)、
  等方入力・異方入力・同時指定エラーの 3 例。
- Hubbard: hop / density-density / onsite の定義(数式)、係数規約
  (−t·hop)、複素 t と反転変換。
- Kondo: 2 ラベル表現、磁場両側適用、サイト倍加(前半=局在)との対応、
  GC 変種(粒子数条件は solver 層、カタログ対象外)。

- [ ] **Step 3: 5 章(wannier90 変換仕様)を書く**

Task 10 Step 1 で確定した規則を模型別に文書化:
- Hubbard: H/U/J 各チャネルの変換(符号、onsite 分離、Hermite 正準対、
  縮退重み、cutoff、λ/α、doublecounting)
- Spin: 超交換 `2|t_mn|^2(1/U_m+1/U_n)` 生成の現行アルゴリズムと、
  直接写像との差異
- example_hubbard.yaml の読み解き

- [ ] **Step 4: Commit**

```bash
git add lattice_catalog/manual.md
git commit -m "Add lattice catalog manual: per-lattice, per-model, and wannier90 chapters"
```

---

### Task 13: manual.md 後半(6–7章)+ README + 最終検証

**Files:**
- Modify: `lattice_catalog/manual.md`(6–7 章を追記)
- Create: `lattice_catalog/README.md`

**Interfaces:**
- Consumes: 全タスクの成果物
- Produces: 完成した一式

- [ ] **Step 1: 6 章(仕様拡張提案)を書く**

スペック §4.5 の 5 項目(fermion site_dof、1 サイト演算子語彙、
名前付き 2 体演算子、`{param}` 参照の許容箇所、catalog ヘッダ)を、
それぞれ tensor_terms への展開形・型・意味論つきで draft への
追記提案として記述。

- [ ] **Step 2: 7 章(既知の制限・未対応事項)を書く**

- oracle 比較(数値同値性)未実施 — 展開エンジン実装時の課題
- ladder は W=2/3 の例のみ(一般 W は生成規則の文書)
- wannier90 は Hubbard 例のみ、Spin/Kondo は仕様記述のみ
- pyrochlore Kondo の `isite+3` 挙動(上流バグ疑い、現行動作を記載)
- box 行列は写像説明のみで YAML 例なし
- 機械的識別子(prime を含む type 名)の分離は仕様確定時の課題
- 多体項(3 体以上)は現行 StdFace に存在しないため対象外

- [ ] **Step 3: README.md を書く**

一式の目的(1 段落)、ディレクトリ構成、リンタの使い方、manual.md への
案内、参照仕様(格子定義仕様 draft 2026/07/27)への言及。

- [ ] **Step 4: 最終検証**

Run: `python3 lattice_catalog/tools/lint_catalog.py`
Expected: `31 files, 0 errors`

Run: `python3 lattice_catalog/tools/keyword_inventory.py | python3 -c "import json,sys; print(len(json.load(sys.stdin)),'keywords')"`
の件数と manual 2 章の表の行数が一致することを確認。

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/
git commit -m "Complete lattice catalog manual and README"
```
