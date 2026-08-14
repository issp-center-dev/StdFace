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
- **向き = ソース順保持**: bonds は参照実装の向きのまま書く(`_BONDS` の `site_i → site_j`、Kondo は `general_j(..., 遍歴, 局在)` の引数順)。正準形への並べ替えは**しない**。各ボンドは一度だけ書き、リンタが反転同値 `(i,j,R)≡(j,i,−R)` の重複を検出する。反転時の係数変換(hopping: 複素共役、交換テンソル: 転置)は消費側規範として CONVENTIONS に記載。
- **value 意味論**: `H = Σ value·operator` の物理ハミルトニアン係数(solver 出力の係数規約とは別物 — スペック §4.3 の検証連鎖)。**符号は必ずデータに持たせる**: param 参照の一般形は `{param: <名>, scale: <実数, 省略時1.0>, default: <最終値, 省略時0>}`、値 = scale × param。符号表(スペック §4.3): t 族 scale −1、mu/h/Gamma/Gamma_y coeff −1、U/V/D/J 族/Kondo J は +。
- J 族は 9 成分 tensor_terms 正準形。**パラメータ解決順**(`input_params.py::_resolve_spin_matrix` と同一): 成分局所 > 成分大域 > スカラー局所(対角) > スカラー大域(対角) > 0。競合規則(スカラー同士・スカラーvs行列・行列vs行列)とプライム系(大域 fallback なし)の解決表を CONVENTIONS に規範として記載。
- twist は度単位、境界 n 回横断で `exp(i·n·π·θ/180)`。`{twist: {param: phase0}}` 等で参照。
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
- Create: `lattice_catalog/tools/test_tools.py`(リンタ・目録の自動テスト)

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
   - 向き = ソース順保持(_BONDS の site_i→site_j、Kondo は遍歴→局在)。
     正準化しない。一意性は反転同値 (i,j,R)≡(j,i,−R) で判定
   - 反転時の係数変換(消費側規範): hopping 複素共役、交換テンソル転置
6. couplings / onsite 規約
   - value 意味論: H = Σ value·operator(物理ハミルトニアン係数)。
     solver 出力係数規約(HPhi trans.def の H=−Σt 暗黙符号)との違いと
     検証連鎖(パラメータ→builder→trans/intr→solver規約→物理符号)
   - 符号表(スペック §4.3 の表を転記): t 族 scale −1、mu/磁場 coeff −1、
     U/V/D/J/Kondo J は +。符号は必ずデータ(scale/coeff)に置く
   - param 参照の一般形 {param, scale(省略時1), default(省略時0)}、
     値 = scale × param。型制約(2S: 正整数 等)
   - J 族 9 成分 tensor_terms 正準形(成分キーワード名一覧、
     ops 対 ↔ 接尾辞対応表 [Sx,Sy]↔xy 等)
   - パラメータ解決順序(実装準拠): 成分局所 > 成分大域 >
     スカラー局所(対角) > スカラー大域(対角) > 0。
     競合規則(スカラー同士/スカラーvs行列/行列vs行列)と
     プライム系(大域 fallback なし)の prefix 別解決表
   - 演算子意味論と端点順序: hop(符号なし)、density-density、
     s_i . S_j(第1端点=遍歴)、onsite 語彙(N, Nup, Ndn, NupNdn,
     Sx, Sy, Sz, Szz)
   - 演算子と site_dof の型整合表(hop/density-density: fermion–fermion、
     s_i . S_j: fermion–spin この順、J テンソル: spin–spin)
   - 模型別 onsite 適用範囲(Global Constraints の表と同内容)
   - Kondo 2 ラベル(_c/_s、同一分率座標、物理的には 1 サイト 2 自由度)
7. 拡張方言一覧(fermion site_dof、演算子語彙、{param,scale,default}参照、
   catalog ヘッダ、4フェルミオン一般項の素描(将来課題))
8. 検算・出典記録の書式(manifest 参照、関数名+コミットハッシュ)
```

- [ ] **Step 2: manifest.yaml の形式を定義し骨格を書く**

```yaml
# lattice_catalog 検算台帳
# 各エントリはリンタが YAML 本体と突合する期待値(全キー必須)。
# coordination: サイトラベル × ボンド type ごとの配位数
#   (min_size_for_check のトーラス上で計数展開により実測比較。
#    全ボンド type を列挙する — 省略不可)
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
#   coordination:
#     A: {J0: 2, "J0'": 2, "J0''": 2}
#   min_size_for_check: [7]
#   source: {file: python/stdface/lattice/chain_lattice.py, func: chain, commit: <hash>}
```

- [ ] **Step 3: lint_catalog.py を書く**

チェック項目(スペック §6.1 と同一。各 ID は診断メッセージの接頭辞):

| ID | 検査内容 |
|----|----------|
| C1 | `catalog.schema == "stdface-catalog/0.1"` **かつ** `catalog.dialect == "experimental"` |
| C2 | schema 検査: 文書が dict、geometry/system/model の必須キーと型、サイトラベル重複なし、dimension は正整数、R/size 要素は整数、site_dof が存在 |
| C3 | R の長さ = dimension = size の長さ |
| C4 | bonds の from/to、**onsite のサイトラベル**が定義済み。`geometry.sites` のラベル集合 = `site_dof` のキー集合 |
| C5 | bonds の type ↔ couplings キー整合(未定義参照・未使用定義) |
| C6 | 反転同値 `(i,j,R) ≡ (j,i,−R)` での重複検出(向きの並べ替え検査はしない — ソース順保持のため) |
| C7 | `{param: ...}` 参照名が目録に存在。`scale` は実数、`default` は数値 |
| C8 | J 型 coupling(tensor_terms が Sx/Sy/Sz の 2 サイト積のもの)は 9 成分を一度ずつ持ち、ops 対と param 接尾辞が対応(`[Sx,Sy]↔…xy` 等) |
| C9 | 演算子と site_dof の型整合: hop/density-density は fermion–fermion、`s_i . S_j` は fermion–spin(この順)、J テンソルは spin–spin。tensor_terms の ops 長 = ボンド 2 / onsite 1 |
| C10 | manifest 全キー(lattice/model/dimension/n_sites_uc/bonds_per_uc/coordination/min_size_for_check/source)の存在と、n_sites_uc・dimension・bonds_per_uc の一致 |
| C11 | **計数展開**: min_size_for_check のトーラス上に bonds を展開し、ラベル×type の配位数を実測して manifest の coordination と比較 |

実装要件:
- 不正 YAML・キー欠落でも例外で落とさず、ファイル単位の診断として報告
  (try/except で C2 診断に変換)。
- 引数のパスは `Path.resolve()` してから ROOT 相対に変換(ROOT 外は
  エラー報告)。引数なしで `lattice_catalog/**/*.yaml` 全件
  (manifest.yaml 除外)。違反があれば exit 1、末尾に
  `{N} files, {M} errors` を出力。
- 目録は `keyword_inventory.py` をサブプロセス実行して取得。

中核アルゴリズム(この通り実装):

```python
def reversal_dup(bonds):
    """C6: 反転同値 (i,j,R)~(j,i,-R) での重複検出。"""
    seen, errs = set(), []
    for b in bonds:
        k = (b["from"], b["to"], tuple(b["R"]))
        rk = (b["to"], b["from"], tuple(-x for x in b["R"]))
        if k in seen or rk in seen:
            errs.append(f"C6: (逆向き)重複ボンド {b}")
        seen.add(k)
    return errs

_J_COMPONENTS = {  # C8: ops 対 ↔ param 接尾辞
    ("Sx", "Sx"): "x",  ("Sy", "Sy"): "y",  ("Sz", "Sz"): "z",
    ("Sx", "Sy"): "xy", ("Sx", "Sz"): "xz", ("Sy", "Sx"): "yx",
    ("Sy", "Sz"): "yz", ("Sz", "Sx"): "zx", ("Sz", "Sy"): "zy",
}

def check_j_coupling(type_name, tensor_terms):
    """C8: 9 成分完全性と ops↔接尾辞対応。"""
    errs, seen = [], set()
    for tt in tensor_terms:
        pair = tuple(tt["ops"])
        suffix = _J_COMPONENTS.get(pair)
        if suffix is None:
            errs.append(f"C8: {type_name}: 不正な ops 対 {pair}")
            continue
        expected = f"{type_name}{suffix}"
        got = tt["coeff"].get("param") if isinstance(tt["coeff"], dict) else None
        if got != expected:
            errs.append(f"C8: {type_name}: param {got} != {expected}")
        seen.add(pair)
    if len(seen) != 9:
        errs.append(f"C8: {type_name}: 成分数 {len(seen)} != 9")
    return errs

def expand_and_count(dim, labels, bonds, size):
    """C11: min_size トーラス上でボンドを展開し、
    ラベル×type の配位数(そのラベルのサイト 1 個に接続する本数)を返す。"""
    import itertools
    ncells = 1
    for s in size:
        ncells *= s
    touch = {lb: {} for lb in labels}          # label -> type -> 接続本数合計
    for cell in itertools.product(*[range(s) for s in size]):
        for b in bonds:
            to_cell = tuple((c + r) % s for c, r, s in zip(cell, b["R"], size))
            t = b["type"]
            touch[b["from"]][t] = touch[b["from"]].get(t, 0) + 1
            touch[b["to"]][t] = touch[b["to"]].get(t, 0) + 1
            if b["from"] == b["to"] and tuple(b["R"]) == (0,) * dim:
                pass  # 自己ループは二重加算しない(実際には存在しない想定)
            _ = to_cell
    # ラベルごとのサイト数 = ncells なので配位数 = 合計 / ncells
    return {lb: {t: n // ncells for t, n in d.items()} for lb, d in touch.items()}
```

(注: `expand_and_count` は from/to 双方の接続を数える。R=0 の同一
ラベル自己ボンドは現行 StdFace に存在しないため考慮不要。)

- [ ] **Step 4: keyword_inventory.py を書く**

```python
#!/usr/bin/env python3
"""StdFace パーサーレジストリ全体からキーワード目録(JSON)を生成する。

母集合: core の _COMMON_KEYWORDS + ソルバーレジストリ経由で列挙した
全プラグインのキーワードテーブル(クラスのハードコード禁止)+
格子名/模型名の alias(lattice/model registry から)。
出力(canonical 単位に集約、安定ソート):
[{"keyword_canonical": "J0x", "keyword_lower": "j0x",
  "sources": ["common"], "kind": "keyword"},
 {"keyword_canonical": "chain", "sources": ["lattice_registry"],
  "kind": "lattice_alias", "canonical_target": "chain"}, ...]
"""
from __future__ import annotations
import json, sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

from stdface.core.keyword_parser import _COMMON_KEYWORDS  # noqa: E402
# ソルバーは plugin レジストリ経由で列挙する。実装時に
# stdface.plugin のレジストリ API(登録済みプラグイン一覧)を確認し、
# 各プラグインのキーワードテーブル属性を取得する。
# 格子 alias は stdface.lattice のレジストリ(LatticePlugin.aliases)から。

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
    agg: dict[tuple, set] = defaultdict(set)   # (canonical, kind) -> sources
    for kw in _COMMON_KEYWORDS:
        agg[(canon(kw), "keyword")].add("common")
    # ソルバーレジストリの全プラグイン: agg[(canon(kw), "keyword")].add(name)
    # 格子レジストリ: agg[(alias, "lattice_alias")].add("lattice_registry")
    # 模型名: agg[(name, "model_alias")].add("model_registry")
    # (実装時にレジストリ API に合わせて記入。既存コードは変更しない)
    entries = [{"keyword_canonical": c, "kind": k, "sources": sorted(s)}
               for (c, k), s in sorted(agg.items())]
    json.dump(entries, sys.stdout, ensure_ascii=False, indent=1)

if __name__ == "__main__":
    main()
```

- [ ] **Step 5: test_tools.py を書く**

tools の自動テスト(pytest 不要、`python3 test_tools.py` で完結する
assert ベースでよい)。最低限:
- inventory: `2S` がソルバー表から取れる / prime 付き J 成分
  (`J0'x` 等)が canonical に含まれる / 出力が安定ソート・重複なし /
  ソルバーレジストリの件数 ≧ 4(将来プラグイン追加の検知)
- リンタ C1–C11: 各チェックの正例・負例(インライン YAML 文字列で
  最小ケースを構成 — 壊れた schema、逆向き重複、9 成分欠落、
  型不整合 fermion–spin、manifest 不一致、計数展開の配位数不一致)
- 不正 YAML(パース不能・null 文書)で例外にならず診断が出ること

Run: `python3 lattice_catalog/tools/test_tools.py`
Expected: `all tools tests passed`

- [ ] **Step 6: 動作確認**

Run: `python3 lattice_catalog/tools/keyword_inventory.py | python3 -c "import json,sys; d=json.load(sys.stdin); print(len(d), 'entries'); assert any(e['keyword_canonical']=='2S' for e in d)"`
Expected: 件数表示、assert 成功

Run: `python3 lattice_catalog/tools/lint_catalog.py`
Expected: `0 files, 0 errors`(カタログ YAML なしの状態で正常終了)

- [ ] **Step 7: Commit**

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
    A: {spin: {param: 2S, scale: 0.5, default: 0.5}}   # S = 0.5 × 2S。既定 S=1/2

  bonds:
    # 参照実装 (_BONDS) のソース順を保持
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
    # H への寄与 = value × operator (物理係数)。-t は scale で表す。
    t0:    {operator: "hop", value: {param: t0, scale: -1.0}}   # -t Σσ(c†c+h.c.)
    "t0'": {operator: "hop", value: {param: "t0'", scale: -1.0}}
    "t0''": {operator: "hop", value: {param: "t0''", scale: -1.0}}
    V0:    {operator: "density-density", value: {param: V0}}    # +V n_i n_j
    "V0'": {operator: "density-density", value: {param: "V0'"}}
    "V0''": {operator: "density-density", value: {param: "V0''"}}

  onsite:
    A:
      hubbard_u:                 # +U n↑n↓
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
    A_s: {spin: {param: 2S, scale: 0.5, default: 0.5}}

  bonds:
    # Kondo 結合はソース順 = general_j(..., 遍歴, 局在) の引数順。
    # 遍歴 A_c が第 1 端点 (s_i . S_j の役割と一致)
    - {from: A_c, to: A_s, R: [0], type: J}   # Kondo 結合 (セル内)
    - {from: A_c, to: A_c, R: [1], type: t0}
    - {from: A_c, to: A_c, R: [2], type: "t0'"}
    - {from: A_c, to: A_c, R: [3], type: "t0''"}
    - {from: A_c, to: A_c, R: [1], type: V0}
    - {from: A_c, to: A_c, R: [2], type: "V0'"}
    - {from: A_c, to: A_c, R: [3], type: "V0''"}

  couplings:
    t0:    {operator: "hop", value: {param: t0, scale: -1.0}}
    "t0'": {operator: "hop", value: {param: "t0'", scale: -1.0}}
    "t0''": {operator: "hop", value: {param: "t0''", scale: -1.0}}
    V0:    {operator: "density-density", value: {param: V0}}
    "V0'": {operator: "density-density", value: {param: "V0'"}}
    "V0''": {operator: "density-density", value: {param: "V0''"}}
    J:     {operator: "s_i . S_j", value: {param: J}}   # +J s·S (第1端点=遍歴)

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

注意: chain_kondo の J は等方スカラーのみ(Kondo の `input_spin` は
スカラー入力)だが、成分キーワード(`Jx` 等)が共通テーブルに存在する
ため `{param: J}` のスカラー参照とし、その旨コメントする。

- [ ] **Step 4: manifest.yaml に chain 3 エントリを追加**

commit ハッシュは `git rev-parse HEAD` の値を使用。coordination は
ラベル×type で**全ボンド type を列挙**(省略不可):

```yaml
files:
  chain/chain_spin.yaml:
    lattice: chain
    model: spin
    dimension: 1
    n_sites_uc: 1
    bonds_per_uc: {J0: 1, "J0'": 1, "J0''": 1}
    coordination:
      A: {J0: 2, "J0'": 2, "J0''": 2}
    min_size_for_check: [7]
    source: {file: python/stdface/lattice/chain_lattice.py, func: chain, commit: <hash>}
  chain/chain_hubbard.yaml:
    lattice: chain
    model: hubbard
    dimension: 1
    n_sites_uc: 1
    bonds_per_uc: {t0: 1, "t0'": 1, "t0''": 1, V0: 1, "V0'": 1, "V0''": 1}
    coordination:
      A: {t0: 2, "t0'": 2, "t0''": 2, V0: 2, "V0'": 2, "V0''": 2}
    min_size_for_check: [7]
    source: {file: python/stdface/lattice/chain_lattice.py, func: chain, commit: <hash>}
  chain/chain_kondo.yaml:
    lattice: chain
    model: kondo
    dimension: 1
    n_sites_uc: 2
    bonds_per_uc: {J: 1, t0: 1, "t0'": 1, "t0''": 1, V0: 1, "V0'": 1, "V0''": 1}
    coordination:
      A_c: {J: 1, t0: 2, "t0'": 2, "t0''": 2, V0: 2, "V0'": 2, "V0''": 2}
      A_s: {J: 1}
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

- [ ] **Step 3: manifest エントリ追加**(3 件。coordination はラベル×type
で全 type 列挙: 最近接は J0+J1 各 2、対角 `J'` 系 4。
min_size_for_check は各方向 L > 2·max|R成分| の奇数)

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
    # Kondo 結合: 現行実装は A3_c (副格子3の遍歴サイト) と全局在スピンを
    # 結合する (他格子の副格子ごと対応と異なる。上流実装のバグの可能性
    # あり — manual 7章参照。上流修正時はカタログの versioning が必要)。
    # ソース順保持: general_j(..., isite+3, jsite+uc_i) の引数順で
    # 遍歴 A3_c が常に第 1 端点 (s_i . S_j の役割と一致)
    - {from: A3_c, to: A0_s, R: [0, 0, 0], type: J}
    - {from: A3_c, to: A1_s, R: [0, 0, 0], type: J}
    - {from: A3_c, to: A2_s, R: [0, 0, 0], type: J}
    - {from: A3_c, to: A3_s, R: [0, 0, 0], type: J}
```

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
データセットを選ぶ。**H(_hr.dat)と U(_ur.dat)チャネルのみ使用**。
J チャネル(_jr.dat: Hund・exchange・pair-hopping)は 4 フェルミオン項の
データモデルが未定義のため対象外(manual 5.5 章で規則を散文記述、
7 章に制限として明記)。

- [ ] **Step 3: example_hubbard.yaml を書く**

- `geometry`: _geom の格子ベクトルと Wannier 中心 → sites(W1, W2, …)
- `bonds` + `couplings`: 選んだ cutoff の**範囲内全要素**を Step 1 の
  規則で変換(上位 N 件抽出ではない)。正準対のみ列挙し、
  排除した Hermite 対の扱いをコメントで明記。R=0 対角は onsite へ。
- 値は実数値を直接記載(外部データ由来。param 参照不可の旨コメント)。
  H チャネルの符号反転など Step 1 で確定した規則の適用は
  ファイル冒頭コメントに規則名で明記する。
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

- **符号表と検証連鎖**(スペック §4.3 の表): StdFace パラメータ →
  builder(`interaction_builder.py`)→ trans/intr 係数 → solver 規約
  (HPhi trans.def の H=−Σt 暗黙符号)→ 物理符号、の導出を明記。
- Spin: 9 成分正準形とパラメータ解決規則(成分キーワード 9 種の表、
  解決順序: 成分局所 > 成分大域 > スカラー局所 > スカラー大域 > 0)、
  等方入力・異方入力・同時指定エラーの 3 例。
- Hubbard: hop / density-density / onsite の定義(数式)、scale による
  符号表現、複素 t と反転同値(共役)。
- Kondo: 2 ラベル表現と端点順序(遍歴第 1)、磁場両側適用、
  サイト倍加(前半=局在)との対応、GC 変種(粒子数条件は solver 層)。

- [ ] **Step 3: 5 章(wannier90 変換仕様)を書く**

Task 10 Step 1 で確定した規則を模型別に文書化:
- Hubbard: H/U チャネルの変換(符号反転、onsite 分離、Hermite 正準対、
  縮退重み、cutoff、λ/α、doublecounting)
- Spin: 超交換 `2|t_mn|^2(1/U_m+1/U_n)` 生成の現行アルゴリズムと、
  直接写像との差異
- **5.5 J チャネル**: Hund・exchange・pair-hopping の変換規則を散文で
  記述し、YAML 例が対象外である理由(4 フェルミオン項スキーマ未定義)と
  一般項スキーマの素描(6 章の拡張提案 6 項)を示す
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

スペック §4.6 の 6 項目(fermion site_dof、1 サイト演算子語彙、
名前付き 2 体演算子と端点順序の意味論、`{param, scale, default}` 参照の
許容箇所、catalog ヘッダ、4 フェルミオン一般項の素描)を、
それぞれ tensor_terms への展開形・型・意味論つきで draft への
追記提案として記述。

- [ ] **Step 2: 7 章(既知の制限・未対応事項)を書く**

- oracle 比較(数値同値性)未実施 — 展開エンジン実装時の課題
- ladder は W=2/3 の例のみ(一般 W は生成規則の文書)
- wannier90 は Hubbard の H/U チャネル例のみ。J チャネル
  (4 フェルミオン項)は散文記述のみ、Spin/Kondo も仕様記述のみ
- pyrochlore Kondo の `isite+3` 挙動(上流バグ疑い、現行動作を記載。
  上流修正時はカタログの versioning が必要)
- box 行列は写像説明のみで YAML 例なし
- 機械的識別子(prime を含む type 名)の分離は仕様確定時の課題
- 多体項(3 体以上)は現行 StdFace に存在しないため対象外
- YAML・manifest・manual の三重記載はリンタと inventory の機械照合で
  緩和しているが、ソース変更時は三箇所の同期更新が必要

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
