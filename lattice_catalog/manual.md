# StdFace 格子カタログ マニュアル

本書は `lattice_catalog/` 配下の YAML カタログの読み方と、StdFace が
定義するキーワード群が新フォーマット(geometry / system / model の
三層構造)のどこに対応するかを説明する。規範文書は
`lattice_catalog/CONVENTIONS.md` であり、本書はその解説・実例・
索引を提供する位置づけである(本書と CONVENTIONS.md が矛盾する場合は
CONVENTIONS.md を優先する)。

---

## 1. 概要と読み方

### 1.1 三層仕様の要約と本カタログの位置づけ

参照仕様「格子定義仕様 (draft)」2026/07/27(以下「draft 仕様」)は、
格子・模型の定義を次の三層に分けて記述する。

| 層 | 内容 |
|---|---|
| `geometry` | 次元、格子ベクトル、単位胞内サイトの分率座標 |
| `system` | 繰り返し数(サイズ)、境界条件、supercell 変換 |
| `model` | サイトの自由度(site_dof)、ボンド(bonds)、相互作用(couplings)、
  オンサイト項(onsite) |

本カタログの各 YAML ファイルはこの三層構造 + `catalog:` ヘッダの
4 トップレベルキーからなる自己完結型ファイルであり(CONVENTIONS.md
§2)、1 ファイル = 1 格子 × 1 模型の対応を持つ。

**dialect 宣言**: 各ファイルは

```yaml
catalog: {schema: stdface-catalog/0.1, dialect: experimental, lattice: <格子名>, model: <模型名>}
```

を先頭に持つ。`dialect: experimental` は、本カタログが draft 仕様
そのものではなく、draft が未規定の事項(fermion site_dof、演算子語彙、
`{param, scale, default}` 参照など、CONVENTIONS.md §7 に列挙)を
補う**実験的方言(experimental dialect)**であることを明示する固定値
である。draft 仕様が規定する部分(geometry / system / model の三層
構造そのもの)は draft に準拠し、draft が触れていない部分は
CONVENTIONS.md を規範とする(CONVENTIONS.md §1)。

**検証水準**: 本カタログが主張する検証水準は以下の 3 点であり、
**実行同値性検証(数値 oracle 比較)ではない**。

1. **ソース突合**: 各 YAML の `bonds` / `couplings` / `onsite` は、
   参照実装(`python/stdface/lattice/*.py` の `_BONDS` テーブル、
   `python/stdface/core/model_plugin.py` の演算子生成、
   `python/stdface/core/input_params.py` のパラメータ解決規則)を
   直接読んで書き写したものであり、各ファイル末尾に
   `source: {file, func, commit}` として出典を記録する
   (CONVENTIONS.md §8)。
2. **リンタ検査**(`tools/lint_catalog.py`、検査 C1–C11): スキーマ・
   ラベル整合・ボンド重複・param 参照整合・J テンソル成分完全性・
   演算子と site_dof の型整合、および manifest 突合を機械的に検査する。
   特に **C11(計数展開)** は `min_size_for_check` で指定したトーラス
   サイズ上に実際にボンドを展開し、サイトラベル×type ごとの配位数を
   実測して `manifest.yaml` の `coordination` と突合する
   (単なる静的スキーマ検査を超えた「計数の実測比較」)。
3. **manifest 検算**: `lattice_catalog/manifest.yaml` に各ファイルの
   期待値(`n_sites_uc`, `bonds_per_uc`, `coordination`,
   `min_size_for_check`, `source`)を機械可読な形で記録し、
   リンタが YAML 本体と突合する(C10)。

すなわち、「stan.in をこの YAML から生成したハミルトニアンで実際に
計算した結果が、StdFace が生成する `trans.def` 等の入力ファイルで
計算した結果と数値的に一致する」という強い意味での実行同値性は、
本カタログの範囲では検証していない。この限界は 7 章(既知の制限)で
改めて明記する。

### 1.2 ディレクトリ構成

```
lattice_catalog/
  CONVENTIONS.md            # 記述規約(規範文書)
  manual.md                 # 本書
  manifest.yaml             # 検算台帳(機械可読)
  tools/
    lint_catalog.py         # 意味検査リンタ(C1–C11、計数展開器を含む)
    keyword_inventory.py    # キーワード目録生成(ソルバーレジストリ走査)
    test_tools.py           # tools の自動テスト(開発時専用)
  chain/       chain_spin.yaml, chain_hubbard.yaml, chain_kondo.yaml
  ladder/      ladder_w2_{spin,hubbard,kondo}.yaml,
               ladder_w3_{spin,hubbard,kondo}.yaml
  square/      square_{spin,hubbard,kondo}.yaml
  triangular/  triangular_{spin,hubbard,kondo}.yaml
  honeycomb/   honeycomb_{spin,hubbard,kondo}.yaml
  kagome/      kagome_{spin,hubbard,kondo}.yaml
  orthorhombic/ orthorhombic_{spin,hubbard,kondo}.yaml
  fc_ortho/    fc_ortho_{spin,hubbard,kondo}.yaml
  pyrochlore/  pyrochlore_{spin,hubbard,kondo}.yaml
  wannier90/   example_hubbard.yaml
```

格子 9 種(chain, ladder(W=2/W=3 の 2 例), square, triangular,
honeycomb, kagome, orthorhombic, fc_ortho, pyrochlore)× 模型 3 種
(Spin / Hubbard / Kondo)で 30 ファイル、これに wannier90 の
Hubbard(H+U チャネル)例 1 ファイルを加えた **31 件**が
`lint_catalog.py` の検査対象である(`manifest.yaml` はデータファイル
であり検査対象ファイル数には含まない)。

### 1.3 CONVENTIONS.md の規約の読み方

以下は CONVENTIONS.md の主要規約を読み解く上での補足である。詳細な
文言・規範性は必ず CONVENTIONS.md 本文を参照すること。

#### (1) geometry / system

- `geometry.sites` は配列であり、配列順が draft 仕様のサイト番号に
  機械的に一致する。サイトラベル(`A`, `A_c` 等)はファイル内で
  一意でなければならない。
- `system.size` は StdFace の `W, L, Height` に、`system.boundary` は
  `phase0`–`phase2` に対応する。境界の位相因子は
  `exp(i · n · π · θ / 180)`(**度単位**、n = 境界を横断する回数)。
- **box → supercell(S) の写像**: StdFace の `box` キーワード群
  (`a0w, a0l, a0h, a1w, a1l, a1h, a2w, a2l, a2h` の 9 つ)は
  3×3 整数行列で、行が超格子ベクトル番号(`a0`→ 0 番目, `a1`→ 1 番目,
  `a2`→ 2 番目)、列が元の単位胞方向成分(`w`→ 0 列, `l`→ 1 列,
  `h`→ 2 列)を表す。これは `system.supercell` の行ベクトル規約
  `A_super = S · A`(`det(S) ≠ 0`)における変換行列 `S` そのものに
  対応する(`box[i][j] = S[i][j]`)。`supercell` と `size`
  (`W/L/Height`)は StdFace 上も排他的に用いられ、本カタログの YAML
  では常に `size` 表現を基本形として採用する(`supercell` 表現は
  YAML では使用していない — CONVENTIONS.md §4)。

#### (2) bonds — ソース順保持とボンド反転同値

- セル差分 `R = cell(to) − cell(from)`、変位
  `δ = (frac_to − frac_from) + R · A` と定義する。
- ボンドの `from → to` の向きは**参照実装のソース順をそのまま保持**
  する(正準形への並べ替えは行わない)。通常ボンドは `_BONDS` テーブル
  の並び、Kondo 結合は `general_j(..., isite, jsite)` の引数順
  (**遍歴サイトが第 1 引数**)。並べ替えを禁止する理由は、並べ替えが
  複素ホッピングの複素共役や交換テンソルの転置を要求し、単一の
  `couplings` キーでは表現できないためである。
- リンタは**反転同値** `(type, i, j, R) ≡ (type, j, i, −R)` での
  重複を検出する(C6)。`type` を同値キーに含めるため、同一の幾何学的
  ボンド上に異なる `type`(例: `t0` と `V0`)が共存するのは正当な
  記述として扱われる。
- カタログを消費する側(resolver/展開エンジン)がボンドを
  `(j, i, −R)` として読み替える場合は、hopping は複素共役
  (`t_ij = conj(t_ji)`)、交換テンソルは転置(`J_ij = J_ji^T`)を
  取る必要がある。**カタログ自体はソース順で一意に記述するため、この
  変換はカタログ制作側の作業ではなく、消費側が反転読み替えを行う場合
  にのみ必要になる規範**である。

#### (3) couplings / onsite — value 意味論と符号表

`model.couplings[<type>].value` および
`model.onsite[<site>][<term>].value` は
`H = Σ_bonds value·operator + Σ_onsite value·operator` という
**物理ハミルトニアンの係数**である。これは HPhi の `trans.def` 等の
ソルバー出力ファイルの係数規約(`H_trans = −Σ t c†c` のように符号が
暗黙に反転している)とは別物であり、両者の対応は
「StdFace パラメータ → builder 呼び出し → trans/intr 係数 →
solver 出力規約 → 物理ハミルトニアンの符号」という検証連鎖で突合する
(導出は 4 章)。

符号は必ずデータ(`scale` / `coeff`)に持たせ、コメントに書いては
ならない(機械的検証可能性のための必須規約)。代表例:

| StdFace | 物理ハミルトニアン寄与 | YAML 表現の要点 |
|---|---|---|
| t 族 | −t Σ_σ (c†c + h.c.) | `value: {param: t0, scale: -1.0}` |
| U | +U n↑n↓ | onsite `coeff: +1.0` |
| V 族 | +V n_i n_j | `value: {param: V0}`(scale 省略時 +1.0) |
| h / Gamma / Gamma_y | −h Sz − Γ Sx − Γy Sy | onsite `coeff: -1.0` |
| Kondo J | +J s·S | `value: {param: J}`(scale 省略時 +1.0) |

(全項目は CONVENTIONS.md §6.2 参照)

#### (4) `{param, scale, default}` の条件分岐

拡張方言の中核をなす一般形は

```yaml
value: {param: <名前>, scale: <実数, 省略時 1.0>, default: <数値, 省略時 0>}
```

であり、**param が指定された場合と指定されない場合で計算式が変わる**
(この分岐そのものが規範である点に注意):

- **param 指定時**: 値 = `scale × param`。
- **param 未指定時**: 値 = `default`(**`scale` は適用しない最終値**)。

例(`chain_spin.yaml` の site_dof): `{param: 2S, scale: 0.5, default: 0.5}`
は、`2S=1` が指定されれば S = 0.5 × 1 = 0.5、`2S` が指定されなければ
S = 0.5(`default` をそのまま採用し `scale` は掛けない)という意味に
なる。同じ一般形が `value`(couplings/onsite)、`coeff`(J 成分)、
`spin`(site_dof の `2S`)、`twist`(境界の `phaseN`)のいずれにも
共通して用いられる。

#### (5) J テンソル 9 成分形とパラメータ解決順序

J 族の交換相互作用(`model.couplings[<type>].operator.tensor_terms`)は
常に **9 成分**を書き切る正準形を取る(等方成分への縮約は行わない)。
ops 対と param 接尾辞の対応:

| ops 対 | 接尾辞 | ops 対 | 接尾辞 |
|---|---|---|---|
| `[Sx,Sx]` | `x` | `[Sy,Sz]` | `yz` |
| `[Sy,Sy]` | `y` | `[Sz,Sx]` | `zx` |
| `[Sz,Sz]` | `z` | `[Sz,Sy]` | `zy` |
| `[Sx,Sy]` | `xy` | `[Sx,Sz]` | `xz` |
| `[Sy,Sx]` | `yx` | | |

各成分の値は `input_params.py::_resolve_spin_matrix` と同一の優先順位
で解決される(カタログはこの解決規則を規範として提供し、実際の解決は
**カタログを消費する側**が行う):

1. 成分局所(例 `J0xy`) 2. 成分大域(例 `Jxy`)
3. スカラー局所・対角のみ(例 `J0`) 4. スカラー大域・対角のみ(例 `J`)
5. 0(既定値)

プライム系(`J0'`, `J0''`, `J1'` 等)には大域 fallback が存在しない
点が例外である(詳細は CONVENTIONS.md §6.5)。

### 1.4 YAML の読み方 — `chain_spin.yaml` / `chain_kondo.yaml` を例に

**`chain/chain_spin.yaml`**(Spin 模型、S=1/2 既定の 1 次元鎖)

- `catalog:` ヘッダの `lattice: chain, model: spin` が
  `manifest.yaml` の `chain/chain_spin.yaml` エントリと突合される
  (C10)。
- `geometry.sites` はサイト 1 つ(`A`)、`lattice_vectors.a1` は
  StdFace の `a`(既定 1.0)に対応。
- `system.size: [16]` は StdFace の `L=16`、`boundary` の
  `{twist: {param: phase0}}` は境界位相を表す(`phase0=0` で
  periodic と等価)。
- `model.site_dof.A` は `{spin: {param: 2S, scale: 0.5, default: 0.5}}`
  — 1.3 節(4)の一般形の実例そのもの。
- `model.bonds` は最近接(`J0`, R=[1])・次近接(`J0'`, R=[2])・
  三次近接(`J0''`, R=[3])の 3 本を `_BONDS` のソース順のまま列挙。
- `model.couplings` の `J0` / `"J0'"` / `"J0''"` はそれぞれ 9 成分
  `tensor_terms` を持ち、`coeff: {param: J0x}` のように成分ごとに
  独立した param 参照を持つ(1.3 節(5)の実例)。プライムを含む
  `type` 名(`"J0'"`)は YAML 上引用符付きで書く。
- `model.onsite.A` には `field_z`(`-h Sz`)、`field_x`
  (`-Gamma Sx`)、`field_y`(`-Gamma_y Sy`)、`aniso_z`(`+D (Sz)^2`)
  の 4 項があり、各項は `operator.tensor_terms` に**符号のみを表す
  リテラル `coeff`**(`-1.0` や `+1.0`)と、外部パラメータを担う
  `value: {param: ...}` の 2 つを持つ(CONVENTIONS.md §6.1 の
  「J 族 / onsite の例外」に対応する実例)。

**`chain/chain_kondo.yaml`**(Kondo 格子模型、2 ラベル規約の実例)

- `geometry.sites` は `A_c`(遍歴電子)と `A_s`(局在スピン)の
  2 ラベルを持つが、両者は**同一の分率座標**(`frac: [0.0]`)を持つ
  ——「1 物理サイト = 2 geometry ラベル」という Kondo の 2 ラベル規約
  (CONVENTIONS.md §6.7)の直接的な表れである。
- `model.site_dof` は `A_c: {fermion: {orbitals: 1}}`(fermion
  site_dof、拡張方言 §7-1)、`A_s: {spin: {param: 2S, ...}}` と、
  ラベルごとに異なる自由度種別を持つ。
- `model.bonds` の先頭行 `{from: A_c, to: A_s, R: [0], type: J}` は
  Kondo 結合であり、**遍歴電子 `A_c` が第 1 引数**という
  `general_j(..., isite, jsite)` の引数順をそのまま反映している
  (ソース順保持規約、1.3 節(2))。続く `t0/t0'/t0''` と
  `V0/V0'/V0''` は `A_c → A_c` のホッピング・サイト間 Coulomb で、
  Hubbard 模型と同じ構造を遍歴自由度にのみ適用する。
- `model.couplings.J` は `operator: "s_i . S_j"` であり、
  端点順序に意味を持つ(第 1 端点 = 遍歴電子スピン、第 2 端点 =
  局在スピン)。
- `model.onsite` は `A_c`(Hubbard 一式 `U, mu` + 磁場 3 成分)と
  `A_s`(磁場 3 成分のみ)に分かれる。磁場(`h/Gamma/Gamma_y`)が
  **両ラベルに適用される**点(CONVENTIONS.md §6.7)がファイル上
  そのまま確認できる。

---

## 2. キーワード対応表

### 2.1 分類方針

母集合は `python3 lattice_catalog/tools/keyword_inventory.py` の出力
(`kind: "keyword"` のエントリ **313 件**)である。このツールは
`stdface.core.keyword_parser._COMMON_KEYWORDS` と、
`stdface.plugin` 経由で列挙した全ソルバープラグイン
(HPhi / HWAVE / UHF / mVMC)のキーワードテーブルを走査して
canonical キーワード単位に集約したものであり、具象クラス名を
ハードコードしていないため、将来プラグインが増えても自動的に反映
される(ツール自体の設計)。

各エントリを以下 7 分類のいずれか 1 つに割り当てる。1 つの分類に
含まれる関連キーワード群(`wx/wy/wz` のような成分 3 つ組、
J 族の 9 成分×全 prefix など)は**族単位で 1 行にまとめる**ことを
許容し、その場合は「StdFace キーワード」列に構成キーワードを列挙し、
「備考」列に命名規則を明記する。

- (a) **geometry**: 格子ベクトルの長さ・方向
- (b) **system**: サイズ・境界・supercell
- (c) **bonds + couplings**: ホッピング/交換/Coulomb 相互作用
- (d) **onsite**: オンサイト項
- (e) **site_dof**: サイト自由度
- (f) **wannier90**: wannier90 変換専用パラメータ(詳細は 5 章)
- (g) **対象外**: 模型定義(三層構造)に現れない計算条件・メタ情報・
  ソルバー固有パラメータ(理由を行ごとに付す)

### 2.2 (a) geometry

| StdFace キーワード | 出典テーブル | 新フォーマットでの対応 | 備考 |
|---|---|---|---|
| `a` | common | `geometry.lattice_vectors` のスケール(既定 1.0) | 格子定数。`chain_spin.yaml` の `a1: [1.0]` 相当 |
| `wx, wy, wz` | common | `geometry.lattice_vectors` の 1 本目のベクトル成分 | `direct[0][0..2]`。W 方向の実空間方向ベクトル |
| `lx, ly, lz` | common | `geometry.lattice_vectors` の 2 本目のベクトル成分 | `direct[1][0..2]`。L 方向 |
| `hx, hy, hz` | common | `geometry.lattice_vectors` の 3 本目のベクトル成分 | `direct[2][0..2]`。Height 方向 |
| `wlength, llength, hlength` | common | `geometry.lattice_vectors` の各方向ベクトルの長さスケール | `length[0..2]`。方向ベクトル `wx/wy/wz` 等と併用しスケールを決める |

### 2.3 (b) system

| StdFace キーワード | 出典テーブル | 新フォーマットでの対応 | 備考 |
|---|---|---|---|
| `W, L, Height` | common | `system.size`(整数配列) | stan.in の代表値。manifest の `min_size_for_check`(最大 \|R\| 成分から独立に導出されるリンタ検査サイズ。各方向 `2*max|R|` を超える最小の奇数)とは独立 |
| `phase0, phase1, phase2` | common | `system.boundary[*].twist`(`{param: phaseN}`) | 度単位。境界を n 回横断する経路の位相因子は `exp(i·n·π·θ/180)` |
| `a0w, a0l, a0h, a1w, a1l, a1h, a2w, a2l, a2h`(box 行列 9 成分) | common | `system.supercell`(`A_super = S·A`、`det(S)≠0`) | 1.3 節(1)参照。`box[i][j] = S[i][j]`(行=超格子ベクトル番号 a0/a1/a2、列=元の方向 w/l/h)。`size` と排他利用。本カタログの YAML は `size` 表現のみ使用 |

### 2.4 (c) bonds + couplings

| StdFace キーワード | 出典テーブル | 新フォーマットでの対応 | 備考 |
|---|---|---|---|
| `t, t', t'', t0, t0', t0'', t1, t1', t1'', t2, t2', t2''`(12 種) | common | `model.bonds[*].type` + `model.couplings[<type>]: {operator: "hop", value: {param: ..., scale: -1.0}}` | ホッピング族。`t0/t1/t2` はボンド type ごとの局所値、無 prefix の `t/t'/t''` は `input_hopp()`(`input_params.py`)による**大域フォールバック**(`t0/t1/t2` 等が未指定のとき `t` を採用。局所・大域の同時指定はエラー)。`'`/`''` は同じ枝分かれを持つ次近接以降の系列を表す接尾辞 |
| `V, V', V'', V0, V0', V0'', V1, V1', V1'', V2, V2', V2''`(12 種) | common | `model.couplings[<type>]: {operator: "density-density", value: {param: ...}}`(scale 省略時 +1.0) | サイト間 Coulomb 族。命名規則・大域フォールバック(`input_coulomb_v()`)の構造は t 族と同一 |
| `J, J', J'', J0, J0', J0'', J1, J1', J1'', J2, J2', J2''`(12 prefix、各「等方スカラー 1 + 異方成分 9」で 120 キーワード) | common | `model.couplings[<type>].operator.tensor_terms`(9 成分正準形。等方版は 1.3 節(5) の解決順序 3–4 でスカラー扱い) | 接尾辞 `x,y,z,xy,xz,yx,yz,zx,zy` は ops 対 `[S?,S?]` に 1 対 1 対応(1.3 節(5)の表、CONVENTIONS.md §6.4)。無 prefix の `J`(スカラー)は Spin 模型の大域等方交換と Kondo 模型の `s_i . S_j` 結合の両方で使われる(模型により演算子が異なる。CONVENTIONS.md §6.2/§6.6) |

族単位の内訳(t 12 + V 12 + J 120 = **144 キーワード**)。

### 2.5 (d) onsite

| StdFace キーワード | 出典テーブル | 新フォーマットでの対応 | 備考 |
|---|---|---|---|
| `U` | common | onsite `{tensor_terms: [{ops: [NupNdn], coeff: +1.0}], value: {param: U}}` | オンサイト Coulomb。Hubbard/Kondo(`_c` ラベル)に適用 |
| `mu` | common | onsite `{tensor_terms: [{ops: [N], coeff: -1.0}], value: {param: mu}}` | 化学ポテンシャル。Hubbard/Kondo(`_c` ラベル)に適用 |
| `D` | common | onsite `{tensor_terms: [{ops: [Szz], coeff: +1.0}], value: {param: D}}` | 単イオン異方性。Spin 模型のみ |
| `h` | common | onsite `{tensor_terms: [{ops: [Sz], coeff: -1.0}], value: {param: h}}` | 磁場 z 成分。全模型に適用(Kondo は両ラベル) |
| `Gamma` | common | onsite `{tensor_terms: [{ops: [Sx], coeff: -1.0}], value: {param: Gamma}}` | 磁場 x 成分。適用範囲は `h` と同じ |
| `Gamma_y` | common | onsite `{tensor_terms: [{ops: [Sy], coeff: -1.0}], value: {param: Gamma_y}}` | 磁場 y 成分。適用範囲は `h` と同じ |

### 2.6 (e) site_dof

| StdFace キーワード | 出典テーブル | 新フォーマットでの対応 | 備考 |
|---|---|---|---|
| `2S` | HPhi | `model.site_dof[<label>]: {spin: {param: 2S, scale: 0.5, default: 0.5}}` | S = 0.5 × 2S。未指定時 S=0.5(1.3 節(4)の条件分岐の実例)。キーワード自体は現行実装では HPhi プラグインのテーブルにのみ登録されているが、Spin 自由度を持つ全模型(Spin/Kondo の `_s` ラベル)に概念上適用される |

### 2.7 (f) wannier90

wannier90 変換専用パラメータは 5 章(wannier90 変換仕様)で詳説する。
ここでは inventory 上の分類のみを示す。

| StdFace キーワード | 出典テーブル | 新フォーマットでの対応 | 備考 |
|---|---|---|---|
| `cutoff_j`, `cutoff_jw/jl/jh`, `cutoff_j_a0w/a0l/a0h/a1w/a1l/a1h/a2w/a2l/a2h`(13 種) | common | wannier90 変換の J チャネル cutoff 指定 → 5 章参照 | J(交換)ホッピング行列要素の打ち切り半径・成分別打ち切り。J チャネルの YAML 例自体は対象外(1.3 節・design §2 参照) |
| `cutoff_t`, `cutoff_tw/tl/th`, `cutoff_t_a0w/a0l/a0h/a1w/a1l/a1h/a2w/a2l/a2h`(13 種) | common | wannier90 変換の hopping(t)チャネル cutoff 指定 → 5 章参照 | `example_hubbard.yaml` の H(hopping)チャネルに対応する打ち切りパラメータ |
| `cutoff_u`, `cutoff_uw/ul/uh`, `cutoff_u_a0w/a0l/a0h/a1w/a1l/a1h/a2w/a2l/a2h`(13 種) | common | wannier90 変換の U(Coulomb)チャネル cutoff 指定 → 5 章参照 | `example_hubbard.yaml` の U チャネルに対応する打ち切りパラメータ |
| `cutoff_length_j, cutoff_length_t, cutoff_length_u`(3 種) | common | wannier90 変換のチャネル別・実距離での打ち切り指定 → 5 章参照 | `cutoff_j/t/u` 系(セル添字ベース)とは別軸の、実空間距離ベースの打ち切り |
| `lambda, lambda_u, lambda_j`(3 種) | common | wannier90 `*_hr.dat` 読み込み時に全行列要素へ一様に乗じるスケール係数 → 5 章参照 | `_hr.dat` 形式自身が持つ縮退重みは読み込み時に破棄される(`wannier90_io.py::_skip_degeneracy_weights`)。`lambda*` はそれとは別に、全行列要素に一様適用されるスケール係数(`wannier90_io.py` の読み込みループで `lam * (dtmp_re + 1j*dtmp_im)` として適用)。`lambda` 未指定時は `lambda_U`/`lambda_J` の共有既定値になる |
| `alpha` | common | wannier90 Hubbard/Hund 二重計数補正の混合重み → 5 章参照 | Hartree-Fock 二重計数補正(`doublecounting` モードで有効化)における混合重み(`wannier90.py` の `_apply_coulomb_terms` 系、範囲 [0,1]・既定 0.5)。Spin チャネルの超交換 `2\|t\|²(1/U_m+1/U_n)` 生成には関与しない(そちらは `StdI.alpha` を参照しない) |
| `doublecounting` | common | wannier90 二重計数補正モード指定 → 5 章参照 | Coulomb 相互作用の二重計数補正方式の選択 |

族単位の内訳(cutoff_j 13 + cutoff_t 13 + cutoff_u 13 + cutoff_length_* 3
+ lambda* 3 + alpha 1 + doublecounting 1 = **47 キーワード**)。

### 2.8 (g) 対象外

三層構造(geometry/system/model)そのものには現れない、計算条件・
メタ情報・ソルバー固有の実行時パラメータ。行ごとに対象外とする理由を
付す。

| StdFace キーワード | 出典テーブル | 新フォーマットでの対応 | 備考(対象外とする理由) |
|---|---|---|---|
| `model, lattice, outputmode`(3 種) | common | なし | 格子/模型の選択とソルバー出力形式の指定というメタ情報。格子/模型は `catalog.lattice` / `catalog.model` ヘッダに既に反映済みであり、`outputmode` は生成器の出力形式選択でありハミルトニアン定義ではない |
| `ncond, nelec, 2Sz`(3 種) | common | なし | 全電子数・全 Sz(磁化セクター)の指定という計算条件(ソルバー入力)であり、模型のハミルトニアン定義そのものではない。`nelec` は `ncond` のエイリアス |
| `K`(1 種) | common | なし | StdFace が予約しているが現行実装の全格子関数で `not_used_d("K", ...)` により明示的に未使用と検査される(例: `chain_lattice.py`, `square_lattice.py` 等)キーワード |
| `Vecpoth, Vecpotl, Vecpotw, calcspec, cdatafilehead, dt, eigenvecio, exct, expandcoef, expecinterval, flgtemp, freq, hamio, initial_iv, initialvectype, lanczos_max, lanczoseps, lanczostarget, largevalue, method, ngpu, nomega, numave, nvec, omegaim, omegamax, omegamin, omegaorg, outputexcitedvec, pumptype, restart, scalapack, spectrumqh, spectrumql, spectrumqw, spectrumtype, tdump, tshift, uquench`(39 種) | HPhi | なし | HPhi 固有の対角化手法(Lanczos 等)・時間発展(pump-probe)・動的相関関数(スペクトル)などの計算条件パラメータ。ハミルトニアンの定義ではなくソルバーの実行制御 |
| `calcmode, exportall, fileprefix, lattice_gp`(4 種) | HWAVE | なし | HWAVE(グリーン関数法)固有の計算モード・出力ファイル設定 |
| `eps, epsslater, iteration_max, mix`(4 種) | HWAVE,UHF | なし | UHF/HWAVE 平均場計算の収束判定閾値・反復回数・混合パラメータ(計算条件) |
| `a0hsub, a0lsub, a0wsub, a1hsub, a1lsub, a1wsub, a2hsub, a2lsub, a2wsub, hsub, lsub, wsub, nmptrans, rndseed`(14 種) | HWAVE,UHF,mVMC | なし | UHF/HWAVE/mVMC が用いる「サブ格子」supercell(`boxsub`)とサイズ(`Hsub/Lsub/Wsub`)、多点並進対称数(`nmptrans`)、乱数種(`rndseed`)。`boxsub` は 2.3 節(b)の `box`(`system.supercell`)とは別物で、平均場/変分モンテカルロの初期対称性(サブ格子秩序パターン)を指定する計算条件であり、物理的な格子構造(geometry/system)の一部ではない |
| `complextype, cparafilehead, dsroptredcut, dsroptstadel, dsroptstepdt, ndataidxstart, ndataqtysmp, nlanczosmode, nspgaussleg, nsplitsize, nspstot, nsrcg, nsroptitrsmp, nsroptitrstep, nstore, nvmccalmode, nvmcinterval, nvmcsample, nvmcwarmup`(19 種) | mVMC | なし | mVMC(変分モンテカルロ)固有の最適化(SR 法)・サンプリング・並列化に関する計算条件パラメータ |

族単位の内訳(3 + 3 + 1 + 39 + 4 + 4 + 14 + 19 = **87 キーワード**)。

### 2.9 格子名別名表(lattice_alias)

`kind: "lattice_alias"` のエントリ(**26 件**)。StdFace の `lattice`
キーワードに指定する別名と、`stdface.lattice` レジストリ上の正規名
(`canonical_target`)の対応。正規名は `lattice_catalog/` の
ディレクトリ名と**必ずしも一致しない**点に注意(`square` →
registry 正規名 `tetragonal`、`fc_ortho` ディレクトリ →
registry 正規名 `fco`。カタログの `catalog.lattice` / ディレクトリ名
はいずれも `manifest.yaml` 側の命名 `fc_ortho` / `square` を用いる
独自の対応であり、StdFace 入力の別名解決とは別レイヤーである)。

| StdFace 別名 | registry 正規名 | 対応するカタログディレクトリ |
|---|---|---|
| `chain`, `chainlattice` | `chain` | `chain/` |
| `square`, `squarelattice`, `tetragonal`, `tetragonallattice` | `tetragonal` | `square/` |
| `ladder`, `ladderlattice` | `ladder` | `ladder/` |
| `triangular`, `triangularlattice` | `triangular` | `triangular/` |
| `honeycomb`, `honeycomblattice` | `honeycomb` | `honeycomb/` |
| `kagome`, `kagomelattice` | `kagome` | `kagome/` |
| `orthorhombic`, `simpleorthorhombic`, `cubic`, `simplecubic` | `orthorhombic` | `orthorhombic/` |
| `fco`, `face-centeredorthorhombic`, `fcorthorhombic`, `face-centeredcubic`, `fccubic`, `fcc` | `fco` | `fc_ortho/` |
| `pyrochlore` | `pyrochlore` | `pyrochlore/` |
| `wannier90` | `wannier90` | `wannier90/` |

### 2.10 模型名表(model_alias)

`kind: "model_alias"` のエントリ(**3 件**)。

| StdFace の `model` 値 | 対応するカタログのファイル名末尾 |
|---|---|
| `spin` | `*_spin.yaml` |
| `hubbard` | `*_hubbard.yaml` |
| `kondo` | `*_kondo.yaml` |

### 2.11 網羅性チェック

`python3 lattice_catalog/tools/keyword_inventory.py` の出力は
合計 **342 件**であり、内訳は `kind: "keyword"` が **313 件**、
`kind: "lattice_alias"` が **26 件**、`kind: "model_alias"` が
**3 件**である(`313 + 26 + 3 = 342`)。

`kind: "keyword"` の 313 件は、2.2–2.8 節の (a)–(g) 分類にすべて
1 回ずつ含まれる。各分類の族単位の内訳を合算すると:

```
(a) geometry           13
(b) system              15
(c) bonds + couplings  144
(d) onsite               6
(e) site_dof             1
(f) wannier90            47
(g) 対象外               87
-----------------------------
合計                    313
```

`13 + 15 + 144 + 6 + 1 + 47 + 87 = 313` となり、`keyword_inventory.py`
が報告する keyword 件数(313)と一致する。`lattice_alias`(26 件、
2.9 節)・`model_alias`(3 件、2.10 節)は別表で全件を網羅しており、
`kind` 別の 3 集合(keyword / lattice_alias / model_alias)を合わせて
inventory の全 342 件が本章の表でちょうど 1 回ずつ分類されている
(欠落・重複なし)。

---

## 3. 格子ごとの解説

本章は 9 格子(chain, ladder(W=2/W=3 の 2 例を 1 節にまとめる), square,
triangular, honeycomb, kagome, orthorhombic, fc_ortho, pyrochlore)の
各 1 節からなる。各節は以下の 4 項目で構成する。

- **幾何**: `geometry.lattice_vectors`(格子ベクトル)と `geometry.sites`
  (副格子表: ラベル・分率座標)。
- **ボンド定義表**: `model.bonds` の `type` / `from`–`to` / `R` / 意味。
  表は代表として Spin 模型のボンド(J 系列)を掲載する。Hubbard/Kondo は
  同じ `(from, to, R)` に対し `type` が `t`/`V` 系列に置き換わるだけの
  並行構造を持つため(1.4 節で確認済みの chain の例と同型)、節ごとの
  差分がある場合のみ本文で個別に注記する。Kondo はこれに加え、遍歴サイト
  (`_c`)を第 1 端点とする Kondo 結合(`type: J`)を持つ(1.3 節(2))。
  ボンドは参照実装の `_BONDS` テーブルのソース順をそのまま列挙し、
  並べ替えは行わない(CONVENTIONS.md §5)。
- **検算**: `manifest.yaml` の `bonds_per_uc` / `coordination` /
  `min_size_for_check` の値がどう導出されるかを、ボンド表からの手計算・
  ソース調査・またはリンタ C11(計数展開)の実行結果で示す。
- **出典**: 参照した Python 実装ファイル・関数・コミットハッシュ
  (`manifest.yaml` の `source` と同一)。3D 格子(orthorhombic 以降)は
  C 実装(`src/*.c`)との突合結果も付す。

### 3.1 chain

**幾何**: `dimension: 1`。`a1: [1.0]`(StdFace の `a`、既定 1.0)。
Spin/Hubbard は単一サイト `A`(`frac: [0.0]`)。Kondo は同一分率座標
`[0.0]` を持つ 2 ラベル `A_c`(遍歴電子)・`A_s`(局在スピン)
(1.3 節(2)の Kondo 2 ラベル規約)。

**内部 2D 表現(`phase0` → 内部 `phase[1]` 転写)**: chain は論理的には
1 次元格子だが、参照実装 `chain_lattice.py::chain` は内部的に `W=1` の
2 次元表現(`StdI.direct` は 2×2 行列、`W` 方向と `L` 方向を持つ)を
用いて構築されている。具体的には(`chain_lattice.py` 92–95 行):

```python
StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)  # ユーザ入力
not_used_d("phase1", StdI.phase[1])
StdI.phase[1] = StdI.phase[0]   # 内部 L 方向(周期)へ転写
StdI.phase[0] = 0.0             # 内部 W 方向(強制 W=1)は常に 0
```

すなわちユーザが指定する `phase0` は内部的には `L` 方向(周期境界を持つ
唯一の方向)の位相 `phase[1]` として扱われ、`W` 方向(強制的に 1 に
固定される非周期方向)の内部 `phase[0]` は常に 0 に固定される。本カタログ
の YAML は `dimension: 1` の論理表現のみを公開し、`system.boundary` には
`phase0` の 1 成分のみを記載する — この内部転写は実装上の詳細であり、
YAML の意味論には影響しない(CONVENTIONS.md §4)。同型の転写は
3.2 節の ladder にも現れる。

**ボンド定義表**(Spin, `_BONDS` 3 行そのまま):

| type | from→to | R | 意味 |
|---|---|---|---|
| `J0` | A→A | `[1]` | 最近接(別名 `J`) |
| `J0'` | A→A | `[2]` | 次近接(別名 `J'`) |
| `J0''` | A→A | `[3]` | 三次近接(別名 `J''`) |

Hubbard/Kondo は `t0/t0'/t0''` + `V0/V0'/V0''` が同じ `R` で並行に存在する。
Kondo はさらに先頭に `{from: A_c, to: A_s, R: [0], type: J}` を持つ
(`general_j(StdI.J, 1, StdI.S2, isite, jsite_kondo, ...)` の引数順 =
遍歴 `A_c` が第 1 端点)。

**検算**: `chain_spin.yaml` の `bonds_per_uc: {J0:1, J0':1, J0'':1}`
(各 1 本/単位胞)、`coordination.A: {J0:2, J0':2, J0'':2}`(±R 両方向で
2)。`min_size_for_check: [7]`(全 `R` のうち最大成分は `J0''` の 3、
`2×3=6` を超える最小の奇数 7)。`chain_kondo.yaml` は `n_sites_uc: 2`、
`coordination.A_c.J: 1` / `coordination.A_s.J: 1`(`A_c≠A_s` の単一方向
ボンドなので二重化されない)。

**出典**: `python/stdface/lattice/chain_lattice.py::chain`(`_BONDS`)、
commit `b96aef2107f1ab200565efd496606a01184e3f46`。

### 3.2 ladder(W=2 / W=3、W 一般化規則)

**幾何**: `dimension: 1`。単位胞 = 1 ラング(rung、W サイト)。`a1: [1.0]`
は脚(leg, L 方向)のみに対応し、ラング(旧 W)方向は非周期のため
`lattice_vectors` に現れない。W=2 はラベル `A0, A1`(脚 0/1、いずれも
frac `[0.0]`)、W=3 は `A0, A1, A2`(脚 0(端)/1(中央)/2(端))。

**横方向座標が落ちる制限**: chain と同様に `system.W` は参照実装内で
強制的に 1 に固定され(`StdI.NsiteUC = StdI.W; StdI.W = 1`、
`ladder.py` 80–81 行)、カタログは `dimension: 1` の論理表現のみを持つ。
そのため脚どうしの空間的な相対位置(旧 W 方向の座標)は `frac` には
一切現れず(全ラベルが同一の `frac: [0.0]` を持つ)、脚は幾何座標では
なく **ラベル名(`A0`/`A1`/…)のみによって区別される**。これは本カタログの
1 次元表現が意図的に持つ制限であり、幾何座標だけを読んでラング方向の
構造を復元することはできない(ボンドの `from`/`to` ラベルの組合せが
構造情報を担う)。

**rung 方向が非周期(open)である根拠**(`ladder.py` 176–188 行):

1. `_BONDS` の各行の `dW` 成分は常に 0 — 脚をまたぐセル方向オフセットは
   一度も生成されない。
2. rung・diagonal ボンドは `if uc_i < NsiteUC - 1:` の条件下でのみ
   生成され、最後の脚(`uc_i = W-1`)から最初の脚(`uc_i = 0`)へ
   "閉じる" 結合は生成されない。

**位相転写**: chain と同型(`ladder.py` 75–78 行)。ユーザ入力 `phase0`
は内部 `phase[1]`(周期方向 = leg 方向)に転写され、内部 `phase[0]`
(強制 W=1 の rung 方向)は常に 0。YAML は `system.boundary` に `phase0`
1 成分のみを持つ。

**W 一般化規則**(`ladder.py` 176–188 行の動的ループを手展開したもの):

```
for uc_i in range(W):
    leg  (脚, 同一脚内):        A{uc_i}-A{uc_i}   R=[1] J1,  R=[2] J1'
    if uc_i < W-1:
        rung (ラング, セル内):    A{uc_i}-A{uc_i+1} R=[0] J0
        diag (斜め, L+1):         A{uc_i}-A{uc_i+1} R=[1] J2
        diag (斜め, L-1):         A{uc_i}-A{uc_i+1} R=[-1] J2'
```

端の脚(`uc_i=0` または `uc_i=W-1`)は rung/diag ボンドを片側しか持たない
ため、`W=2` は両脚とも端脚で対称、`W=3` は中央脚(`A1`)のみが
rung/diag を両側(`A0` 側・`A2` 側)に持ち、端脚(`A0`/`A2`)より配位数が
2 倍になる。任意の `W` へは、この `for uc_i in range(W)` ループを
そのまま伸ばせばよい。

**ボンド定義表**(W=2, Spin):

| type | from→to | R | 意味 |
|---|---|---|---|
| `J1` | A0→A0 | `[1]` | leg 最近接(`uc_i=0`) |
| `J1'` | A0→A0 | `[2]` | leg 次近接(`uc_i=0`) |
| `J0` | A0→A1 | `[0]` | rung, セル内(`uc_i=0<W-1`) |
| `J2` | A0→A1 | `[1]` | diag L+1(`uc_i=0<W-1`) |
| `J2'` | A0→A1 | `[-1]` | diag L-1(`uc_i=0<W-1`) |
| `J1` | A1→A1 | `[1]` | leg 最近接(`uc_i=1`、端) |
| `J1'` | A1→A1 | `[2]` | leg 次近接(`uc_i=1`、端) |

W=3 は上記に `A1-A2` の rung/diag ブロック(`J0/J2/J2'`、`uc_i=1<W-1=2`)
が追加され、`A2-A2` の leg ブロック(`J1/J1'`)で終わる(合計 12 行)。

**検算**: W=2: `n_sites_uc: 2`、`bonds_per_uc: {J1:2, J1':2, J0:1, J2:1,
J2':1}`、`coordination.A0`/`coordination.A1` はいずれも
`{J1:2, J1':2, J0:1, J2:1, J2':1}`(両脚とも端で対称)。W=3:
`n_sites_uc: 3`、`bonds_per_uc: {J1:3, J1':3, J0:2, J2:2, J2':2}`、
中央脚 `coordination.A1: {J0:2, J2:2, J2':2}` は端脚
`coordination.A0/A2: {J0:1, J2:1, J2':1}` の 2 倍。leg 系列(`J1`/`J1'`)は
全脚とも配位数 2(並進で折り返す)。`min_size_for_check: [5]`(最大
`|R|` 成分は `J1'` の 2、`2×2=4` を超える最小の奇数 5)。W=2/W=3 いずれも
`ladder()` 参照実装を実際に実行し `exchange.def` 出力と手展開結果が
一致することを確認済み(Task 6)。

**出典**: `python/stdface/lattice/ladder.py::ladder`(`_BONDS`,
176–188 行)、commit `e733ff893a52b81addd27ef7e347f77273c0d129`。

### 3.3 square

**幾何**: `dimension: 2`。`a1: [1.0, 0.0]`, `a2: [0.0, 1.0]`(単位行列、
StdFace の `Wx/Wy/Lx/Ly` 既定値)。単一サイト `A`(`frac: [0.0, 0.0]`)。

**ボンド定義表**(Spin, `_BONDS` 6 行):

| type | from→to | R | 意味 |
|---|---|---|---|
| `J0` | A→A | `[1, 0]` | 最近接 W 方向 |
| `J1` | A→A | `[0, 1]` | 最近接 L 方向 |
| `J0'` | A→A | `[1, 1]` | 次近接 W+L 対角 |
| `J1'` | A→A | `[1, -1]` | 次近接 W-L 対角 |
| `J0''` | A→A | `[2, 0]` | 三次近接 2W |
| `J1''` | A→A | `[0, 2]` | 三次近接 2L |

**検算**: `n_sites_uc: 1`、全 type とも `coordination.A: 2`。最近接系列
合計(`J0+J1`)= 4(正方格子の物理配位数と一致)、対角系列・軸三次近接
系列もそれぞれ合計 4。`min_size_for_check: [5, 5]`(最大 `|R|` 成分 2、
`2×2=4` を超える最小の奇数 5)。

**出典**: `python/stdface/lattice/square_lattice.py::tetragonal`
(`_BONDS`)、commit `fe8d43f4722a880072c3ac0f3cafd99d741dff77`。関数名の
注: brief 上の想定名は `square` だったが実装上の関数名は `tetragonal`
であり(`SquarePlugin` が `["tetragonal", "square", ...]` のエイリアス
経由で `tetragonal()` に委譲、2.9 節のレジストリ正規名表を参照)、
manifest には実際のソースに合わせて `func: tetragonal` を記録している。

### 3.4 triangular

**幾何**: `dimension: 2`。`a1: [1.0, 0.0]`,
`a2: [0.5, 0.8660254037844386]`(= `(1/2, √3/2)`、非直交)。単一サイト
`A`(`frac: [0.0, 0.0]`)。

**ボンド定義表**(Spin, `_BONDS` 9 行。次近接ブロックはソース上
`J1', J2', J0'` の順(数値順ではない)であり、CONVENTIONS.md §5 の
ソース順保持規約によりそのまま列挙する):

| type | from→to | R | 意味 |
|---|---|---|---|
| `J0` | A→A | `[1, 0]` | 最近接 W |
| `J1` | A→A | `[0, 1]` | 最近接 L |
| `J2` | A→A | `[1, -1]` | 最近接 W-L |
| `J1'` | A→A | `[2, -1]` | 次近接 2W-L |
| `J2'` | A→A | `[1, 1]` | 次近接 W+L |
| `J0'` | A→A | `[-1, 2]` | 次近接 -W+2L |
| `J0''` | A→A | `[2, 0]` | 三次近接 2W |
| `J1''` | A→A | `[0, 2]` | 三次近接 2L |
| `J2''` | A→A | `[2, -2]` | 三次近接 2W-2L |

**検算**: 最近接配位数 = `J0+J1+J2` = 2+2+2 = **6**(三角格子の物理配位数
と一致)。`min_size_for_check: [5, 5]`(最大 `|R|` 成分 2)。

**出典**: `python/stdface/lattice/triangular_lattice.py::triangular`
(`_BONDS`)、commit `c9ef7c9df93fc73afcf9076cd8bacddbf1d9bf24`。

### 3.5 honeycomb

**幾何**: `dimension: 2`, `n_sites_uc: 2`。`a1: [1.0, 0.0]`,
`a2: [0.5, 0.8660254037844386]`(triangular と同じ格子ベクトル式)。
`A`(`frac: [0.0, 0.0]`)、`B`(`frac: [0.3333333333333333,
0.3333333333333333]` = `(1/3, 1/3)`)。

**ボンド定義表**(Spin, `_BONDS` 12 行。最近接の Kitaev 対応:
`J0` = セル内ボンド(z 型)、`J1` = W 方向(x 型)、`J2` = L 方向(y 型)):

| type | from→to | R | 意味 |
|---|---|---|---|
| `J0` | A→B | `[0, 0]` | 最近接 セル内(Kitaev z) |
| `J1` | B→A | `[1, 0]` | 最近接 W 方向(Kitaev x) |
| `J2` | B→A | `[0, 1]` | 最近接 L 方向(Kitaev y) |
| `J2'` | A→A | `[1, 0]` | 次近接 W, 副格子内(0→0) |
| `J2'` | B→B | `[1, 0]` | 次近接 W, 副格子内(1→1) |
| `J1'` | A→A | `[0, 1]` | 次近接 L, 副格子内(0→0) |
| `J1'` | B→B | `[0, 1]` | 次近接 L, 副格子内(1→1) |
| `J0'` | A→A | `[1, -1]` | 次近接 W-L, 副格子内(0→0) |
| `J0'` | B→B | `[1, -1]` | 次近接 W-L, 副格子内(1→1) |
| `J1''` | A→B | `[1, -1]` | 三次近接 |
| `J0''` | A→B | `[-1, -1]` | 三次近接 |
| `J2''` | A→B | `[-1, 1]` | 三次近接 |

**検算**: `A`/`B` は対称。最近接(`J0/J1/J2` 各 1 タイプ)配位数 1 ずつ
(合計 3 = ハニカム格子の物理配位数)、次近接(同一副格子内、`J0'/J1'/J2'`)
配位数 2 ずつ、三次近接(`J0''/J1''/J2''`)配位数 1 ずつ。
`min_size_for_check: [3, 3]`(最大 `|R|` 成分 1)。Kondo は 2 物理サイト
×2 ラベル = 4 ラベル(`A_c, A_s, B_c, B_s`)、Kondo 結合 2 本
(`A_c→A_s`, `B_c→B_s`、いずれも `R=[0,0]`)。

**出典**: `python/stdface/lattice/honeycomb_lattice.py::honeycomb`
(`_BONDS`)、commit `e42ef547c5fb073d78200752c5c5ff01375d4619`。

### 3.6 kagome

**幾何**: `dimension: 2`, `n_sites_uc: 3`。`a1: [1.0, 0.0]`,
`a2: [0.5, 0.8660254037844386]`(honeycomb と同じ格子ベクトル式)。
`A`(`frac: [0.0, 0.0]`)、`B`(`frac: [0.5, 0.0]`)、
`C`(`frac: [0.0, 0.5]`)。三次近接(`J''` 系列)は参照実装
(`kagome.py`)に一切現れず、本カタログにも存在しない。

**ボンド定義表**(Spin, `_BONDS` 12 行: 最近接 6 + 次近接 6):

| type | from→to | R | 意味 |
|---|---|---|---|
| `J2` | A→B | `[0, 0]` | 最近接 intra 0→1 |
| `J1` | A→C | `[0, 0]` | 最近接 intra 0→2 |
| `J0` | B→C | `[0, 0]` | 最近接 intra 1→2 |
| `J2` | B→A | `[1, 0]` | 最近接 W 方向 |
| `J1` | C→A | `[0, 1]` | 最近接 L 方向 |
| `J0` | B→C | `[1, -1]` | 最近接 W-L 方向 |
| `J1'` | C→A | `[1, 0]` | 次近接 W, 2→0 |
| `J0'` | B→C | `[1, 0]` | 次近接 W, 1→2 |
| `J2'` | B→A | `[0, 1]` | 次近接 L, 1→0 |
| `J0'` | C→B | `[0, 1]` | 次近接 L, 2→1 |
| `J1'` | A→C | `[1, -1]` | 次近接 W-L, 0→2 |
| `J2'` | A→B | `[-1, 1]` | 次近接 L-W, 0→1 |

**検算**: 最近接配位数 4/サイト(各ラベルは 3 タイプ中 2 タイプに触れ、
それぞれ配位数 2: A は J1/J2、B は J0/J2、C は J0/J1)。次近接系列も
同じパターンで配位数 4/サイト。`min_size_for_check: [3, 3]`(最大
`|R|` 成分 1)。

**出典**: `python/stdface/lattice/kagome.py::kagome`(`_BONDS`)、
commit `20b2c3cda3a489b5dcd4c9180a19ca23e9280c73`。

### 3.7 orthorhombic(単純立方、3D 基本テンプレート)

**幾何**: `dimension: 3`。`a1: [1,0,0], a2: [0,1,0], a3: [0,0,1]`
(単位行列)。単一サイト `A`(`frac: [0,0,0]`)。**境界位相の扱いが
chain/ladder と異なる**: 内部転写は行われず、`phase0/1/2` はそのまま
`system.boundary[0..2]` に 1:1 対応する(`orthorhombic.py` は
`StdI.phase[0..2]` を直接 `phase0/1/2` から設定する)。

**ボンド定義表**(Spin, `_BONDS` 13 行: 最近接 3 + 次近接(面対角)6 +
三次近接(体対角)4):

| type | from→to | R | 意味 |
|---|---|---|---|
| `J0` | A→A | `[1,0,0]` | 最近接 W 方向 |
| `J1` | A→A | `[0,1,0]` | 最近接 L 方向 |
| `J2` | A→A | `[0,0,1]` | 最近接 H 方向 |
| `J0'` | A→A | `[0,1,1]` | 次近接 面対角 +L+H |
| `J0'` | A→A | `[0,1,-1]` | 次近接 面対角 +L-H |
| `J1'` | A→A | `[1,0,1]` | 次近接 面対角 +H+W |
| `J1'` | A→A | `[-1,0,1]` | 次近接 面対角 +H-W |
| `J2'` | A→A | `[1,1,0]` | 次近接 面対角 +W+L |
| `J2'` | A→A | `[1,-1,0]` | 次近接 面対角 +W-L |
| `J''` | A→A | `[1,1,1]` | 三次近接 体対角(単一大域 type) |
| `J''` | A→A | `[-1,1,1]` | 三次近接 体対角 |
| `J''` | A→A | `[1,-1,1]` | 三次近接 体対角 |
| `J''` | A→A | `[1,1,-1]` | 三次近接 体対角 |

三次近接だけは `J0''/J1''/J2''` の系列別ではなく、単一の大域 `type: "J''"`
にまとめられる点が他の格子(chain/square/triangular)の命名規則と異なる。

**検算**: 最近接配位数 2/type × 3type = 6、次近接 4/type × 3type = 12、
三次近接 8(単一 type)— 単純立方格子の物理配位数(6/12/8)と一致。
`min_size_for_check: [3, 3, 3]`(最大 `|R|` 成分 1)。C 実装
`src/Orthorhombic.c` の 13 回の `StdFace_FindSite` 呼出しと 1 対 1 で
突合済み(Task 7)。この 3D 構造(`system.boundary` の 3 成分、13 行の
`_BONDS`)が fc_ortho・pyrochlore の 3D テンプレートとして踏襲される。

**出典**: `python/stdface/lattice/orthorhombic.py::orthorhombic`
(`_BONDS`)、commit `1cf9f86a394d01f9efe180031e14bfa5d344bfcc`
(`src/Orthorhombic.c` と突合済み)。

### 3.8 fc_ortho(面心直方格子)

**幾何**: `dimension: 3`。orthorhombic の単位行列とは異なり真に
面心構造: `a1: [0.0, 0.5, 0.5]`, `a2: [0.5, 0.0, 0.5]`,
`a3: [0.5, 0.5, 0.0]`。単一サイト `A`(`frac: [0,0,0]`)。境界は
orthorhombic と同じ 3 成分 `phase0/1/2`。

**ボンド定義表**(Spin, `_BONDS` 9 行: 最近接 6(3 系列×等価方向 2) +
次近接 3(系列毎 1 方向)):

| type | from→to | R | 意味 |
|---|---|---|---|
| `J0` | A→A | `[1,0,0]` | 最近接 W(等価対 1/2) |
| `J0` | A→A | `[0,1,-1]` | 最近接 W(等価対 2/2) |
| `J1` | A→A | `[0,1,0]` | 最近接 L(等価対 1/2) |
| `J1` | A→A | `[-1,0,1]` | 最近接 L(等価対 2/2) |
| `J2` | A→A | `[0,0,1]` | 最近接 H(等価対 1/2) |
| `J2` | A→A | `[1,-1,0]` | 最近接 H(等価対 2/2) |
| `J0'` | A→A | `[-1,1,1]` | 次近接 -W+L+H |
| `J1'` | A→A | `[1,-1,1]` | 次近接 -L+H+W |
| `J2'` | A→A | `[1,1,-1]` | 次近接 -H+W+L |

**`J''`/`t''`/`V''` は stan.in で指定してもエラーにならず、かつ物理
ハミルトニアンには一切反映されない**という結果だけを見ると 3 者は
同じに見えるが、実装機構は異なる。本カタログのコードベースには
「未使用キーワード検出器」(未参照の入力キーワードを走査して警告する
仕組み)は存在しないため、最終的な挙動はいずれも「黙って受理・黙って
破棄・エラーなし」に帰着する — 以下はその内部で何が起きているかを
区別する。

1. **`J0''`/`J1''`/`J2''`(意図的な accept-then-drop)**: Spin 分岐は
   `input_spin_nn(StdI.Jpp, StdI.JppAll, StdI.J0pp, ..., "J0''")` 等
   (`fc_ortho.py` 89–91 行)によって明示的に読み取り、解決済みの値を
   `StdI.J0pp` 等に格納する。しかしこの値は `_BONDS` テーブルにも
   `general_j` 呼出しにも一切現れない — **読み取ってから捨てる**、
   意図的な設計判断であることが C 実装 `src/FCOrtho.c` の相互作用
   配列サイジング式のコメントから確認できる(210 行):
   ```c
   nintrMax = StdI->NCell * (StdI->NsiteUC/*D*/ + 6/*J*/ + 3/*J'*/ + 0/*J''*/) ...
   ```
   `+ 0/*J''*/` という記述自体が、C 実装作者が `J''` に配列スロットを
   意図的に一切割り当てていない(= 対応するボンドが存在しないことを
   見越している)ことを示す。
2. **`t''`(`StdI.tpp`)(未参照)**: Hubbard/Kondo 分岐(`else`,
   108–121 行)は `t0/t1/t2` を `input_hopp(StdI.t, ...)` から、
   `t0'/t1'/t2'` を `input_hopp(StdI.tp, ...)` から解決するが、
   単一プライムの系列で打ち切られており、`t0''/t1''/t2''` を
   `input_hopp(StdI.tpp, ...)` から解決する呼出しは**存在しない**
   — すなわち `StdI.tpp` は Hubbard/Kondo 分岐のロジックから一度も
   参照されない(共通キーワードパーサが stan.in の `t''` の値を
   `StdI.tpp` に格納するだけで、fc_ortho 側は読みも捨てもしない)。
   J'' のような「読み取ってから捨てる」明示的な処理とは異なる。
   (Spin 分岐にある `not_used_d("t''", StdI.tpp)`(103 行)は、
   t 族全体を spin 模型では使わないという定型的な拒否リストの一部
   であり、Hubbard/Kondo 分岐での `t''` の扱いとは無関係。)
   C 実装(`src/FCOrtho.c` 147–152 行)も同型で、`InputHopp` 呼出しは
   `t0'/t1'/t2'` までで打ち切られている。
3. **`V''`(`StdI.Vpp`)(未参照、かつ拒否リストからも欠落)**:
   `t''` と同様、Hubbard/Kondo 分岐に `V0''/V1''/V2''` を
   `input_coulomb_v(StdI.Vpp, ...)` から解決する呼出しは存在しない
   (153–158 行は `V0'/V1'/V2'` までで打ち切り)。**さらに**、Spin
   分岐の `not_used_d` 拒否リスト(93–107 行)も `V'`(`StdI.Vp`、
   107 行)までで打ち切られており、`not_used_d("V''", StdI.Vpp)` に
   相当する行が**存在しない**。これは兄弟格子である
   `square_lattice.py`(Spin 分岐、110 行:
   `not_used_d("V''", StdI.Vpp)`)や `triangular_lattice.py`
   (Spin 分岐、113 行:同様の呼出し)が `V''` を明示的に拒否リストへ
   含めているのと対照的であり、C 実装 `src/FCOrtho.c` を確認しても
   同じ欠落が存在する(139–142 行の `StdFace_NotUsed_d` 系列は
   `V'` で打ち切られ、`V''` の行がない)。C/Python 両実装で一貫して
   `V''` だけが拒否リストから欠けているため、`J''` のような意図的な
   設計判断というよりは**上流(C 実装)側の実装上の空隙(不整合)が
   Python 移植にもそのまま引き継がれた可能性が高い** — 7 章
   (既知の制限)で改めて記載する。

以上のとおり、`J''` は「読んで捨てる」、`t''` は「(Hubbard/Kondo
分岐からは)読まれもしない」、`V''` は「読まれもせず、かつ Spin 分岐の
拒否リストからも漏れている」という 3 通りの異なる内部機構でありながら、
ユーザから見た挙動(エラーなく受理され、ハミルトニアンには反映されない)
だけが共通している。本カタログはこの現行実装のとおり `J''`/`t''`/
`V''` を `bonds`/`couplings` から除外している。

**検算**: 最近接配位数 4/type(2 ソース行×2)、次近接配位数 2/type
(1 ソース行×2)。`min_size_for_check: [3, 3, 3]`(最大 `|R|` 成分 1)。

**出典**: `python/stdface/lattice/fc_ortho.py::fc_ortho`(`_BONDS`)、
commit `8e3601167d1254807f58ce365107333f1d505fa0`
(`src/FCOrtho.c::StdFace_FCOrtho` と byte-for-byte の呼出し列一致を
確認済み)。

### 3.9 pyrochlore

**幾何**: `dimension: 3`, `n_sites_uc: 4`(4 面体の 4 頂点)。格子ベクトルは
fc_ortho と同じ FCC: `a1: [0,0.5,0.5], a2: [0.5,0,0.5], a3: [0.5,0.5,0]`。
`A0`(`frac: [0,0,0]`)、`A1`(`frac: [0.5,0,0]`)、`A2`(`frac: [0,0.5,0]`)、
`A3`(`frac: [0,0,0.5]`)。

**ボンド定義表**(Spin, `_BONDS` 12 行: 四面体内(`R=[0,0,0]`)6 +
四面体間 6):

| type | from→to | R | 意味 |
|---|---|---|---|
| `J0` | A0→A1 | `[0,0,0]` | 四面体内 W |
| `J1` | A0→A2 | `[0,0,0]` | 四面体内 L |
| `J2` | A0→A3 | `[0,0,0]` | 四面体内 H |
| `J0'` | A2→A3 | `[0,0,0]` | 四面体内 L-H |
| `J1'` | A3→A1 | `[0,0,0]` | 四面体内 H-W |
| `J2'` | A1→A2 | `[0,0,0]` | 四面体内 W-L |
| `J0` | A1→A0 | `[1,0,0]` | 四面体間 W |
| `J1` | A2→A0 | `[0,1,0]` | 四面体間 L |
| `J2` | A3→A0 | `[0,0,1]` | 四面体間 H |
| `J0'` | A3→A2 | `[0,-1,1]` | 四面体間 L-H |
| `J1'` | A1→A3 | `[1,0,-1]` | 四面体間 H-W |
| `J2'` | A2→A1 | `[-1,1,0]` | 四面体間 W-L |

四面体内の `(A2,A3)`/`(A3,A1)`/`(A1,A2)` という一見エンドポイントの
向きが揃っていない組合せは `_BONDS` の生の順序をそのまま保持したもので、
並べ替えは行っていない(1.3 節(2))。

**Kondo の非対称 J 結合(pyrochlore 固有の例外)**: 他の 8 格子はいずれも
`X_c → X_s`(自分自身の副格子内で閉じる)という Kondo 結合パターンを
持つが、pyrochlore のみ現行実装(C/Python 共通)が非対称である。
`_local` 関数(`pyrochlore.py` 169–173 行)/`Pyrochlore.c`(246–251 行)は
`general_j(StdI, StdI.J, 1, StdI.S2, isite + 3, jsite + uc_i)` を
`uc_i = 0..3` についてループしており、**副格子 3 の遍歴サイト `A3_c`
のみ**が全 4 副格子の局在スピン(`A0_s..A3_s`)全てと結合する:

```yaml
- {from: A3_c, to: A0_s, R: [0, 0, 0], type: J}
- {from: A3_c, to: A1_s, R: [0, 0, 0], type: J}
- {from: A3_c, to: A2_s, R: [0, 0, 0], type: J}
- {from: A3_c, to: A3_s, R: [0, 0, 0], type: J}
```

`A0_c`/`A1_c`/`A2_c` は Kondo 結合(`J`)を一切持たない(`coordination`
に `J` キーが現れない)。この非対称性は上流実装のバグである可能性が
あるが、本カタログの検証水準は「現行実装ソースとの突合」であるため
上流の挙動をそのまま採用する(CONVENTIONS.md §6.7)。上流が将来
修正された場合はカタログの schema バージョン上げが必要になる旨は
7 章(既知の制限)で改めて記載する。

**検算**: Spin/Hubbard の最近接配位数 6/サイト(3 type × 2)。
`min_size_for_check: [3, 3, 3]`(最大 `|R|` 成分 1)。Kondo:
`coordination.A3_c.J: 4`、`coordination.A{0,1,2,3}_s.J: 1`(各 1)、
`A0_c`/`A1_c`/`A2_c` に `J` キーなし(計数 0 のため省略、manifest.yaml
の記載規則どおり)。

**出典**: `python/stdface/lattice/pyrochlore.py::pyrochlore`(`_local`,
`_BONDS`)、commit `70fbecd02aeb60df04625843217a08b545c2c8b0`
(`src/Pyrochlore.c` と byte-for-byte の一致を確認済み、Task 9)。

---

## 4. 模型ごとの演算子対応

### 4.1 符号表と検証連鎖

1.3 節(3)/CONVENTIONS.md §6.2 の符号表は「規範」として提示済みだが、
本節ではその**導出**を示す。検証連鎖は

```
StdFace パラメータ → builder 呼び出し(interaction_builder.py)
→ trans/intr 係数(StdI.trans_list / StdI.intr_list / StdI.Cintra_list / StdI.Cinter_list)
→ solver 出力規約(HPhi trans.def の H = −Σ(trans 値) c†c という暗黙符号)
→ 物理ハミルトニアンの符号
```

の 4 段階からなる。**HPhi の `trans.def` の係数規約(書かれた値の符号を
反転して物理ハミルトニアンに使う)はソルバー自身の入力ファイル仕様であり、
本リポジトリの Python 実装から導出されるものではなく、既知の事実として
検証連鎖の最終段に用いる**(CONVENTIONS.md §6.1)。

**チャネルごとに反転の有無が異なる**点が本節の要点である。

- **`trans_list` 系チャネル**(`hop`、onsite の `mu`/`h`/`Gamma`/
  `Gamma_y`): これらはすべて `interaction_builder.py` の `trans()`
  (115–144 行)を経由して `StdI.trans_list` に積まれ、`trans.def` の
  `H = −Σ(trans 値) c†c` という暗黙反転を**1 回だけ**受ける。
  代表例として `hop`(t 族)を追跡する:
  1. StdFace パラメータ: `input_hopp(StdI.t, StdI.t0, "t0")` が
     stan.in から `t0`(複素数可)を解決する。
  2. builder 呼び出し: 標準格子の `_dispatch_bond_interaction`
     (`interaction_builder.py` 630–663 行、661 行)は
     `hopping(StdI, Cphase * t, isite, jsite, dR)` を呼ぶ —
     **符号は反転されず**、境界位相 `Cphase`(`|Cphase|=1`)を掛けた
     生の `t0` がそのまま渡される。
  3. trans 係数: `hopping()`(147–191 行)は両スピンについて
     `trans(StdI, t0, jsite, σ, isite, σ)` と
     `trans(StdI, conj(t0), isite, σ, jsite, σ)` を追加する —
     ここでも符号反転はない(`trans_list` の値 = 生の `t0`)。
  4. solver 規約: `trans.def` はこの値をそのまま出力し(`writer/
     common_writer.py`)、HPhi は `H_trans = −Σ(trans.def 値) c†c` と
     解釈する。
  5. 物理符号: `H = −t0 Σσ(c†c + h.c.)`。したがって
     `couplings[t0].value = {param: t0, scale: -1.0}` が
     `scale × param = −t0` を与え、この物理係数と一致する。

  onsite の `mu`/`h`/`Gamma`/`Gamma_y` も同じ経路(`hubbard_local_terms`/
  `mag_field_terms`、いずれも `_trans_term` 経由で `trans_list` へ)を通り、
  builder 側の生値がそのまま `trans_list` に入る(例:
  `HubbardModel.build_local_terms` は `hubbard_local_terms(StdI.mu,
  -StdI.h, -StdI.Gamma, -StdI.Gamma_y, StdI.U, isite)` を呼び、
  `mu0 = StdI.mu` は無反転で渡る一方、`h0 = -StdI.h` のように**呼び出し側
  で先に符号を作ってから渡す**項もある — いずれにせよ `trans_list` に
  積まれた後は同じ 1 回の solver 側反転を受けるので、最終的な物理符号は
  1.3 節(3)の表(`mu`: `−mu N`、`h`: `−h Sz` 等)と一致する)。

- **`intr_list`/`Cintra_list`/`Cinter_list` 系チャネル**(`density-density`
  の `V`、onsite の `U`、J 族の交換相互作用): これらは `trans.def` を
  経由しないため、solver 側の暗黙反転を**受けない**。
  - `coulomb(StdI, V, i, j)`(546–560 行)は `StdI.Cinter_list.append((V,
    i, j))` と無反転で追加する → 物理係数 = `+V`(`couplings[V0].value
    = {param: V0}`、`scale` 省略時 `+1.0`)。
  - `hubbard_local_terms`(225–242 行)は `terms.Cintra.append((U0,
    isite))` と無反転で追加する → 物理係数 = `+U`(onsite `coeff:
    +1.0`)。
  - `general_j_terms`(456–543 行)は `J[2,2]*Siz*Sjz` 等を直接
    `terms.intr` に積む(無反転) → 物理係数 = `+J_ab S^a S^b`
    (J 族 `coeff` はそのまま `+1.0`/`-1.0` の符号のみを表すリテラル)。
    Kondo の `J`(`s_i . S_j`)も `general_j(StdI.J, 1, StdI.S2, isite,
    jsite_kondo)` として同じ経路(`intr_list`)を通るため無反転
    (`couplings[J].value = {param: J}`、`scale` 省略時 `+1.0`)。

- **wannier90 の hop チャネルは例外的に 2 回反転する**:
  `wannier90.py::_apply_hopping_terms` は非局所項を
  `hopping(StdI, -Cphase * tUJ[0][it], jsite, isite, dR)` として渡す
  — 標準格子の `_dispatch_bond_interaction` とは異なり、**builder
  呼び出し自体に符号反転(`-Cphase`)が入っている**。この builder 側の
  反転と solver 側の反転(上記 trans_list 系と同じ 1 回)が相殺し、
  結果として `物理ホッピング係数 = +H_mn(R)`(符号反転なし)になる。
  これが chain/square 系の t 族(`scale: -1.0`)と wannier90 の hop
  (`scale: +1.0`)が逆に見える理由であり、詳細な検算は 5.1 節で示す。

以上から、1.3 節(3)の符号表の各行は「`trans_list` 経由か
`intr_list`/`Cintra_list`/`Cinter_list` 経由か」という**チャネルの違い**
だけで一意に説明できる: `trans_list` 経由の項(t, mu, h, Gamma, Gamma_y)
は solver 側の 1 回反転を打ち消すために `scale`/`coeff` に `-1.0`
を持たせ、それ以外(V, U, J)は無反転のため `+1.0`(または省略)で
物理係数と直接一致させる。

### 4.2 Spin: 9 成分正準形とパラメータ解決規則

キーワード成分表(1.3 節(5)/CONVENTIONS.md §6.4 の再掲。ops 対と
param 接尾辞は 1 対 1 対応):

| ops 対 | 接尾辞 | ops 対 | 接尾辞 | ops 対 | 接尾辞 |
|---|---|---|---|---|---|
| `[Sx,Sx]` | `x` | `[Sz,Sz]` | `z` | `[Sy,Sx]` | `yx` |
| `[Sy,Sy]` | `y` | `[Sx,Sy]` | `xy` | `[Sy,Sz]` | `yz` |
| | | `[Sx,Sz]` | `xz` | `[Sz,Sx]` | `zx` |
| | | | | `[Sz,Sy]` | `zy` |

**解決順序**(`input_params.py::_resolve_spin_matrix`、135–183 行。
prefix `J0` を例に取ると各成分 `(i,j)` について):

1. 成分局所(`J0xy` 等。`J0[i,j]` が既に設定済み)
2. 成分大域(`Jxy` 等。フォールバック行列 `J[i,j]`)
3. スカラー局所・対角のみ(`J0`。`i==j` の場合のみ `J0All`)
4. スカラー大域・対角のみ(`J`。`i==j` の場合のみ `JAll`)
5. 0(既定値)

プライム系(`J0'`, `J0''`, `J1'` 等、`input_spin`)には大域 fallback が
存在しない例外がある(1.3 節(5)、CONVENTIONS.md §6.5)。

**等方入力・異方入力・同時指定エラーの 3 例**(いずれも
`input_spin_nn`、`input_params.py` 186–227 行の実装に基づく。`J0`
系列を例に取る):

- **例 1(等方入力)**: stan.in で大域スカラー `J = 1.0` のみを指定し、
  `J0` や `J0x` 等は無指定のまま。解決順序 1–2 は不成立(局所・大域とも
  成分は未設定)、3 も不成立(局所スカラー `J0All` 未設定)、4 で大域
  スカラー `JAll = 1.0` が対角成分(`x, y, z`)にのみ採用される。結果:
  `J0x = J0y = J0z = 1.0`、非対角成分はすべて 0(等方 Heisenberg 交換)。
- **例 2(異方入力)**: stan.in で成分局所値 `J0x = 1.0, J0y = 1.0,
  J0z = 2.0` を指定。解決順序 1 でそれぞれの `(x,x)`/`(y,y)`/`(z,z)`
  成分が直接採用される。結果: `J0x = J0y = 1.0 ≠ J0z = 2.0`(XXZ 型
  異方交換)。
- **例 3(同時指定エラー)**: stan.in で局所スカラー `J0 = 1.0` と
  成分局所 `J0z = 2.0` を同時に指定。`input_spin_nn` は
  `_check_scalar_vs_matrix(J0All, J0(行列), J0name, J0name)`
  (`input_params.py` 220 行)を呼び、`J0All` が非 NaN かつ `J0` 行列内に
  既に非 NaN 要素(`z,z`)があることを検出して
  `ValueError("J0 and J0z conflict !")` を送出する。すなわち
  **等方(スカラー)指定と異方(成分)指定は同一 prefix 内で同時に
  行えない**(スカラー同士・スカラー対行列・行列対行列の全組合せが
  同様にエラーとなる。CONVENTIONS.md §6.5 の競合規則)。

### 4.3 Hubbard: hop / density-density / onsite

演算子の数式(CONVENTIONS.md §6.6):

- `hop` = `Σ_σ (c†_iσ c_jσ + h.c.)`(演算子自体は符号を持たない。符号は
  4.1 節のとおり `value.scale = -1.0` に持たせ、物理寄与は
  `−t Σσ(c†c + h.c.)`)。
- `density-density` = `n_i n_j`(物理寄与 `+V n_i n_j`、`scale` 省略時
  `+1.0`)。
- onsite `U`: `tensor_terms: [{ops: [NupNdn], coeff: +1.0}]`
  (物理寄与 `+U n↑n↓`)。

**複素 `t` と反転同値(共役)**: `hopping()`(interaction_builder.py
147–191 行)は 1 つの `(isite, jsite, R, t)` から
`trans(StdI, t, jsite, σ, isite, σ)` と `trans(StdI, conj(t), isite, σ,
jsite, σ)` の両方を生成する — すなわち `hop` 演算子の `+ h.c.` は
実装上リテラルに複素共役を取ることで実現されている。カタログの
ボンドは `(from, to, R)` の**片方向のみ**を記述するソース順保持規約
(CONVENTIONS.md §5)を取るため、もし読み手が `(from, to, R)` を
`(to, from, −R)` として反転読み替えする場合(通常は行わない)、その
値は元の複素共役 `conj(t)` でなければならない — これは 1.3 節(2)の
「反転時の係数変換(消費側規範)」がまさにこの `hopping()` の実装から
導かれた規則であることを示す。実数の `t`(位相のない通常のホッピング)
では `conj(t) = t` となり反転同値は自明になるが、境界位相(`phase0`
等によるツイスト境界)や磁束(Peierls 位相)が乗ると `t` は複素数
になり得るため、共役を取る規則は一般に必要である。

### 4.4 Kondo: 2 ラベル・端点順序・サイト倍加・GC 変種

**2 ラベルと端点順序**: 1.3 節(2)/CONVENTIONS.md §6.7 のとおり、Kondo
は 1 物理サイトを `<X>_c`(遍歴電子、fermion)・`<X>_s`(局在スピン、
spin)の 2 geometry ラベルとして表現する。Kondo 結合(`s_i . S_j`)は
`general_j(StdI.J, 1, StdI.S2, isite, jsite_kondo, ...)` の引数順が
そのままボンドの `from`(`isite` = 遍歴 `_c`)/`to`(`jsite_kondo` =
局在 `_s`)に写るため、**第 1 端点(`from`)= 遍歴電子スピン、第 2 端点
(`to`)= 局在スピン**という順序に意味を持つ(3 章の各格子のボンド表
参照。pyrochlore のみ非対称な例外を持つ、3.9 節)。

**磁場の両側適用**: `KondoModel.build_local_terms`
(`model_plugin.py` 87–101 行)は `hubbard_local_terms(...)` で遍歴
サイト(`isite`)に Hubbard 一式(`mu`, `U`, 磁場)を適用したのち、
`mag_field_terms(StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y,
jsite_kondo)` で局在サイト(`jsite_kondo`)にも磁場(`h`/`Gamma`/
`Gamma_y`)を**改めて**適用する。したがって `A_c`/`A_s` の両方の
onsite ブロックに `field_z`/`field_x`/`field_y` が現れる(chain_kondo
の 1.4 節の例で確認済み)一方、`U`/`mu`/`D` は `_c` にのみ現れる。

**サイト倍加(前半 = 局在)との対応**: 参照実装は Kondo 模型の内部
サイトインデックスを **倍加**し(`set_local_spin_flags`、
`site_util.py` 704–736 行: `StdI.nsite *= 2`)、**前半**
(`locspinflag[:half] = StdI.S2`)を局在スピン、**後半**
(`locspinflag[half:] = 0`)を遍歴電子(fermion)に割り当てる。この
対応は `expand_bonds_2d`(`interaction_builder.py` 807–848 行、838–845
行)の呼び出しにそのまま現れる:

```python
kondo_off = StdI.NsiteUC * StdI.NCell if StdI.model == ModelType.KONDO else 0
...
add_local_terms(StdI, base + uc + kondo_off, base + uc)
#                      ^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^
#                      isite = 遍歴(後半、offset あり)   jsite_kondo = 局在(前半、offset なし)
```

すなわち内部インデックス空間では「前半(offset なし)= 局在スピン、
後半(`kondo_off` 分オフセット)= 遍歴電子」という対応になっており、
`add_local_terms(StdI, isite, jsite_kondo)` の引数順(遍歴が第 1
引数)とは向きが逆であることに注意(引数順は「意味上の役割」、
サイトインデックスの前半/後半は「メモリレイアウト」であり、両者は
独立した規約である)。カタログの `geometry.sites` はこの内部
インデックスの前半/後半を直接は反映せず、`_c`/`_s` ラベルという
意味論のみを公開する。

**GC 変種(粒子数条件は solver 層)**: HPhi は `CalcModel` の値
(`0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC` —
この説明文字列は `calcmod.def` に出力される、`solvers/hphi/writer.py`
306 行。対応表そのものは同ファイル `MODEL_GC_TO_CALC_MODEL` 辞書、
80–87 行)により通常(canonical)模型とグランドカノニカル(GC)模型を
切り替える。この
切り替えは `(ModelType, lGC)` の組で決まり、`lGC` は全電子数
(`ncond`/`nelec`)・全 `Sz`(`2Sz`)という**粒子数セクターの指定**から
solver 出力層で決定される — ハミルトニアンの `bonds`/`couplings`/
`onsite` 定義そのもの(本カタログが記述する対象)には一切依存しない。
これは 2.8 節で `ncond, nelec, 2Sz` を「対象外」(計算条件であり
模型のハミルトニアン定義ではない)に分類した理由そのものであり、
Kondo の GC 変種(`KondoGC`)であってもカタログ上の `bonds`/
`couplings`/`onsite` は通常の Kondo と同一である。

---

## 5. wannier90 変換仕様

### 5.1 概要

wannier90 格子(`lattice_catalog/wannier90/`)は、他の 9 格子のように
StdFace の解析的な格子ベクトル・ボンド式から `bonds`/`couplings`/
`onsite` を導出するのではなく、RESPACK/Wannier90 が出力する
`*_hr.dat`(hopping)・`*_ur.dat`(オンサイト/サイト間 Coulomb)・
`*_jr.dat`(Hund/交換/ペアホッピング)という**外部データファイル**を
読み込んで変換する。対応する変換ロジックは
`python/stdface/lattice/wannier90.py`(呼出しエントリ `wannier90()`、
842–933 行)・`wannier90_io.py::_read_w90`(データ読込み)・
`interaction_builder.py::hopping`/`trans`/`coulomb`(4 章で確認した
標準の trans/intr 経路)の 3 ファイルにまたがる。本章の規則(W1–W6)は
`lattice_catalog/wannier90/example_hubbard.yaml` のヘッダコメントに
逐語的な導出として既に記載されており(Task 10 で `trans.def`/
`coulombintra.def` の実出力との数値突合により検証済み)、本章はそれを
manual 側の記述として集約する。

**対応模型は Spin・Hubbard のみ**であり、**Kondo は非対応**である
(`wannier90.py::_validate_wannier_params` 618–641 行が
`ModelType.KONDO` の場合に `ValueError("wannier + Kondo is not
available !")` を送出する)。したがって本章に Kondo 節はない。

### 5.2 Hubbard: H/U チャネル

**規則 W1(`H_mn(R)` の符号 — 検証済みの正しい規則)**:
`wannier90.py::_apply_hopping_terms`(389–406 行)は非局所項を

```python
hopping(StdI, -Cphase * tUJ[0][it], jsite, isite, dR)
```

として渡す。4.1 節で確認したとおり標準格子の `_dispatch_bond_interaction`
は `hopping(StdI, Cphase * t, ...)` と**無反転**で呼ぶのに対し、
wannier90 経路はここで明示的に `-Cphase`(符号反転)を掛けている。この
builder 側の反転と、trans.def の solver 側反転(4.1 節で確認した
`H = −Σ(trans 値) c†c`)が相殺し、

```
物理ホッピング係数 = +H_mn(R)   (符号反転なし。Cphase は境界位相であり別枠)
```

となる。標準格子の t 族が `scale: -1.0` を持つのに対し、wannier90 の
`hop` 係数は**符号反転なし**(生データをそのまま物理係数として使う)
という違いが生じるのはこのためである。Task 10 では `W=4, L=4,
Height=1` で実際に `stdface_main` を実行し、`trans.def` の出力値
(`+1.0`)が `H_00(R) = -1.0` の符号反転(`-(+1.0) = -1.0`)である
ことを確認しており、`example_hubbard.yaml` のヘッダコメント
「検算 1/3」に手順が記載されている。

**onsite 分離(`R=0`, `m=n` の対角要素)**: `_apply_hopping_terms`
380–388 行は `R=(0,0,0)` かつ `m=n` の項を `bonds` ではなく
`isite` ごとの onsite 一体項として `trans_list` に分離する
(`(-tUJ[0][it], isite, spin, isite, spin)` を両スピンに追加)。
trans.def 規約と合わせると物理係数は `+H_mm(0)`(両スピンの数演算子和
に掛かる)。`example_hubbard.yaml` のデータでは `H_11(0) = 0.0` であり、
振幅 cutoff(既定 `1e-8`)未満のため実際には onsite 項として現れない
(ゼロなので省略。規則自体は非ゼロなら `{ops:[N], coeff:1.0},
value:<H_mm(0)>` という形になる)。

**Hermite 正準対の選択(first-in-file-wins)**: `wannier90_io.py::
_read_w90`(355–361 行)は WSC(R ベクトル)ブロック単位で「ファイル
出現順で先着優先」の規則を持つ: あるブロックの `R` が、既に読んだ
(より前に出現した)ブロックの `−R` と一致する場合、そのブロック
全体(全 `m,n` 組)を丸ごとゼロ化する。`R=(0,0,0)` のブロックのみ
別途、`m>n`(下三角)要素を追加でゼロ化する(`m<=n` が正準)。これは
「`R>0`」のような固定的な数学規則ではなく、**ファイル中の出現順**に
依存する規則である点に注意。排除された反転対は `hop` 演算子の定義
(`Σσ(c†c + h.c.)`)により自動的に再構成されるため、カタログの
`bonds` には正準対のみを列挙する(4.3 節の複素共役規則と同じ機構)。

**`_hr.dat`/`_ur.dat` 自身の縮退重みは読み捨てる**:
`wannier90_io.py::_skip_degeneracy_weights`(56–72 行、呼び出しは
315 行)はファイルヘッダの `ndegen` 列を読み込んだ後**使用せず捨てる**。
実際に乗算される重みは `Weight_tot`(初期値 1.0)のみで、これは
`_apply_boundary_weights`(129–169 行、呼び出しは 364 行)が `W/L/
Height` から計算する**有限クラスタ境界での 0.5 halving**であり、
ファイルの `ndegen` 値とは無関係。本カタログの `bonds` は無限格子の
単位胞相対表現であるためこの境界重みも適用しない。

**cutoff の適用順**(`_read_w90`、255–370 行): (1) 実空間長の cutoff
(`cutoff_length_t/u/j`、既定は `t` が無効・`U`/`J` が 0.3)、
(2) `cutoff_Vec` 設定時はそれ、未設定なら整数 box cutoff
(`cutoff_R` 系。`U` チャネルの既定 `cutoff_UR` は `(0,0,0)` — 既定では
`R=0` のみ通過)、(3) 正準対選択(上記)、(4) 境界重み(有限クラスタ
固有、カタログでは不適用)、(5) 振幅 cutoff(`cutoff_t`/`cutoff_u`/
`cutoff_j`、既定 `1.0e-8`)未満の要素を EFFECTIVE term から除外。
`example_hubbard.yaml` は既定の振幅 cutoff のみを適用し、長さ/R/Vec
による追加制限は行わない(データの `R` は全て `|R 各成分| <= 1`)。

**`lambda`/`alpha`/`doublecounting`**: `lambda_U`/`lambda_J`(既定
1.0)は各チャネルの生データ読込み時に一様に乗算されるスケール
(`wannier90_io.py` 353 行: `Mat_tot = lam * (re + i*im)`)。`alpha`
(既定 0.5)と `doublecounting_mode` は `_dr.dat` 密度行列ファイルが
存在し `doublecounting_mode != "none"` の場合のみ、onsite 二重計上
補正項(および FULL モードでは追加の hopping 補正項)を
`trans_list` に加える(`_apply_coulomb_terms` 458–509 行)。
`example_hubbard.yaml` には `_dr.dat` が存在せず既定モード `"none"`
のため一切発動しない。

### 5.3 Spin: 超交換アルゴリズムと直接写像との差異

**現行アルゴリズム**: Spin 模型の wannier90 変換は `_jr.dat`
(Hund/交換チャネル)を Spin 用の交換相互作用として直接読み込む
のではなく、`_hr.dat`(hopping)と `_ur.dat`(オンサイト U)から
**2 次摂動の超交換公式**によって等方 Heisenberg 交換を**自動生成**
する。`_apply_hopping_terms`(397–404 行):

```python
if StdI.model == ModelType.SPIN:
    diag_val = (2.0 * tUJ[0][it] * np.conj(tUJ[0][it])
                * (1.0/Uspin[m] + 1.0/Uspin[n])).real
    Jtmp = np.diag([diag_val, diag_val, diag_val])
    general_j(StdI, Jtmp, StdI.S2, StdI.S2, isite, jsite)
```

すなわち各非局所ホッピング要素 `t_mn = H_mn(R)` と、同じ `_ur.dat`
から抽出したオンサイト U(`Uspin[m] = U_mm(0)`, `Uspin[n] = U_nn(0)`、
684–694 行で `_apply_coulomb_terms` と同じ「局所項」判定ロジックを
使って別途抽出)から、`J = 2|t_mn|²(1/U_m + 1/U_n)` という単一の
等方対角テンソル(`J_x = J_y = J_z = J`, 非対角 0)を計算し、
**Hubbard 模型なら `hop` ボンドになる同じ `(isite, jsite, R)` の位置に**
J テンソルの `couplings` として書き込む。

**直接写像との差異**: 「直接写像」— すなわち `_jr.dat` の交換
チャネルの値をそのまま Spin 模型の交換係数として使う経路 — も
別途存在するが(`_apply_hund_terms` 604–610 行、`StdI.Ex_list` に
`ex_val`(mVMC ソルバーは `+tUJ[2][it].real`、それ以外は
`-tUJ[2][it].real` — ソルバーによって符号が異なる)を追加する)、
これは `_jr.dat` が存在する場合に**上記の超交換生成に追加で**働く
副次的な経路である。つまり:

- 超交換生成(`_hr.dat` + `_ur.dat` から)は `_jr.dat` の有無に関わらず
  **常に**実行され、Spin 模型の主要な近接交換はこの経路で決まる。
- `_jr.dat` に基づく直接写像(`Ex_list`)は `_jr.dat` が存在する場合
  のみ**追加的に**発生し、超交換の代わりにはならない(両者は加算的)。
- 超交換で生成される J テンソルは常に等方(対角のみ、非対角 0)である
  のに対し、`_jr.dat` 直接写像は交換相互作用専用の別チャネルの生値を
  そのまま使う。

この非対称性(自動生成が既定・データ由来の直接値が例外的な追加項)は
StdFace 実装が「多軌道 Hubbard/Kondo 模型を強相関極限で有効スピン
模型へ落とし込む」という物理的な近似(超交換描像)を組み込んだ
コンバータであることを反映しており、Wannier 軌道間ホッピングを字面
どおり `hop` ボンドとして書き写す Hubbard 変換(5.2 節)とは根本的に
異なるアルゴリズムである点に注意が必要である。本カタログには Spin
用 wannier90 の YAML 例は用意していない(§7 既知の制限参照)。

### 5.4 `example_hubbard.yaml` の読み解き

`lattice_catalog/wannier90/example_hubbard.yaml` は
`test/wannier90_data/{zvo_geom.dat, zvo_hr.dat, zvo_ur.dat}`
(2D 正方格子、単一 Wannier 軌道、最近接 `t = -1.0 eV`、オンサイト
`U = 8.0 eV`、`zvo_jr.dat` は存在しない)を Hubbard 模型として変換した
最小例である。

- `geometry.lattice_vectors`: `a1=[3,0,0], a2=[0,3,0], a3=[0,0,10]`
  (`zvo_geom.dat` 1–3 行目そのまま。`a3` の 10 は真空層 — 物理的には
  2D スラブ)。`sites: [{label: W1, frac: [0,0,0]}]`
  (`zvo_geom.dat` 5 行目の Wannier 中心、`NsiteUC=1`)。
- `model.site_dof.W1: {fermion: {orbitals: 1}}`(拡張方言、
  CONVENTIONS.md §7-1)。
- `model.bonds`: 規則 W3 により正準対のみを 2 本列挙する
  (`R=[-1,0,0]` type `t0`、`R=[0,-1,0]` type `t1`。`zvo_hr.dat` の
  WSC 出現順 `(-1,0,0)→(0,-1,0)→(0,0,0)→(0,1,0)→(1,0,0)` のうち
  後半 2 つは前半 2 つの反転と一致するため除外される)。
- `model.couplings.t0`/`t1`: `{operator: "hop", value: -1.0}` —
  **`{param: ...}` ではなくリテラル数値**を直接記載している。これは
  wannier90 由来の値が外部データそのものであり stan.in の
  パラメータ参照に対応しないためで、`lint_catalog.py` の C7(param
  参照の目録整合検査)は `value` が dict でない場合はチェックを
  スキップする設計になっており、この直値記載でリンタを通過する。
  値は規則 W1 により符号反転なし(`+H_00(R) = -1.0`、そのまま検算 1
  に対応)。
- `model.onsite.W1.hubbard_u`: `value: 8.0`(`U_11(0) = 8.0` の直値、
  検算 2 に対応)。同じ `onsite.W1` ブロックには `chemical_potential`
  (`{param: mu}`)・`field_z`/`field_x`/`field_y`
  (`{param: h}`/`{param: Gamma}`/`{param: Gamma_y}`)という、
  wannier データに依存しない標準 StdFace パラメータへの `{param:
  ...}` 参照も**混在**している — wannier90 格子でも `mu`/`h`/`Gamma`/
  `Gamma_y` は通常の stan.in パラメータとして受理されるため(4 章の
  onsite 適用範囲の表と同じ扱い)。ただし `U` だけは例外で
  `{param: U}` という参照形は使えない: `wannier90.py::
  _validate_wannier_params` は `lattice="wannier90"` に対して
  `not_used_d("U", StdI.U)` を呼ぶため、`U` は `_ur.dat` から供給
  されなければならず、stan.in の `U` キーワードで与えることは
  許されない(この YAML では `{param: U}` の代わりに直値 `8.0` を
  使っている理由)。
- ヘッダコメント末尾の「規則 W2」注記のとおり、`R=0, m=n` の wannier
  由来 onsite 項(`H_11(0)=0.0`)は振幅 cutoff 未満のため省略されて
  おり、`onsite.W1` に wannier 由来の追加項(`{ops:[N], coeff:1.0}`
  型)は現れない。

### 5.5 J チャネル(Hund・交換・ペアホッピング)

`_jr.dat` から読み込まれる J チャネル(Hund coupling)は
`wannier90.py::_apply_hund_terms`(512–611 行)で処理される。局所項
(`R=0, m=n`)は計算対象外とし(553–555 行)、非局所項(または `R=0`
かつ `m≠n` の軌道間項)についてのみ、単一の係数 `J_mn = tUJ[2][it]`
(実数部)から以下の**3 つの異なる 2 体演算子**を**同時に**生成する
(いずれも符号反転なし、無反転で `StdI.Hund_list`/`Ex_list`/
`PairHopp_list` に積まれ、それぞれ HPhi の `hund.def`/`exchange.def`/
`pairhopp.def` に対応する — `writer/interaction_writer.py` の
`_InteractionMeta` テーブルで確認済み):

- **Hund 項**(`Hund_list`、Hubbard/Spin 共通で常に追加): 軌道 `m,n`
  間の Hund 型密度相互作用。
- **交換項**(`Ex_list`): Hubbard では Hund と同じ係数で追加。Spin
  では 5.3 節の直接写像経路として、ソルバーにより符号が異なる
  (mVMC: `+J_mn`、それ以外: `−J_mn`)値で追加される。
- **ペアホッピング項**(`PairHopp_list`、Hubbard のみ): 軌道間の
  ペア(2 電子)ホッピング `c†_{i↑}c†_{i↓}c_{j↓}c_{j↑}` 型。

すなわち **1 個の J チャネル行列要素が、Hund 密度項・スピン交換項・
ペアホッピング項という 3 つの異なる 4 フェルミオン演算子を同時に
規定する**(多軌道 Hubbard-Kanamori 型相互作用に共通する構造)。
二重計数補正は U チャネル(5.2 節、重み `alpha`)と非対称で、J チャネル
側は重み `(1 - alpha)` を使う(`_apply_hund_terms` 577, 585 行:
`-(1.0 - StdI.alpha) * tUJ[2][it].real * DenMat0`)。

**YAML 例が対象外である理由**: 本カタログの拡張方言(1.3 節(3)/
CONVENTIONS.md §7)が定義する演算子語彙は、2 端点の名前付き演算子
(`hop`, `density-density`, `s_i . S_j`)と 9 成分 J テンソル
(spin–spin)の 2 種類のみであり、いずれも「1 つのボンド `(i,j,R)`
に対して 1 つの物理量」という前提を持つ。J チャネルは 1 つの係数
から Hund/交換/ペアホッピングという**独立した 3 つの 4 フェルミオン
演算子**を同時に生成するため、現行の演算子語彙では表現できない —
これを表現するには、順序付き生成消滅演算子列と site/orbital 束縛を
持つ**一般的な 4 フェルミオン項スキーマ**が必要になる。この
スキーマは draft 仕様にも本カタログの拡張方言にも未定義であるため
(CONVENTIONS.md §7 項目 6)、`lattice_catalog/wannier90/` には J
チャネルの YAML 例を用意していない。

一般項スキーマの素描(6 章の拡張提案 6 項目のうち第 6 項として
記載予定 — 順序付き生成消滅演算子列 `[c†_{i,σ,orb}, c_{j,σ',orb'},
...]` と、各演算子への site/orbital 束縛、係数 1 つを持つ表現。
名前付き演算子(`hop` 等)はこの一般形の省略記法という位置づけに
なる)は、本タスクの範囲外である 6 章(仕様拡張提案)に委ねる。

---

## 6. 仕様拡張提案(draft への追記提案)

### 6.0 本章の位置づけ

CONVENTIONS.md §7(拡張方言一覧)は、本カタログが `dialect:
experimental` の下で採用している draft 仕様の未規定事項を 6 項目
列挙している。本章はその 6 項目それぞれを、**draft 仕様本体への
追記提案**として、(1) `tensor_terms`/YAML への展開形、(2) 各
フィールドの型、(3) 意味論、の 3 点セットで記述し直したものである。
提案の実体は CONVENTIONS.md §7 および design spec §4.6 に既にある
判断の**書き起こし**であり、新たな設計判断を追加するものではない。
1–5 は本カタログの 31 ファイルで実装・リンタ検査済みの確定仕様、
6 のみは素描(実装対象外)であることを明示する。

### 6.1 提案 1 — フェルミオン site_dof

**展開形**(`model.site_dof.<label>`):

```yaml
site_dof:
  <label>: {fermion: {orbitals: <n>}}
```

**型**: `orbitals` は正整数(`n ≥ 1`)。本カタログの全実例は
`orbitals: 1`(単一軌道)であり、複数軌道の実例は未収録
(§7.7 参照)。

**意味論**: ラベル `<label>` の自由度がフェルミオン(生成消滅演算子
`c†, c` を持つ)であり、軌道数が `n` であることを表す。draft 仕様が
既に規定する `spin: {...}` 系(1.3 節(4)、`{param, scale, default}`
一般形を伴う)と対になる site_dof の型の 1 つとして、`fermion` を
draft の site_dof 語彙に追加することを提案する。CONVENTIONS.md §6.6
の演算子型整合表(リンタ C9)は、`fermion` site_dof を持つラベルにのみ
`hop`/`density-density`/`N`/`Nup`/`Ndn`/`NupNdn` を許可し、`spin`
site_dof を持つラベルにのみ J テンソル・onsite スピン演算子を許可する
——すなわち本提案は site_dof の型を演算子の適用可否と結び付ける
機械的検査の前提でもある。Kondo 模型の `<X>_c` ラベル(1.4 節、
CONVENTIONS.md §6.7)と wannier90 の `W1` ラベル(5.4 節)が実例。

### 6.2 提案 2 — 1 サイト演算子語彙と tensor_terms 展開形

**展開形**(`model.onsite.<label>.<term>`、CONVENTIONS.md §6.1 の
onsite 例外形。要素 1 個の `tensor_terms` 配列):

```yaml
<term>:
  operator: {tensor_terms: [{ops: [<OP>], coeff: <±1.0 リテラル>}]}
  value: {param: <名前>, scale: <実数>, default: <数値>}
```

**型**: `<OP>` は以下 8 語の enum。

| `<OP>` | 意味 | 要求される site_dof |
|---|---|---|
| `N` | 数演算子 `n↑+n↓` | fermion |
| `Nup` | 上向きスピン数演算子 `n↑` | fermion |
| `Ndn` | 下向きスピン数演算子 `n↓` | fermion |
| `NupNdn` | 二重占有 `n↑n↓` | fermion |
| `Sx` | スピン x 成分 | spin(または電子スピンとして fermion) |
| `Sy` | スピン y 成分 | spin(または電子スピンとして fermion) |
| `Sz` | スピン z 成分 | spin(または電子スピンとして fermion) |
| `Szz` | 単イオン異方性 `(Sz)^2` | spin(または電子スピンとして fermion) |

**意味論**: `tensor_terms[0].coeff` は符号のみを表すリテラル数値
(`+1.0`/`-1.0`)であり、外部パラメータへの参照は同じ項の `value`
(提案 4、`{param, scale, default}`)が担う——すなわち
`H_onsite = value × coeff × <OP>` という 2 段構成が本提案の核心
である(1.3 節(3)、CONVENTIONS.md §6.2 の符号表がこの 8 語×`±1.0`
の具体例)。この分離により、同じ物理演算子(例: `Sz`)に対して模型
ごとに異なる符号規約(`h`: `-Sz`)を、演算子語彙自体を増やすことなく
`coeff` の値だけで表現できる。

### 6.3 提案 3 — 名前付き 2 体演算子と端点順序の意味論

**展開形**(`model.couplings.<type>.operator`、文字列形):

```yaml
couplings:
  <type>: {operator: "hop" | "density-density" | "s_i . S_j", value: {...}}
```

**型**: `operator` は 3 値 enum の文字列(J 族の `tensor_terms`
辞書形とは排他)。

**意味論**:

| `operator` | 数式 | 端点(site_dof)要件 | 端点順序の意味 |
|---|---|---|---|
| `hop` | `Σ_σ (c†_iσ c_jσ + h.c.)` | fermion – fermion | 対称(`+ h.c.` により順序非依存。ただし複素 `t` では反転時に共役、4.3 節) |
| `density-density` | `n_i n_j` | fermion – fermion | 対称 |
| `s_i . S_j` | `s_i・S_j`(Kondo 結合) | fermion – spin(この順) | **非対称**: 第 1 端点(`from`, i)= 遍歴電子スピン、第 2 端点(`to`, j)= 局在スピン(CONVENTIONS.md §6.6、1.3 節(2)) |

端点順序が意味を持つ(`s_i . S_j`)場合、`model.bonds` の `from`/`to`
順序がそのまま演算子の第 1/第 2 引数に対応する——これは 1.3 節(2)の
ソース順保持規約(参照実装の `general_j(..., isite, jsite)` 引数順を
そのまま書き写す)が、`operator` フィールドの意味論と整合する
ように設計されている点を明示する提案である。`hop`/`density-density`
は演算子自体が対称なため、端点順序はソース順保持のためだけに保たれ、
物理的な非対称性は持たない。

### 6.4 提案 4 — `{param, scale, default}` 参照の許容箇所

**展開形**(4 箇所すべてで共通の一般形):

```yaml
value: {param: <名前>, scale: <実数, 省略時 1.0>, default: <数値, 省略時 0>}
```

**型**: `param` は文字列(参照先パラメータ名)。`scale` は実数。
`default` は数値(整数・実数いずれも可)。

**意味論**: 1.3 節(4)/CONVENTIONS.md §6.3 のとおり、**param が
指定された場合と指定されない場合で計算式が変わる**条件分岐そのもの
が規範である:

- param 指定時: 値 = `scale × param`。
- param 未指定時: 値 = `default`(**scale は適用しない最終値**)。

この一般形は draft 仕様上、以下の**4 箇所**に共通して現れることを
提案する(いずれも本カタログで実例あり):

| 出現箇所 | 例 | 対応する物理量 |
|---|---|---|
| `couplings[<type>].value` | `{param: t0, scale: -1.0}` | ホッピング・Coulomb・Kondo 結合の係数 |
| J テンソル `tensor_terms[*].coeff` | `{param: J0xy}` | 交換相互作用の各成分 |
| `site_dof.<label>.spin` | `{param: 2S, scale: 0.5, default: 0.5}` | スピン量子数 `S`(既定 `S=0.5`) |
| `system.boundary[*].twist` | `{param: phase0}` | 境界位相(度単位) |

draft 仕様は現状これらの各箇所に対して個別の型を規定していないため、
本提案はこの 1 つの一般形を draft の共通プリミティブとして採用する
ことを求めるものである。

### 6.5 提案 5 — `catalog:` ヘッダ

**展開形**(各 YAML ファイル先頭、トップレベルキーの 1 つ):

```yaml
catalog: {schema: stdface-catalog/0.1, dialect: experimental, lattice: <格子名>, model: <模型名>}
```

**型**: 4 フィールドいずれも文字列。`schema` は
`<名前>/<メジャー>.<マイナー>` 形式のバージョン文字列、`dialect` は
`experimental`(現状唯一の値)、`lattice`/`model` は
`manifest.yaml` のエントリキー(`<lattice>/<lattice>_<model>.yaml`)
と突合される(C10、1.4 節)。

**意味論**: ファイル単体を自己記述的にする**メタデータヘッダ**として
draft 仕様のトップレベルキーに `catalog:` を追加することを提案する。
`schema` は将来の破壊的変更に備えたバージョン管理の起点、`dialect`
は draft 本体からの逸脱(本提案 1–4, 6)を宣言する場であり、
`lattice`/`model` は geometry/system/model の三層構造だけからは
機械的に復元できない「このファイルが何の格子・模型のインスタンス
であるか」という索引情報を提供する(ディレクトリ構成・ファイル名
からの推測に依存しないための冗長化)。

### 6.6 提案 6(素描のみ)— 4 フェルミオン一般項

**位置づけ**: 本項目のみ、他の 5 項目と異なり**実装対象外の素描**
である(5.5 節・CONVENTIONS.md §7 項目 6 で予告済み)。

**動機**: 5.5 節で確認したとおり、wannier90 の J チャネル
(`_jr.dat`)は 1 つの係数行列要素 `J_mn` から Hund 密度項・交換項
(`Ex_list`)・ペアホッピング項(`PairHopp_list`)という**独立した
3 つの 4 フェルミオン演算子**を同時に生成する。提案 3(名前付き
2 体演算子)も J テンソル(spin–spin)も「1 ボンドに 1 物理量」を
前提とするため、この構造を表現できない。

**展開形の素描**(未確定。方向性のみ):

```yaml
# 素描 — draft 仕様にも本カタログの拡張方言にも未確定
couplings:
  <type>:
    operator:
      general_terms:
        - ops:
            - {op: "c+", site: i, orbital: m, spin: up}
            - {op: "c",  site: j, orbital: n, spin: up}
            - {op: "c+", site: j, orbital: n, spin: dn}
            - {op: "c",  site: i, orbital: m, spin: dn}
          coeff: {param: J_mn}
```

**型**: 未確定。少なくとも (a) 順序付き生成消滅演算子列
(`c+`/`c` の並び)、(b) 各演算子への site/orbital/spin 束縛、
(c) 演算子列全体に対する係数 1 つ、の 3 要素を持つ必要がある
(4 フェルミオン項なので演算子列の長さは 4)。

**意味論**: 名前付き演算子(`hop`, `density-density`, `s_i . S_j`)
は、この一般形において特定の演算子列パターン(例: `hop` =
`c†_{i↑}c_{j↑} + c†_{i↓}c_{j↓} + h.c.` を `general_terms` で書いた
もの)に展開できる**省略記法**として再定義できる、という位置づけを
提案する。すなわち提案 3 は本提案の特殊ケースになる。

**残された課題**: 演算子列の順序(反交換関係との整合)、site/orbital
束縛の型、複数演算子列の線形結合としての 1 つの `couplings` エントリ
の表現(Hund/交換/ペアホッピングを 1 つの `type` にまとめるか複数
`type` に分けるか)はいずれも未確定であり、draft 仕様側の設計判断を
要する。本カタログは 3 体以上の多体項(現行 StdFace に存在しない、
7.7 節)とは異なり、4 フェルミオン項自体は wannier90 の J チャネルと
いう**既に存在するデータソース**を持つため、上記が確定すれば
`lattice_catalog/wannier90/` に J チャネルの YAML 例を追加できる
見込みである。

---

## 7. 既知の制限・未対応事項

本章は本カタログの検証範囲・記述範囲の外側にある事項を列挙する。
1.1 節で述べたとおり、本カタログの検証水準は「ソース突合 + リンタ +
manifest 検算」であり、数値的な実行同値性検証ではない——この限界を
筆頭に、以下の 10 項目を既知の制限として明記する。

### 7.1 oracle 比較(数値同値性)未実施

本カタログは「YAML の `bonds`/`couplings`/`onsite` が参照実装の
ソースコードと一致すること」(ソース突合)を検証しているが、
「YAML からハミルトニアンを実際に展開し、その数値が StdFace の
生成する `trans.def`/`interall.def` 等を使った計算結果と一致する
こと」(oracle 比較、実行同値性)は検証していない(1.1 節)。この
検証には YAML → ハミルトニアン行列(または `trans.def` 相当の
出力)への**展開エンジン**の実装が前提となり、現時点では未着手
である。展開エンジン実装時の最初の検証課題として引き継ぐ。

### 7.2 ladder は W=2/W=3 の例のみ、および内部 1 次元表現の制限

3.2 節で導出した「W 一般化規則」(`for uc_i in range(W)` ループの
手展開)は任意の `W` に対する `_BONDS` の**生成規則を文書化した
もの**であり、`W=2`/`W=3` 以外の `W` に対する YAML ファイルは
用意していない(検証済みなのはこの 2 例のみ、Task 6)。加えて、
ladder は(chain と同様)参照実装内部で `W` が強制的に 1 に固定される
1 次元表現を取るため、本カタログの `geometry.frac` には脚どうしの
横方向の相対位置が一切現れず、脚は `frac` 座標ではなく**ラベル名
のみ**によって区別される(3.2 節)。これは意図的な設計上の制限で
あり、幾何座標だけを読んでラング方向の構造を復元することはできない。

### 7.3 wannier90 は Hubbard の H/U チャネル例のみ

`lattice_catalog/wannier90/example_hubbard.yaml` は Hubbard 模型の
hopping(H)チャネルと onsite Coulomb(U)チャネルのみを対象とする
(5.2 節・5.4 節)。以下は本カタログの範囲外である:

- **J チャネル**(Hund・交換・ペアホッピングの 4 フェルミオン項、
  5.5 節): 変換規則は散文で記述したのみで、YAML 例は存在しない。
  4 フェルミオン一般項スキーマが未確定であるため(6.6 節)。
- **Spin 模型の wannier90 変換**(5.3 節、超交換 `2|t|²(1/U_m+1/U_n)`
  の自動生成): アルゴリズムは manual 上に記述したが、YAML 例は
  用意していない。
- **Kondo 模型の wannier90 変換**: `wannier90.py::
  _validate_wannier_params` が `ModelType.KONDO` を明示的に拒否する
  (5.1 節)ため、そもそも実装上非対応である。この点は仕様記述
  (5.1 節)のみで、YAML 例が存在しないのは当然の帰結である。

### 7.4 pyrochlore Kondo の `isite+3` 挙動(上流バグ疑い)

3.9 節・CONVENTIONS.md §6.7 で確認したとおり、pyrochlore の Kondo
結合は他の 8 格子と異なり非対称である: 現行実装(C/Python 共通)は
`general_j(..., isite + 3, jsite + uc_i)` という形で**副格子 3 の
遍歴サイト `A3_c` のみ**が全 4 副格子の局在スピンと結合する構造を
生成し、`A0_c`/`A1_c`/`A2_c` は Kondo 結合を一切持たない。本カタログ
は検証水準の方針(現行実装ソースとの突合)に従いこの挙動をそのまま
`pyrochlore_kondo.yaml` に記載している。この非対称性が意図した設計
なのか上流(C 実装)のバグなのかは本カタログの範囲では判定できない。
**上流が将来この挙動を修正した場合、本カタログは `catalog.schema`
のバージョンを上げた上で該当ファイルを更新する必要がある**(現行
`schema: stdface-catalog/0.1` は「現行実装の挙動」を記述している
という前提に立つため)。

### 7.5 fc_ortho の `V''` が `not_used` 拒否リストからも欠落している上流の不整合

3.8 節で詳述したとおり、`fc_ortho` の Spin 分岐は `J''` を「読んで
から捨てる」、Hubbard/Kondo 分岐の `t''` は「読まれもしない」という
挙動を持つのに対し、`V''`(`StdI.Vpp`)は Hubbard/Kondo 分岐から
読まれないことに加え、**Spin 分岐の `not_used_d` 拒否リストからも
`V''` の行が欠落している**(兄弟格子の square/triangular は
`not_used_d("V''", StdI.Vpp)` を明示的に持つのに対し、fc_ortho の
C 実装(`src/FCOrtho.c`)・Python 実装(`fc_ortho.py`)双方でこの
1 行が存在しない)。本カタログはこの現行実装の挙動どおり `V''` を
`bonds`/`couplings` から除外しているが、この欠落自体が上流実装の
意図的な設計なのか単なる記述漏れなのかは不明であり、上流での修正を
待つ必要がある未解決事項として記載する。

### 7.6 box 行列は写像説明のみで YAML 例なし

1.3 節(1)・2.3 節で `box`(`a0w, a0l, ..., a2h` の 9 成分)から
`system.supercell`(`S`、`A_super = S・A`)への写像規則を文書化した
が、本カタログの 31 ファイルはいずれも `supercell` 表現ではなく
`size`(`W/L/Height`)表現のみを採用しており(CONVENTIONS.md §4)、
`system.supercell` を実際に使う YAML 例は 1 件も存在しない。写像
規則自体は検証済みだが、`supercell` フィールドのスキーマ的な実例
(リンタが `supercell` を検査するケース)は本カタログの範囲外である。

### 7.7 機械的識別子(prime を含む type 名)の分離は仕様確定時の課題

`"J0'"`, `"J0''"`, `"t0''"` のような prime を含む `type` 名は、YAML
上引用符付き文字列として扱わざるを得ず(1.4 節)、`'` 自体は YAML の
構文上特別扱いされないが人間の可読性・機械パーサ双方にとって
曖昧さの余地を残す(例: `J0'` と `J0''` を文字列比較以外の方法で
「同じ系列の異なる次数」と機械的に認識する手段が draft 仕様上
規定されていない)。本カタログは現行 StdFace の命名規則
(`t/t'/t''`, `J/J'/J''` 等)をそのまま `type` 名に転写しているが、
`type` 名から「系列」(prefix)と「次数」(prime の数)を分離した
機械可読な識別子体系(例: `{series: J0, order: 2}` のような構造化
表現)を draft 仕様側で規定するかどうかは、本カタログの範囲を超える
仕様確定時の課題として残す。

### 7.8 多体項(3 体以上)は対象外

本カタログが扱う `couplings`/`onsite` はいずれも 1 体項(onsite)
または 2 体項(bonds 上の couplings、6.6 節の素描を含めても
4 フェルミオン=2 体の相互作用)までであり、3 サイト以上にまたがる
真の多体相互作用(3 スピン交換等)は扱っていない。これは現行の
StdFace(C/Python 実装いずれも)がそもそも 3 体以上の相互作用を
生成する機構を持たないためであり、本カタログにとっての制限という
よりは**参照実装の機能範囲そのものの反映**である。将来 StdFace が
多体項生成機能を追加した場合でも、draft 仕様側の site_dof/演算子
機構(1 体・2 体の一般化)はそのまま拡張の土台として使える見込み
である。

### 7.9 YAML・manifest・manual の三重記載と同期コスト

各カタログエントリの情報(ボンド定義・検算値・出典)は、
(a) YAML 本体(`bonds`/`couplings`/`onsite`)、(b)
`manifest.yaml`(`bonds_per_uc`/`coordination`/`min_size_for_check`/
`source`)、(c) manual.md 3 章の解説文、の**3 箇所**に重複して
記載されている。この三重記載は `lint_catalog.py` の manifest 突合
(C10)・計数展開(C11)、および `keyword_inventory.py` による inventory
の機械照合(2 章)によって「(a) と (b) の不一致」は検出できるが、
**(c) manual.md の解説文と (a)/(b) の不一致はリンタでは検出できない**
(manual.md は自然文であり機械検査の対象外)。したがって、いずれかの
格子・模型のソース実装(`python/stdface/lattice/*.py`)が変更された
場合、変更者は (a) YAML、(b) manifest.yaml、(c) manual.md の**3 箇所
すべて**を手動で同期更新する必要がある。この同期漏れを防ぐ機構
(例: manual.md からボンド表を自動生成する等)は本カタログの範囲外
である。

### 7.10 manifest の `source` が単一 file/func 形式であることの限界

`manifest.yaml` の各エントリが持つ `source: {file, func, commit}`
(8 章)は「1 ファイル・1 関数」を前提とした形式である。ほとんどの
格子(chain〜pyrochlore の 9 格子)はこの形式で十分に出典を表現できる
一方、wannier90 の変換ロジックは 5.1 節で確認したとおり
`wannier90.py`(呼出しエントリ・hopping/coulomb 適用)、
`wannier90_io.py`(`_read_w90` によるデータ読込み)、
`interaction_builder.py`(`hopping`/`trans`/`coulomb` という標準
trans/intr 経路)の**3 ファイルにまたがる**導出である。現状は
`manifest.yaml` の単一 `source` エントリに加え、
`example_hubbard.yaml` のヘッダコメント(規則 W1–W6、5.1 節・5.4 節)
に補足の出典・導出をコメント併記することでこの限界を補っている。
複数ファイルに跨る出典を構造化して表現できる `source` の schema
拡張(例: `source: [{file, func, commit}, ...]` のような配列化)は
本カタログの範囲を超える将来課題として残す。
