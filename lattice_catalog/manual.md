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
| `W, L, Height` | common | `system.size`(整数配列) | 検査で用いる代表値。manifest の `min_size_for_check` と対応 |
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
