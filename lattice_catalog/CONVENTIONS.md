# lattice_catalog 記述規約 (stdface-catalog/0.1)

本書は `lattice_catalog/` 配下の YAML カタログが従うべき規約を定める規範文書である。
参照仕様は「格子定義仕様 (draft)」2026/07/27 (以下「draft 仕様」)。
本書と draft 仕様が矛盾する場合、draft 仕様準拠部分は draft を優先し、
draft が触れていない事項(実験的拡張)は本書を規範とする。
設計上の根拠は `docs/superpowers/specs/2026-08-15-lattice-catalog-design.md` §4/§6 を参照。

---

## 1. 位置づけと dialect 宣言

各 YAML ファイルは先頭に必ず以下のヘッダを持つ。

```yaml
catalog:
  schema: stdface-catalog/0.1
  dialect: experimental
  lattice: <格子名>   # 例: chain, square, kagome, ...
  model: <模型名>     # 例: spin, hubbard, kondo
```

- `schema` は本規約のバージョン識別子。固定値 `stdface-catalog/0.1`。
- `dialect` は draft 仕様に対する拡張方言であることを示す固定値
  `experimental`(§7 参照)。
- `lattice` / `model` は manifest 突合(リンタ C10)の照合元となる。
  `lattice_catalog/manifest.yaml` の対応エントリの `lattice` / `model` と
  一致しなければならない。
- draft 仕様準拠部分(geometry / system / model の三層構造そのもの)と、
  本書 §7 に列挙する拡張方言(fermion site_dof、演算子語彙、
  `{param, scale, default}` 参照、`catalog:` ヘッダなど)は区別して扱う。
  draft 未規定の事項は本書の規約が唯一の規範となる。

## 2. ファイル構成

各 YAML ファイルは自己完結型とし、以下の四要素をトップレベルキーに持つ。

```yaml
catalog:  { ... }   # §1
geometry: { ... }   # §3
system:   { ... }   # §4
model:    { ... }   # §5, §6
```

他ファイルへの参照(include 等)は行わない。1 ファイル = 1 格子 × 1 模型。

## 3. geometry 規約

- `geometry.sites` は配列であり、**配列順が draft 仕様 §4.2 の
  「サイト番号」に一致する**。すなわち `sites[0]` が仕様上のサイト 0、
  `sites[1]` がサイト 1、という対応を機械的に維持する。
- サイトラベル(`sites[*].label`)は仕様書やソース実装の慣例名
  (`A`, `B`, `A3` 等)を用い、ファイル内で重複してはならない(C2)。
- `geometry.lattice_vectors`(単位胞の格子ベクトル。次元ぶんの named
  ベクトル `a1`, `a2`, `a3`, ... からなる辞書。§5 の行列 `A` は
  これらを行ベクトルとして積み上げたもの)と各サイトの分率座標
  `frac` を持つ。

## 4. system 規約

- `W, L, Height` などの繰り返し数は `system.size`(整数配列)に写す。
  検査で用いる代表値は `system.size` の各成分そのもの
  (manifest の `min_size_for_check` と対応)。
- `phase0`–`phase2` は `system.boundary`(配列。各要素は
  `{twist: {param: phaseN}}`)に写し、境界を n 回横断する経路の
  位相因子は `exp(i · n · π · θ / 180)` として消費側が解釈する
  (**度単位**)。`{param: phaseN}` 自体は §6.3 の一般形で表現する。
- chain 格子は論理的には 1 次元だが、現行実装は内部的に
  `W=1` の 2 次元表現を用いる。`phase0` → 内部 `phase[1]` への転写は
  実装上の詳細であり、manual 3 章で説明する(YAML の意味論には影響しない)。
- `box`(supercell 変換行列)は `system.supercell`(行ベクトル規約
  `A_super = S · A`、`det(S) ≠ 0`)に対応する。`supercell` と
  `W/L/Height`(= `size`)は排他的に用いる。この対応関係の詳細は
  manual にのみ記載し、本カタログの YAML では `size` 表現を基本とする。

## 5. bonds 規約

- **セル差分 R の定義**: `R = cell(to) − cell(from)`(整数ベクトル、
  次元は `dimension` と一致)。
- **変位 δ の定義**: `δ = (frac_to − frac_from) + R · A`
  (`A` は `geometry.lattice_vectors` の行ベクトル `a1`, `a2`, ... を
  積み上げた行列)。
- **type 名**: StdFace のキーワード名をそのまま用いる
  (`J0`, `J0'`, `t0` 等)。プライムを含む名前は YAML 上
  引用符付きで記述する(`"J0'"`)。プライムの機械的識別子への分離は
  将来課題(manual 7 章)。
- **向き(ソース順保持)**: ボンドの `from` → `to` の向きは
  参照実装のソース順をそのまま保持する。
  - 通常のボンド: `_BONDS` テーブルの `site_i → site_j` の並び。
  - Kondo 結合: `general_j(..., isite, jsite)` の引数順
    (**遍歴サイトが第 1 引数**)。
  正準形への並べ替えは行わない。並べ替えは複素ホッピングの複素共役、
  交換テンソルの転置を要求し、単一の `couplings` キーでは表現できないため
  (design §4.2, Round 2 判断)。
- **一意性判定**: 各ボンドは一度だけ記述する。リンタは**反転同値**
  `(type, i, j, R) ≡ (type, j, i, −R)` での重複を検出する(C6)。
  同一の幾何学的ボンド上に異なる `type`(例: `t0` と `V0`)が
  共存するのは正当であり、`type` を同値キーに含めることでこれを
  区別する。
- **反転時の係数変換(消費側規範)**: `(i,j,R)` を `(j,i,−R)` として
  読み替える場合、
  - hopping(複素数): 複素共役を取る(`t_ij = conj(t_ji)`)。
  - 交換テンソル(J 族 3×3): 転置を取る(`J_ij = J_ji^T`)。
  カタログ自体はソース順で一意に記述するため、この変換は**カタログを
  消費する側**(resolver/展開エンジン)が反転読み替えを行う場合にのみ
  必要になる規範であり、カタログ制作側はソース順のまま書けばよい。

## 6. couplings / onsite 規約

### 6.1 value 意味論

`model.couplings[<type>].value`(スカラー/ベクトル系 couplings。`hop` /
`density-density` / `s_i . S_j`)および
`model.onsite[<site>][<term>].value` は

```
H = Σ_bonds value · operator + Σ_onsite value · operator
```

という**物理ハミルトニアンの係数**である。これは HPhi の
`trans.def` 等ソルバー出力ファイルの係数規約
(`H_trans = −Σ t c†c` のように符号が暗黙に反転している)とは**別物**
であることに注意する。両者の対応は以下の検証連鎖で突合する
(manual 4 章に導出付きで記載):

```
StdFace パラメータ → builder 呼び出し (interaction_builder.py)
→ trans/intr 係数 → solver 出力規約 → 物理ハミルトニアンの符号
```

**J 族 / onsite の例外**: J 族の交換相互作用
(`model.couplings[<type>].operator.tensor_terms`)には共有の `value` は
なく、各成分項が個別に `coeff: {param: ...}` を持つ(§6.4)。onsite の
各項(`model.onsite[<site>][<term>]`)も同じく
`operator: {tensor_terms: [{ops: [...], coeff: <数値>}]}` の形を取るが、
こちらの `tensor_terms` は要素 1 個の配列で、その `coeff` は符号のみを
表す**リテラル数値**(`+1.0` / `-1.0` 等)である。実際の外部パラメータ
参照は同じ項の `value`(`{param, scale, default}`、§6.3)が担う。

### 6.2 符号表(規範)

| StdFace | 物理ハミルトニアン寄与 | YAML 表現 |
|---|---|---|
| t 族(ホッピング) | −t Σ_σ (c†c + h.c.) | `couplings[t0]: {operator: hop, value: {param: t0, scale: -1.0}}` |
| mu(化学ポテンシャル) | −mu N | onsite `operator.tensor_terms: [{ops: [N], coeff: -1.0}]`, `value: {param: mu}` |
| U(オンサイト Coulomb) | +U n↑n↓ | onsite `operator.tensor_terms: [{ops: [NupNdn], coeff: 1.0}]`, `value: {param: U}` |
| V 族(サイト間 Coulomb) | +V n_i n_j | `couplings[V0]: {operator: density-density, value: {param: V0}}`(scale 省略可、既定 +1.0) |
| J 族(交換相互作用) | +Σ_ab J_ab S^a S^b | `couplings[J0]: {operator: {tensor_terms: [...]}}`(§6.4) |
| h / Gamma / Gamma_y(磁場) | −h Sz − Γ Sx − Γy Sy | onsite `operator.tensor_terms: [{ops: [Sz/Sx/Sy], coeff: -1.0}]`, `value: {param: h/Gamma/Gamma_y}` |
| D(単イオン異方性) | +D (Sz)² | onsite `operator.tensor_terms: [{ops: [Szz], coeff: 1.0}]`, `value: {param: D}` |
| Kondo J(s·S 結合) | +J s·S | `couplings[J]: {operator: "s_i . S_j", value: {param: J}}`(scale 省略時 +1.0) |

**符号は必ずデータ(`scale` / `coeff`)に持たせ、コメントとして記述しては
ならない。** これは検証可能性(リンタ・後続処理での機械的突合)を
担保するための必須規約である。

### 6.3 param 参照の一般形

`{param, scale, default}` は以下の一般形を持つ拡張方言(§7)である。

```yaml
value: {param: <名前>, scale: <実数, 省略時 1.0>, default: <数値, 省略時 0>}
```

- **param が(消費側の入力で)指定された場合**: 値 = `scale × param`。
- **param が指定されない場合**: 値 = `default`(**scale は適用しない
  最終値として扱う**)。

すなわち `scale` は「param が与えられたときにだけ効く倍率」であり、
`default` は「param 不在時の最終値」であって `scale` との積ではない。
この条件分岐は分岐そのものが規範であり、C7 検査およびツールテスト
(`test_tools.py`)で固定する。

例: `{param: 2S, scale: 0.5, default: 0.5}`
- `2S=1` が指定された場合 → S = 0.5 × 1 = 0.5。
- `2S` が指定されない場合 → S = 0.5(`default` をそのまま採用、
  `scale` を掛けない)。

型制約: `scale` は実数、`default` は数値(整数・実数いずれも可)。
個別パラメータの型制約は各キーワードの意味論に従う
(例: `2S` は正整数)。

wannier90 の hop チャネルは上記とは別の規則を持つ: `wannier90.py::
_apply_hopping_terms` は非局所項を `hopping(StdI, -Cphase*tUJ[0][it],
jsite, isite, dR)` として渡しており(builder 呼び出し自体が符号反転
を持つ)、これと trans.def の solver 側符号反転(§6.1)が相殺した
結果、**物理ホッピング係数は `+H_mn(R)` そのもの(符号反転なし)**と
なる(`H_mn → −H_mn` という反転は生じない)。検証済みの導出は
§5.5 / manual 5 章に別途明記する(本カタログの `{param, scale,
default}` 一般形の対象外)。

### 6.4 J 族 9 成分 tensor_terms 正準形

J 族の交換相互作用は `model.couplings[<type>].operator.tensor_terms` の
9 成分として表現する(等方成分 `J0` 等の 1 パラメータへの縮約は行わない
— 常に 9 項を書き切る)。各項は 2 スピン演算子の積 `ops: [S?, S?]` と
その係数 `coeff: {param: <name>}` からなる。

ops 対と param 接尾辞の対応表:

| ops 対 | 接尾辞 | 例(prefix=J0) |
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

対角成分(`x`, `y`, `z`)は `ops` の両端が同じ演算子、非対角成分
(`xy` 等)は異なる演算子の順序付き積であることに注意する
(`xy` と `yx` は独立パラメータ)。

### 6.5 パラメータ解決順序(実装準拠)

各成分の値は `input_params.py::_resolve_spin_matrix` と同一の
優先順位で解決される(カタログはこの解決規則を規範として提供し、
実際の解決は**カタログ消費側**(resolver)が行う。カタログの出力契約は
「ボンドごとの解決済み数値係数」):

1. 成分局所(例: `J0xy`)
2. 成分大域(例: `Jxy`)
3. スカラー局所・対角のみ(例: `J0`、`x=y=z` 成分にのみ適用)
4. スカラー大域・対角のみ(例: `J`、`x=y=z` 成分にのみ適用)
5. 0(既定値)

**競合規則**(`input_spin_nn` と同一。以下はすべてエラーとして
検出される組合せ):

- スカラー同士: `J`(大域スカラー)と `J0`(局所スカラー)の同時指定。
- スカラー vs 行列(4 組合せすべて): `J` vs `J`(行列)、
  `J` vs `J0`(行列)、`J0`(スカラー)vs `J`(行列)、
  `J0`(スカラー)vs `J0`(行列)。
- 行列 vs 行列: `J0` の成分と `J` の成分の同時指定。

**プライム系の例外**(`input_spin` と同一): `J0'`, `J0''`, `J1'` 等の
プライム付き系列には**大域 fallback が存在しない**。すなわち
`J0'xy` のスカラー版 `J0'` と成分版 `J0'xy` の競合のみが検査対象で、
対応する大域変数(`J'` 等)は fallback チェーンに現れない。

prefix ごとの解決表(概要):

| prefix 系列 | 局所成分 | 局所スカラー | 大域成分 fallback | 大域スカラー fallback |
|---|---|---|---|---|
| `J0`, `J1`, `J2`(隣接系) | あり | あり | `J` の成分 | `J`(スカラー) |
| `J0'`, `J0''`, `J1'`, ...(プライム系) | あり | あり | なし | なし |
| `J`(大域そのもの) | — | — | — | — |

### 6.6 演算子意味論と端点順序

- `hop` = Σ_σ (c†_iσ c_jσ + h.c.)(符号なし。符号は `value` の
  `scale`/`coeff` に持たせる — §6.2)。
- `density-density` = n_i n_j。
- `s_i . S_j`(Kondo 結合): **第 1 端点(i)が遍歴電子スピン、
  第 2 端点(j)が局在スピン**。端点順序に意味があり、
  §5 のソース順保持規約により順序が保証される。
- onsite 演算子語彙: `N`, `Nup`, `Ndn`, `NupNdn`, `Sx`, `Sy`, `Sz`, `Szz`。
  これは `model.onsite[<site>][<term>].operator.tensor_terms[0].ops` の
  要素(単一演算子。§6.1 の例外参照)として現れる。

演算子と `site_dof` の型整合表(リンタ C9 が検査):

| 演算子 | 要求される site_dof の型(端点順) |
|---|---|
| `hop`(couplings, `operator` 文字列) | fermion – fermion |
| `density-density`(couplings, `operator` 文字列) | fermion – fermion |
| `s_i . S_j`(couplings, `operator` 文字列) | fermion – spin(この順。第 1 端点が fermion) |
| J テンソル(couplings, `operator.tensor_terms`、ops が `S?`) | spin – spin |
| onsite スピン演算子(`Sx`/`Sy`/`Sz`/`Szz`) | spin、または fermion(電子スピンとして) |
| onsite `N`/`Nup`/`Ndn`/`NupNdn` | fermion |

模型別 onsite 適用範囲(Global Constraints の表と同内容):

| 模型 | onsite に現れうる項 |
|---|---|
| Spin | 磁場(`h`, `Gamma`, `Gamma_y`)、単イオン異方性 `D` |
| Hubbard | `mu`(化学ポテンシャル)、`U`、磁場(電子スピンに対する `h`/`Gamma`/`Gamma_y`) |
| Kondo | `_c` ラベルに Hubbard 一式(`mu`, `U`, 磁場)、`_s` ラベルに磁場、両ラベルとも磁場を持つ(§6.7) |

### 6.7 Kondo の 2 ラベル規約

Kondo 模型は 1 物理サイトを **2 つの geometry ラベル**
(`<X>_c` = 遍歴電子の fermion 自由度、`<X>_s` = 局在スピンの
spin 自由度)として表現する。両ラベルは**同一の分率座標**を持つが、
物理的には 1 サイト上の 2 自由度である。

- Hubbard 一式(`t`, `U`, `mu`)は `_c` ラベルにのみ適用。
- 磁場(`h`, `Gamma`, `Gamma_y`)は `_c`, `_s` **両ラベル**に適用
  (§6.2 の符号表通り)。
- Kondo `J`(`s_i . S_j`)は `_c → _s` の順(遍歴が第 1 端点、§6.6)。
- **pyrochlore の例外**: 現行実装(C/Python 共通)は J を
  「副格子 3 の遍歴サイト ↔ 全副格子の局在スピン」という非対称な
  形で生成する(`src/Pyrochlore.c` の `GeneralJ` /
  `python/stdface/lattice/pyrochlore.py::_local`)。本カタログは
  現行実装の動作をそのまま忠実に再現する(`from` は遍歴サイト
  `A3_c` 固定)。これは上流実装のバグである可能性があるが、
  本カタログの検証水準は「現行実装ソースとの突合」であるため、
  上流の挙動を優先する。上流が将来修正された場合はカタログの
  versioning(schema バージョン上げ)が必要になる旨を manual に
  注記する。

## 7. 拡張方言一覧(experimental extensions)

`dialect: experimental` の下で本カタログが用いる、draft 仕様が
明示的には規定していない拡張は以下の通り。draft への提案としても
manual 6 章に記載する。

1. **fermion site_dof**: `{fermion: {orbitals: n}}` 形式でフェルミオン
   自由度を表現する(spin 自由度との区別、§6.6 の型整合表参照)。
2. **1 サイト演算子語彙**: `N`, `Nup`, `Ndn`, `NupNdn`, `Sx`, `Sy`, `Sz`,
   `Szz` とその展開形(§6.6)。
3. **名前付き 2 体演算子**: `hop`, `density-density`, `s_i . S_j` と、
   端点順序に意味を持たせる規約(§6.6)。
4. **param 参照** `{param, scale, default}`(§6.3)。`value` / `coeff` /
   `spin`(`2S`)/ `twist`(`phaseN`)のいずれにも同一の一般形を用いる。
5. **`catalog:` ヘッダ**(`schema` / `dialect` / `lattice` / `model`、§1)。
6. **(素描のみ・将来課題)4 フェルミオン一般項**: 順序付き生成消滅
   演算子列と site/orbital 束縛を持つ項表現。名前付き演算子
   (`hop` 等)はその省略記法という位置づけになる。wannier90 の
   J チャネル(Hund・exchange・pair-hopping)の記述に将来必要となるが、
   今回のカタログでは対象外(design §2 対象外, §5.5)。

## 8. 検算・出典記録の書式

- 各 YAML エントリの検算値(`n_sites_uc`, `bonds_per_uc`,
  `coordination`, `min_size_for_check` 等)は
  `lattice_catalog/manifest.yaml` に記録し、リンタ(C10, C11)が
  自動突合する。
- 出典は「ファイルパス + 関数名 + 参照コミットハッシュ」の組で記録する
  (行番号は補助情報として付記してもよいが、コミットハッシュが
  一意な参照点となる)。

```yaml
source:
  file: python/stdface/lattice/chain_lattice.py
  func: chain
  commit: <hash>
```

- manifest 側のフォーマット定義とキー一覧は `manifest.yaml` 冒頭の
  コメントを参照。
