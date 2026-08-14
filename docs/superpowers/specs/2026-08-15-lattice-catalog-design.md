# StdFace 格子カタログ(新フォーマット対応一式)設計書

- 日付: 2026-08-15(Codex デザインレビュー Round 1・Round 2 反映済み)
- ブランチ: `lattice-catalog`
- 参照仕様: 「格子定義仕様 (draft)」2026/07/27
  (PASUMS 共有ドライブ `19. StdFace/1. 打ち合わせ/20260727/stdface_lattice_spec.pdf`)

## 1. 目的

現行 StdFace が定義する全格子・全相互作用タイプを、新しい三層格子定義仕様
(geometry / system / model)へ網羅的に対応させ、

1. 格子×模型ごとの新フォーマット YAML 定義ファイル(静的カタログ)
2. 対応関係を説明する日本語マニュアル
3. カタログの構造的正しさを機械検査する意味検査リンタ + manifest

の一式を作成する。変換ツール(stan.in → YAML)と展開エンジンはスコープ外。

### 位置づけ

- YAML は draft 仕様に**拡張方言**(experimental extensions、§4.6)を加えた
  形式。各ファイル先頭で
  `catalog: {schema: stdface-catalog/0.1, dialect: experimental}` を宣言。
- 主張する検証水準: (a) 現行実装ソースとの突合、(b) リンタによる構造・
  整合検査(計数展開による配位数実測を含む)、(c) manifest 検算。
  数値同値性(oracle 比較)は主張しない。

## 2. スコープ

### 対象

- **格子 9 種**: chain, ladder(W=2 / W=3 の 2 例), square, triangular,
  honeycomb, kagome, orthorhombic, fc_ortho, pyrochlore
- **wannier90**: 模型別変換仕様の章 + **H(hopping)・U(Coulomb)チャネル
  のみの Hubbard 小規模例** 1 件(J チャネルは §5.5 参照)
- **模型 3 種**: Spin / Hubbard / Kondo(GC 変種は粒子数条件のみ →
  solver 層の扱いとしてマニュアルに明記)
- **相互作用の語彙**(現行実装から抽出・検証済み):
  - ホッピング t 族(複素)、`U`、`mu`、V 族、J 族(3×3 テンソル、成分
    キーワード `<prefix>{x,y,z,xy,xz,yx,yz,zx,zy}`、`x`=xx 対角)
  - `D`、`h, Gamma, Gamma_y`、`phase0–2`(度)、`box`
- **YAML 件数**: (8 格子 + ladder×2) × 3 模型 = 30 件 + wannier90 1 件 = 31 件

### 対象外

- 展開エンジン、stan.in 変換ツール、数値 oracle 比較
- wannier90 J チャネル(Hund・exchange・pair-hopping)の YAML 例
  (4 フェルミオン項のデータモデルが未定義のため。§5.5 で規則を散文
  記述し、一般項スキーマの素描を将来課題として示す)
- 計算条件系キーワード → 対応表で「対象外」と分類
- draft 本体の改訂(拡張は提案として記載)

## 3. 成果物の構成

```
lattice_catalog/
  README.md / manual.md / CONVENTIONS.md
  manifest.yaml                  # 機械可読の検算台帳
  tools/
    lint_catalog.py              # 意味検査リンタ(計数展開器を含む)
    keyword_inventory.py         # キーワード目録生成(レジストリ走査)
    test_tools.py                # tools の自動テスト(開発時専用)
  chain/ ladder/ square/ triangular/ honeycomb/ kagome/
  orthorhombic/ fc_ortho/ pyrochlore/
  wannier90/example_hubbard.yaml
```

## 4. 対応規約

### 4.1 geometry / system

- 幾何は現行実装から抽出。サイト順序は `geometry.sites` の**配列順** =
  仕様 §4.2 の「サイト番号」。
- `W, L, Height` → `system.size`。`phase0–2` → `{twist: {param: phaseN}}`
  (**度単位**、境界 n 回横断で `exp(i·n·π·θ/180)`)。
- chain は論理 1 次元(内部 W=1 2D 表現、`phase0`→`phase[1]` 転写は
  manual 3 章)。
- `box` → `supercell(S)` 写像(`A_super = S·A` 行ベクトル規約、
  det(S)≠0、`W/L/Height` との排他)はマニュアル説明のみ。

### 4.2 ボンド(model.bonds)— ソース順保持

- **R の定義**: `R = cell(to) − cell(from)`。δ = (frac_to − frac_from) + R·A。
- **向き**: **参照実装のソース順をそのまま保持する**
  (`_BONDS` テーブルの `site_i → site_j`、Kondo 結合は
  `general_j(..., isite, jsite)` の引数順 = 遍歴が第 1)。
  正準形への並べ替えは**行わない**(並べ替えは複素 hopping の共役・
  交換テンソルの転置を要求し、単一の couplings キーでは表現できない
  ため。Round 2 指摘)。
- **一意性**: 各ボンドは一度だけ書く。リンタが**反転同値**
  `(i,j,R) ≡ (j,i,−R)` での重複を検出する(反転時の係数変換
  — hopping: 複素共役、交換テンソル: 転置 — は消費側の規範として
  CONVENTIONS に記載)。
- type 名 = StdFace キーワード名(プライムは引用符付き)。
  機械的識別子の分離は将来課題(manual 7 章)。
- 出典: ファイル + 関数名 + 参照コミットハッシュ(行番号は補助)。

### 4.3 係数の意味論と符号(規範)

- **value 意味論**: `H = Σ_bonds value·(operator) + Σ_onsite value·(operator)`
  の**物理ハミルトニアン係数**。solver 出力ファイル(HPhi trans.def 等)の
  係数規約(H_trans = −Σ t c†c の暗黙符号)とは**別物**であることを
  manual に明記し、以下の**検証連鎖**で突合する:
  `StdFace パラメータ → builder 呼び出し(interaction_builder.py)
  → trans/intr 係数 → solver 規約 → 物理符号`。
- **符号表**(検証済み。manual 4 章に導出つきで記載):

  | StdFace | 物理ハミルトニアン寄与 | YAML 表現 |
  |---|---|---|
  | t 族 | −t Σσ (c†c + h.c.) | `value: {param: t0, scale: -1.0}` |
  | mu | −mu N | onsite `coeff: -1.0`(ops [N]) |
  | U | +U n↑n↓ | `coeff: +1.0` |
  | V 族 | +V n_i n_j | `scale: +1.0`(省略可) |
  | J 族 | +Σ_ab J_ab S^a S^b | tensor_terms(§4.4) |
  | h/Gamma/Gamma_y | −h Sz −Γ Sx −Γy Sy | `coeff: -1.0` |
  | D | +D (Sz)² | `coeff: +1.0` |
  | Kondo J | +J s·S | `scale: +1.0` |

- **符号は必ずデータ(scale / coeff)に持たせ、コメントに置かない。**
- **param 参照の一般形**(拡張方言):
  `{param: <名>, scale: <実数, 省略時 1.0>, default: <最終値, 省略時 0>}`
  → 値 = scale × param(param 未指定時は default)。
- wannier90 の `H_mn → −H_mn` 反転は別規則として §5.5 / manual 5 章に明記。

### 4.4 J テンソルとパラメータ解決(規範)

- J 族は 9 成分 tensor_terms 正準形(`coeff: {param: J0xy}` 等)。
- **解決順序**(`input_params.py::_resolve_spin_matrix` と同一):
  1. 局所成分(`J0xy` 等)
  2. **大域成分**(`Jxy` 等)
  3. 局所スカラー `J0`(対角のみ)
  4. 大域スカラー `J`(対角のみ)
  5. 0
- **競合規則**(`input_spin_nn` と同一): スカラー同士(J vs J0)、
  スカラー vs 行列(全 4 組合せ)、行列 vs 行列(J0 成分と J 成分の
  同時指定)はエラー。プライム系(`input_spin`)は大域 fallback なし
  (`J0'` スカラー vs `J0'` 成分のみ)。prefix ごとの解決表・競合表を
  CONVENTIONS に規範として列挙する。
- 解決の実行者は**カタログ消費側**(resolver)。カタログは規範表を
  提供する(入力契約)。出力契約は「ボンドごとの解決済み数値係数」。

### 4.5 模型ごとの演算子

- 演算子語彙: `hop` = Σσ (c†_iσ c_jσ + h.c.)(符号なし)、
  `density-density` = n_i n_j、`s_i . S_j`(**第 1 端点 = 遍歴**。
  端点順序に意味があり、ソース順保持(§4.2)により保証)、
  onsite: `N, Nup, Ndn, NupNdn, Sx, Sy, Sz, Szz`。
- 演算子と site_dof の型整合(リンタ検査): `hop`/`density-density` は
  fermion–fermion、`s_i . S_j` は fermion–spin(この順)、J テンソルは
  spin–spin、onsite スピン演算子は spin または fermion(電子スピン)。
- **Spin**: J 族 §4.4、D、磁場(符号表)。
- **Hubbard**: t/U/mu/V + 磁場(電子スピン)。
- **Kondo**: `<X>_c`(fermion)+ `<X>_s`(spin、同一分率座標、
  1 物理サイト 2 自由度)。Hubbard 一式は `_c`、磁場は**両側**、
  J は `_c → _s`(遍歴第 1)。
  - **pyrochlore の例外**: 現行実装(C/Python 共通)は J を
    「副格子 3 の遍歴サイト ↔ 全副格子の局在スピン」に生成
    (`src/Pyrochlore.c` GeneralJ / `pyrochlore.py::_local`)。
    カタログは現行動作を忠実に再現し(from は遍歴 `A3_c`)、
    上流バグの可能性と、上流修正時のカタログ versioning が必要に
    なる旨を注記する。
- **2S**: `spin: {param: 2S, scale: 0.5, default: 0.5}`
  (spin 値 = S = 0.5 × 2S。param 未指定時 S=0.5。型: 2S は正整数)。

### 4.6 仕様への拡張提案(マニュアル 6 章)

1. `site_dof` フェルミオン `{fermion: {orbitals: n}}`
2. 1 サイト演算子語彙(N 系・S 系・Szz)と展開形
3. 名前付き 2 体演算子(`hop`, `density-density`, `s_i . S_j`)と
   端点順序の意味論
4. param 参照 `{param, scale, default}`(value / coeff / spin / twist)
5. `catalog:` ヘッダ(schema/dialect)
6. (素描のみ)4 フェルミオン一般項: 順序付き生成消滅演算子列と
   site/orbital 束縛を持つ項表現。named operator はその省略記法。
   wannier90 J チャネルの記述に将来必要(今回は対象外)。

## 5. マニュアル構成(manual.md)

1. 概要と読み方(位置づけ、CONVENTIONS 詳説、param 参照と解決規則)
2. キーワード対応表(母集合 = レジストリ走査による機械的 inventory。
   canonical キーワード単位に集約し出典を配列で保持。格子名・模型名の
   alias は registry から別表で生成)
3. 格子ごとの解説(幾何、ボンド表、manifest 検算の根拠、出典)
4. 模型ごとの解説(**符号表と検証連鎖の導出**、解決規則、Kondo 対応、
   GC 変種)
5. wannier90 変換仕様(模型別):
   - Hubbard: H_mn(R) の符号反転、R=0 対角の onsite 分離、Hermite
     正準対と縮退重み、cutoff、λ/α、doublecounting
   - Spin: 超交換 `2|t|²(1/U_m+1/U_n)` 生成
   - **5.5 J チャネル**: Hund・exchange・pair-hopping の変換規則を散文で
     記述し、YAML 例は対象外である理由(4 フェルミオン項スキーマ未定義)
     と §4.6-6 の素描を示す
6. 仕様拡張提案(§4.6)
7. 既知の制限(oracle 比較なし、ladder W=2/3 のみ、wannier90 J チャネル
   例なし、pyrochlore Kondo 現行動作固定と versioning、機械的識別子)

## 6. 検証方法

### 6.1 リンタ(`tools/lint_catalog.py`)

構造検査(C1–C9)+ manifest 突合(C10)+ **計数展開**(C11):

- C1 catalog ヘッダ(schema **と dialect** の値)
- C2 schema 検査(null 文書、必須キー、型、サイトラベル重複、
  dimension/R/size の整数・次元整合)
- C3 R 次元 = dimension = size 次元
- C4 ラベル整合(bonds の from/to、**onsite のラベル**、
  geometry.sites と site_dof の集合一致)
- C5 type ↔ couplings 整合(未定義参照・未使用定義)
- C6 反転同値 `(i,j,R)≡(j,i,−R)` での重複検出
- C7 param 参照の目録整合(scale/default の型検査を含む)
- C8 J coupling の 9 成分完全性と **ops 対 ↔ param 接尾辞の対応**
  (`[Sx,Sy] ↔ …xy` 等、過不足・重複なし)
- C9 演算子と site_dof の型整合(§4.5 の表)、tensor_terms の ops 長
- C10 manifest 全項目突合(lattice/model/dimension/n_sites_uc/
  bonds_per_uc/source の存在)
- C11 **計数展開**: min_size トーラス上に bonds を展開し、
  サイトラベル×type ごとの配位数を実測して manifest の
  `coordination: {label: {type: count}}` と比較(ladder W=3 の
  非一様配位、Kondo 重複座標もラベル単位で扱える)
- 実装要件: 不正 YAML でも例外で落ちずファイル単位で診断、
  相対パス入力の解決、`test_tools.py` に C1–C11 の正例・負例テスト

### 6.2 その他

- ソース突合(`_BONDS`・`model_plugin.py`・`input_params.py`)
- `keyword_inventory.py` はソルバー**レジストリ経由**で列挙
  (クラスのハードコード禁止)、canonical 集約、格子/模型 alias 出力
- 数値 oracle 比較は対象外(§2)

## 7. 進め方

Task 0(規約・schema・リンタ・manifest・目録 + tools 自動テスト)→
Task 1 chain 雛形 → 各格子 → wannier90 → manual。コミットは Task 単位。

## 8. 決定事項の記録

初回(ユーザー承認済み): 静的カタログ+マニュアル / 全 3 模型 /
wannier90 例あり / `lattice_catalog/` / 日本語 / 自己完結型 /
type 名 = StdFace 名 / Kondo 2 ラベル

Round 1 対応(ユーザー判断): リンタ+manifest / pyrochlore 現行動作忠実 /
wannier90 模型別仕様+Hubbard 例 / ladder W=2・W=3

Round 2 対応(設計判断、検証に基づく):
- 係数 = 物理ハミルトニアン規約と確定(Codex の符号指摘は HPhi の
  H_trans = −Σt 規約の見落としで、物理符号としては −t/−mu/−h が正。
  ただし意味論が未規範だった点は正当 → 符号をデータに持たせ scale を導入)
- ボンドはソース順保持(正準化強制を撤回)
- 解決順序を実装通り(成分局所 > 成分大域 > スカラー局所 > スカラー大域)に修正
- 2S は scale 付き param 参照
- wannier90 例は H+U チャネルのみ(J チャネルはデータモデル未定義のため
  散文記述+将来スキーマ素描)
- manifest coordination をラベル×type に、リンタに計数展開(C11)を追加
