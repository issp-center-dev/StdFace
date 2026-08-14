# StdFace 格子カタログ(新フォーマット対応一式)設計書

- 日付: 2026-08-15(Codex デザインレビュー Round 1 反映済み)
- ブランチ: `lattice-catalog`
- 参照仕様: 「格子定義仕様 (draft)」2026/07/27
  (PASUMS 共有ドライブ `19. StdFace/1. 打ち合わせ/20260727/stdface_lattice_spec.pdf`)

## 1. 目的

現行 StdFace が定義する全格子・全相互作用タイプを、新しい三層格子定義仕様
(geometry / system / model)へ網羅的に対応させ、

1. 格子×模型ごとの新フォーマット YAML 定義ファイル(静的カタログ)
2. 対応関係を説明する日本語マニュアル
3. カタログの構造的正しさを機械検査する意味検査リンタ + manifest

の一式を作成する。カタログは新仕様 §8「カタログ層」の実体となることを想定した
参照データであり、変換ツール(stan.in → YAML コンバータ)と展開エンジンは
今回のスコープ外。

### 位置づけ(重要)

- 本カタログの YAML は draft 仕様に**拡張方言**(experimental extensions、
  §4.5)を加えた形式で書く。各ファイル先頭で
  `catalog: {schema: stdface-catalog/0.1, dialect: experimental}` を宣言し、
  draft 準拠部分と拡張部分を区別できるようにする。
- 成果物は「展開エンジンで実行検証済み」とは主張しない。主張するのは
  (a) 現行実装ソースとの突合、(b) リンタによる構造検査、(c) manifest に
  基づく検算、の 3 点である。

## 2. スコープ

### 対象

- **格子 9 種(固定幾何)**: chain, ladder(W=2 と W=3 の 2 例), square,
  triangular, honeycomb, kagome, orthorhombic, fc_ortho, pyrochlore
- **wannier90**: 模型別変換仕様の章 + Hubbard 小規模例 YAML 1 件
- **模型 3 種**: Spin / Hubbard / Kondo(GC 変種は粒子数条件の違いのみ。
  粒子数条件はカタログ外(solver 層)に属することをマニュアルに明記)
- **相互作用の語彙**(現行実装から抽出済み):
  - ホッピング: `t, t', t'', t0–t2` とそのプライム変種(複素数)
  - オンサイト斥力: `U`、化学ポテンシャル: `mu`
  - オフサイトクーロン: `V` 族(`V0–V2` とプライム変種)
  - スピン交換: `J` 族(`J0–J2` とプライム変種)、各々 **3×3 テンソル**。
    成分キーワードは `<prefix>{x,y,z,xy,xz,yx,yz,zx,zy}`
    (`x`=xx 対角、`keyword_parser.py::_j_matrix_keywords`)
  - 単イオン異方性 `D`、磁場 `h, Gamma, Gamma_y`
  - 位相(twist)`phase0–2`(**度単位**、境界横断ごとに
    `exp(i·π/180·phase)` を乗算; `site_util.py::ExpPhase`)
  - 一般化セル `box`(a0W 等の行列指定)

- **YAML 件数**: (8 格子 + ladder×2) × 3 模型 = 30 件 + wannier90 1 件 = 31 件

### 対象外

- 展開エンジンの実装、stan.in → 新フォーマット変換ツール
- StdFace 生成の相互作用リストとの数値同値性比較(oracle 比較)。
  リンタは構造検査までとし、同値性検証は展開エンジン実装時の課題とする
- 計算条件系キーワード(`method` 等)→ 対応表で「対象外」と分類
- 新仕様 draft 本体の改訂(拡張は「提案」としてマニュアルに記載するのみ)

## 3. 成果物の構成

```
lattice_catalog/
  README.md                      # 一式の入口
  manual.md                      # 対応マニュアル本体(日本語)
  CONVENTIONS.md                 # 記述規約(Task 0 で確定、manual 1章に転記)
  manifest.yaml                  # 機械可読の検算台帳(全ファイル分)
  tools/
    lint_catalog.py              # 意味検査リンタ(開発時専用、PyYAML 使用可)
    keyword_inventory.py         # パーサーテーブル走査によるキーワード目録生成
  chain/
    chain_spin.yaml
    chain_hubbard.yaml
    chain_kondo.yaml
  ladder/
    ladder_w2_spin.yaml …        # W=2 / W=3 各 3 模型(6 件)
  square/ … triangular/ … honeycomb/ … kagome/ …
  orthorhombic/ … fc_ortho/ … pyrochlore/
  wannier90/
    example_hubbard.yaml         # Hubbard の小規模厳密例
```

## 4. 対応規約

### 4.1 geometry / system

- 格子ベクトル・分率座標・サイトラベルは現行実装から抽出。
- サイトの順序は `geometry.sites` の**配列順**を正とし、仕様 §4.2 の
  「サイト番号」に対応させる(ラベル辞書順ではない)。
- `W, L, Height` → `system.size`。境界は既定 `periodic`。
- `phase0–2` → 方向別 `{twist: θ}`。**単位は度**、境界を n 回横切る
  ボンドには位相 `exp(i·n·π·θ/180)` が乗る(現行実装と同一規約)。
  twist の値も `{param: phase0}` 形式で参照する(コメント扱いにしない)。
- chain は論理 1 次元で記述する。現行実装の内部表現(W=1 の 2D、
  `phase0` 入力が内部で L 方向 `phase[1]` に転写される)との対応は
  マニュアル 3 章に明記する。
- 一般化セル(`box` 行列)は仕様 §5.1 `supercell(S)` への写像として
  マニュアルで説明する。行列の作用方向(`A_super = S·A`、行ベクトル規約)、
  det(S) ≠ 0 条件、`W/L/Height` 指定との排他関係を明記する。
  カタログ YAML 自体は標準セルで記述する。

### 4.2 ボンド(model.bonds)

- **R の定義**: `R = cell(to) − cell(from)`(整数ベクトル)。
  実変位は δ = (frac_to − frac_from) + R·A。
- **ボンド type 名 = StdFace パラメータ名**(`J0, J0', t1, V2, …`)。
  プライムを含む名前は YAML では引用符付き文字列とする。
  機械的識別子との分離は仕様確定時の課題として manual 7 章に記載
  (本カタログでは StdFace 原名を正とする)。
- **向きの一意化**: 各ボンドは一度だけ書く。正準形は
  from のサイト番号(sites 配列順)≦ to、同一サイトなら R 辞書順で正
  (最初の非零成分が正)。
- **反転変換**: `(i, j, R)` を `(j, i, −R)` に反転する場合、
  ホッピングは係数を複素共役、スピン交換は 3×3 テンソルを転置
  (`J_ab → J_ba`)する。この規則を CONVENTIONS.md に明記し、
  リンタで「逆向き重複」を検出する。
- 各格子プラグインの `_BONDS` テーブルから代表ボンドを全数抽出。
  出典は**関数名 + テーブル内容 + 参照コミットハッシュ**で記録する
  (行番号のみの参照は陳腐化するため補助情報とする)。

### 4.3 係数参照とパラメータ解決(拡張方言の核)

- 係数は `{param: <StdFaceキーワード名>}` で参照する。数値は書かない
  (wannier90 例のみ外部データ由来の実数値を書く)。
- **スピン交換の正準形**: J 族はすべて 9 成分の tensor_terms で書く:

  ```yaml
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
  ```

  (tensor_terms の `coeff` に `{param: ...}` を許すのは拡張方言。)
- **パラメータ解決規則**(StdFace の `input_spin_nn` / `input_hopp` /
  `input_coulomb_v` と同一。CONVENTIONS.md に規範として記載):
  1. 成分キーワード(`J0x` 等)が指定されればその値。
  2. 未指定成分は、等方スカラー `J0` が指定されていれば対角に `J0`、
     非対角に 0。
  3. `J0` も未指定なら大域 `J`(および `J'→J0'` 等の別名)へフォールバック。
  4. スカラーと成分の同時指定は StdFace 同様エラー。
  5. どこにも指定がなければ 0。
- ホッピング t 族・クーロン V 族の別名フォールバック(`t→t0` 等)も
  同じ形式で規範化する。

### 4.4 模型ごとの演算子(model.couplings / onsite)

- **演算子の意味論**(CONVENTIONS.md に規範として記載):
  - `hop` = Σσ (c†_iσ c_jσ + h.c.)。**符号は演算子に埋め込まず**、
    ハミルトニアンへの寄与は `−t · hop`(係数規約として文書化)。
    複素 t の h.c. は向き規約 + 反転変換(共役)で処理。
  - `density-density` = n_i n_j。
  - `s_i . S_j` = 遍歴スピン s と局在スピン S の交換(Kondo)。
  - onsite 語彙: `N, Nup, Ndn, NupNdn, Sx, Sy, Sz, Szz(=(Sz)^2)`。
- **Spin**: J 族 → §4.3 正準形。`D` → onsite `Szz`。
  磁場 → onsite `field_z/x/y`(係数 −h, −Gamma, −Gamma_y)。
- **Hubbard**: t 族 → `hop`、`U` → onsite `NupNdn`、`mu` → onsite `−N`、
  V 族 → `density-density`。磁場(h/Gamma/Gamma_y)は電子スピンに適用
  (`model_plugin.py::HubbardModel` と同一)。
- **Kondo**: 遍歴 `<X>_c`(fermion)+ 局在 `<X>_s`(spin、同一分率座標)
  の 2 ラベル。物理的には 1 サイト上の 2 種の自由度であり、ラベルが
  自由度を区別する(仕様 §2.1「サイトラベルは元素種とは別概念」の応用)。
  - Hubbard 項一式(U, mu, t, V, 磁場)は `<X>_c` に適用。
  - **磁場(h/Gamma/Gamma_y)は `<X>_s` にも適用**(両側適用。
    `model_plugin.py::KondoModel` で確認済み)。
  - Kondo 結合 `J` は R=0 の `<X>_c–<X>_s` ボンド。
  - **pyrochlore の例外**: 現行実装(C/Python 共通)は Kondo 結合を
    `isite + 3`(副格子 3 の遍歴サイト)と全副格子の局在スピンの間に
    生成する(`src/Pyrochlore.c:249`)。カタログは**現行動作を忠実に
    再現**し、上流実装のバグの可能性がある旨を YAML コメントと
    マニュアルに注記する。
- **2S の表現**: `site_dof` の spin 値も `{param: 2S}` で参照する
  (拡張方言)。例: `A: {spin: {param: 2S, default: 0.5}}`。
  既定値は StdFace の既定(2S=1 → S=1/2)。

### 4.5 仕様への拡張提案(マニュアル 6 章)

draft 仕様への追記提案として以下を明記(いずれも本カタログで使用):

1. `site_dof` のフェルミオン自由度 `{fermion: {orbitals: n}}`
2. 1 サイト演算子語彙: `Cdag_up, Cdag_dn, C_up, C_dn, N, Nup, Ndn,
   NupNdn, Szz`(tensor_terms への展開形を明記)
3. 名前付き 2 体演算子: `hop`, `density-density`, `s_i . S_j`
4. `{param: ...}` 参照: coupling の value、tensor_terms の coeff、
   site_dof の spin、boundary の twist で許す
5. `catalog:` ヘッダ(schema/dialect 宣言)

## 5. マニュアル構成(manual.md)

1. 概要と読み方(三層仕様の要約、位置づけ、CONVENTIONS の詳説)
2. StdFace キーワード → 新フォーマット対応表。
   **母集合はパーサーレジストリ全体を機械走査して確定**
   (`_COMMON_KEYWORDS` + 全ソルバープラグインのテーブル(`2s` 等)+
   格子名/模型名の別名 + wannier90 系キーワード)。
   `tools/keyword_inventory.py` の出力を基に、
   geometry / system / bonds+couplings / onsite / site_dof /
   wannier90 / 対象外 に分類。
3. 格子ごとの解説(幾何、ボンド定義表、検算、出典(関数名+コミット))
4. 模型ごとの解説(演算子対応、パラメータ解決規則、Kondo 2 ラベルと
   サイト倍加の対応、GC 変種の扱い)
5. wannier90 の変換仕様(**模型別**):
   - Hubbard: H_mn(R) → hop(符号反転あり)、R=0 & m=n は onsite
     一体項として分離、U_mn(R) → density-density、J_mn(R) → Hund・
     exchange・pair-hopping。Hermite 対 `(m,n,R)/(n,m,−R)` の正準対
     選択と二重計数防止、縮退重み、cutoff(値・R 範囲・長さ)、
     `lambda*/alpha/doublecounting` の扱い。
   - Spin: H から超交換 `2|t_mn|^2(1/U_m + 1/U_n)` を生成する現行
     アルゴリズムの記述と、直接写像との差異。
   - 例 YAML は Hubbard 1 件のみ(小規模データ、cutoff 内全要素)。
6. 仕様拡張提案(§4.5)
7. 既知の制限・未対応事項(ladder は W=2/3 の例のみ、oracle 比較未実施、
   機械的識別子の分離は将来課題、等)

## 6. 検証方法

1. **リンタ**(`tools/lint_catalog.py`、各タスクで実行):
   - YAML 構文 + カタログ schema 検証(必須キーの存在・型)
   - R の次元 = geometry.dimension
   - bonds の from/to ラベルと site_dof のラベル整合
   - bonds の type ↔ couplings のキー整合(未定義参照・未使用定義)
   - 向き規約違反・逆向き重複(反転変換を考慮した同値判定)
   - `{param: ...}` 参照名がキーワード目録に存在すること
   - tensor_terms の ops 長(arity)と bond/term のサイト数の整合
   - manifest.yaml との突合: ファイルごとの期待ボンド本数・サイト数・
     配位数(**折り返し縮退のない最小サイズ条件**を manifest に併記)
2. **ソース突合**: `_BONDS` テーブルと 1 行ずつ突合(出典記録付き)。
3. **onsite 項の突合**: `model_plugin.py` の項生成と YAML の onsite を照合。
4. 数値同値性(oracle 比較)は対象外(§2)。

## 7. 進め方

- Task 0(規約・schema・リンタ・manifest 形式・キーワード目録)を先行
  させ、chain 着手前にレビュー可能にする(規約欠陥の 30 ファイル波及を防ぐ)。
- Task 1: chain 3 模型で雛形を確定 → 残り格子へ展開(ladder は W=2/W=3)。
- 後続: wannier90(Hubbard 例)、manual、README。
- コミットは Task 単位。

## 8. 決定事項の記録

初回(ユーザー承認済み):
- 静的カタログ + マニュアル / 全 3 模型 / wannier90 例あり /
  `lattice_catalog/` / 日本語 / 自己完結型ファイル構成 /
  type 名 = StdFace パラメータ名 / Kondo 2 ラベル

Codex デザインレビュー Round 1 対応(ユーザー判断):
- 検証: 意味検査リンタ + manifest(oracle 比較はスコープ外と明記)
- pyrochlore Kondo: 現行動作(`isite+3`)を忠実に再現、バグ疑いを注記
- wannier90: 模型別変換仕様章 + Hubbard 小例 1 件
- ladder: W=2 / W=3 の 2 例 + 一般化規則の文書化
