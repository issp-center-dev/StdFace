# StdFace 格子カタログ(新フォーマット対応一式)設計書

- 日付: 2026-08-15
- ブランチ: `lattice-catalog`
- 参照仕様: 「格子定義仕様 (draft)」2026/07/27
  (PASUMS 共有ドライブ `19. StdFace/1. 打ち合わせ/20260727/stdface_lattice_spec.pdf`)

## 1. 目的

現行 StdFace が定義する全格子・全相互作用タイプを、新しい三層格子定義仕様
(geometry / system / model)へ網羅的に対応させ、

1. 格子×模型ごとの新フォーマット YAML 定義ファイル(静的カタログ)
2. 対応関係を説明する日本語マニュアル

の一式を作成する。カタログは新仕様 §8「カタログ層」の実体となることを想定した
参照データであり、変換ツール(stan.in → YAML コンバータ)は今回のスコープ外。

## 2. スコープ

### 対象

- **格子 9 種(固定幾何)**: chain, ladder, square, triangular, honeycomb,
  kagome, orthorhombic, fc_ortho, pyrochlore
- **wannier90**: 写像方針の章 + 小規模な具体例 YAML 1 件
- **模型 3 種**: Spin / Hubbard / Kondo(GC 変種は粒子数条件の違いのみで
  格子定義としては同一 → マニュアルで注記)
- **相互作用の語彙**(現行実装から抽出済み):
  - ホッピング: `t, t', t'', t0–t2` とそのプライム変種(複素数)
  - オンサイト斥力: `U`、化学ポテンシャル: `mu`
  - オフサイトクーロン: `V` 族(`V0–V2` とプライム変種)
  - スピン交換: `J` 族(`J0–J2` とプライム変種)、各々 **3×3 テンソル**
    (`jx, jy, jz, jxy, …` 成分キーワード)
  - 単イオン異方性 `D`、磁場 `h, Gamma, Gamma_y`
  - 位相(twist)`phase0–2`、一般化セル `box`(a0W 等の行列指定)

各格子が実際に受け付けるパラメータ集合は実装
(`python/stdface/lattice/*.py`、必要に応じ `src/*.c`)から抽出する。
確認済みの概要(属性アクセスベース、実装ステップで精査):

| 格子 | 特記事項 |
|---|---|
| chain | t/J/V の 0–2 + p/pp まで広く使用 |
| ladder | 脚方向 `J1p, J2p` 等、梯子特有の割当 |
| square, triangular, honeycomb | pp(第三近接)まで使用 |
| kagome | p まで(pp なし) |
| orthorhombic, fc_ortho, pyrochlore | 3D、p/pp 一部使用 |

### 対象外

- 展開エンジン(YAML → サイト・ボンドリスト)の実装
- stan.in → 新フォーマットの変換ツール
- 計算条件系キーワード(`method`, `Lanczos_max` 等、格子定義でないもの)
  → マニュアルの対応表に「対象外」と明記する
- 新仕様 draft 本体の改訂(拡張は「提案」としてマニュアルに記載するのみ)

## 3. 成果物の構成

```
lattice_catalog/
  README.md                      # 一式の入口(マニュアルへの案内)
  manual.md                      # 対応マニュアル本体(日本語)
  chain/
    chain_spin.yaml
    chain_hubbard.yaml
    chain_kondo.yaml
  square/ …                      # 9 格子 × 3 模型 = 27 ファイル
  triangular/ …
  honeycomb/ …
  kagome/ …
  ladder/ …
  orthorhombic/ …
  fc_ortho/ …
  pyrochlore/ …
  wannier90/
    example.yaml                 # wannier90 入力の写像例
```

各 YAML は仕様 §7.3(Shastry-Sutherland 例)と同形式の
`geometry / system / model` 完結の「動く見本」。

## 4. 対応規約

### 4.1 geometry / system

- 格子ベクトル・分率座標・サイトラベルは現行実装から抽出。
  ラベル慣例: honeycomb A/B、kagome A/B/C、pyrochlore 4 副格子、
  ladder は leg 別ラベル。
- `W, L, Height` → `system.size`。境界は既定 `periodic`。
  `phase0–2` → 方向別 `{twist: θ}`。
- 一般化セル(`box` 行列、`wx/wy/…` 直接指定)は仕様 §5.1 `supercell(S)`
  への写像としてマニュアルで説明。カタログ YAML 自体は標準セル
  (`size: [W, L(, H)]`)で記述する。

### 4.2 ボンド(model.bonds)

- **ボンド type 名 = StdFace パラメータ名**(`J0, J0', J0'', t1, V2, …`)。
- 各格子プラグインのボンド生成ループから代表ボンド (from, to, R) を
  **全数抽出**し、仕様 §4.2 の向きの一意化規約
  (from のサイト番号 → to、同じなら R 辞書順)に従い記載。
- 係数は `value: {param: <StdFaceキーワード名>}` 形式で参照し、
  StdFace 入力との対応が一目で分かるようにする。

### 4.3 模型ごとの演算子(model.couplings / onsite)

- **Spin**:
  - J 族 3×3 テンソル → `tensor_terms`(対角 `[Sx,Sx]` 等 + 非対角
    `[Sx,Sy]` 等)。`jx/jxy` 等の成分キーワードとの対応表を付す。
  - `D` → `onsite` の `Szz`(複合シンボル)。
  - `h/Gamma/Gamma_y` → `onsite` の `field_z/field_x/field_y`。
- **Hubbard**:
  - t 族 → ホッピング(複素対応。エルミート共役は向き規約で暗黙に処理)。
  - `U` → onsite `n↑n↓`、`mu` → onsite `N`、V 族 → ボンド `n_i n_j`。
- **Kondo**:
  - 各サイトを遍歴電子ラベル(例 `A_c`)+ 局在スピンラベル(例 `A_s`、
    同一分率座標)の 2 ラベルで表現。
  - Kondo 結合 `J` は R=0 の `A_c–A_s` ボンド(`s_i · S_j` 型演算子)。
  - 現行実装のサイト倍加(`jsite_kondo = isite + nsite/2`)の自然な写像
    であることをマニュアルで説明。

### 4.4 仕様への拡張提案(マニュアル 1 章として)

draft 仕様の演算子語彙はスピン想定のため、以下を**追記提案**として明記:

- `site_dof` へのフェルミオン自由度: `A: {fermion: {orbitals: 1}}`
  (Kondo の局在側は従来どおり `{spin: 0.5}`)。
- 1 サイト演算子語彙の追加: `Cdag_up, Cdag_dn, C_up, C_dn, N, Nup, Ndn`。
- 頻出 2 体項の名前付き短縮形: `hop`(= Σσ c†_iσ c_jσ + h.c.)、
  `density-density`(= n_i n_j)、`s_i . S_j`(Kondo 結合)等。
  いずれも `tensor_terms`(1 サイト演算子のテンソル積の線形結合)へ
  機械的に展開できることを示す。

## 5. マニュアル構成(manual.md)

1. 概要と読み方(三層仕様の要約、カタログの使い方)
2. StdFace キーワード → 新フォーマット対応表(全キーワード網羅。
   計算条件系は「対象外」と明記)
3. 格子ごとの解説(幾何、ボンド定義、検算、実装ソースとの突合記録)
4. 模型ごとの解説(演算子対応: Spin / Hubbard / Kondo)
5. wannier90 の写像(H_mn(R)・U_mn(R)・J_mn(R) → bonds/couplings、
   cutoff 系キーワードの扱い)
6. 仕様拡張提案(§4.4 の内容)
7. 既知の制限・未対応事項

## 6. 検証方法

- 各 YAML に仕様 §7.3 流の**検算**(配位数・単位胞あたりボンド本数・
  対称性チェック)をコメントで付記。
- ボンド定義の正しさは現行 Python 実装のボンド生成コード(必要に応じ
  C 実装)との突合で担保し、突合結果(対応する関数・ループ)を
  マニュアル 3 章に記録。
- 展開エンジンは未実装のため実行検証はスコープ外。ただし YAML の
  構文妥当性は `python3 -c "import yaml; yaml.safe_load(...)"` 相当で
  機械チェックする(PyYAML はリポジトリ依存に追加しない。
  チェックは開発時のみの補助とする)。

## 7. 進め方

- 実装計画(writing-plans)で格子ごとに 1 ステップに分割:
  各ステップ = ボンド抽出 → YAML 3 件(spin/hubbard/kondo)作成 → 検算。
- 先行ステップ: カタログ規約の雛形(chain で 3 模型分を最初に作り、
  規約を固めてから残り 8 格子へ展開)。
- 後続ステップ: wannier90 例、manual.md、README.md。
- コミットは格子単位。マニュアルは章単位。

## 8. 決定事項の記録(ユーザー承認済み)

- 成果物: 静的カタログ + マニュアル(コンバータは作らない)
- 模型範囲: Spin + Hubbard + Kondo 全部(フェルミオン拡張は提案として記載)
- wannier90: 例ファイルまで作る
- 置き場所: リポジトリ内 `lattice_catalog/`
- 言語: 日本語
- ファイル構成: 自己完結型(格子×模型で 1 ファイル)
- ボンド type 名: StdFace パラメータ名をそのまま使用
- Kondo: 2 ラベル(遍歴 + 局在)表現
