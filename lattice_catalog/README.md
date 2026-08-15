# StdFace 格子カタログ

`lattice_catalog/` は、StdFace が組み込みでサポートする格子(chain,
ladder, square, triangular, honeycomb, kagome, orthorhombic, fc_ortho,
pyrochlore)と模型(Spin, Hubbard, Kondo)の全組み合わせ、および
wannier90 変換の一例を、参照仕様「格子定義仕様 (draft)」2026/07/27
(以下「draft 仕様」)が定める geometry / system / model の三層構造
YAML として書き起こしたものである。各 YAML は StdFace の C/Python
実装(`_BONDS` テーブル・演算子生成ロジック・パラメータ解決規則)を
直接読んで書き写した**ソース突合済みの記述**であり、`tools/
lint_catalog.py` によるリンタ検査(スキーマ・ボンド重複・param 参照
整合・演算子型整合・`manifest.yaml` との計数突合)を通過している。
draft 仕様が未規定の事項(fermion site_dof、1 サイト演算子語彙、
`{param, scale, default}` 参照、`catalog:` ヘッダ等)は
`dialect: experimental` として明示された実験的拡張方言で補っている。
ただし本カタログが保証するのは「参照実装ソースとの静的な一致」で
あり、生成したハミルトニアンを実際に数値計算した結果が StdFace の
出力と一致するという実行同値性(oracle 比較)までは検証していない
(詳細は `manual.md` 1.1 節・7 章)。

## ディレクトリ構成

```
lattice_catalog/
  README.md                 # 本ファイル
  CONVENTIONS.md             # 記述規約(規範文書)
  manual.md                  # 解説書(読み方・実例・全キーワード対応表)
  manifest.yaml               # 検算台帳(bonds_per_uc/coordination 等、機械可読)
  tools/
    lint_catalog.py           # 意味検査リンタ(C1–C12)
    keyword_inventory.py      # StdFace 全キーワードの目録生成
    test_tools.py              # tools 自体の自動テスト(開発時専用)
  chain/        chain_{spin,hubbard,kondo}.yaml
  ladder/       ladder_w2_{spin,hubbard,kondo}.yaml, ladder_w3_{spin,hubbard,kondo}.yaml
  square/       square_{spin,hubbard,kondo}.yaml
  triangular/   triangular_{spin,hubbard,kondo}.yaml
  honeycomb/    honeycomb_{spin,hubbard,kondo}.yaml
  kagome/       kagome_{spin,hubbard,kondo}.yaml
  orthorhombic/ orthorhombic_{spin,hubbard,kondo}.yaml
  fc_ortho/     fc_ortho_{spin,hubbard,kondo}.yaml
  pyrochlore/   pyrochlore_{spin,hubbard,kondo}.yaml
  wannier90/    example_hubbard.yaml
```

格子 9 種 × 模型 3 種で 30 ファイル、これに wannier90 の Hubbard 例
1 ファイルを加えた合計 **31 ファイル**が `lint_catalog.py` の検査
対象である(`manifest.yaml` はデータファイルであり検査対象ファイル数
には含まない)。

## リンタの使い方

開発時依存: 本リンタは PyYAML (`pyyaml`) を必要とする(開発時専用
ツールであり、`python/pyproject.toml` の実行時依存には含めていない。
未インストールの場合は `pip install pyyaml` を促す明確なメッセージで
終了する)。

リポジトリルートから実行します:

```bash
python3 lattice_catalog/tools/lint_catalog.py
```

全 31 ファイルに対して、スキーマ検査(ジオメトリ/サイト/ボンド/
site_dof の深部型検査を含む)・ラベル整合・ボンド反転同値の重複検出・
param 参照整合・J テンソル成分完全性・演算子と site_dof の型整合・
符号規約(CONVENTIONS.md §6.2)、および `manifest.yaml` との突合
(期待される `bonds_per_uc`/`coordination` を実際に指定サイズの
トーラス上にボンドを展開して実測比較する「計数展開」検査、および
`min_size_for_check` 自体の厳密性検査を含む)を行う。エラーがなければ

```
31 files, 0 errors
```

と出力される。個々の検査項目(C1–C12)の詳細は `CONVENTIONS.md` と
`manual.md` 1.1 節を参照。

`tools/` には他に、StdFace 全ソルバー(HPhi/HWAVE/UHF/mVMC)の
キーワードを走査して一覧化する `keyword_inventory.py`、および
リンタ・inventory ツール自体の振る舞いを検査する `test_tools.py`
(開発時専用、pytest 不使用の自前テストランナー)がある。

## manual.md への案内

本カタログの読み方・実例・全キーワードとの対応表は `manual.md` に
まとめている。構成は以下のとおり:

1. **概要と読み方** — 三層仕様の要約、`dialect: experimental` の
   位置づけ、本カタログの検証水準(ソース突合・リンタ・manifest
   検算であって数値 oracle 比較ではない点)、ディレクトリ構成、
   CONVENTIONS.md の主要規約の補足説明。
2. **キーワード対応表** — `keyword_inventory.py` が報告する
   StdFace 全キーワード(351 件: keyword 313 + lattice_alias 26 +
   model_alias 12)を、geometry / system / bonds+couplings / onsite /
   site_dof / wannier90 / 対象外の 7 分類に**全数**割り当てた表。
3. **格子ごとの解説** — 9 格子それぞれの幾何・ボンド定義表・
   manifest 検算根拠・出典。
4. **模型ごとの演算子対応** — 符号表の導出(検証連鎖)、Spin の
   J テンソル解決規則、Hubbard/Kondo の演算子と符号、GC 変種。
5. **wannier90 変換仕様** — RESPACK/Wannier90 出力からの変換規則
   (H/U チャネルの符号・cutoff・Hermite 正準対選択、Spin の超交換
   自動生成、J チャネルの散文記述)。
6. **仕様拡張提案** — `dialect: experimental` が導入する 6 項目の
   拡張(fermion site_dof、1 サイト演算子語彙、名前付き 2 体演算子、
   `{param, scale, default}` 参照、`catalog:` ヘッダ、4 フェルミオン
   一般項の素描)を、draft 仕様への追記提案として整理。
7. **既知の制限・未対応事項** — oracle 比較未実施、ladder の
   W=2/3 限定、wannier90 の対応範囲、上流実装の既知の不整合など。

規範文書(本書と矛盾する場合に優先される文書)は `CONVENTIONS.md`
であり、`manual.md` はその解説・実例・索引という位置づけである。

## 参照仕様

本カタログの geometry / system / model 三層構造は、参照仕様
「格子定義仕様 (draft)」2026/07/27(`docs/superpowers/specs/
2026-08-15-lattice-catalog-design.md` に整理されている設計判断の
拠り所)に準拠する。draft が明示的には規定していない事項は、本
README・`CONVENTIONS.md` §7・`manual.md` 6 章で「実験的拡張方言」
として文書化した上で採用している。
