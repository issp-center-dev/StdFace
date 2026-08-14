# 格子カタログ(新フォーマット対応一式)実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 現行 StdFace の全格子・全相互作用を新三層仕様(geometry/system/model)の YAML カタログ + 日本語マニュアルとして `lattice_catalog/` に作成する。

**Architecture:** 格子×模型で自己完結 YAML(9格子×3模型=27件 + wannier90 例)。ボンド type 名 = StdFace パラメータ名、係数は `{param: <キーワード>}` 参照。Task 1 (chain) で規約を確定し、残り格子へ機械的に展開する。

**Tech Stack:** YAML(手書き)、Markdown。検証は PyYAML(ローカル開発時のみ、依存追加しない)。

**Spec:** `docs/superpowers/specs/2026-08-15-lattice-catalog-design.md`

## Global Constraints

- 成果物はすべて `lattice_catalog/` 配下。リポジトリの既存コードは一切変更しない。
- マニュアル・YAML コメントは日本語。
- ボンド type 名は StdFace キーワード名そのまま(`J0'` のようにプライム含む)。
- 向きの一意化規約(仕様 §4.2): 各ボンドは一度だけ書く。from のサイトラベル辞書順 → to、同じなら R が辞書順で正(最初の非零成分が正)。
- 係数は `value: {param: <StdFaceキーワード>}` で参照。数値は書かない(system.size のみ代表値+コメント)。
- 各 YAML 冒頭コメントに: 対応する stan.in 記述例、検算(配位数・単位胞あたりボンド本数)。
- コミットは Task ごと。コミットメッセージは中立的な表現(内部ツール名を出さない)。
- 依存パッケージは追加しない(PyYAML はチェック専用でローカル利用のみ)。
- 参照実装: `python/stdface/lattice/<lattice>.py` の `_BONDS` テーブル。形式は
  `(dW, dL, site_i, site_j, nn_level, J系, t系, V系)`(3D は `(dW, dL, dH, ...)`)。
  `dW, dL(, dH)` → R、`site_i/site_j` → from/to ラベル、J系/t系/V系の変数名 → type 名。

---

### Task 1: カタログ規約の確定 + chain(3模型)

**Files:**
- Create: `lattice_catalog/chain/chain_spin.yaml`
- Create: `lattice_catalog/chain/chain_hubbard.yaml`
- Create: `lattice_catalog/chain/chain_kondo.yaml`
- Create: `lattice_catalog/CONVENTIONS.md`(規約の1ページ要約。manual.md 執筆時に取り込む)

**Interfaces:**
- Consumes: `python/stdface/lattice/chain_lattice.py:77-184`(幾何と `_BONDS`)
- Produces: 後続タスク全てが従う YAML 雛形と規約(CONVENTIONS.md)

chain の抽出済み事実(`chain_lattice.py` で確認済み):
- 1 サイト/単位胞、tau=(0,0,0)。1 次元(内部では W=1 の 2D 表現だがカタログでは dimension: 1)
- ボンド: `J0/t0/V0` = R=[1]、`J0'/t0'/V0'` = R=[2]、`J0''/t0''/V0''` = R=[3]
- 別名: `J→J0, J'→J0', J''→J0'', t→t0, t'→t0', t''→t0'', V→V0, V'→V0', V''→V0''`
- Spin: `2S, D, h, Gamma, Gamma_y`。Hubbard: `mu, U` + t/V 族 + 磁場。
  Kondo: Hubbard 項 + `J`(遍歴-局在結合、2S)。磁場は局在スピン側。
- `phase0` → 周期方向の `{twist: phase0}`

- [ ] **Step 1: CONVENTIONS.md を書く**

内容(全文):

```markdown
# lattice_catalog 記述規約

1. ファイル = 格子×模型で自己完結(geometry / system / model)。
2. ボンド type 名 = StdFace キーワード名(`J0`, `J0'`, `t1`, `V2` など)。
3. 係数は `value: {param: <StdFaceキーワード>}`。数値を直接書かない。
4. 向きの一意化: from ラベル辞書順 → to、同ラベルなら R 辞書順で正
   (最初の非零成分が正)。各ボンドは一度だけ書く。
5. `system.size` のみ代表値(コメントで StdFace キーワードを注記)。
   境界は periodic を既定、`phase*` は `{twist: ...}` で注記。
6. スピン交換の既定は等方 `operator: "S_i . S_j"`。異方成分
   (`J0x`, `J0xy` 等)は tensor_terms 展開(manual 4章・6章参照)。
7. フェルミオン演算子(拡張提案): `hop` = Σσ (c†_iσ c_jσ + h.c.)、
   `density-density` = n_i n_j、onsite `Nup*Ndn`, `N`, `Sz` 等。
8. Kondo は 2 ラベル表現: `<X>_c`(遍歴, fermion)と `<X>_s`(局在, spin)
   を同一分率座標に置き、`J` は R=0 の `<X>_c–<X>_s` ボンド。
9. 各ファイル冒頭に stan.in 対応例と検算(配位数・ボンド本数)をコメント。
10. 別名キーワード(`J→J0` 等)はコメントで注記。
```

- [ ] **Step 2: chain_spin.yaml を書く**

```yaml
# StdFace chain / Spin 模型
# stan.in 対応例:
#   model = "Spin"
#   lattice = "chain"
#   L = 16
#   J = 1.0        # J は J0 の別名
# 検算: 配位数 2 (J0)。単位胞あたり J0, J0', J0'' 各 1 本。

geometry:
  dimension: 1
  lattice_vectors:
    a1: [1.0]                    # StdFace: a (格子定数, 既定 1.0)
  sites:
    - {label: A, frac: [0.0]}

system:
  finite: true
  size: [16]                     # StdFace: L
  boundary:
    - periodic                   # StdFace: phase0 → {twist: phase0}

model:
  site_dof:
    A: {spin: 0.5}               # StdFace: 2S = 1 (spin = 2S/2)

  bonds:
    - {from: A, to: A, R: [1], type: J0}    # 最近接   (別名 J)
    - {from: A, to: A, R: [2], type: J0'}   # 次近接   (別名 J')
    - {from: A, to: A, R: [3], type: J0''}  # 三次近接 (別名 J'')

  couplings:
    J0:   {operator: "S_i . S_j", value: {param: J0}}
    J0':  {operator: "S_i . S_j", value: {param: J0'}}
    J0'': {operator: "S_i . S_j", value: {param: J0''}}
    # 異方交換 (J0x, J0xy 等) は tensor_terms 展開を用いる (manual 4章)

  onsite:
    A:
      field_z:                   # Zeeman: -h Sz
        operator: {tensor_terms: [{ops: [Sz], coeff: -1.0}]}
        value: {param: h}
      field_x:                   # 横磁場: -Gamma Sx
        operator: {tensor_terms: [{ops: [Sx], coeff: -1.0}]}
        value: {param: Gamma}
      field_y:                   # 横磁場: -Gamma_y Sy
        operator: {tensor_terms: [{ops: [Sy], coeff: -1.0}]}
        value: {param: Gamma_y}
      aniso_z:                   # 単イオン異方性: D (Sz)^2
        operator: {tensor_terms: [{ops: [Szz], coeff: 1.0}]}
        value: {param: D}
```

- [ ] **Step 3: chain_hubbard.yaml を書く**

```yaml
# StdFace chain / Hubbard 模型
# stan.in 対応例:
#   model = "Hubbard"
#   lattice = "chain"
#   L = 16
#   t = 1.0        # t は t0 の別名
#   U = 4.0
# 検算: 配位数 2 (t0)。単位胞あたり t0, t0', t0'' / V0, V0', V0'' 各 1 本。

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
    - periodic                   # StdFace: phase0 → {twist: phase0}

model:
  site_dof:
    A: {fermion: {orbitals: 1}}  # 拡張提案 (manual 6章)

  bonds:
    - {from: A, to: A, R: [1], type: t0}    # 最近接ホッピング (別名 t)
    - {from: A, to: A, R: [2], type: t0'}   # 次近接 (別名 t')
    - {from: A, to: A, R: [3], type: t0''}  # 三次近接 (別名 t'')
    - {from: A, to: A, R: [1], type: V0}    # 最近接クーロン (別名 V)
    - {from: A, to: A, R: [2], type: V0'}   # 次近接 (別名 V')
    - {from: A, to: A, R: [3], type: V0''}  # 三次近接 (別名 V'')

  couplings:
    # hop = -t Σσ (c†_iσ c_jσ + h.c.)。t は複素可、h.c. は向き規約で暗黙。
    t0:   {operator: "hop", value: {param: t0}}
    t0':  {operator: "hop", value: {param: t0'}}
    t0'': {operator: "hop", value: {param: t0''}}
    V0:   {operator: "density-density", value: {param: V0}}   # V n_i n_j
    V0':  {operator: "density-density", value: {param: V0'}}
    V0'': {operator: "density-density", value: {param: V0''}}

  onsite:
    A:
      hubbard_u:                 # U n↑ n↓
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

- [ ] **Step 4: chain_kondo.yaml を書く**

```yaml
# StdFace chain / Kondo 格子模型
# stan.in 対応例:
#   model = "Kondo"
#   lattice = "chain"
#   L = 16
#   t = 1.0
#   J = 1.0        # Kondo 結合 (遍歴-局在)
# 検算: 遍歴鎖の配位数 2 (t0)。J は各セル内 A_c–A_s の 1 本 (R=0)。
# 現行実装ではサイト倍加 (jsite_kondo = isite + nsite/2) で表現される。

geometry:
  dimension: 1
  lattice_vectors:
    a1: [1.0]                    # StdFace: a
  sites:
    - {label: A_c, frac: [0.0]}  # 遍歴電子
    - {label: A_s, frac: [0.0]}  # 局在スピン (同一座標)

system:
  finite: true
  size: [16]                     # StdFace: L
  boundary:
    - periodic                   # StdFace: phase0 → {twist: phase0}

model:
  site_dof:
    A_c: {fermion: {orbitals: 1}}   # 拡張提案 (manual 6章)
    A_s: {spin: 0.5}                # StdFace: 2S = 1

  bonds:
    - {from: A_c, to: A_c, R: [1], type: t0}
    - {from: A_c, to: A_c, R: [2], type: t0'}
    - {from: A_c, to: A_c, R: [3], type: t0''}
    - {from: A_c, to: A_c, R: [1], type: V0}
    - {from: A_c, to: A_c, R: [2], type: V0'}
    - {from: A_c, to: A_c, R: [3], type: V0''}
    - {from: A_c, to: A_s, R: [0], type: J}    # Kondo 結合

  couplings:
    t0:   {operator: "hop", value: {param: t0}}
    t0':  {operator: "hop", value: {param: t0'}}
    t0'': {operator: "hop", value: {param: t0''}}
    V0:   {operator: "density-density", value: {param: V0}}
    V0':  {operator: "density-density", value: {param: V0'}}
    V0'': {operator: "density-density", value: {param: V0''}}
    J:    {operator: "s_i . S_j", value: {param: J}}   # 遍歴スピン s と局在 S

  onsite:
    A_c:
      hubbard_u:
        operator: {tensor_terms: [{ops: [NupNdn], coeff: 1.0}]}
        value: {param: U}
      chemical_potential:
        operator: {tensor_terms: [{ops: [N], coeff: -1.0}]}
        value: {param: mu}
    A_s:
      field_z:                   # 磁場は局在スピン側
        operator: {tensor_terms: [{ops: [Sz], coeff: -1.0}]}
        value: {param: h}
      field_x:
        operator: {tensor_terms: [{ops: [Sx], coeff: -1.0}]}
        value: {param: Gamma}
```

- [ ] **Step 5: 構文チェック**

Run: `python3 -c "import yaml,glob; [yaml.safe_load(open(f)) for f in sorted(glob.glob('lattice_catalog/chain/*.yaml'))]; print('OK')"`
Expected: `OK`

- [ ] **Step 6: 検算(目視+ソース突合)**

`python/stdface/lattice/chain_lattice.py:178-182` の `_BONDS` と YAML の bonds を 1 行ずつ突合。
(0,1,0,0,1,J0,t0,V0) → R=[1] の J0/t0/V0、(0,2,...) → R=[2]、(0,3,...) → R=[3]。
不一致があれば YAML を修正。

さらに onsite 項(h/Gamma/Gamma_y/D/U/mu)の模型ごとの適用範囲を
`python/stdface/core/model_plugin.py`(SpinModel/HubbardModel/KondoModel の
項生成)と突合し、各模型で実際に生成されない項は YAML から削除する
(例: 非 Spin 模型で D は `not_used`)。確定した適用範囲は
CONVENTIONS.md に追記し、後続の格子タスクでも同じ範囲を使う。

- [ ] **Step 7: Commit**

```bash
git add lattice_catalog/CONVENTIONS.md lattice_catalog/chain/
git commit -m "Add lattice catalog conventions and chain definitions (spin/Hubbard/Kondo)"
```

---

### Task 2: square(3模型)

**Files:**
- Create: `lattice_catalog/square/square_spin.yaml`
- Create: `lattice_catalog/square/square_hubbard.yaml`
- Create: `lattice_catalog/square/square_kondo.yaml`

**Interfaces:**
- Consumes: Task 1 の CONVENTIONS.md と chain の 3 YAML(雛形として全面的に踏襲)
- Produces: なし(独立成果物)

- [ ] **Step 1: ソースから幾何とボンドを抽出**

`python/stdface/lattice/square_lattice.py` を読み、以下を確定:
- `NsiteUC`、`tau`(サイト数・分率座標)、`direct` の既定(格子ベクトル)
- `_BONDS` テーブル全行 → `(dW,dL)` を R、変数名を type 名に変換
- 使用パラメータ(確認済み概要: t/J/V の 0,1 系 + p/pp。`J0`=W方向nn, `J1`=L方向nn 等の方向割当をソースコメントと `input_spin_nn` の名前引数で確認)
- Spin/Hubbard/Kondo 各分岐で `not_used_*` になるパラメータ(YAML に含めない)

- [ ] **Step 2: 3 つの YAML を書く**

chain の雛形(Task 1 Step 2–4)と同一構造で作成。相違点のみ:
- `geometry.dimension: 2`、`lattice_vectors: a1: [1.0, 0.0], a2: [0.0, 1.0]`(ソースの既定値を確認して記載)
- `system.size: [4, 4]  # StdFace: W, L`、`boundary: [periodic, periodic]  # phase0, phase1`
- bonds/couplings は Step 1 の抽出結果を全数記載
- 冒頭コメントの検算: 最近接配位数 4(J0+J1 で W/L 方向各 2)、
  次近接(対角 J')配位数 4

- [ ] **Step 3: 構文チェック**

Run: `python3 -c "import yaml,glob; [yaml.safe_load(open(f)) for f in sorted(glob.glob('lattice_catalog/square/*.yaml'))]; print('OK')"`
Expected: `OK`

- [ ] **Step 4: 検算(ソース突合)**

`_BONDS` テーブルと YAML bonds を 1 行ずつ突合。本数一致・向き規約準拠を確認。

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/square/
git commit -m "Add square lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 3: triangular(3模型)

**Files:**
- Create: `lattice_catalog/triangular/triangular_spin.yaml`
- Create: `lattice_catalog/triangular/triangular_hubbard.yaml`
- Create: `lattice_catalog/triangular/triangular_kondo.yaml`

**Interfaces:**
- Consumes: Task 1 の規約・雛形
- Produces: なし

- [ ] **Step 1: `python/stdface/lattice/triangular_lattice.py` から抽出**

Task 2 Step 1 と同じ手順。注意点:
- 格子ベクトル既定は a1=(1,0), a2=(1/2,√3/2)(ソースの `direct` 既定値で確認)
- t/J/V は 0,1,2 系 + p/pp まで使用(3 方向の nn を J0/J1/J2 で区別)
- 検算: 最近接配位数 6(J0/J1/J2 各 2)

- [ ] **Step 2: 3 つの YAML を書く**(chain 雛形踏襲、抽出結果を全数記載)

- [ ] **Step 3: 構文チェック**

Run: `python3 -c "import yaml,glob; [yaml.safe_load(open(f)) for f in sorted(glob.glob('lattice_catalog/triangular/*.yaml'))]; print('OK')"`
Expected: `OK`

- [ ] **Step 4: 検算(ソース突合)** — `_BONDS` と 1 行ずつ突合

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/triangular/
git commit -m "Add triangular lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 4: honeycomb(3模型)

**Files:**
- Create: `lattice_catalog/honeycomb/honeycomb_spin.yaml`
- Create: `lattice_catalog/honeycomb/honeycomb_hubbard.yaml`
- Create: `lattice_catalog/honeycomb/honeycomb_kondo.yaml`

**Interfaces:**
- Consumes: Task 1 の規約・雛形
- Produces: なし

- [ ] **Step 1: `python/stdface/lattice/honeycomb_lattice.py` から抽出**

注意点:
- 2 サイト/単位胞(A/B)。`tau` から分率座標を確定。
- Kitaev 的異方性の慣例: J0/J1/J2 が 3 種のボンド方向に対応
  (x/y/z ボンド)。この対応を YAML コメントに明記。
- 検算: 最近接配位数 3(A から J0/J1/J2 各 1 本)
- Kondo は 4 ラベル(A_c, A_s, B_c, B_s)

- [ ] **Step 2: 3 つの YAML を書く**

- [ ] **Step 3: 構文チェック**

Run: `python3 -c "import yaml,glob; [yaml.safe_load(open(f)) for f in sorted(glob.glob('lattice_catalog/honeycomb/*.yaml'))]; print('OK')"`
Expected: `OK`

- [ ] **Step 4: 検算(ソース突合)**

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/honeycomb/
git commit -m "Add honeycomb lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 5: kagome(3模型)

**Files:**
- Create: `lattice_catalog/kagome/kagome_spin.yaml`
- Create: `lattice_catalog/kagome/kagome_hubbard.yaml`
- Create: `lattice_catalog/kagome/kagome_kondo.yaml`

**Interfaces:**
- Consumes: Task 1 の規約・雛形
- Produces: なし

- [ ] **Step 1: `python/stdface/lattice/kagome.py` から抽出**

確認済み事実: 3 サイト/単位胞、tau = (0,0), (1/2,0), (0,1/2)。
`_BONDS` は `kagome.py:164-179`。p 系まで使用(pp なし)。
検算: 最近接配位数 4。

- [ ] **Step 2: 3 つの YAML を書く**(Kondo は 6 ラベル)

- [ ] **Step 3: 構文チェック**

Run: `python3 -c "import yaml,glob; [yaml.safe_load(open(f)) for f in sorted(glob.glob('lattice_catalog/kagome/*.yaml'))]; print('OK')"`
Expected: `OK`

- [ ] **Step 4: 検算(ソース突合)**

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/kagome/
git commit -m "Add kagome lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 6: ladder(3模型)

**Files:**
- Create: `lattice_catalog/ladder/ladder_spin.yaml`
- Create: `lattice_catalog/ladder/ladder_hubbard.yaml`
- Create: `lattice_catalog/ladder/ladder_kondo.yaml`

**Interfaces:**
- Consumes: Task 1 の規約・雛形
- Produces: なし

- [ ] **Step 1: `python/stdface/lattice/ladder.py` から抽出**

注意点:
- W = 脚(leg)本数。単位胞内に W サイト(ラベル例: A0..A{W-1})か、
  W をサイズ側に持つかをソースで確認し、**単位胞 = 1 ラング(W サイト)、
  size = [L]** の 1 次元表現を採用(梯子は rung 方向非周期のため)。
- `_BONDS` は動的構築(`ladder.py:176-187`)なので、ループを読んで
  脚内 (J1, J1')・ラング (J0)・斜め (J2, J2') の割当を確定。
- 検算: 2 本脚のとき脚方向配位数 2、ラング 1。
- 代表値: `system.size: [16]`、W=2 の 2 ラベル(脚 2 本)を既定例とし、
  W 一般化はコメントで説明。

- [ ] **Step 2: 3 つの YAML を書く**

- [ ] **Step 3: 構文チェック**

Run: `python3 -c "import yaml,glob; [yaml.safe_load(open(f)) for f in sorted(glob.glob('lattice_catalog/ladder/*.yaml'))]; print('OK')"`
Expected: `OK`

- [ ] **Step 4: 検算(ソース突合)**

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/ladder/
git commit -m "Add ladder lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 7: orthorhombic(3模型)

**Files:**
- Create: `lattice_catalog/orthorhombic/orthorhombic_spin.yaml`
- Create: `lattice_catalog/orthorhombic/orthorhombic_hubbard.yaml`
- Create: `lattice_catalog/orthorhombic/orthorhombic_kondo.yaml`

**Interfaces:**
- Consumes: Task 1 の規約・雛形
- Produces: 3D 格子の雛形(Task 8, 9 が参照)

- [ ] **Step 1: `python/stdface/lattice/orthorhombic.py` から抽出**

注意点:
- 3 次元: `dimension: 3`、`lattice_vectors: a1/a2/a3`、
  `system.size: [4, 4, 4]  # StdFace: W, L, Height`、boundary 3 方向。
- `_BONDS` は 3D 形式 `(dW, dL, dH, ...)` → R = [dW, dL, dH]。
- J0/J1/J2 = 3 軸方向 nn。pp(体対角など)の割当をソースで確認。
- 検算: 最近接配位数 6。

- [ ] **Step 2: 3 つの YAML を書く**

- [ ] **Step 3: 構文チェック**

Run: `python3 -c "import yaml,glob; [yaml.safe_load(open(f)) for f in sorted(glob.glob('lattice_catalog/orthorhombic/*.yaml'))]; print('OK')"`
Expected: `OK`

- [ ] **Step 4: 検算(ソース突合)**

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/orthorhombic/
git commit -m "Add orthorhombic lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 8: fc_ortho(3模型)

**Files:**
- Create: `lattice_catalog/fc_ortho/fc_ortho_spin.yaml`
- Create: `lattice_catalog/fc_ortho/fc_ortho_hubbard.yaml`
- Create: `lattice_catalog/fc_ortho/fc_ortho_kondo.yaml`

**Interfaces:**
- Consumes: Task 1 の規約、Task 7 の 3D 雛形
- Produces: なし

- [ ] **Step 1: `python/stdface/lattice/fc_ortho.py` から抽出**

注意点: 面心系の格子ベクトル既定(`direct`)をソースで確認して記載。
J0/J1/J2 の方向割当と p/pp の使用範囲を `_BONDS` から確定。

- [ ] **Step 2: 3 つの YAML を書く**

- [ ] **Step 3: 構文チェック**

Run: `python3 -c "import yaml,glob; [yaml.safe_load(open(f)) for f in sorted(glob.glob('lattice_catalog/fc_ortho/*.yaml'))]; print('OK')"`
Expected: `OK`

- [ ] **Step 4: 検算(ソース突合)**

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/fc_ortho/
git commit -m "Add face-centered orthorhombic lattice catalog definitions"
```

---

### Task 9: pyrochlore(3模型)

**Files:**
- Create: `lattice_catalog/pyrochlore/pyrochlore_spin.yaml`
- Create: `lattice_catalog/pyrochlore/pyrochlore_hubbard.yaml`
- Create: `lattice_catalog/pyrochlore/pyrochlore_kondo.yaml`

**Interfaces:**
- Consumes: Task 1 の規約、Task 7 の 3D 雛形
- Produces: なし

- [ ] **Step 1: `python/stdface/lattice/pyrochlore.py` から抽出**

注意点: 4 副格子(tau 4 点)。FCC 格子ベクトル既定を確認。
四面体内ボンド(R=0)と四面体間ボンドの区別をコメントに明記。
検算: 最近接配位数 6。Kondo は 8 ラベル。

- [ ] **Step 2: 3 つの YAML を書く**

- [ ] **Step 3: 構文チェック**

Run: `python3 -c "import yaml,glob; [yaml.safe_load(open(f)) for f in sorted(glob.glob('lattice_catalog/pyrochlore/*.yaml'))]; print('OK')"`
Expected: `OK`

- [ ] **Step 4: 検算(ソース突合)**

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/pyrochlore/
git commit -m "Add pyrochlore lattice catalog definitions (spin/Hubbard/Kondo)"
```

---

### Task 10: wannier90 写像例

**Files:**
- Create: `lattice_catalog/wannier90/example.yaml`

**Interfaces:**
- Consumes: Task 1 の規約。`python/stdface/lattice/wannier90.py`(読み込み処理)、
  `test/wannier90_data/` または `samples/`(実在の小規模データ)
- Produces: manual 5 章(Task 12)が参照する具体例

- [ ] **Step 1: 入力データを選ぶ**

`test/wannier90_data/` と `samples/` を調べ、軌道数が最小の
wannier90 テストケース(_hr.dat / _geom)を 1 つ選ぶ。

- [ ] **Step 2: example.yaml を書く**

選んだデータから:
- `geometry`: _geom の格子ベクトルと Wannier 中心 → sites(軌道 = ラベル W1, W2, …)
- `bonds`: _hr.dat の H_mn(R) の代表要素(絶対値上位 + cutoff コメント)
  → `{from: Wm, to: Wn, R: [...], type: t_mn_R...}`。
  全要素は書かず、**写像規則の実演として上位 10 件程度**に留め、
  「全要素は同規則で機械変換」とコメントで明記。
- `couplings`: hop 演算子 + 数値(この例のみ `value:` に実数値を書く。
  外部データ由来のため param 参照できない旨コメント)
- U_mn(R)/J_mn(R)(`cutoff_u`, `cutoff_j`)がテストデータにあれば
  density-density / S_i.S_j の例を各 2–3 件追加。

- [ ] **Step 3: 構文チェック**

Run: `python3 -c "import yaml; yaml.safe_load(open('lattice_catalog/wannier90/example.yaml')); print('OK')"`
Expected: `OK`

- [ ] **Step 4: 検算** — 記載した H_mn(R) 要素を元ファイルの行と突合(行番号をコメントに記録)

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/wannier90/
git commit -m "Add wannier90 mapping example for the lattice catalog"
```

---

### Task 11: manual.md 前半(1–2章: 概要・キーワード対応表)

**Files:**
- Create: `lattice_catalog/manual.md`(1–2 章)

**Interfaces:**
- Consumes: `python/stdface/core/keyword_parser.py:337-449`(`_COMMON_KEYWORDS` 全表)、
  Task 1 の CONVENTIONS.md
- Produces: manual.md の章立て(Task 12, 13 が追記)

- [ ] **Step 1: 1 章(概要と読み方)を書く**

内容: 三層仕様の要約(geometry/system/model、本カタログが仕様 §8 カタログ層の
実体であること)、ファイル構成一覧、CONVENTIONS.md の規約 10 項目の転記・詳説、
`value: {param: ...}` の読み方。

- [ ] **Step 2: 2 章(キーワード対応表)を書く**

`_COMMON_KEYWORDS` の全キーワードを列挙した表:

| StdFace キーワード | 新フォーマットでの対応 | 備考 |

分類: (a) geometry へ(`a`, `wx/wy/...`, `wlength/llength`)、
(b) system へ(`W/L/Height`, `phase0-2`, `box` → supercell(S))、
(c) model.bonds/couplings へ(t/J/V 族全変種)、
(d) model.onsite へ(`U`, `mu`, `D`, `h`, `Gamma`, `Gamma_y`)、
(e) site_dof へ(`2S`)、
(f) wannier90 専用(`cutoff_*`, `lambda*`, `alpha`, `doublecounting` → 5 章参照)、
(g) 対象外(`model`, `lattice`, `outputmode`, `ncond/nelec`, `2Sz`, `K` —
計算条件・粒子数条件であり格子定義でない旨明記)

- [ ] **Step 3: 構文・整合チェック** — 表のキーワード数と `_COMMON_KEYWORDS` のエントリ数を照合(欠落なし)

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
- Consumes: Task 1–10 の全 YAML と各タスクの検算結果
- Produces: なし

- [ ] **Step 1: 3 章(格子ごとの解説)を書く**

9 格子それぞれ 1 節: 幾何(格子ベクトル・副格子表)、ボンド定義表
(type / from–to / R / 意味)、検算(配位数・本数)、
実装ソース突合記録(`<file>.py` の `_BONDS` 行番号)。
各 YAML の冒頭コメントと同内容をまとめ直し、格子間の比較ができる形にする。

- [ ] **Step 2: 4 章(模型ごとの演算子対応)を書く**

- Spin: `S_i . S_j` と 3×3 テンソルの tensor_terms 展開規則。
  成分キーワード対応表(`jx→[Sx,Sx]`, `jxy→[Sx,Sy]`, … 9 成分)と
  ワーク例 1 件(J0 異方の完全展開)。
- Hubbard: `hop` / `density-density` / onsite 項の定義(数式付き)。
- Kondo: 2 ラベル表現の説明と、現行実装のサイト倍加
  (`site_util.py` の該当処理)との対応。GC 変種は粒子数条件のみの
  違いで格子定義は同一である旨注記。

- [ ] **Step 3: 5 章(wannier90 の写像)を書く**

H_mn(R)/U_mn(R)/J_mn(R) → bonds+couplings の一般規則、
`cutoff_*`(値・R 範囲・長さ)の扱い、`lambda*/alpha/doublecounting` の位置づけ、
Task 10 の example.yaml の読み解き。

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

draft 仕様への追記提案として:
- `site_dof` のフェルミオン自由度 `{fermion: {orbitals: n}}`
- 1 サイト演算子語彙: `Cdag_up, Cdag_dn, C_up, C_dn, N, Nup, Ndn, NupNdn`
- 名前付き 2 体演算子: `hop`, `density-density`, `s_i . S_j`
  (それぞれ tensor_terms への展開形を明記)
- 異方交換の係数参照: tensor_terms の `coeff` に `{param: ...}` を
  許す拡張(2 章・4 章で使用した形式)

- [ ] **Step 2: 7 章(既知の制限・未対応事項)を書く**

- 展開エンジン未実装のため実行検証はしていない(ソース突合のみ)
- 一般化セル(`box` 行列)は supercell(S) 写像の説明のみで YAML 例なし
- wannier90 は全要素変換でなく規則の実演
- 多体項(3 体以上)は現行 StdFace に存在しないため対象外
  (仕様側の機構がそのまま使える旨注記)

- [ ] **Step 3: README.md を書く**

内容: 一式の目的(1 段落)、ディレクトリ構成、manual.md への案内、
参照仕様(格子定義仕様 draft 2026/07/27)への言及。

- [ ] **Step 4: 最終検証**

Run: `python3 -c "import yaml,glob; fs=sorted(glob.glob('lattice_catalog/**/*.yaml', recursive=True)); [yaml.safe_load(open(f)) for f in fs]; print(len(fs), 'files OK')"`
Expected: `28 files OK`(27 + wannier90)

manual.md 2 章の表と全 YAML の `{param: ...}` 参照名を突合し、
表にないパラメータ参照・参照されない表エントリがないか確認。

- [ ] **Step 5: Commit**

```bash
git add lattice_catalog/
git commit -m "Complete lattice catalog manual and README"
```
