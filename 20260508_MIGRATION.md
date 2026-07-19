# CUDA 12.2 → 13.2 Migration Summary

## 新規ファイル

### `cuda13_compat.h`
CUDA 12.2 / 13.x 両方でコンパイル可能な互換性ヘッダ。

| マクロ / 関数 | 説明 |
|---|---|
| `ENABLE_SMEM_SPILLING()` | CUDA 13.0+ でレジスタスピルを共有メモリへ誘導。12.x では no-op。 |
| `safe_prefetch()` | `cudaMemPrefetchAsync` の API 変更（CUDA 12.2+）に対応するラッパー |
| `cuda_select_best_gpu()` | `cudaDeviceProp` 経由の GPU 選択ヘルパー |

## 各ファイルの変更内容

### cbp.cu（CT 再構成 FBP カーネル）
1. `cuda13_compat.h` のインクルード追加
2. **テクスチャリファレンスコードの完全削除**（CUDA 13.0 で `texture<>` ヘッダが削除済み）
   - コメントアウトされていた `texture<float2,1>` 宣言、`tex1Dfetch` カーネル（`BP`）、`cudaBindTexture`/`cudaUnbindTexture` を削除
   - `BP_GMF`（グローバルメモリフェッチ版）を唯一のバックプロジェクションカーネルとして使用
3. 未使用の `TEXTURE_LIMIT` マクロと `GMF` 変数を削除
4. `BP_GMF` カーネルに `ENABLE_SMEM_SPILLING()` を追加

### tif_adf_g.cu（異方性拡散フィルタ）
1. `cuda13_compat.h` のインクルード追加
2. `aniso_diffusion_kernel` に `ENABLE_SMEM_SPILLING()` 追加
3. `uint32` → `uint32_t`、`uint16` → `uint16_t`（libtiff 4.4+ 互換）

### tif_blf_g.cu（バイラテラルフィルタ）
1. `cuda13_compat.h` のインクルード追加
2. `bilateral_filter_3d_kernel` に `ENABLE_SMEM_SPILLING()` 追加
3. `uint32` → `uint32_t`、`uint16` → `uint16_t`

### tif_bm4d_g.cu（BM4D デノイズ）
1. `cuda13_compat.h` のインクルード追加
2. `bm4d_filter_kernel` に `ENABLE_SMEM_SPILLING()` 追加
   - **最も効果が期待されるカーネル**（`match_dist/sx/sy/sz[16]` 配列によるレジスタ圧力が極めて高い）
3. `uint32` → `uint32_t`、`uint16` → `uint16_t`

### tif_gsf_g.cu（ガウシアンフィルタ）
1. `cuda13_compat.h` のインクルード追加
2. `gaussian_filter_3d_kernel` に `ENABLE_SMEM_SPILLING()` 追加
3. `uint32` → `uint32_t`、`uint16` → `uint16_t`

### tif_mdf_g.cu（メディアンフィルタ）
1. `cuda13_compat.h` のインクルード追加
2. `median_filter_3d_kernel` に `ENABLE_SMEM_SPILLING()` 追加
   - `neighbors[]` ソート配列がローカルメモリにスピルするため効果大
3. `uint32` → `uint32_t`、`uint16` → `uint16_t`

### tif_nlm_g.cu（非局所平均フィルタ）
1. `cuda13_compat.h` のインクルード追加
2. `nlm_filter_3d_kernel` に `ENABLE_SMEM_SPILLING()` 追加
   - 6 重ネストループ + 多数のローカル変数でレジスタ圧が高い
3. `uint32` → `uint32_t`、`uint16` → `uint16_t`

### tif_tvd_g.cu（Total Variation デノイズ）
1. `cuda13_compat.h` のインクルード追加
2. `dual_update_kernel` / `primal_update_kernel` の両方に `ENABLE_SMEM_SPILLING()` 追加
3. `uint32` → `uint32_t`、`uint16` → `uint16_t`

### tif_wvd_g.cu（ウェーブレットデノイズ）
1. `cuda13_compat.h` のインクルード追加
2. Haar DWT カーネルはレジスタ使用量が少ないため `ENABLE_SMEM_SPILLING()` は未追加
3. `uint32` → `uint32_t`、`uint16` → `uint16_t`

## コンパイル方法

### CUDA 13.2 でのビルド（推奨）
```bash
# Linux
nvcc -O3 -o tif_nlm_g tif_nlm_g.cu -use_fast_math -ltiff

# Windows
nvcc -O3 -I%CUDAINCL% -o tif_nlm_g.exe tif_nlm_g.cu -use_fast_math -Xcompiler "/wd 4819" libtiff.lib
```

### CUDA 12.2 でのビルド（後方互換あり）
変更なしでそのままビルド可能。`ENABLE_SMEM_SPILLING()` は CUDA 12.x では no-op に展開されます。

## 期待される性能改善

| カーネル | レジスタ圧力 | smem spilling 効果（推定） |
|---|---|---|
| `bm4d_filter_kernel` | 極めて高 | 5〜10% |
| `nlm_filter_3d_kernel` | 高 | 5〜10% |
| `median_filter_3d_kernel` | 高（ソート配列） | 5〜8% |
| `bilateral_filter_3d_kernel` | 中〜高 | 3〜7% |
| `BP_GMF` | 中 | 2〜5% |
| `aniso_diffusion_kernel` | 中 | 2〜5% |
| `dual/primal_update_kernel` | 中 | 1〜3% |

## 注意事項

- `cuda13_compat.h` は各 `.cu` ファイルと同じディレクトリに配置してください
- `cbp.cu` は `cu.h` / `cbp.h` に依存しており、それらは変更不要です
- CUDA 13.x では **Maxwell/Pascal/Volta** GPU はサポート外です（Turing 以降が必要）
- 最低ドライバ: r580 系（580.65.06+）
