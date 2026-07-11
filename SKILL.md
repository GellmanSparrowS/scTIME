# ONC scRNA-seq Analysis Skill

## Overview
Patterns, pitfalls, and reusable strategies for multi-timepoint single-cell RNA-seq analysis, developed during GSE137398 ONC RGC project.

---

## Data Loading

### SCP metadata quirks
- Cell barcodes are in the **DataFrame index**, NOT in a named column
- `NAME` column = cluster/subtype annotation, not cell ID
- First row = `TYPE` header artifact → always `df.iloc[1:]`
- Read with `pd.read_csv(..., sep="\t", index_col=0)`

### Count matrices
- GEO matrices: gene×cell, first column unnamed → `"Unnamed: 0"`
- Cell columns: `{SampleID}_{10xBarcode}` format
- All time-point matrices share identical gene sets (verify early)
- Chunked reading: `chunksize=5000`, `usecols=["Unnamed: 0"] + keep_barcodes`

---

## Memory Management

### Large sparse matrices (40K genes × 100K cells)
- Convert to `scipy.sparse.csr_matrix` with `dtype=np.float32`
- A single `.copy()` on int64 sparse can OOM (>1.5GB) — avoid deep copies
- Use `adata = adata[mask]` without `.copy()`
- Call `adata.X = adata.X.astype(np.float32)` before heavy ops

### Checkpoint strategy
- Save `.h5ad` after each major step: merge → QC → normalize → preprocess
- Check `Path.exists()` before expensive computations
- Saves hours of recomputation on crash

---

## Preprocessing Pipeline

### QC (ONC RGC-specific)
```
min_genes=200, max_genes=6000, min_counts=500, max_mito=20%
```
- Mouse MT genes: `mt-` prefix (lowercase)
- ONC RGCs: low mito% (<5%), ~8% cells removed

### Normalize → HVG → PCA → UMAP
- `sc.pp.normalize_total(target_sum=10000)` + `sc.pp.log1p`
- HVG: `seurat` flavor (avoid `seurat_v3` → needs scikit-misc → numpy conflict)
- 3000 HVGs from 40790 genes
- PCA 50 comps, neighbors k=15 with 30 PCs

---

## DEG Strategy

### Two complementary approaches
1. **Single-cell wilcoxon** (`sc.tl.rank_genes_groups`): High sensitivity, ranks all genes. Useful for exploration but treat p-values cautiously (cells aren't independent replicates).
2. **Pseudobulk t-test**: Aggregate per sample, then `scipy.stats.ttest_ind` + BH correction. Lower sensitivity but statistically valid biological replication.

### Pseudobulk construction
- Sum raw counts per sample: `adata[cells].X.sum(axis=0)`
- Normalize to CPM: `X / lib_size * 1e6`
- 32 GSM → ~21 samples after QC

### Key DEG genes (confirmed by both methods)
- Up: `Ddit3`, `Atf3`, `Gadd45a`, `Hmox1`, `Cd63`, `Tubb6`, `Phgdh`
- Down: `Chrna6`, `Shroom2`, `Nptx1`, `Cacng5`, `Snap25`

---

## Time-course Analysis

### Gene clustering
- Compute mean expression per gene × time point
- Z-score across time, then KMeans (n=6)
- Sort clusters by peak expression time

### Transition genes
- Largest delta between consecutive time points
- Key transitions: 1d→2d (`Tuba1a`), 2d→4d (`Ecel1`)

---

## Subtype Vulnerability

### Analysis approach
- `pd.crosstab(subtype, time_point)` → proportions
- Vulnerability score = late_pct / ctrl_pct
- < 0.5 = vulnerable, > 1.5 = resilient
- Chi-squared test for composition shift significance

### Key finding
- "Unassigned" cells balloon from 0% → 16.9% — loss of subtype identity under stress

---

## Pathway Enrichment

### gseapy / Enrichr
- Organism: `"mouse"` (lowercase)
- Gene sets: `"GO_Biological_Process_2023"`
- Network required; wrap in try/except

### Alternatively
- Rank genes by trend (late - early), save as `.rnk` for GSEA

---

## Environment

### Version constraints
- scanpy 1.12 needs `numpy < 2.5`
- scikit-misc needs `numpy >= 2.5` → avoid, use `seurat` HVG flavor
- Python 3.14 works but prefer 3.10-3.12 for stability
- Install to local `.venv/` with `--target` flag

### Output structure
```
output/
  merged_onc.h5ad      # 76646 cells × 40790 genes
  qc_filtered.h5ad     # 70049 cells after QC
  normalized_hvg.h5ad  # 70049 cells × 40790 genes, normalized
  preprocessed.h5ad    # 70049 cells × 3000 HVGs
  deg/                 # DEG CSVs + pseudobulk
  trajectory/          # gene clusters, markers, transitions
  subtype/             # composition, vulnerability
  pathway/             # GO enrichment
