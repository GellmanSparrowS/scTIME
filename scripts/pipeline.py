"""
GSE137398 ONC RGC Analysis Pipeline
Memory-efficient loading, QC, preprocessing, checkpoint support via AnnData.
"""
import os
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
from tqdm import tqdm

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


class ONCPipeline:
    """Main pipeline for GSE137398 ONC RGC time-course analysis."""

    # --- Paths ---
    BASE = Path(r"F:\ONC\GSE137398_ONC_RGC_project")
    MATRIX_DIR = BASE / "02_processed_counts" / "original_matrix_extracted"
    METADATA_DIR = BASE / "01_metadata" / "SCP509"
    OUTPUT_DIR = Path(r"F:\ONC\output")

    # --- Time points ---
    TIME_POINTS = {
        "Ctrl": "control", "12h": "12h_afterCrush", "1d": "1d_afterCrush",
        "2d": "2d_afterCrush", "4d": "4d_afterCrush", "1w": "1w_afterCrush",
        "2w": "2w_afterCrush",
    }
    TIME_ORDER = ["Ctrl", "12h", "1d", "2d", "4d", "1w", "2w"]
    TIME_DAYS = {"Ctrl": 0, "12h": 0.5, "1d": 1, "2d": 2, "4d": 4, "1w": 7, "2w": 14}
    TIME_STAGE = {
        "Ctrl": "Control", "12h": "Early", "1d": "Early", "2d": "Early",
        "4d": "Middle", "1w": "Middle", "2w": "Late",
    }

    # --- QC thresholds (adjustable) ---
    MIN_GENES = 200
    MAX_GENES = 6000
    MIN_COUNTS = 500
    MAX_MITO_PCT = 20.0
    MITO_PREFIX = "mt-"

    # --- Checkpoints ---
    CACHE_MERGED = OUTPUT_DIR / "merged_onc.h5ad"
    CACHE_QC = OUTPUT_DIR / "qc_filtered.h5ad"
    CACHE_NORM = OUTPUT_DIR / "normalized_hvg.h5ad"
    CACHE_PREPROC = OUTPUT_DIR / "preprocessed.h5ad"

    def __init__(self, output_dir=None):
        if output_dir:
            self.OUTPUT_DIR = Path(output_dir)
        os.makedirs(self.OUTPUT_DIR, exist_ok=True)
        self.adata = None

    # ==================================================================
    # Step 1: Load metadata
    # ==================================================================
    def load_metadata(self):
        print("Loading metadata...")
        meta_path = self.METADATA_DIR / "RGC_ONC_metadata.txt"
        self.meta = pd.read_csv(meta_path, sep="\t", index_col=0, low_memory=False)
        if self.meta.index[0] == "TYPE":
            self.meta = self.meta.iloc[1:]
        print(f"  ONC+Atlas metadata: {self.meta.shape[0]} cells")

        cluster_path = self.METADATA_DIR / "cluster_annotation_original.tsv"
        self.cluster = pd.read_csv(cluster_path, sep="\t", index_col=0, low_memory=False)
        if self.cluster.index[0] == "TYPE":
            self.cluster = self.cluster.iloc[1:]
        for c in ["X", "Y"]:
            if c in self.cluster.columns:
                self.cluster[c] = pd.to_numeric(self.cluster[c], errors="coerce")
        print(f"  ONC cluster annotation: {self.cluster.shape[0]} cells with UMAP")

        subtype_path = self.METADATA_DIR / "rgc_subtype_annotation_original.tsv"
        self.subtype = pd.read_csv(subtype_path, sep="\t", index_col=0, low_memory=False)
        if self.subtype.index[0] == "TYPE":
            self.subtype = self.subtype.iloc[1:]
        for c in ["X", "Y"]:
            if c in self.subtype.columns:
                self.subtype[c] = pd.to_numeric(self.subtype[c], errors="coerce")
        print(f"  Atlas subtype annotation: {self.subtype.shape[0]} cells")

        self._build_cell_lookup()
        return self

    def _build_cell_lookup(self):
        self.onc_cell_set = set(str(x) for x in self.meta.index)
        idx = self.meta.index.tolist()
        self.bc_to_subtype = dict(zip(idx, self.meta["NAME"]))
        self.bc_to_tp = dict(zip(idx, self.meta["Cluster"]))
        self.bc_to_sample = dict(zip(idx, self.meta["SampleID"]))
        self.bc_to_batch = dict(zip(idx, self.meta["BatchID"]))
        tp_counts = self.meta["Cluster"].value_counts()
        print(f"  Cells per group: {dict(tp_counts)}")

    # ==================================================================
    # Step 2: Load count matrices -> merged AnnData
    # ==================================================================
    def load_and_merge(self, chunk_size=5000):
        if self.adata is not None:
            print("AnnData already loaded, skipping.")
            return self
        print("Loading count matrices (chunked -> sparse -> merge)...")
        all_adatas = []
        for tp_label in self.TIME_ORDER:
            tp_suffix = self.TIME_POINTS[tp_label]
            fpath = self.MATRIX_DIR / f"GSE137398_ONCRGCs_{tp_suffix}_count_mat.csv"
            if not fpath.exists():
                print(f"  [SKIP] {fpath} not found")
                continue
            print(f"  Processing {tp_label}...")
            adata_tp = self._load_single_matrix(fpath, tp_label, chunk_size)
            if adata_tp is not None:
                all_adatas.append(adata_tp)
        if not all_adatas:
            raise RuntimeError("No matrices loaded.")
        self.adata = ad.concat(all_adatas, join="inner", index_unique="-")
        del all_adatas
        print(f"\nMerged AnnData: {self.adata.n_obs} cells x {self.adata.n_vars} genes")
        self._add_derived_obs()
        return self

    def _load_single_matrix(self, fpath, tp_label, chunk_size):
        header = pd.read_csv(fpath, nrows=0, index_col=0)
        all_barcodes = list(header.columns)
        keep_barcodes = [b for b in all_barcodes if str(b) in self.onc_cell_set]
        n_total, n_keep = len(all_barcodes), len(keep_barcodes)
        print(f"    {n_total} -> {n_keep} annotated cells")
        if n_keep == 0:
            print(f"    [WARN] No ONC cells in {tp_label}, skipping")
            return None
        genes = []
        rows_list = []
        for chunk in tqdm(
            pd.read_csv(fpath, index_col="Unnamed: 0",
                        usecols=["Unnamed: 0"] + keep_barcodes, chunksize=chunk_size),
            desc=f"    {tp_label} chunks",
        ):
            genes.extend(chunk.index.tolist())
            rows_list.append(sp.csr_matrix(chunk.values, dtype=np.float32))
        X_sparse = sp.vstack(rows_list, format="csr").T.tocsr()
        obs = pd.DataFrame(index=keep_barcodes)
        obs["time_point"] = tp_label
        obs["time_days"] = self.TIME_DAYS[tp_label]
        obs["stage"] = self.TIME_STAGE[tp_label]
        obs["rna_subtype"] = [self.bc_to_subtype.get(b, "unknown") for b in keep_barcodes]
        obs["sample_id"] = [self.bc_to_sample.get(b, "unknown") for b in keep_barcodes]
        obs["batch_id"] = [self.bc_to_batch.get(b, "unknown") for b in keep_barcodes]
        var = pd.DataFrame(index=genes)
        return ad.AnnData(X=X_sparse, obs=obs, var=var)

    def _add_derived_obs(self):
        if self.adata is None:
            return
        self.adata.obs["time_point"] = pd.Categorical(
            self.adata.obs["time_point"], categories=self.TIME_ORDER, ordered=True)
        self.adata.obs["stage"] = pd.Categorical(
            self.adata.obs["stage"], categories=["Control", "Early", "Middle", "Late"], ordered=True)
        if hasattr(self, "cluster") and self.cluster is not None and len(self.cluster) > 0:
            coord_map = self.cluster[["X", "Y"]]
            for c in ["X", "Y"]:
                col_name = f"UMAP_{c.lower()}"
                self.adata.obs[col_name] = self.adata.obs.index.map(
                    lambda x, cm=coord_map, cname=c: cm.loc[x, cname] if x in cm.index else np.nan)

    # ==================================================================
    # Step 3: QC
    # ==================================================================
    def qc_filter(self):
        if self.adata is None:
            raise ValueError("No AnnData. Run load_and_merge() first.")
        print("Computing QC metrics...")
        self.adata.X = self.adata.X.astype(np.float32)
        self.adata.raw = self.adata
        self.adata.obs["n_genes"] = (self.adata.X > 0).sum(axis=1).A1
        self.adata.obs["n_counts"] = self.adata.X.sum(axis=1).A1
        mito_genes = [g for g in self.adata.var_names if g.lower().startswith(self.MITO_PREFIX)]
        if mito_genes:
            mito_idx = [i for i, g in enumerate(self.adata.var_names) if g in mito_genes]
            mito_counts = np.array(self.adata.X[:, mito_idx].sum(axis=1)).flatten()
            self.adata.obs["pct_mito"] = mito_counts / self.adata.obs["n_counts"] * 100
        else:
            self.adata.obs["pct_mito"] = 0.0
        keep = (
            (self.adata.obs["n_genes"] >= self.MIN_GENES)
            & (self.adata.obs["n_genes"] <= self.MAX_GENES)
            & (self.adata.obs["n_counts"] >= self.MIN_COUNTS)
            & (self.adata.obs["pct_mito"] <= self.MAX_MITO_PCT)
        )
        n_before = self.adata.n_obs
        self.adata = self.adata[keep]
        n_after = self.adata.n_obs
        pct = (1 - n_after / max(n_before, 1)) * 100
        print(f"  QC: {n_before} -> {n_after} cells (removed {n_before - n_after}, {pct:.1f}%)")
        self.qc_report()
        return self

    def qc_report(self):
        """Print QC stats per time point."""
        if self.adata is None:
            return
        print("  QC stats by time point:")
        print(f"  {'TP':5s} {'cells':>6s} {'genes':>7s} {'counts':>8s} {'mito%':>6s}")
        for tp in self.TIME_ORDER:
            sub = self.adata[self.adata.obs["time_point"] == tp]
            n = sub.n_obs
            if n == 0:
                continue
            print(f"  {tp:5s} {n:6d} {sub.obs['n_genes'].median():7.0f} "
                  f"{sub.obs['n_counts'].median():8.0f} {sub.obs['pct_mito'].median():6.1f}")

    # ==================================================================
    # Step 4: Normalize
    # ==================================================================
    def normalize(self, target_sum=10000):
        print(f"Normalizing (target_sum={target_sum})...")
        sc.pp.normalize_total(self.adata, target_sum=target_sum)
        sc.pp.log1p(self.adata)
        self.adata.raw = self.adata
        print(f"  Done: mean expr = {self.adata.X.mean():.2f}")
        return self

    # ==================================================================
    # Step 5: HVG
    # ==================================================================
    def find_hvg(self, n_top_genes=3000, flavor="seurat"):
        print(f"Selecting HVGs (n={n_top_genes}, flavor={flavor})...")
        if flavor == "seurat_v3" and "counts" not in self.adata.layers:
            if hasattr(self.adata, "raw") and self.adata.raw is not None:
                self.adata.layers["counts"] = self.adata.raw.X.copy()
            else:
                print("  Warning: no counts, falling back to seurat")
                flavor = "seurat"
        sc.pp.highly_variable_genes(self.adata, n_top_genes=n_top_genes, flavor=flavor)
        n_hvg = self.adata.var["highly_variable"].sum()
        print(f"  {n_hvg} HVGs selected")
        return self

    # ==================================================================
    # Step 6: PCA
    # ==================================================================
    def run_pca(self, n_comps=50):
        print(f"Scaling + PCA (n_comps={n_comps})...")
        sc.pp.scale(self.adata, max_value=10)
        sc.tl.pca(self.adata, n_comps=n_comps, svd_solver="arpack")
        print(f"  PCA: {self.adata.obsm['X_pca'].shape}")
        return self

    # ==================================================================
    # Step 7: Neighbors + UMAP
    # ==================================================================
    def run_neighbors_umap(self, n_neighbors=15, n_pcs=30):
        print(f"Computing neighbors (k={n_neighbors}, PCs={n_pcs})...")
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        print("Computing UMAP...")
        sc.tl.umap(self.adata)
        print("  Done")
        return self

    # ==================================================================
    # Checkpoint I/O
    # ==================================================================
    def save_checkpoint(self, path=None):
        if path is None:
            path = self.CACHE_MERGED
        os.makedirs(Path(path).parent, exist_ok=True)
        self.adata.write_h5ad(path)
        print(f"Checkpoint saved: {path}")

    def load_checkpoint(self, path=None):
        if path is None:
            path = self.CACHE_MERGED
        if not Path(path).exists():
            print(f"No checkpoint at {path}")
            return False
        self.adata = ad.read_h5ad(path)
        print(f"Checkpoint loaded: {self.adata.n_obs} cells x {self.adata.n_vars} genes")
        return True

    # ==================================================================
    # Run helpers
    # ==================================================================
    def run_to_merged(self):
        if self.load_checkpoint(self.CACHE_MERGED):
            return self
        self.load_metadata()
        self.load_and_merge()
        return self

    def run_qc(self):
        self.qc_filter()
        return self

    def run_preprocessing(self, n_hvg=3000, n_pcs=50, n_neighbors=15):
        """Full preprocessing chain with checkpoints at each stage."""
        if self.adata is None:
            if not self.load_checkpoint(self.CACHE_MERGED):
                self.run_to_merged()
                self.save_checkpoint(self.CACHE_MERGED)

        # QC
        if self.CACHE_QC.exists():
            print("QC checkpoint found, loading...")
            self.load_checkpoint(self.CACHE_QC)
        else:
            self.qc_filter()
            self.save_checkpoint(self.CACHE_QC)

        # Normalize + HVG + PCA + UMAP
        if self.CACHE_NORM.exists():
            print("Norm checkpoint found, loading...")
            self.load_checkpoint(self.CACHE_NORM)
        else:
            self.normalize()
            self.find_hvg(n_top_genes=n_hvg)
            self.run_pca(n_comps=n_pcs)
            self.run_neighbors_umap(n_neighbors=n_neighbors)
            self.save_checkpoint(self.CACHE_NORM)

        # Subset to HVGs for final preprocessed
        if "highly_variable" in self.adata.var:
            self.adata = self.adata[:, self.adata.var["highly_variable"]].copy()
            print(f"Subset to HVGs: {self.adata.n_obs} cells x {self.adata.n_vars} genes")

        self.save_checkpoint(self.CACHE_PREPROC)
        print("Preprocessing complete.")
        return self


if __name__ == "__main__":
    pipe = ONCPipeline()
    pipe.run_preprocessing()
