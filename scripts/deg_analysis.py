"""
GSE137398 ONC RGC - Differential Expression Analysis
Uses preprocessed AnnData checkpoint. Pseudobulk + wilcoxon per time point.
"""
import os
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from tqdm import tqdm

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


class DEGAnalysis:
    """Differential expression for ONC time-course."""

    TIME_POINTS = ["Ctrl", "12h", "1d", "2d", "4d", "1w", "2w"]
    OUTPUT_DIR = Path(r"F:\ONC\output\deg")
    CACHE_PREPROC = Path(r"F:\ONC\output\preprocessed.h5ad")

    def __init__(self):
        os.makedirs(self.OUTPUT_DIR, exist_ok=True)
        self.adata = None

    def load(self, path=None):
        if path is None:
            path = self.CACHE_PREPROC
        self.adata = ad.read_h5ad(path)
        print(f"Loaded: {self.adata.n_obs} cells x {self.adata.n_vars} genes")
        return self

    def run_single_de(self, case_tp, ctrl_tp="Ctrl", method="wilcoxon"):
        """Run DEG: case vs control for a single time point comparison."""
        key = f"{case_tp}_vs_{ctrl_tp}"
        print(f"\n=== {key} DEG ===")

        sub = self.adata[self.adata.obs["time_point"].isin([ctrl_tp, case_tp])].copy()
        sub.obs["group"] = sub.obs["time_point"].astype(str)
        n_ctrl = (sub.obs["group"] == ctrl_tp).sum()
        n_case = (sub.obs["group"] == case_tp).sum()
        print(f"  {ctrl_tp}: {n_ctrl} cells, {case_tp}: {n_case} cells")

        sc.tl.rank_genes_groups(
            sub, groupby="group", groups=[case_tp], reference=ctrl_tp,
            method=method, n_genes=sub.n_vars, corr_method="benjamini-hochberg",
        )

        result = self._extract_result(sub, key)
        out_path = self.OUTPUT_DIR / f"{key}_deg.csv"
        result.to_csv(out_path)
        n_sig = (result["pvals_adj"] < 0.05).sum()
        n_up = (result["logfoldchanges"] > 0).sum()
        n_down = (result["logfoldchanges"] < 0).sum()
        print(f"  DEGs: {len(result)} total, {n_sig} significant (FDR<0.05)")
        print(f"  Up: {n_up}, Down: {n_down}")
        print(f"  Top 10 up: {', '.join(result.head(10).index.tolist())}")
        print(f"  Saved: {out_path}")
        return result

    def _extract_result(self, adata, key):
        """Extract DEG results from scanpy rank_genes_groups output."""
        result = adata.uns["rank_genes_groups"]
        groups = result["names"].dtype.names
        dfs = []
        for group in groups:
            df = pd.DataFrame({
                "gene": result["names"][group],
                "logfoldchanges": result["logfoldchanges"][group],
                "pvals": result["pvals"][group],
                "pvals_adj": result["pvals_adj"][group],
                "scores": result["scores"][group],
            })
            df = df.set_index("gene")
            dfs.append(df)
        return pd.concat(dfs) if len(dfs) > 1 else dfs[0]

    def run_all_degs(self):
        """Run DEG for all 6 crush time points vs Ctrl."""
        all_results = {}
        for tp in self.TIME_POINTS[1:]:  # skip Ctrl
            result = self.run_single_de(tp, "Ctrl")
            all_results[tp] = result
        return all_results

    def build_pseudobulk(self):
        """Aggregate to pseudobulk per sample for DESeq2-style analysis."""
        print("Building pseudobulk matrix...")
        # Use raw counts from original data
        raw = ad.read_h5ad(Path(r"F:\ONC\output\merged_onc.h5ad"))
        raw = raw[self.adata.obs_names].copy()

        samples = raw.obs["sample_id"].unique()
        pb_list = []
        pb_obs = []
        for s in tqdm(samples):
            cells = raw.obs["sample_id"] == s
            if cells.sum() == 0:
                continue
            pb_vec = np.array(raw[cells].X.sum(axis=0)).flatten()
            pb_list.append(pb_vec)
            pb_obs.append({
                "sample_id": s,
                "time_point": raw.obs.loc[cells, "time_point"].iloc[0],
            })
        pb_mat = np.vstack(pb_list)
        pb_adata = ad.AnnData(
            X=pb_mat, obs=pd.DataFrame(pb_obs), var=raw.var.copy()
        )
        pb_adata.obs["time_point"] = pd.Categorical(
            pb_adata.obs["time_point"], categories=self.TIME_POINTS, ordered=True
        )
        out_path = self.OUTPUT_DIR / "pseudobulk.h5ad"
        pb_adata.write_h5ad(out_path)
        print(f"Pseudobulk: {pb_adata.n_obs} samples x {pb_adata.n_vars} genes")
        print(f"Saved: {out_path}")
        return pb_adata


if __name__ == "__main__":
    deg = DEGAnalysis()
    deg.load()
    # Pilot: A02 1d vs Ctrl
    deg.run_single_de("1d", "Ctrl")
    # Run all
    deg.run_all_degs()
    # Pseudobulk matrix for later use
    deg.build_pseudobulk()
    print("\nAll DEG analyses complete.")
