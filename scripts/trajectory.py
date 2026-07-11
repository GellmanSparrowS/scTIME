"""
GSE137398 ONC RGC - Time-course trajectory and gene clustering.
"""
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from tqdm import tqdm
from sklearn.cluster import KMeans

warnings.filterwarnings("ignore")

CACHE = Path(r"F:\ONC\output\normalized_hvg.h5ad")
OUTPUT = Path(r"F:\ONC\output\trajectory")


class TimeCourseAnalysis:
    TIME_ORDER = ["Ctrl", "12h", "1d", "2d", "4d", "1w", "2w"]
    TIME_DAYS = {"Ctrl": 0, "12h": 0.5, "1d": 1, "2d": 2, "4d": 4, "1w": 7, "2w": 14}

    def __init__(self):
        OUTPUT.mkdir(parents=True, exist_ok=True)
        self.adata = None
        self.tc_mean = None

    def load(self):
        self.adata = ad.read_h5ad(CACHE)
        print(f"Loaded: {self.adata.n_obs} cells x {self.adata.n_vars} genes")
        return self

    def compute_time_means(self):
        print("Computing time-course mean expression...")
        tp_list = []
        for tp in self.TIME_ORDER:
            cells = self.adata.obs["time_point"] == tp
            if cells.sum() == 0:
                tp_list.append(np.zeros(self.adata.n_vars))
            else:
                tp_list.append(np.array(self.adata[cells].X.mean(axis=0)).flatten())
        self.tc_mean = pd.DataFrame(
            np.column_stack(tp_list),
            index=self.adata.var_names,
            columns=self.TIME_ORDER,
        )
        print(f"  Shape: {self.tc_mean.shape}")
        return self

    def cluster_genes(self, n_clusters=6):
        print(f"Clustering genes (n={n_clusters})...")
        tc_arr = self.tc_mean.values.astype(np.float64)
        tc_z = tc_arr - tc_arr.mean(axis=1, keepdims=True)
        tc_std = tc_arr.std(axis=1, keepdims=True)
        tc_std[tc_std == 0] = 1.0
        tc_z = tc_z / tc_std

        labels = KMeans(n_clusters=n_clusters, random_state=42, n_init=10).fit_predict(tc_z)
        self.tc_mean["cluster"] = labels.astype(str)

        # Reorder by peak time
        peak_order = {}
        for c in range(n_clusters):
            c_genes = self.tc_mean[self.tc_mean["cluster"] == str(c)]
            peak_tp = c_genes[self.TIME_ORDER].mean(axis=0).idxmax()
            peak_idx = self.TIME_ORDER.index(peak_tp)
            peak_order[c] = peak_idx
        sorted_clusters = sorted(peak_order.items(), key=lambda x: x[1])
        cluster_map = {str(old): str(new) for new, (old, _) in enumerate(sorted_clusters)}
        self.tc_mean["cluster"] = self.tc_mean["cluster"].map(cluster_map)

        for c in range(n_clusters):
            n = (self.tc_mean["cluster"] == str(c)).sum()
            c_genes = self.tc_mean[self.tc_mean["cluster"] == str(c)]
            peak_tp = c_genes[self.TIME_ORDER].mean(axis=0).idxmax()
            genes_top = c_genes.nlargest(3, peak_tp).index.tolist()
            print(f"  Cluster {c}: {n} genes, peak at {peak_tp} | top: {', '.join(genes_top)}")

        out = OUTPUT / "gene_clusters.csv"
        self.tc_mean.to_csv(out)
        print(f"  Saved: {out}")
        return self

    def find_markers_per_time(self, min_lfc=1.0, top_n=20):
        print("Finding time-point marker genes...")
        markers = {}
        for tp in self.TIME_ORDER[1:]:
            lfc = self.tc_mean[tp] - self.tc_mean["Ctrl"]
            sig = lfc[lfc.abs() >= min_lfc].sort_values(ascending=False)
            markers[tp] = sig.head(top_n)
            print(f"  {tp}: {len(sig)} genes with |LFC|>={min_lfc}")
        max_len = max(len(markers[tp]) for tp in self.TIME_ORDER[1:])
        marker_data = {}
        for tp in self.TIME_ORDER[1:]:
            genes = markers[tp].index.tolist()
            genes += [""] * (max_len - len(genes))
            marker_data[f"{tp}_gene"] = genes
        marker_df = pd.DataFrame(marker_data)
        out = OUTPUT / "time_markers.csv"
        marker_df.to_csv(out)
        return markers

    def compute_transition_genes(self):
        print("Computing transition genes...")
        transitions = {}
        for i in range(len(self.TIME_ORDER) - 1):
            t1, t2 = self.TIME_ORDER[i], self.TIME_ORDER[i + 1]
            delta = self.tc_mean[t2] - self.tc_mean[t1]
            top = delta.abs().nlargest(20)
            transitions[f"{t1}_to_{t2}"] = top
            print(f"  {t1}->{t2}: {top.index[0]} (delta={top.iloc[0]:.3f})")
        trans_df = pd.DataFrame({
            f"{k}_gene": transitions[k].index.tolist() for k in transitions
        })
        out = OUTPUT / "transition_genes.csv"
        trans_df.to_csv(out)
        return transitions


if __name__ == "__main__":
    tc = TimeCourseAnalysis()
    tc.load()
    tc.compute_time_means()
    tc.cluster_genes(n_clusters=6)
    tc.find_markers_per_time(min_lfc=1.0, top_n=20)
    tc.compute_transition_genes()
    print("\nTime-course analysis complete.")
