# scTIME: 单细胞时间序列分析流程

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

一套模块化、支持断点续跑的单细胞 RNA-seq 时间序列分析流程。适用于多条件、多时间点的实验设计，内置质量控制、差异表达、基因轨迹聚类、细胞类型组成分析和通路富集功能。

## 特性

- **断点续跑**：每个主要步骤保存中间 .h5ad 文件，重复运行自动跳过已完成步骤
- **内存高效**：分块读取 CSV 并自动转换为稀疏矩阵
- **双重差异策略**：单细胞级 Wilcoxon（灵敏度高）+ 伪Bulk 级 t 检验（统计严谨）
- **时间序列感知**：按表达轨迹聚类基因，识别时间点间转换基因
- **细胞类型脆弱性**：组成比例偏移分析，含卡方检验
- **通路富集**：通过 gseapy 集成 Enrichr
- **出版级图表**：Nature 风格 PNG + PDF，附带源数据 CSV

## 分析流程

`
原始计数 → 质控过滤 → 归一化 → 高变基因 → PCA → UMAP
    │
    ├── 单细胞差异表达 (Wilcoxon)
    ├── 伪Bulk差异表达 (t检验)
    ├── 基因时间轨迹聚类 (KMeans)
    ├── 细胞类型比例分析
    └── GO通路富集 (Enrichr)
`

## 快速开始

### 环境依赖
`ash
pip install pandas numpy scipy anndata scanpy statsmodels gseapy scikit-learn seaborn matplotlib tqdm
`

### 运行
`ash
python scripts/pipeline.py          # 质控+预处理（支持续跑）
python scripts/deg_analysis.py      # 单细胞差异+伪Bulk
python scripts/pb_deg.py            # 伪Bulk t检验差异
python scripts/trajectory.py        # 基因聚类+转换分析
python scripts/subtype_analysis.py  # 细胞类型组成
python scripts/pathway_analysis.py  # GO富集
python scripts/make_figures.py      # 出版级图表
`

所有脚本均包含 `if __name__ == "__main__"` 入口，可独立运行。

## 仓库结构

`
├── scripts/                  # 分析脚本（可独立运行）
├── manuscript/               # 手稿（本地）
├── SKILL.md                  # 方法参考与最佳实践
├── PROGRESS.md               # 分析进度日志
└── README.md / README_zh.md  # 说明文档
`

## 输出结构

`
output/
├── merged.h5ad               # 原始合并计数
├── qc_filtered.h5ad          # 质控后
├── normalized_hvg.h5ad       # 归一化+高变基因
├── preprocessed.h5ad         # PCA+UMAP，高变基因子集
├── deg/                      # 差异表达表
├── trajectory/               # 基因聚类与转换
├── subtype/                  # 细胞类型组成分析
├── pathway/                  # GO富集结果
├── figures/                  # PNG+PDF图表
└── figure_data/              # 图表源数据
`

## 参数配置

关键参数在 `ONCPipeline` 类中作为类属性定义，运行前可调整：

- `MIN_GENES` / `MAX_GENES` — 基因数阈值
- `MAX_MITO_PCT` — 线粒体比例上限
- `n_top_genes` — 高变基因数量
- `n_clusters` — 基因轨迹聚类数

## 许可证

MIT
