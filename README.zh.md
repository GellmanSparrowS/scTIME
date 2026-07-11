# scTIME

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

一套面向单细胞 RNA-seq 时间序列实验的模块化、断点续跑分析框架。

`mermaid
flowchart TD
    A[原始计数] --> B[质控 + 归一化]
    B --> C[高变基因 + PCA + UMAP]
    C --> D1[单细胞差异]
    C --> D2[伪Bulk差异]
    C --> D3[基因聚类]
    C --> D4[亚型组成]
    C --> D5[通路富集]
    D1 --> E[图表输出]
    D2 --> E
    D3 --> E
    D4 --> E
    D5 --> E
`

## 特性

与大多数将每个时间点独立处理的流程不同，scTIME 专为时间序列设计，核心差异在于：

- **双重差异策略** — 单细胞级 Wilcoxon（灵敏度）+ 伪Bulk级 t 检验与 BH 校正（统计严谨性），消除伪重复问题
- **时间序列感知** — 全基因组表达轨迹 k 均值聚类；相邻时间点间转换基因鉴定
- **细胞类型脆弱性评分** — 基于倍变化的定量分类（脆弱/稳定/韧性），结合卡方检验
- **时间分辨 GO 富集** — 每个时间点独立进行 Enrichr 富集，而非仅比较汇总条件
- **断点续跑** — 每个步骤保存中间 .h5ad 文件，重复运行自动跳过已完成步骤
- **内存高效** — 分块 CSV 读取、自动稀疏矩阵转换、float32 存储
- **出版级图表** — Nature 风格 PNG + PDF 图表，附带对应 CSV 源数据

## 快速开始

`ash
pip install pandas numpy scipy anndata scanpy statsmodels gseapy scikit-learn seaborn matplotlib tqdm

python scripts/pipeline.py          # 合并 + 质控 + 归一化 + HVG + PCA + UMAP
python scripts/deg_analysis.py      # 单细胞差异 + 伪Bulk构建
python scripts/pb_deg.py            # 伪Bulk t检验差异
python scripts/trajectory.py        # 基因共表达聚类 + 转换分析
python scripts/subtype_analysis.py  # 细胞类型组成 + 脆弱性评分
python scripts/pathway_analysis.py  # 时间点 GO 富集
python scripts/make_figures.py      # 出版图表
`

## 参数配置

所有参数均为类属性，运行前可调整：

`python
MIN_GENES = 200        # 每细胞最小基因数
MAX_GENES = 6000       # 每细胞最大基因数
MIN_COUNTS = 500       # 最小UMI计数
MAX_MITO_PCT = 20.0    # 最大线粒体比例
n_top_genes = 3000     # 高变基因数
n_comps = 50           # PCA成分数
n_clusters = 6         # 基因轨迹聚类数
`

## 许可证

MIT
