#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/06/29
# Author        : William GoGo
"""按 cluster 计算基因表达均值 → Pearson 相关性矩阵 → 热图

输入:
    --h5ad        AnnData 文件 (h5ad)
    --cluster-key obs 中的聚类列名 (例: leiden_r025, bm_predLables)
    --output-dir  输出目录,产出 *_corr_matrix.xlsx + *_corr_heatmap.png

通常在 Pipeline/03 注释之后调用,既可基于 leiden cluster 也可基于 celltypist 注释。
"""

import os
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from loguru import logger


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--h5ad', required=True, help='输入 h5ad 文件')
    p.add_argument('--cluster-key', required=True, help='obs 中作为分组的列名')
    p.add_argument('--use-raw', action='store_true', default=True,
                   help='使用 adata.raw 中的表达矩阵 (默认 True)')
    p.add_argument('--no-raw', dest='use_raw', action='store_false',
                   help='不使用 raw, 改用 adata.X')
    p.add_argument('--output-dir', required=True, help='输出目录')
    p.add_argument('--prefix', default='cluster_corr', help='输出文件前缀')
    p.add_argument('--annot-fontsize', type=int, default=8, help='热图数值字体大小')
    return p.parse_args()


def calculate_cluster_correlation(adata, cluster_key, use_raw=True):
    """按 cluster 计算基因表达均值, 返回相关性矩阵与均值矩阵。"""
    if use_raw and adata.raw is not None:
        expr_matrix = adata.raw.X.T
        var_names = adata.raw.var_names
    else:
        expr_matrix = adata.X.T
        var_names = adata.var_names

    expr_df = pd.DataFrame(
        expr_matrix.toarray() if hasattr(expr_matrix, 'toarray') else expr_matrix,
        index=var_names,
        columns=adata.obs_names,
    )

    cluster_means = []
    clusters = adata.obs[cluster_key].unique()
    for cluster in clusters:
        cluster_cells = adata.obs[adata.obs[cluster_key] == cluster].index
        cluster_means.append(expr_df[cluster_cells].mean(axis=1))

    cluster_expr_df = pd.concat(cluster_means, axis=1)
    cluster_expr_df.columns = clusters

    corr_matrix = cluster_expr_df.corr(method='pearson')
    return corr_matrix, cluster_expr_df


def plot_correlation_heatmap(corr_matrix, out_file, title='Cluster Correlation Heatmap',
                             figsize=(16, 12), cmap='coolwarm', annot_fontsize=8):
    plt.figure(figsize=figsize)
    sns.heatmap(corr_matrix, cmap=cmap, annot=True,
                annot_kws={'size': annot_fontsize}, fmt='.2f',
                square=True, linewidths=0.5,
                cbar_kws={'label': 'Pearson correlation coefficient'})
    plt.title(title, fontsize=16, pad=20)
    plt.tight_layout()
    plt.savefig(out_file, bbox_inches='tight', dpi=300)
    plt.close()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    logger.info(f'读取 {args.h5ad}')
    adata = sc.read_h5ad(args.h5ad)
    if args.cluster_key not in adata.obs.columns:
        raise KeyError(f'obs 中不存在列 {args.cluster_key}')

    corr_matrix, cluster_expr_df = calculate_cluster_correlation(
        adata, args.cluster_key, use_raw=args.use_raw)

    matrix_file = os.path.join(args.output_dir, f'{args.prefix}_matrix.xlsx')
    corr_matrix.to_excel(matrix_file)
    logger.success(f'相关性矩阵 → {matrix_file}')

    png_file = os.path.join(args.output_dir, f'{args.prefix}_heatmap.png')
    plot_correlation_heatmap(corr_matrix, png_file,
                             title=f'Cluster Correlation ({args.cluster_key})',
                             annot_fontsize=args.annot_fontsize)
    logger.success(f'相关性热图 → {png_file}')


if __name__ == '__main__':
    main()
