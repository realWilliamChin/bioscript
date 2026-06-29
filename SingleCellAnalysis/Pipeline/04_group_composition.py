#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/06/29
# Author        : William GoGo
"""Pipeline 第 04 步: 样本/组 × cluster 组成统计 + 组间比例检验

输入: 03 步注释后的 h5ad
输出:
    <out-dir>/count/group_<cluster_key>_clusterCellsCountAnno.xlsx
    <out-dir>/count/sample_<cluster_key>_clusterCellsStatAnno.csv
    <out-dir>/figures/cellsamples_bars.png   样本×cluster 堆叠柱状图
    <out-dir>/stat/cluster_kw_stat.xlsx + Cluster<id>_dunn.xlsx
        (调用 utils/prop_kw_dunn.py 的逻辑做组间 KW + Dunn)
"""

import os
import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from loguru import logger

# 复用 utils
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, '..', 'utils'))
from prop_kw_dunn import kw_dunn_one_cluster  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input-h5ad', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--cluster-key', default='leiden_r025',
                   help='obs 中的 cluster 列名')
    p.add_argument('--group-key', default='group')
    p.add_argument('--sample-key', default='sampleID')
    p.add_argument('--marker-xlsx', default=None,
                   help='可选,02 步导出的 marker 表(用于添加 Core_Markers 列)')
    p.add_argument('--top-markers', type=int, default=5)
    p.add_argument('--sig-cutoff', type=float, default=0.05)
    return p.parse_args()


def load_top_markers(marker_xlsx, top_n):
    if not marker_xlsx or not os.path.exists(marker_xlsx):
        return {}
    df = pd.read_excel(marker_xlsx)
    cluster_col = 'clusterID' if 'clusterID' in df.columns else df.columns[0]
    gene_col = 'genes' if 'genes' in df.columns else 'names'
    out = {}
    for cluster, sub in df.groupby(cluster_col):
        out[str(cluster)] = sub.head(top_n)[gene_col].tolist()
    return out


def group_x_cluster(adata, group_key, cluster_key, marker_map):
    cnt = pd.crosstab(adata.obs[group_key], adata.obs[cluster_key]).T
    total = cnt.sum().sum()
    prop = cnt.div(cnt.sum(axis=0), axis=1).round(4) * 100

    cnt = cnt.reset_index().rename(columns={cluster_key: 'clusterID'}).set_index('clusterID')
    prop = prop.reset_index().rename(columns={cluster_key: 'clusterID'}).set_index('clusterID')

    cnt['Cell_Count'] = cnt.sum(axis=1)
    cnt['Cell_Percentage'] = (100 * cnt['Cell_Count'] / total).round(2)
    if marker_map:
        cnt['Core_Markers'] = cnt.index.map(
            lambda c: ','.join(marker_map.get(str(c), [])))

    prop.columns = [f'{c}(%)' for c in prop.columns]
    return pd.concat([cnt, prop], axis=1)


def sample_x_cluster(adata, sample_key, cluster_key, marker_map, out_dir):
    cnt = pd.crosstab(adata.obs[sample_key], adata.obs[cluster_key]).T
    total = cnt.sum().sum()
    prop = cnt.div(cnt.sum(axis=0), axis=1)

    # 堆叠柱状图
    fig_dir = os.path.join(out_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    ax = prop.T.plot(kind='bar', stacked=True, colormap='tab20',
                     figsize=(10, 6))
    ax.set_title('Cluster composition of samples', fontweight='bold')
    ax.set_xlabel('sample')
    ax.set_ylabel('Proportion')
    ax.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45)
    plt.savefig(os.path.join(fig_dir, 'cellsamples_bars.png'), bbox_inches='tight', dpi=300)
    plt.close()

    cnt = cnt.reset_index().rename(columns={cluster_key: 'clusterID'}).set_index('clusterID')
    cnt['Cell_Count'] = cnt.sum(axis=1)
    cnt['Cell_Percentage'] = (100 * cnt['Cell_Count'] / total).round(2)
    if marker_map:
        cnt['Core_Markers'] = cnt.index.map(
            lambda c: ','.join(marker_map.get(str(c), [])))
    prop_pct = (prop * 100).round(2)
    prop_pct.columns = [f'{c}(%)' for c in prop_pct.columns]
    prop_pct = prop_pct.reset_index().rename(columns={cluster_key: 'clusterID'}).set_index('clusterID')
    return pd.concat([cnt, prop_pct], axis=1), prop_pct


def kw_dunn_across_groups(prop_pct, sample2group, out_dir, sig_cutoff):
    stat_dir = os.path.join(out_dir, 'stat')
    os.makedirs(stat_dir, exist_ok=True)
    records = []
    # prop_pct 列形如 SLE1(%), 行为 cluster
    prop = prop_pct.copy()
    prop.columns = [c.split('(')[0] for c in prop.columns]
    for cluster in prop.index:
        kw_stat, kw_p, dunn_df = kw_dunn_one_cluster(prop.loc[cluster], sample2group)
        if kw_p is None:
            continue
        records.append({'clusterID': cluster, 'stat': kw_stat, 'p_value': kw_p})
        if kw_p < sig_cutoff and dunn_df is not None:
            dunn_df.to_excel(os.path.join(stat_dir, f'Cluster{cluster}_dunn.xlsx'))
    pd.DataFrame(records).to_excel(
        os.path.join(stat_dir, 'cluster_kw_stat.xlsx'), index=False)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    count_dir = os.path.join(args.output_dir, 'count')
    os.makedirs(count_dir, exist_ok=True)

    logger.info(f'读取 {args.input_h5ad}')
    adata = sc.read_h5ad(args.input_h5ad)
    marker_map = load_top_markers(args.marker_xlsx, args.top_markers)

    grp_df = group_x_cluster(adata, args.group_key, args.cluster_key, marker_map)
    grp_df.to_excel(os.path.join(count_dir,
                                 f'group_{args.cluster_key}_clusterCellsCountAnno.xlsx'))
    logger.success('组 × cluster 统计完成')

    sample_df, prop_pct = sample_x_cluster(adata, args.sample_key,
                                           args.cluster_key, marker_map,
                                           args.output_dir)
    sample_df.to_csv(os.path.join(count_dir,
                                  f'sample_{args.cluster_key}_clusterCellsStatAnno.csv'))
    logger.success('样本 × cluster 统计完成')

    sample_meta = adata.obs[[args.sample_key, args.group_key]].drop_duplicates()
    sample2group = dict(zip(sample_meta[args.sample_key].astype(str),
                            sample_meta[args.group_key]))
    kw_dunn_across_groups(prop_pct, sample2group, args.output_dir, args.sig_cutoff)
    logger.success('组间比例 KW + Dunn 完成')


if __name__ == '__main__':
    main()
