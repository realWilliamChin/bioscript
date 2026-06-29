#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/06/29
# Author        : William GoGo
"""组间细胞比例的 Kruskal-Wallis + Dunn post-hoc 检验

输入:
    --prop-tsv   行=cluster, 列=sample(%) 的占比表 (Pipeline/04 输出)
    --meta-tsv   样本→组 映射, 两列 sampleID/group
    --output-dir 输出目录
                  ├── cluster_kw_stat.xlsx        每个 cluster 的 KW 统计 + p
                  └── cluster<id>_dunn.xlsx      显著 cluster 的 Dunn 矩阵

KW < 0.05 的 cluster 才写 Dunn 矩阵。
"""

import os
import argparse

import pandas as pd
from scipy import stats
import scikit_posthocs as sp
from loguru import logger


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--prop-tsv', required=True,
                   help='cluster×sample 占比表, 第一列为 clusterID, 其余为样本占比(%)')
    p.add_argument('--meta-tsv', required=True,
                   help='样本→组 映射, 两列 sampleID/group, tab 分隔')
    p.add_argument('--output-dir', required=True)
    p.add_argument('--sig-cutoff', type=float, default=0.05, help='KW p 显著阈值')
    return p.parse_args()


def load_prop_table(prop_tsv):
    df = pd.read_csv(prop_tsv, sep=None, engine='python', index_col=0)
    # 去掉列名中的 (%) 后缀, 统一对齐 sampleID
    df.columns = [c.split('(')[0].strip() for c in df.columns]
    return df


def load_sample_meta(meta_tsv):
    df = pd.read_csv(meta_tsv, sep=None, engine='python')
    if 'sampleID' not in df.columns or 'group' not in df.columns:
        raise ValueError('meta-tsv 必须包含 sampleID 与 group 两列')
    return dict(zip(df['sampleID'].astype(str), df['group']))


def kw_dunn_one_cluster(prop_series, sample2group):
    """对单个 cluster 的 sample 占比做 KW 检验; 显著则返回 Dunn 矩阵, 否则 None"""
    grp2vals = {}
    for sample, val in prop_series.items():
        grp = sample2group.get(str(sample))
        if grp is None:
            continue
        grp2vals.setdefault(grp, []).append(val)
    if len(grp2vals) < 2:
        return None, None, None

    kw_stat, kw_p = stats.kruskal(*grp2vals.values())

    dfs = []
    for grp, vals in grp2vals.items():
        dfs.append(pd.DataFrame({'value': vals, 'group': grp}))
    long_df = pd.concat(dfs, ignore_index=True)
    dunn_df = sp.posthoc_dunn(long_df, val_col='value', group_col='group',
                              p_adjust='holm')
    return kw_stat, kw_p, dunn_df


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    prop_df = load_prop_table(args.prop_tsv)
    sample2group = load_sample_meta(args.meta_tsv)

    kw_records = []
    for cluster in prop_df.index:
        kw_stat, kw_p, dunn_df = kw_dunn_one_cluster(prop_df.loc[cluster], sample2group)
        if kw_p is None:
            continue
        kw_records.append({'clusterID': cluster, 'stat': kw_stat, 'p_value': kw_p})
        if kw_p < args.sig_cutoff and dunn_df is not None:
            dunn_file = os.path.join(args.output_dir, f'Cluster{cluster}_dunn.xlsx')
            dunn_df.to_excel(dunn_file)
            logger.info(f'cluster {cluster}: KW p={kw_p:.3e} 显著, Dunn → {dunn_file}')

    kw_df = pd.DataFrame(kw_records)
    kw_file = os.path.join(args.output_dir, 'cluster_kw_stat.xlsx')
    kw_df.to_excel(kw_file, index=False)
    logger.success(f'KW 统计 → {kw_file}')


if __name__ == '__main__':
    main()
