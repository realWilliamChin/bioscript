#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/06/29
# Author        : William GoGo
"""Pipeline 第 05 步: scanpy 组对差异基因 (rank_genes_groups, wilcoxon)

输入: 03/04 步后的 h5ad
输出:
    <out-dir>/<group1>_vs_<group2>_DEG.xlsx     每个组对的 DEG 列表
    <out-dir>/heatmap_<group1>_vs_<group2>.png  top50 基因热图(可选)

可指定组对列表(--pairs "A,B;C,D"), 默认对 group 列做两两组合。
"""

import os
import argparse
from itertools import combinations

import pandas as pd
import scanpy as sc
from loguru import logger


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input-h5ad', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--group-key', default='group')
    p.add_argument('--pairs', default=None,
                   help='指定组对, 形如 "HD,Class_III;HD,NoLN"; 不传则两两组合')
    p.add_argument('--padj-cutoff', type=float, default=0.05)
    p.add_argument('--logfc-cutoff', type=float, default=1.0)
    p.add_argument('--top-n-heatmap', type=int, default=50,
                   help='热图取每组对的 topN 显著基因')
    p.add_argument('--skip-heatmap', action='store_true')
    return p.parse_args()


def build_pairs(adata, group_key, pairs_arg):
    if pairs_arg:
        pairs = []
        for seg in pairs_arg.split(';'):
            a, b = seg.split(',')
            pairs.append((a.strip(), b.strip()))
        return pairs
    groups = list(adata.obs[group_key].unique())
    return list(combinations(groups, 2))


def run_one_pair(adata, group_key, ref, target, args, out_dir):
    logger.info(f'DEG: {ref} vs {target}')
    sc.tl.rank_genes_groups(adata, groupby=group_key,
                            reference=ref, groups=[target],
                            use_raw=True, pts=True, rankby_abs=True,
                            method='wilcoxon')
    df = sc.get.rank_genes_groups_df(adata, group=None)
    df['abs_logFC'] = df['logfoldchanges'].abs()
    pair_name = f'{ref.replace(" ", "")}_vs_{target.replace(" ", "")}'
    df.to_excel(os.path.join(out_dir, f'{pair_name}_DEG.xlsx'), index=False)

    if args.skip_heatmap:
        return df
    top = df.query(
        f'pvals_adj<{args.padj_cutoff} & logfoldchanges.abs()>={args.logfc_cutoff}'
    ).sort_values('abs_logFC', ascending=False)['names'].head(args.top_n_heatmap).tolist()
    if not top:
        logger.warning(f'{pair_name}: 无显著基因, 跳过热图')
        return df

    sub_cells = (adata.obs[group_key] == ref) | (adata.obs[group_key] == target)
    sub = adata[sub_cells, :].copy()
    sc.pl.heatmap(sub, top, use_raw=False, groupby=group_key, swap_axes=True,
                  vmin=-5, vmax=5, show=False,
                  save=f'_{pair_name}_top{args.top_n_heatmap}.png')
    return df


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir

    logger.info(f'读取 {args.input_h5ad}')
    adata = sc.read_h5ad(args.input_h5ad)
    pairs = build_pairs(adata, args.group_key, args.pairs)
    logger.info(f'共 {len(pairs)} 组对')

    for ref, target in pairs:
        run_one_pair(adata, args.group_key, ref, target, args, args.output_dir)

    logger.success(f'DEG 全部完成 → {args.output_dir}')


if __name__ == '__main__':
    main()
