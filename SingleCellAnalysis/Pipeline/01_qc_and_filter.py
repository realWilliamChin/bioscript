#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/06/29
# Author        : William GoGo
"""Pipeline 第 01 步: QC 指标计算 + 可视化 + 阈值过滤 + scrublet 去 doublet

输入: 00 步合并后的 h5ad
输出:
    <out-dir>/00_Raw_data_after_quality_control/
        QC_metrics_beforeQC_violin.png
        QC_metrics_beforeQC_scatter.png
        QC_metrics_afterQC_violin.png
        QC_metrics_afterQC_scatter.png
        QC_bySample_stat_df.xlsx
    <out-h5ad> : 过滤 + doublet 去除后的 h5ad

阈值通过命令行可配置, 不传则使用默认 (200<=n_genes<=7500, 500<=counts<=50000,
mt<=20%, hb<=0.5%)。
"""

import os
import argparse
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from loguru import logger

HB_GENES_DEFAULT = ['HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBZ', 'HBQ1']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input-h5ad', required=True, help='00 步合并的 h5ad')
    p.add_argument('--output-h5ad', required=True, help='过滤后的 h5ad')
    p.add_argument('--output-dir', required=True, help='QC 图与统计表输出目录')
    p.add_argument('--mt-prefix', default='MT-', help='线粒体基因前缀 (人=MT-, 鼠=mt-)')
    p.add_argument('--ribo-prefix', default='RPS,RPL', help='核糖体基因前缀, 逗号分隔')
    p.add_argument('--hb-genes', default=','.join(HB_GENES_DEFAULT),
                   help='血红蛋白基因列表, 逗号分隔')
    p.add_argument('--min-genes', type=int, default=200)
    p.add_argument('--max-genes', type=int, default=7500)
    p.add_argument('--min-counts', type=int, default=500)
    p.add_argument('--max-counts', type=int, default=50000)
    p.add_argument('--max-mt', type=float, default=20.0, help='线粒体比例上限 (%)')
    p.add_argument('--max-hb', type=float, default=0.5, help='血红蛋白比例上限 (%)')
    p.add_argument('--doublet-rate', type=float, default=0.06,
                   help='scrublet expected_doublet_rate')
    p.add_argument('--skip-doublet', action='store_true', help='跳过 scrublet')
    return p.parse_args()


def annotate_qc_genes(adata, mt_prefix, ribo_prefixes, hb_genes):
    adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
    adata.var['ribo'] = adata.var_names.str.startswith(tuple(ribo_prefixes))
    adata.var['hb'] = False
    hb_found = [g for g in hb_genes if g in adata.var_names]
    if hb_found:
        adata.var.loc[hb_found, 'hb'] = True
        logger.info(f'检测到血红蛋白基因: {hb_found}')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], inplace=True)


def plot_qc_violin_scatter(adata, out_dir, tag):
    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    for i, metric in enumerate(['n_genes_by_counts', 'total_counts']):
        sc.pl.violin(adata, metric, groupby='sample', rotation=45,
                     stripplot=False, ax=axes[0, i], show=False)
    for i, metric in enumerate(['pct_counts_mt', 'pct_counts_hb']):
        sc.pl.violin(adata, metric, groupby='sample', rotation=45,
                     stripplot=False, ax=axes[1, i], show=False)
    if tag == 'afterQC':
        axes[1, 0].set_ylim(0, 100)
        axes[1, 1].set_ylim(0, 100)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'QC_metrics_{tag}_violin.png'),
                bbox_inches='tight', dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(21, 5))
    sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt',
                  color='sampleID', ax=axes[0], show=False)
    axes[0].set_title('Genes vs MT%')
    sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_hb',
                  color='sampleID', ax=axes[1], show=False)
    axes[1].set_title('Genes vs HB%')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
                  color='sampleID', ax=axes[2], show=False)
    axes[2].set_title('UMI vs Genes')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'QC_metrics_{tag}_scatter.png'),
                bbox_inches='tight', dpi=300)
    plt.close(fig)


def filter_cells(adata, args):
    mask = (
        (adata.obs.n_genes_by_counts >= args.min_genes) &
        (adata.obs.n_genes_by_counts <= args.max_genes) &
        (adata.obs.total_counts >= args.min_counts) &
        (adata.obs.total_counts <= args.max_counts) &
        (adata.obs.pct_counts_mt <= args.max_mt) &
        (adata.obs.pct_counts_hb <= args.max_hb)
    )
    return adata[mask, :].copy()


def run_scrublet_per_sample(adata, doublet_rate):
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata_hvg = adata[:, adata.var.highly_variable].copy()

    adata_hvg.obs['doublet_score'] = np.nan
    adata_hvg.obs['predicted_doublet'] = False
    for s in adata_hvg.obs['sample'].unique():
        idx = adata_hvg.obs['sample'] == s
        tmp = adata_hvg[idx].copy()
        sc.pp.scrublet(tmp, expected_doublet_rate=doublet_rate)
        adata_hvg.obs.loc[idx, 'doublet_score'] = tmp.obs['doublet_score']
        adata_hvg.obs.loc[idx, 'predicted_doublet'] = tmp.obs['predicted_doublet']

    adata.obs['doublet_score'] = adata_hvg.obs['doublet_score']
    adata.obs['predicted_doublet'] = adata_hvg.obs['predicted_doublet']
    keep = ~adata.obs['predicted_doublet'].fillna(False).astype(bool).values
    return adata[keep].copy()


def summarize_qc(adata_merged, adata_filtered, adata_doublets_rmv, args, out_dir):
    df = pd.DataFrame(adata_merged.obs.groupby('sampleID')['sampleID'].count())
    df.columns = ['rawCellsCount']
    df[f'pct_counts_mt>={args.max_mt}'] = adata_merged.obs.query(
        f'pct_counts_mt>={args.max_mt}').groupby('sampleID')['sampleID'].count()
    df[f'pct_counts_hb>={args.max_hb}'] = adata_merged.obs.query(
        f'pct_counts_hb>={args.max_hb}').groupby('sampleID')['sampleID'].count()
    if 'predicted_doublet' in adata_filtered.obs.columns:
        df['doublets'] = adata_filtered.obs.query(
            'predicted_doublet==True').groupby('sampleID')['sampleID'].count()
    df[f'total_counts<{args.min_counts} or >{args.max_counts}'] = adata_merged.obs.query(
        f'total_counts<{args.min_counts} | total_counts>{args.max_counts}'
    ).groupby('sampleID')['sampleID'].count()
    df[f'n_genes_by_counts<{args.min_genes} or >{args.max_genes}'] = adata_merged.obs.query(
        f'n_genes_by_counts<{args.min_genes} | n_genes_by_counts>{args.max_genes}'
    ).groupby('sampleID')['sampleID'].count()
    df['filteredCellsCount'] = adata_doublets_rmv.obs.groupby('sampleID')['sampleID'].count()
    df['filteredCellsFraction(%)'] = (100 * df['filteredCellsCount']
                                      / df['rawCellsCount']).round(2)
    df.to_excel(os.path.join(out_dir, 'QC_bySample_stat_df.xlsx'))
    return df


def main():
    args = parse_args()
    qc_dir = os.path.join(args.output_dir, '00_Raw_data_after_quality_control')
    os.makedirs(qc_dir, exist_ok=True)

    logger.info(f'读取 {args.input_h5ad}')
    adata = sc.read_h5ad(args.input_h5ad)

    ribo_prefixes = tuple(p.strip() for p in args.ribo_prefix.split(','))
    hb_genes = [g.strip() for g in args.hb_genes.split(',')]
    annotate_qc_genes(adata, args.mt_prefix, ribo_prefixes, hb_genes)

    plot_qc_violin_scatter(adata, qc_dir, tag='beforeQC')

    logger.info('阈值过滤 ...')
    adata_filtered = filter_cells(adata, args)
    logger.info(f'阈值过滤后: {adata_filtered.n_obs} cells')

    if args.skip_doublet:
        adata_final = adata_filtered
    else:
        logger.info('scrublet 去 doublet ...')
        adata_final = run_scrublet_per_sample(adata_filtered, args.doublet_rate)
        logger.info(f'去 doublet 后: {adata_final.n_obs} cells')

    plot_qc_violin_scatter(adata_final, qc_dir, tag='afterQC')

    summarize_qc(adata, adata_filtered, adata_final, args, qc_dir)

    adata_final.write_h5ad(args.output_h5ad)
    logger.success(f'过滤完成 → {args.output_h5ad}')


if __name__ == '__main__':
    main()
