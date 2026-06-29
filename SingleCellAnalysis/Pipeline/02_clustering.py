#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/06/29
# Author        : William GoGo
"""Pipeline 第 02 步: HVG + PCA + (可选 Harmony) + 多分辨率 leiden 聚类

输入: 01 步 QC 后的 h5ad
输出:
    <out-dir>/01_批次校正/                批次校正前后 UMAP
    <out-dir>/02_细胞聚类/
        count/                              每个分辨率的 cluster 细胞数
        figures/                            每个分辨率的 UMAP + DotPlot + rank_genes
        marker_gene/                        选定分辨率的 marker 表
    <out-h5ad>                              带 X_pca/X_pca_harmony/leiden_r0XX 的 h5ad
"""

import os
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
from loguru import logger

try:
    from harmonypy import run_harmony
    HARMONY_OK = True
except ImportError:
    HARMONY_OK = False


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input-h5ad', required=True)
    p.add_argument('--output-h5ad', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--n-top-genes', type=int, default=2000)
    p.add_argument('--n-pcs', type=int, default=50)
    p.add_argument('--batch-key', default='group',
                   help='Harmony 校正变量, 不做校正传空字符串')
    p.add_argument('--resolutions', default='0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50',
                   help='leiden 分辨率列表 (逗号分隔)')
    p.add_argument('--force-markers', default=None,
                   help='强制标记为 highly_variable 的基因列表(逗号分隔), 例: marker rescue')
    p.add_argument('--selected-resolution', default=None,
                   help='可选,指定一个分辨率写出 marker 表(例 0.25 → leiden_r025)')
    return p.parse_args()


def prepare_hvg(adata, n_top_genes, force_markers=None):
    if 'log1p' not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes,
                                min_mean=0.01, max_mean=5, min_disp=0.5,
                                flavor='seurat')
    if force_markers:
        for g in force_markers:
            if g in adata.var_names:
                adata.var.loc[g, 'highly_variable'] = True


def run_pca_and_neighbors(adata, n_pcs):
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack', use_highly_variable=True)


def harmony_correct(adata, batch_key):
    if not HARMONY_OK:
        raise ImportError('harmonypy 未安装, 无法做 Harmony 批次校正')
    pca_matrix = adata.obsm['X_pca']
    harmony_obj = run_harmony(data_mat=pca_matrix.T,
                              meta_data=adata.obs,
                              vars_use=[batch_key],
                              nclust=10, max_iter_harmony=10, verbose=True)
    adata.obsm['X_pca_harmony'] = harmony_obj.Z_corr.T


def umap_noBC_vs_BC(adata, out_dir, batch_key):
    os.makedirs(out_dir, exist_ok=True)
    sc.pp.neighbors(adata, key_added='noBC')
    sc.tl.leiden(adata, neighbors_key='noBC', key_added='noBC_leiden')
    sc.tl.umap(adata, neighbors_key='noBC')
    adata.obsm['X_umap_noBC'] = adata.obsm['X_umap'].copy()

    sc.pl.embedding(adata, basis='umap_noBC',
                    color=['noBC_leiden', batch_key],
                    title='UMAP (no Batch Correction)',
                    save=False, show=False)
    import matplotlib.pyplot as plt
    plt.savefig(os.path.join(out_dir, 'umap_noBC.png'), bbox_inches='tight', dpi=300)
    plt.close()

    if 'X_pca_harmony' in adata.obsm:
        sc.pp.neighbors(adata, key_added='BC', use_rep='X_pca_harmony')
        sc.tl.leiden(adata, neighbors_key='BC', key_added='BC_leiden')
        sc.tl.umap(adata, neighbors_key='BC')
        adata.obsm['X_umap_BC'] = adata.obsm['X_umap'].copy()
        sc.pl.embedding(adata, basis='umap_BC',
                        color=['BC_leiden', batch_key],
                        title='UMAP (Batch Correction)', save=False, show=False)
        plt.savefig(os.path.join(out_dir, 'umap_BC.png'), bbox_inches='tight', dpi=300)
        plt.close()


def multi_resolution_cluster(adata, resolutions, out_dir):
    count_dir = os.path.join(out_dir, 'count')
    fig_dir = os.path.join(out_dir, 'figures')
    os.makedirs(count_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    # 用 (可能存在的) BC 邻居图; 没有就重新跑一次
    use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.umap(adata)

    res2key = {}
    for r in resolutions:
        key = f'leiden_r0{int(r*100):02d}'
        sc.tl.leiden(adata, resolution=r, key_added=key)
        res2key[r] = key

    summary = []
    for r, key in res2key.items():
        cnt = adata.obs.groupby(key)[key].count()
        cnt.to_frame('cellsCount').reset_index().rename(
            columns={key: 'clusterID'}).to_csv(
            os.path.join(count_dir, f'resolution{int(r*100):02d}_clusterCellCount.csv'),
            index=False)
        summary.append({'resolution': f'{r:.2f}', 'clusterCount': cnt.shape[0]})
    pd.DataFrame(summary).to_csv(
        os.path.join(count_dir, 'resolution_allCount.csv'), index=False)

    import matplotlib.pyplot as plt
    sc.settings.figdir = fig_dir
    for key in res2key.values():
        sc.pl.umap(adata, color=[key], title='UMAP (Clusters)',
                   save=f'_{key}.png', show=False)
    return res2key


def export_markers(adata, resolution_key, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    sc.tl.rank_genes_groups(adata, groupby=resolution_key,
                            use_raw=False, pts=True, method='wilcoxon')
    df = sc.get.rank_genes_groups_df(adata, group=None)
    df.rename(columns={'group': 'clusterID', 'names': 'genes'}, inplace=True)
    df.to_excel(os.path.join(out_dir, f'{resolution_key}_all.xlsx'), index=False)
    for cluster in df['clusterID'].unique():
        sub = df[df['clusterID'] == cluster]
        sub.head(5).to_csv(os.path.join(out_dir, f'{resolution_key}_Cluster{cluster}_top5.csv'),
                           index=False)
        sub['genes'].head(200).to_csv(
            os.path.join(out_dir, f'{resolution_key}_Cluster{cluster}_top200.csv'),
            index=False)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    logger.info(f'读取 {args.input_h5ad}')
    adata = sc.read_h5ad(args.input_h5ad)

    force_markers = ([g.strip() for g in args.force_markers.split(',')]
                     if args.force_markers else None)
    prepare_hvg(adata, n_top_genes=args.n_top_genes, force_markers=force_markers)

    run_pca_and_neighbors(adata, n_pcs=args.n_pcs)

    batch_dir = os.path.join(args.output_dir, '01_批次校正')
    if args.batch_key:
        logger.info(f'Harmony 校正 batch_key={args.batch_key}')
        harmony_correct(adata, args.batch_key)
    umap_noBC_vs_BC(adata, batch_dir, batch_key=args.batch_key or 'sample')

    cluster_dir = os.path.join(args.output_dir, '02_细胞聚类')
    resolutions = [float(r) for r in args.resolutions.split(',')]
    res2key = multi_resolution_cluster(adata, resolutions, cluster_dir)

    if args.selected_resolution:
        r = float(args.selected_resolution)
        key = res2key.get(r) or f'leiden_r0{int(r*100):02d}'
        logger.info(f'写出 {key} 的 marker 表 ...')
        export_markers(adata, key, os.path.join(cluster_dir, 'marker_gene'))

    adata.write_h5ad(args.output_h5ad)
    logger.success(f'聚类完成 → {args.output_h5ad}')


if __name__ == '__main__':
    main()
