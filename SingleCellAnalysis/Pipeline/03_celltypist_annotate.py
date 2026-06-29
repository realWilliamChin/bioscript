#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/06/29
# Author        : William GoGo
"""Pipeline 第 03 步: celltypist 细胞类型注释 (best-match + majority-voting 双模式)

输入: 02 步聚类后的 h5ad
输出:
    <out-dir>/celltypist_cellAnnotStat.xlsx     注释结果汇总
    <out-h5ad>                                   obs 中新增
        - bm_predLables, bm_conf_score          best-match
        - mv_predicted_labels, mv_majority_voting, mv_over_clustering, mv_conf_score
                                                 majority-voting
"""

import os
import argparse

import scanpy as sc
import celltypist
from celltypist import models
from loguru import logger


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input-h5ad', required=True)
    p.add_argument('--output-h5ad', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--model', default='Immune_All_Low.pkl',
                   help='celltypist 模型名 (例: Immune_All_Low.pkl, Immune_All_High.pkl)')
    p.add_argument('--mode', default='both', choices=['bm', 'mv', 'both'],
                   help='bm=best match, mv=majority voting, both=两者都跑')
    return p.parse_args()


def annotate_majority_voting(adata, model):
    logger.info('majority voting 注释 ...')
    pred = celltypist.annotate(adata, model=model,
                               majority_voting=True, mode='prob match')
    cols = ['predicted_labels', 'over_clustering', 'majority_voting']
    adata.obs[['mv_predicted_labels', 'mv_over_clustering', 'mv_majority_voting']] = \
        pred.predicted_labels[cols].values
    # conf_score 在 obs 中, 单独取
    pred_adata = pred.to_adata()
    if 'conf_score' in pred_adata.obs.columns:
        adata.obs['mv_conf_score'] = pred_adata.obs['conf_score'].values


def annotate_best_match(adata, model):
    logger.info('best match 注释 ...')
    pred = celltypist.annotate(adata, model=model, mode='best match')
    pred_adata = pred.to_adata()
    adata.obs[['bm_predLables', 'bm_conf_score']] = \
        pred_adata.obs[['predicted_labels', 'conf_score']].values


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    logger.info(f'读取 {args.input_h5ad}')
    adata = sc.read_h5ad(args.input_h5ad)

    logger.info(f'加载 celltypist 模型: {args.model}')
    model = models.Model.load(model=args.model)

    if args.mode in ('mv', 'both'):
        annotate_majority_voting(adata, model)
    if args.mode in ('bm', 'both'):
        annotate_best_match(adata, model)

    adata.write_h5ad(args.output_h5ad)

    if 'bm_predLables' in adata.obs.columns:
        stat_file = os.path.join(args.output_dir, 'celltypist_bm_cellAnnotStat.xlsx')
        adata.obs.groupby('bm_predLables')['bm_predLables'].count() \
            .sort_values(ascending=False).to_excel(stat_file)
        logger.success(f'best-match 统计 → {stat_file}')
    if 'mv_majority_voting' in adata.obs.columns:
        stat_file = os.path.join(args.output_dir, 'celltypist_mv_cellAnnotStat.xlsx')
        adata.obs.groupby('mv_majority_voting')['mv_majority_voting'].count() \
            .sort_values(ascending=False).to_excel(stat_file)
        logger.success(f'majority-voting 统计 → {stat_file}')

    logger.success(f'注释完成 → {args.output_h5ad}')


if __name__ == '__main__':
    main()
