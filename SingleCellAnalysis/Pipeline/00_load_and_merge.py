#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/06/29
# Author        : William GoGo
"""Pipeline 第 00 步: 读入多个 10x mtx 样本目录, 加样本元信息, 合并为单个 AnnData。

输入:
    --raw-dir       存放每个样本子目录的根路径, 每个子目录含 barcodes/features/matrix
    --sample-meta   样本元信息表(xlsx/tsv/csv), 必含 Sample_ID/Group_ID/Data_ID
                    Data_ID 与 --raw-dir 下的子目录名匹配
    --sample-sheet  仅当 sample-meta 是 xlsx 时使用, 指定 sheet 名
    --target-samples 可选,要保留的样本 ID 列表(逗号分隔, 对应 Data_ID)
    --output-h5ad   合并后的 h5ad 输出路径

步骤:
    1. 对每个目标样本读取 10x mtx -> AnnData
    2. 重建唯一 cellID (sample_cell_<i>), 基因名大写、去版本号、去重
    3. 把 sampleID/groupID/dataID 写入 obs
    4. ad.concat 合并所有样本, label='sample'
"""

import os
import re
import argparse

import pandas as pd
import scanpy as sc
import anndata as ad
from loguru import logger


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--raw-dir', required=True, help='样本子目录根路径')
    p.add_argument('--sample-meta', required=True, help='样本元信息表(xlsx/tsv/csv)')
    p.add_argument('--sample-sheet', default=None, help='xlsx 的 sheet 名')
    p.add_argument('--header-row', type=int, default=0, help='元信息表头行号 (默认 0)')
    p.add_argument('--target-samples', default=None,
                   help='可选, 要保留的 Data_ID 列表 (逗号分隔)')
    p.add_argument('--output-h5ad', required=True, help='合并后的 h5ad 输出路径')
    return p.parse_args()


def load_sample_meta(path, sheet=None, header=0):
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.xlsx', '.xls'):
        df = pd.read_excel(path, sheet_name=sheet, header=header)
    else:
        df = pd.read_csv(path, sep=None, engine='python', header=header)
    required = {'Sample_ID', 'Group_ID', 'Data_ID'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f'元信息表缺少列: {missing}')
    df = df[['Sample_ID', 'Group_ID', 'Data_ID']].dropna()
    df['Group_ID'] = df['Group_ID'].apply(lambda x: str(x).replace(' ', ''))
    return df.set_index('Data_ID')


def load_one_sample(sample_id, mtx_dir, sample_meta_row):
    adata = sc.read_10x_mtx(mtx_dir)
    # 1. 重建唯一 cellID
    new_cell_ids = [f'{sample_id}_cell_{i}' for i in range(adata.n_obs)]
    # 2. 修复基因名 (大写 + 去版本号 + 去重)
    adata.var_names = adata.var_names.str.upper()
    adata.var_names = adata.var_names.str.replace(r'\.\d+$', '', regex=True)
    adata = adata[:, ~adata.var_names.duplicated()].copy()
    # 3. 写 obs 元信息
    adata.obs['sample'] = sample_id
    adata.obs['sampleID'] = sample_meta_row['Sample_ID']
    adata.obs['group'] = sample_meta_row['Group_ID']
    adata.obs['cellID'] = new_cell_ids
    adata.obs_names = new_cell_ids
    assert len(adata.obs_names.unique()) == adata.n_obs, f'{sample_id} 细胞 ID 重复'
    assert len(adata.var_names.unique()) == adata.n_vars, f'{sample_id} 基因名重复'
    return adata


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.output_h5ad)) or '.', exist_ok=True)

    sample_meta = load_sample_meta(args.sample_meta,
                                   sheet=args.sample_sheet,
                                   header=args.header_row)
    target_samples = (set(s.strip() for s in args.target_samples.split(','))
                      if args.target_samples else None)

    adatas = {}
    for entry in sorted(os.listdir(args.raw_dir)):
        sample_id = entry.split('_')[0]
        if target_samples is not None and sample_id not in target_samples:
            continue
        if sample_id not in sample_meta.index:
            logger.warning(f'{sample_id} 不在元信息表中, 跳过')
            continue
        sample_dir = os.path.join(args.raw_dir, sample_id)
        if not os.path.isdir(sample_dir):
            logger.warning(f'{sample_dir} 不存在, 跳过')
            continue
        logger.info(f'读取 {sample_id} <- {sample_dir}')
        adatas[sample_id] = load_one_sample(sample_id, sample_dir,
                                            sample_meta.loc[sample_id])

    if not adatas:
        raise RuntimeError('未读取到任何样本, 检查 --raw-dir / --target-samples')

    logger.info(f'合并 {len(adatas)} 个样本 ...')
    adata_merged = ad.concat(list(adatas.values()),
                             label='sample', keys=list(adatas.keys()),
                             join='outer', fill_value=0)

    adata_merged.write_h5ad(args.output_h5ad)
    logger.success(f'合并完成 ({adata_merged.n_obs} cells × {adata_merged.n_vars} genes) → {args.output_h5ad}')


if __name__ == '__main__':
    main()
