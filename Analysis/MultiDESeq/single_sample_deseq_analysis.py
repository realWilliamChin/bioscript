#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2025/04/15 14:46
# Author        : WilliamGoGo

import os, sys
import argparse
import pandas as pd
import numpy as np
import math
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='输入 log2_fpkm 文件')
    parser.add_argument('-o', type=str, help='输出目录')
    parser.add_argument('--fpkm-reads-merged', type=str, help='基因所有表达量和注释的文件')
    parser.add_argument('--deg-value', type=float, help='输入差异倍数阈值')
    args = parser.parse_args()
    return args

        
def main():
    args = parse_input()
    log2_matrix = load_table(args.i)
    fpkm_reads_merged_df = load_table(args.fpkm_reads_merged)

    log2_matrix['GeneID'] = log2_matrix['GeneID'].astype(str)

    numeric_cols = log2_matrix.columns.difference(['GeneID'])
    for col in numeric_cols:
        log2_matrix[col] = pd.to_numeric(log2_matrix[col], errors='coerce')
    
    bs_pos = math.log2(args.deg_value)

    DEG_analysis_results = os.path.join(args.o, 'DEG_analysis_results')
    Expression_data = os.path.join(DEG_analysis_results, 'Expression_data')
    for dir in [args.o, DEG_analysis_results, Expression_data]:
        os.makedirs(dir, exist_ok=True)
    
    deg_summary_df = pd.DataFrame(columns=['Comparisons', 'Total DEGs', 'Up regulated', 'Down regulated'])
    for col in numeric_cols:
        comparsion_name = col.replace('log2_', '')
        treat = comparsion_name.split('-vs-')[0]
        control = comparsion_name.split('-vs-')[1]
        first_kns_cols = [col for col in fpkm_reads_merged_df.columns if col == 'NR_ID']
        if first_kns_cols:
            first_kns_col_idx = fpkm_reads_merged_df.columns.get_loc(first_kns_cols[0])
            kns_cols = fpkm_reads_merged_df.columns[first_kns_col_idx:].tolist()
        else:
            kns_cols = []
        treat_control_kns_cols = ['GeneID', f'{treat}_fpkm', f'{control}_fpkm', f'{treat}_reads', f'{control}_reads'] + kns_cols
        up_genes = log2_matrix[(log2_matrix[col] >= bs_pos) | np.isinf(log2_matrix[col])]
        up_genes = up_genes.rename(columns={col: f'{col}_FC'})
        up_genes = up_genes[['GeneID', f'{col}_FC']]
        up_genes = pd.merge(up_genes, fpkm_reads_merged_df[treat_control_kns_cols], on='GeneID', how='left')
        down_genes = log2_matrix[(log2_matrix[col] <= -bs_pos) | np.isneginf(log2_matrix[col])]
        down_genes = down_genes.rename(columns={col: f'{col}_FC'})
        down_genes = down_genes[['GeneID', f'{col}_FC']]
        down_genes = pd.merge(down_genes, fpkm_reads_merged_df[treat_control_kns_cols], on='GeneID', how='left')
        
        deg_summary_df = pd.concat([
            deg_summary_df,
            pd.DataFrame({
                'Comparisons': [comparsion_name],
                'Total DEGs': [len(up_genes) + len(down_genes)],
                'Up regulated': [len(up_genes)],
                'Down regulated': [len(down_genes)]
            })
        ], ignore_index=True)
        
        write_output_df(up_genes['GeneID'], os.path.join(DEG_analysis_results, f'{comparsion_name}_Up_ID.txt'), index=False)
        write_output_df(down_genes['GeneID'], os.path.join(DEG_analysis_results, f'{comparsion_name}_Down_ID.txt'), index=False)
        write_output_df(up_genes, os.path.join(Expression_data, f'{comparsion_name}_Up_DEG_data.txt'), index=False)
        write_output_df(down_genes, os.path.join(Expression_data, f'{comparsion_name}_Down_DEG_data.txt'), index=False)

    with open(os.path.join(DEG_analysis_results, 'DEG_summary.txt'), 'w') as f:
        f.write(f'# 筛选条件：FoldChange > {args.deg_value}\n')
        deg_summary_df.to_csv(f, sep='\t', index=False)

    logger.success(f'{__file__} Done!')


if __name__ == '__main__':
    main()