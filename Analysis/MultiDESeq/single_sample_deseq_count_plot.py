#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2025/04/15 14:46
# Author        : WilliamGoGo

import os, sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from loguru import logger

# 忽略除零警告
warnings.filterwarnings('ignore', category=RuntimeWarning)

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--fpkm', help='输入 FPKM 文件')
    parser.add_argument('-s', '--samples', help='输入样本信息文件')
    parser.add_argument('-c', '--compare', help='输入比较信息文件')
    parser.add_argument('-o', '--output-prefix', dest='output_prefix', default='comparisons_all_deg_counts',
                        help='输出文件, 默认是 comparisons_all_deg_counts.jpg 和 comparisons_all_deg_counts.txt')
    args = parser.parse_args()
    return args


def plot_fold_change_distribution(log2_matrix, output_jpg, output_table):
    """
    对每个比对组的差异倍数分布进行绘图
    
    Args:
        log2_matrix: 包含比对结果的DataFrame
    """
    fold_changes = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4]
    log2_fold_changes = np.log2(fold_changes)
    
    sns.set_style("whitegrid")
    plt.figure(figsize=(12, 8))
    
    # 初始化空的DataFrame
    all_deg_counts = pd.DataFrame(columns=['Comparison'] + [f'deg{fc}' for fc in fold_changes])
    all_counts = {}  # 存储所有比较组的counts

    for col in log2_matrix.columns:
        if col.startswith('log2_'):
            # 统计每个差异倍数对应的基因数量
            counts = []
            for i, fc in enumerate(log2_fold_changes):
                up_count = sum((log2_matrix[col] >= fc) & (log2_matrix[col] != np.inf))
                down_count = sum((log2_matrix[col] <= -fc) & (log2_matrix[col] != -np.inf))
                total_count = up_count + down_count
                counts.append(total_count)
            
            # 创建新的行并添加到DataFrame
            new_row = {'Comparison': col}
            for i, fc in enumerate(fold_changes):
                new_row[f'deg{fc}'] = counts[i]
            all_deg_counts = pd.concat([all_deg_counts, pd.DataFrame([new_row])], ignore_index=True)
            
            all_counts[col] = counts  # 保存当前比较组的counts
            
            # 折线图
            plt.plot(log2_fold_changes, counts, marker='o', label=col)
    
    write_output_df(all_deg_counts, output_table, index=False)
    
    # 设置图表属性
    plt.xlabel('Fold Change')
    plt.ylabel('Number of Differentially Expressed Genes')
    plt.title('Distribution of Differentially Expressed Genes by Fold Change')
    plt.legend(title='Comparison Groups', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    
    # 设置 x 轴刻度
    plt.xticks(log2_fold_changes, [f"{fc:.1f}" for fc in fold_changes])
    
    # 调整布局并显示
    plt.tight_layout()
    plt.savefig(output_jpg, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    args = parse_input()
    compare_df = load_table(args.compare)
    fpkm_df = load_table(args.fpkm, dtype={'GeneID': str})
    numeric_cols = fpkm_df.columns.difference(['GeneID'])
    fpkm_df[numeric_cols] = fpkm_df[numeric_cols].replace(0, 0.001)

    log2_matrix = fpkm_df[['GeneID']].copy()
    for i, compare in compare_df.iterrows():
        sample_treat = compare['Treat']
        sample_control = compare['Control']
        compare_name = f'{sample_treat}-vs-{sample_control}'
        log2_matrix[compare_name] = np.log2(fpkm_df[sample_treat] / fpkm_df[sample_control])
    
    write_output_df(log2_matrix, 'log2_fpkm_matrix.txt', index=False)
 
    plot_fold_change_distribution(
        log2_matrix,
        args.output_prefix+ '.jpg',
        args.output_prefix + '.txt'
    )
    
    logger.success(f'{__file__} Done!')


if __name__ == '__main__':
    main()
