#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2025/04/16 16:37
# Author        : WilliamGoGo
# Description   : 绘制差异表达基因(DEGs)比较统计图

import os, sys
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, default='DEG_summary.txt', help='输入文件 DEG_summary.txt')
    parser.add_argument('-o', type=str, default='DEGs_comparison_count_barplot.jpeg',
                        help='输出文件')
    args = parser.parse_args()
    return args


def deg_summary_plot(deg_summary_df, output_file):
    deg_summary_df.columns = ['Comparisons', 'Total DEGs', 'Up regulated', 'Down regulated']
    numeric_cols = ['Total DEGs', 'Up regulated', 'Down regulated']
    for col in numeric_cols:
        deg_summary_df[col] = pd.to_numeric(deg_summary_df[col], errors='coerce').fillna(0).astype(int)

    # 转为长表
    plot_df = pd.melt(
        deg_summary_df,
        id_vars=['Comparisons'],
        value_vars=['Up regulated', 'Down regulated'],
        var_name='Regulation',
        value_name='Count'
    )

    n_groups = deg_summary_df.shape[0]
    fig_width = max(6, n_groups * 1.5)
    fig_height = 6

    plt.figure(figsize=(fig_width, fig_height))
    ax = sns.barplot(
        data=plot_df,
        x='Comparisons',
        y='Count',
        hue='Regulation',
        palette={'Up regulated': '#1f77b4', 'Down regulated': '#d62728'},
    )

    # 添加数值标签
    for p in ax.patches:
        height = int(p.get_height())
        if height > 0:
            ax.annotate(f'{height}',
                        (p.get_x() + p.get_width() / 2, height),
                        ha='center', va='bottom', fontsize=10, xytext=(0, 3), textcoords='offset points')

    ax.set_ylabel('Number of Genes', fontsize=12)
    ax.set_title('Differentially Expressed Genes (DEGs)', fontsize=14)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    ax.legend(title=None)
    ax.set_ylim(0, plot_df['Count'].max() * 1.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')


def main():
    args = parse_input()
    deg_summary_df = load_table(args.i, comment='#', skipinitialspace=True)
    deg_summary_plot(deg_summary_df, args.o)


if __name__ == '__main__':
    main()
