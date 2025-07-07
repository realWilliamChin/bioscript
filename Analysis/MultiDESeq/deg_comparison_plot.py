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
    parser.add_argument('-i', type=str, help='输入文件 DEG_summary.txt')
    parser.add_argument('-o', type=str, default='DEGs_comparison_count_barplot.jpeg',
                        help='输出文件')
    args = parser.parse_args()
    return args


def deg_summary_plot(deg_summary_df, output_file):
    deg_summary_df.columns = ['Comparisons', 'Total DEGs', 'Up regulated', 'Down regulated']
    numeric_cols = ['Total DEGs', 'Up regulated', 'Down regulated']
    for col in numeric_cols:
        deg_summary_df[col] = pd.to_numeric(deg_summary_df[col], errors='coerce').fillna(0).astype(int)
    
    comparisons = deg_summary_df['Comparisons']
    deg_summary_df = deg_summary_df.set_index('Comparisons')

    x = np.arange(len(comparisons))  # 生成组的位置
    width = 0.2  # 柱子的宽度

    # 创建图表
    fig, ax = plt.subplots(figsize=(10, 6))

    # 绘制柱状图
    rects1 = ax.bar(x - width/2, deg_summary_df["Up regulated"], width, 
                label='Up regulated genes', color='#1f77b4')  # 蓝色
    rects2 = ax.bar(x + width/2, deg_summary_df["Down regulated"], width, 
                label='Down regulated genes', color='#d62728')  # 红色

    # 添加标签和标题
    ax.set_ylabel('Number of Genes', fontsize=12)
    ax.set_title('Differentially Expressed Genes (DEGs)', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(comparisons, rotation=45, ha='right', fontsize=10)
    ax.legend()

    # 自动调整y轴范围
    max_value = max(max(deg_summary_df["Up regulated"]), max(deg_summary_df["Down regulated"]))
    ax.set_ylim(0, max_value * 1.15)

    # 添加数值标签
    for rects in [rects1, rects2]:
        for rect in rects:
            height = rect.get_height()
            ax.annotate(f'{height}',
                        xy=(rect.get_x() + rect.get_width()/2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    # 调整布局并保存
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')


def main():
    args = parse_input()
    deg_summary_df = load_table(args.i, comment='#', skipinitialspace=True)
    deg_summary_plot(deg_summary_df, args.o)


if __name__ == '__main__':
    main()
