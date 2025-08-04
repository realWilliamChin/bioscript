#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/08/01 16:00
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='输入table，index=True')
    p.add_argument('-o', '--output-prefix', default='', help='输出前缀 + _')
    
    args = p.parse_args()
    
    if not args.output_prefix.endswith('_'):
        args.output_prefix += '_'
    
    return args


def plot_sample_level_percentages(stats_df, output_prefix):
    """
    绘制非零值百分比分布图
    
    Args:
        stats_df (pd.DataFrame): 包含非零值统计的DataFrame (行是分类级别，列是样本)
        output_prefix (str): 输出前缀
    """
    # 设置seaborn样式
    sns.set_theme(style="whitegrid")
    
    # 确保数据是数值类型
    stats_df = stats_df.set_index(stats_df.columns[0])
    stats_df = stats_df.apply(pd.to_numeric, errors='coerce')
    
    # 计算每个样本的总非零数
    sample_total_counts = stats_df.sum(axis=0)
    
    # 计算每个分类级别在每个样本中的百分比
    # 为了防止除以零，对于总数为零的样本，百分比设为0
    percentage_df = stats_df.div(sample_total_counts, axis=1) * 100
    percentage_df = percentage_df.fillna(0) # 填充NaN值为0
    
    # 转置 DataFrame 以便样本作为行索引，分类级别作为列，方便绘制堆叠柱状图
    percentage_df_transposed = percentage_df.T
    
    write_output_df(percentage_df_transposed, f'{output_prefix}level_percentages.txt', index=True)
    
    # 设置颜色调色板
    colors = sns.color_palette("viridis", n_colors=len(stats_df.index))
    
    # 绘制百分比图
    plt.figure(figsize=(14, 8))
    ax = percentage_df_transposed.plot(kind='bar', stacked=True, ax=plt.gca(), color=colors)
    plt.title('Percentages', fontsize=16, fontweight='bold')
    plt.xlabel('Sample', fontsize=14)
    plt.ylabel('Percentage (%)', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    
    # 添加网格线
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # 美化图例
    plt.legend(title='Levels', bbox_to_anchor=(1.05, 1), loc='upper left', 
               frameon=True, fancybox=True, shadow=True) # 添加图例并调整位置
    
    # 添加数值标签（如果样本数量不多）
    if len(percentage_df_transposed) < 10:
        for container in ax.containers:
            ax.bar_label(container, fmt='%.1f%%', label_type='center', fontsize=8)
    
    plt.tight_layout() # 调整布局以防止标签重叠
    plt.savefig(f'{output_prefix}level_percentages.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 绘制 Counts 图
    plt.figure(figsize=(14, 8))
    ax = stats_df.T.plot(kind='bar', stacked=True, ax=plt.gca(), color=colors)
    plt.title('Counts', fontsize=16, fontweight='bold')
    plt.xlabel('Sample', fontsize=14)
    plt.ylabel('Counts', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    
    # 添加网格线
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # 美化图例
    plt.legend(title='Levels', bbox_to_anchor=(1.05, 1), loc='upper left', 
               frameon=True, fancybox=True, shadow=True) # 添加图例并调整位置
    
    # 添加数值标签（如果样本数量不多且数值不会太拥挤）
    if len(stats_df.T) < 10 and stats_df.max().max() < 1000:
        for container in ax.containers:
            ax.bar_label(container, fmt='%d', label_type='center', fontsize=8)
    
    plt.tight_layout() # 调整布局以防止标签重叠
    plt.savefig(f'{output_prefix}level_counts.png', dpi=300, bbox_inches='tight')
    plt.close()


def main():
    args = parse_input()
    stats_df = load_table(args.input)
    plot_sample_level_percentages(stats_df, args.output_prefix)


if __name__ == '__main__':
    main()
