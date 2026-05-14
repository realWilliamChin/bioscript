#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/25 12:41
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df
from data_check import convert_numeric_columns


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', required=True,
                   help='输入table文件，包含 Func.refGene 和深度列')
    p.add_argument('-o', '--output', default='Gene_regions_average_depth_stat.jpg', help='输出图片名称')
    p.add_argument('--width', type=float, default=12, help='图片宽度')
    p.add_argument('--height', type=float, default=8, help='图片高度')
    p.add_argument('--dpi', type=int, default=300, help='图片分辨率')
    p.add_argument('--depth-col', default='Otherinfo3', help='深度列的列名，默认是Otherinfo3，恢复列名后可以指定为INFO、样本名等实际列名')
    
    args = p.parse_args()
    return args


def main():
    args = parse_input()
    data = load_table(args.input, header=0, usecols=['Func.refGene', args.depth_col])
    data.columns = ['Category', 'depth']
    data = convert_numeric_columns(data, exclude_columns=['Category'])
    
    # 强制转换深度列为数值类型，处理'.'等无效值
    data['depth'] = pd.to_numeric(data['depth'], errors='coerce')
    # 过滤掉深度列为空的行
    data = data.dropna(subset=['depth'])
    
    refGene_mapping = pd.DataFrame({
        'Category': ['exonic', 'exonic;splicing', 'splicing', 'intronic', 'intron', 'UTR5',
        'UTR3', 'UTR5;UTR3', 'downstream', 'upstream', 'upstream;downstream', 'intergenic',
        'ncRNA_exonic', 'ncRNA_exonic;splicing', 'ncRNA_splicing', 'ncRNA_intronic']
    })
    
    mean_values = data.groupby('Category')['depth'].mean().reset_index()
    mean_values = pd.merge(refGene_mapping, mean_values, on='Category', how='outer')

    # 设置图形样式
    plt.figure(figsize=(args.width, args.height))
    sns.set(style="whitegrid")
    
    # 创建条形图
    ax = sns.barplot(
        y='Category', 
        x='depth', 
        data=mean_values,
        hue='Category',
        palette="viridis",
        width=0.7,
        legend=False
    )
    
    current_xmax = ax.get_xlim()[1]
    new_xmax = current_xmax * 1.2
    ax.set_xlim(0, new_xmax)
    
    # 在每个条形上添加数值标签
    for i, (_, row) in enumerate(mean_values.iterrows()):
        if pd.notna(row.depth):
            ax.text(
                row.depth + 0.01 * mean_values['depth'].max(),
                i,
                '{:.2f}'.format(row.depth),  # 格式化数值
                ha='left',  # 水平对齐
                va='center',  # 垂直对齐
                fontsize=9
            )
    
    plt.title('Average sequencing depth statistics for gene regions', fontsize=14)
    plt.xlabel('Sequencing Depth', fontsize=12)
    plt.ylabel('Category')
    
    plt.tight_layout()
    plt.rc('font',family='SimHei')
    
    # 保存图片
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
    logger.info(f"图表已保存至: {args.output}")


if __name__ == '__main__':
    main()