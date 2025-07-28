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


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', required=True,
                   help='输入table文件，包含 Func.refGene 和 Otherinfo3（Sequencing depth） 列')
    p.add_argument('-o', '--output', required=True, default='Dene_regions_average_depth_stat.jpg',
                   help='输出图片名称')
    p.add_argument('--width', type=float, default=12, help='图片宽度')
    p.add_argument('--height', type=float, default=8, help='图片高度')
    p.add_argument('--dpi', type=int, default=300, help='图片分辨率')
    
    args = p.parse_args()
    return args


def main():
    args = parse_input()
    data = load_table(args.input, header=0, usecols=['Func.refGene', 'Otherinfo3'])
    data.columns = ['Category', 'Otherinfo3']
    
    refGene_mapping = pd.DataFrame({
        'Category': ['exonic', 'exonic;splicing', 'splicing', 'intronic', 'intron', 'UTR5',
        'UTR3', 'UTR5;UTR3', 'downstream', 'upstream', 'upstream;downstream', 'intergenic',
        'ncRNA_exonic', 'ncRNA_exonic;splicing', 'ncRNA_splicing', 'ncRNA_intronic']
    })
    
    mean_values = data.groupby('Category')['Otherinfo3'].mean().reset_index()
    mean_values = pd.merge(refGene_mapping, mean_values, on='Category', how='outer')

    # 设置图形样式
    plt.figure(figsize=(args.width, args.height))
    sns.set(style="whitegrid")
    
    # 创建条形图
    ax = sns.barplot(
        y='Category', 
        x='Otherinfo3', 
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
        ax.text(
            row.Otherinfo3 + 0.01 * mean_values['Otherinfo3'].max(),
            i,  # Y位置
            f'{row.Otherinfo3:.2f}',  # 格式化数值
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