#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/04/28 17:40
# Author        : William GoGo
import os, sys
import re
import argparse
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import write_output_df

def parse_args():
    parser = argparse.ArgumentParser(description='统计VCF文件中SUPP值分布并自动生成柱状图')
    parser.add_argument('-i', '--input', required=True, help='输入带SUPP字段的VCF文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出统计结果TSV文件路径')
    parser.add_argument('--img-output', help='输出图片路径，默认和TSV同目录同名.png，支持png/jpg/pdf/svg')
    parser.add_argument('--color', default='#1f77b4', help='柱子颜色，默认蓝色: #1f77b4')
    parser.add_argument('--figsize', nargs=2, type=int, default=(10, 6), help='图片尺寸，默认: 10 6')
    parser.add_argument('--dpi', type=int, default=300, help='输出图片DPI，默认: 300')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    logger.info(f'开始统计VCF文件: {args.input}')
    
    supp_counts = defaultdict(int)
    total_variants = 0
    
    with open(args.input, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 8:
                continue
            
            info = parts[7]
            # 提取SUPP值
            supp_match = re.search(r'SUPP=(\d+)', info)
            if not supp_match:
                logger.warning('该记录未找到SUPP字段，跳过')
                continue
            
            supp = int(supp_match.group(1))
            supp_counts[supp] += 1
            total_variants += 1
    
    if total_variants == 0:
        logger.error('未找到任何有效变异记录')
        sys.exit(1)
    
    logger.info(f'总共有 {total_variants} 条变异记录')
    
    # 按SUPP从小到大排序
    sorted_supp = sorted(supp_counts.keys())
    max_supp = max(sorted_supp) if sorted_supp else 0
    
    # 计算累计剩余：最少需要min_supp个样本支持的变异数量（SUPP≥min_supp）
    min_supp_list = list(range(1, max_supp + 1))
    remaining_counts = []
    for min_supp in min_supp_list:
        count = sum(supp_counts[s] for s in supp_counts if s >= min_supp)
        remaining_counts.append(count)
    
    # 输出两种统计结果
    logger.info('=== 1. 每个SUPP对应的变异数量 ===')
    cumulative = 0
    for supp in sorted_supp:
        count = supp_counts[supp]
        cumulative += count
        remaining = total_variants - cumulative
        logger.info(f'SUPP={supp}: {count} 条, 去掉SUPP≤{supp}后剩余 {remaining} 条')
    
    logger.info('\n=== 2. 不同最低支持数对应的剩余变异数量 ===')
    for min_supp, count in zip(min_supp_list, remaining_counts):
        logger.info(f'至少需要{min_supp}个样本支持: {count} 条变异')
    
    # 生成DataFrame并保存
    df_remaining = pd.DataFrame({
        '最少需要N个样本支持': min_supp_list,
        '满足条件的变异数': remaining_counts
    })
    
    write_output_df(df_remaining, args.output, index=False)
    logger.info(f'过滤后剩余变异数统计已保存到: {args.output}')
    
    # 自动生成图片
    if not args.img_output:
        args.img_output = os.path.splitext(args.output)[0] + '.png'
    
    logger.info(f'开始绘制柱状图，保存到: {args.img_output}')
    
    # 绘制柱状图（用剩余数量数据）
    fig, ax = plt.subplots(figsize=args.figsize)
    bars = ax.bar(df_remaining['最少需要N个样本支持'].astype(str), df_remaining['满足条件的变异数'], color=args.color, edgecolor='black')
    
    # 设置标题和标签
    ax.set_title('Distribution of SV', fontsize=14, pad=20)
    ax.set_xlabel('SV sample number', fontsize=12)
    ax.set_ylabel('Number', fontsize=12)
    
    # 在柱子上显示数值
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height + max(df_remaining['满足条件的变异数'])*0.01,
               f'{int(height)}',
               ha='center', fontsize=10)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig(args.img_output, dpi=args.dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f'柱状图绘制完成! 图片已保存到: {args.img_output}')

if __name__ == '__main__':
    main()