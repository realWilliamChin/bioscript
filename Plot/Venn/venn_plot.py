#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/03/10 12:41
# Author        : William GoGo

import argparse
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from venn import venn

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def read_single_file(file_path):
    """读取单文件（TSV格式），每列为一个集合"""
    df = load_table(file_path, dtype=str)
    data = {col: set(df[col].dropna()) for col in df.columns}
    return data


def read_multiple_files(file_paths):
    """读取多文件，每个文件内容为一个集合"""
    data = {}
    for path in file_paths:
        name = os.path.splitext(os.path.basename(path))[0]
        with open(path, 'r') as f:
            elements = {line.strip() for line in f if line.strip()}
        data[name] = elements
    return data


def draw_venn(data, pic_name):
    """绘制 Venn 图"""
    plt.figure(figsize=(8, 6), dpi=150)
    # ========== 全局样式设置 ==========
    plt.style.use('seaborn-v0_8-whitegrid')
    # rcParams['font.family'] = 'Arial'
    rcParams['axes.labelcolor'] = '#2d3436'
    rcParams['axes.edgecolor'] = '#dfe6e9'

    venn(
        data,
        cmap="plasma",
        alpha=0.7,
        fontsize=12,
        legend_loc="upper left"
    )

    plt.title("Venn Diagram", fontsize=16, pad=20)
    plt.gca().set_facecolor('#f8f9fa')
    plt.grid(False)
    plt.savefig(pic_name, bbox_inches='tight', dpi=300)
    plt.close()


def generate_intersection_summary(sets, labels, output_name):
    """生成所有交集的汇总表"""
    from itertools import combinations
    
    # 创建结果列表
    results = []
    
    # 处理每个集合的独有元素
    for i, (label, s) in enumerate(zip(labels, sets)):
        others = sets[:i] + sets[i+1:]
        unique = s - set.union(*others) if others else s
        results.append({
            '组合': label,
            '数量': len(unique)
        })
    
    # 处理所有可能的交集
    n = len(sets)
    for r in range(2, n+1):
        for indices in combinations(range(n), r):
            current_sets = [sets[i] for i in indices]
            current_labels = [labels[i] for i in indices]
            intersection = set.intersection(*current_sets)
            
            # 获取其他集合的并集
            other_indices = set(range(n)) - set(indices)
            other_sets = [sets[i] for i in other_indices]
            other_union = set.union(*other_sets) if other_sets else set()
            
            # 从交集中排除其他集合的元素
            result = intersection - other_union
            
            results.append({
                'Fraction_ID': '_'.join(current_labels),
                'Count': len(result)
            })
    
    # 将结果转换为DataFrame并保存
    df = pd.DataFrame(results)
    write_output_df(df, output_name, index=False)


def analyze_sets(sets, labels, output_dir):
    """分析集合并输出结果"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 计算所有集合的交集
    if len(sets) >= 2:
        common = set.intersection(*sets)
        with open(os.path.join(output_dir, '_'.join(labels)+'_common.txt'), 'w') as f:
            f.write('\n'.join(sorted(common)))
    else:
        print("Warning: 至少需要2个集合计算共同元素")

    # 计算每个集合的独有元素
    for i, (label, s) in enumerate(zip(labels, sets)):
        others = sets[:i] + sets[i+1:]
        unique = s - set.union(*others) if others else s
        with open(os.path.join(output_dir, f'{label}_only.txt'), 'w') as f:
            f.write('\n'.join(sorted(unique)))

    # 计算所有可能的复杂交集
    from itertools import combinations
    n = len(sets)
    
    # 生成所有可能的组合（从2个到n-1个集合的组合）
    for r in range(2, n):
        for indices in combinations(range(n), r):
            # 获取当前组合的集合和标签
            current_sets = [sets[i] for i in indices]
            current_labels = [labels[i] for i in indices]
            
            # 计算当前组合的交集
            intersection = set.intersection(*current_sets)
            
            # 获取其他集合的并集
            other_indices = set(range(n)) - set(indices)
            other_sets = [sets[i] for i in other_indices]
            other_union = set.union(*other_sets) if other_sets else set()
            
            # 从交集中排除其他集合的元素
            result = intersection - other_union
            
            # 生成文件名（使用标签的组合）
            filename = '_'.join(current_labels) + '_common_only.txt'
            
            # 写入结果
            with open(os.path.join(output_dir, filename), 'w') as f:
                f.write('\n'.join(sorted(result)))
    


def parse_input():
    parser = argparse.ArgumentParser(description="Venn 图分析工具")
    parser.add_argument('-i', '--input', nargs='+', required=True, help="输入文件路径")
    parser.add_argument('-o', '--output', help="common_only 文件目录")
    parser.add_argument('--pic-name', dest='pic_name', default='venn_diagram.jpeg', help='输出图片名称')
    parser.add_argument('--summary-name', dest='summary_name', default='intersection_summary.csv', help='输出汇总表名称')
    args = parser.parse_args()
    return args


def main():
    args = parse_input()

    if len(args.input) == 1:
        data = read_single_file(args.input[0])
    else:
        data = read_multiple_files(args.input)

    sets = list(data.values())
    labels = list(data.keys())

    draw_venn(data, args.pic_name)
    if args.output:
        analyze_sets(sets, labels, args.output)
        generate_intersection_summary(sets, labels, args.summary_name)


if __name__ == "__main__":
    main()