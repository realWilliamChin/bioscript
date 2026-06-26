#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/03/10 12:41
# Author        : William GoGo

import argparse
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from venn import venn

# 添加CommonTools到路径（根据脚本位置自动计算）
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../'))
from CommonTools.load_input import load_table, write_output_df


def read_single_file(file_path, remove_header=False, first_column_only=False):
    """读取单文件（TSV格式），每列为一个集合"""
    if remove_header:
        df = load_table(file_path, dtype=str, skiprows=1)
        # 如果没有列名，使用数字作为列名
        if df.columns.tolist() == list(range(len(df.columns))):
            df.columns = [f'Column_{i+1}' for i in range(len(df.columns))]
    else:
        df = load_table(file_path, dtype=str)
    
    # 如果设置了只读取第一列，则只保留第一列
    if first_column_only:
        first_col = df.columns[0]
        df = df[[first_col]]
    
    data = {col: set(df[col].dropna()) for col in df.columns}
    return data


def read_multiple_files(file_paths, remove_header=False, first_column_only=False):
    """读取多文件，每个文件内容为一个集合"""
    data = {}
    for path in file_paths:
        name = os.path.splitext(os.path.basename(path))[0]
        if first_column_only:
            # 如果是表格格式，只读取第一列
            df = load_table(path, dtype=str)
            if remove_header:
                df = df.iloc[1:]
            first_col = df.columns[0]
            elements = set(df[first_col].dropna())
        else:
            # 原始方式：按行读取
            with open(path, 'r') as f:
                lines = f.readlines()
                if remove_header and lines:
                    lines = lines[1:]  # 跳过第一行
                elements = {line.strip() for line in lines if line.strip()}
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
            'Fraction_ID': f'{label}_only',
            'Count': len(unique)
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
                'Fraction_ID': f'{"_".join(current_labels)}_common_only',
                'Count': len(result)
            })
    
    # 将结果转换为DataFrame并保存
    df = pd.DataFrame(results)
    write_output_df(df, output_name, index=False)


def analyze_sets(sets, labels, output_dir, definition_file=None):
    """分析集合并输出结果"""
    os.makedirs(output_dir, exist_ok=True)
    
    if definition_file:
        def_df = load_table(definition_file)
    
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
        label_only_file = os.path.join(output_dir, f'{label}_only.txt')
        with open(label_only_file, 'w') as f:
            f.write('\n'.join(sorted(unique)))
        if definition_file:
            label_only_df = load_table(label_only_file, header=None, names=[def_df.columns[0], ])
            label_only_df = pd.merge(label_only_df, def_df, on=def_df.columns[0], how='left')
            label_only_df.drop_duplicates(subset=[def_df.columns[0]], inplace=True)
            write_output_df(label_only_df, label_only_file, index=False)

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
            set_file = os.path.join(output_dir, '_'.join(current_labels) + '_common_only.txt')
            
            # 写入结果
            with open(set_file, 'w') as f:
                f.write('\n'.join(sorted(result)))
                
            if definition_file:
                set_df = load_table(set_file, header=None, names=[def_df.columns[0], ])
                set_df = pd.merge(set_df, def_df, on=def_df.columns[0], how='left')
                set_df.drop_duplicates(subset=[def_df.columns[0]], inplace=True)
                write_output_df(set_df, set_file, index=False)


def parse_input():
    parser = argparse.ArgumentParser(description="Venn 图分析工具")
    parser.add_argument('-i', '--input', nargs='+', required=True, help="输入文件路径")
    parser.add_argument('-o', '--output', help="common_only 文件目录")
    parser.add_argument('--pic-name', dest='pic_name', default='venn_diagram.jpeg', help='输出图片名称')
    parser.add_argument('--summary-name', dest='summary_name', default='intersection_summary.csv', help='输出汇总表名称')
    parser.add_argument('--remove-header', dest='remove_header', action='store_true', help='去掉读取表的第一行（header）')
    parser.add_argument('--first-column-only', dest='first_column_only', action='store_true', help='只读取表的第一列')
    parser.add_argument('-d', '--definition', help='输入定义文件，输出的每个文件添加定义')
    args = parser.parse_args()
    return args


def main():
    args = parse_input()

    if len(args.input) == 1:
        data = read_single_file(args.input[0], args.remove_header, args.first_column_only)
    else:
        data = read_multiple_files(args.input, args.remove_header, args.first_column_only)

    sets = list(data.values())
    labels = list(data.keys())

    draw_venn(data, args.pic_name)
    if args.output:
        analyze_sets(sets, labels, args.output, args.definition)
        generate_intersection_summary(sets, labels, args.summary_name)


if __name__ == "__main__":
    main()