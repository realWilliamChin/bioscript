#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/04/15 11:03
# Author        : William GoGo
import os, sys
import argparse
import ast
import pandas as pd
from loguru import logger
sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df
from data_check import convert_numeric_columns

def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='输入 taxonomy.txt 文件')
    p.add_argument('-o', '--output', default='taxonomy_formatted.txt', help='输出 taxonomy_formatted.txt 文件')
    p.add_argument('--output-dir', help='输出各个分类的文件')
    p.add_argument('-r', '--raw-taxonomy', help='原始 taxonomy.txt 文件（未过滤，用于计算 Raw_read_count）')
    p.add_argument('-s', '--samples-described', help='样本描述文件路径（包含 group 和 sample 列）')
    p.add_argument('--summary-output', help='输出总结表文件路径')
    
    return p.parse_args()


def split_taxonomy(taxonomy_str):
    """将分类学字符串拆分成各个分类级别，识别 d__/k__, p__, c__, o__, f__, g__, s__ 等前缀"""
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    # 支持 d__ (Domain) 和 k__ (Kingdom) 两种格式
    prefixes = [['d__', 'k__'], 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    parts = taxonomy_str.split(';')
    result = {}
    
    # 初始化所有级别为默认值（使用 k__ 作为默认，microeco 兼容）
    result['Kingdom'] = 'k__'
    for level, prefix in zip(levels[1:], prefixes[1:]):
        result[level] = prefix
    
    # 处理每个部分
    for part in parts:
        part = part.strip()
        if not part:
            continue
        # 处理 Kingdom/Domain 级别（支持 d__ 和 k__）
        if part.startswith('d__') or part.startswith('k__'):
            result['Kingdom'] = part if part.startswith('k__') else part.replace('d__', 'k__', 1)
        # 处理其他级别
        else:
            for i, prefix in enumerate(prefixes[1:], 1):
                if part.startswith(prefix):
                    result[levels[i]] = part
                    break
    
    return result


def normalize_taxonomy(raw_taxonomy):
    """
    将 taxonomy 字段统一为 dict:
    - 若已是 dict 直接返回
    - 若是形如 "{'Kingdom': 'k__Bacteria', ...}" 的字符串，使用 literal_eval 解析
    - 其他字符串则按分号分隔格式交给 split_taxonomy
    """
    if isinstance(raw_taxonomy, dict):
        return raw_taxonomy
    if isinstance(raw_taxonomy, str):
        text = raw_taxonomy.strip()
        if text.startswith('{') and text.endswith('}'):
            try:
                return ast.literal_eval(text)
            except Exception:
                logger.warning(f"literal_eval 解析失败，按分号格式处理: {text}")
        return split_taxonomy(text)
    # 兜底返回空 dict
    return {}


def generate_summary_table(input_file, raw_taxonomy_file, output_dir, samples_described_file, output_summary_file):
    """
    生成总结表，包含各个样本和分组的统计信息
    
    参数:
        input_file: 输入的 ASV 表文件路径
        output_dir: 各个分类层级表的目录
        samples_described_file: 样本描述文件路径（包含 group 和 sample 列）
        output_summary_file: 输出总结表文件路径
    """
    logger.info("开始生成总结表...")
    
    # 读取输入 ASV 表（过滤后的，用于 ASV_count 和各分类层级统计）
    input_df = load_table(input_file)
    # 读取原始 ASV 表（未过滤，用于 Raw_read_count）
    raw_df = load_table(raw_taxonomy_file)
    
    # 读取样本描述文件
    samples_df = load_table(samples_described_file)
    if 'group' not in samples_df.columns or 'sample' not in samples_df.columns:
        raise ValueError("samples_described.txt 必须包含 'group' 和 'sample' 列")
    
    # 获取样本列（来自 samples_described.txt 的 sample 列）
    sample_cols = samples_df['sample'].tolist()
    
    # 读取各个分类层级的表
    taxonomy_levels = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    level_dfs = {}
    for level in taxonomy_levels:
        level_file = os.path.join(output_dir, f'{level}.txt')
        if os.path.exists(level_file):
            level_dfs[level] = load_table(level_file)
        else:
            logger.error(f"未找到 {level_file}，请检查输出目录是否正确，跳过该层级")
    
    # 初始化结果列表
    results = []
    
    # 1. 计算 Total（所有样本的总和）
    # Raw_read_count 使用原始表计算
    total_raw_read_count = raw_df[sample_cols].sum().sum()
    total_asv_count = (input_df[sample_cols].sum(axis=1) > 0).sum()
    total_stats = {
        'Type': 'Total',
        'Raw_read_count': total_raw_read_count,
        'ASV_count': total_asv_count
    }
    # 计算各个分类层级的数量
    for level in taxonomy_levels:
        if level in level_dfs:
            level_df = level_dfs[level]
            level_sample_cols = [col for col in level_df.columns if col in sample_cols]
            if level_sample_cols:
                # 统计非零的分类单元数
                total_stats[level] = (level_df[level_sample_cols].sum(axis=1) > 0).sum()
            else:
                total_stats[level] = 0
        else:
            total_stats[level] = 0
    results.append(total_stats)
    
    # 2. 计算各个样本的统计
    for sample in sample_cols:
        # Raw_read_count 使用原始表
        sample_read_count = raw_df[sample].sum()
        sample_asv_count = (input_df[sample] > 0).sum()
        sample_stats = {
            'Type': sample,
            'Raw_read_count': sample_read_count,
            'ASV_count': sample_asv_count
        }
        # 计算各个分类层级的数量
        for level in taxonomy_levels:
            if level in level_dfs:
                level_df = level_dfs[level]
                if sample in level_df.columns:
                    sample_stats[level] = (level_df[sample] > 0).sum()
                else:
                    sample_stats[level] = 0
            else:
                sample_stats[level] = 0
        results.append(sample_stats)
    
    # 3. 计算分组的统计
    groups = samples_df['group'].unique()
    group_stats_dict = {}
    # 控制输出顺序，分别保存 total / ave / coverage 到不同的列表
    group_total_rows = []
    group_ave_rows = []
    group_coverage_rows = []
    
    for group in groups:
        # 获取该组的样本列表
        group_samples = samples_df[samples_df['group'] == group]['sample'].tolist()
        # 过滤出实际存在的样本列
        group_sample_cols = [col for col in group_samples if col in sample_cols]
        
        if not group_sample_cols:
            logger.warning(f"分组 {group} 没有找到对应的样本列")
            continue
        
        # 计算分组总和
        # Raw_read_count 使用原始表
        group_read_count = raw_df[group_sample_cols].sum().sum()
        group_asv_count = (input_df[group_sample_cols].sum(axis=1) > 0).sum()
        
        group_total_stats = {
            'Type': f'{group}_total',
            'Raw_read_count': group_read_count,
            'ASV_count': group_asv_count
        }
        
        # 计算各个分类层级的数量
        for level in taxonomy_levels:
            if level in level_dfs:
                level_df = level_dfs[level]
                level_group_cols = [col for col in group_sample_cols if col in level_df.columns]
                if level_group_cols:
                    group_total_stats[level] = (level_df[level_group_cols].sum(axis=1) > 0).sum()
                else:
                    group_total_stats[level] = 0
            else:
                group_total_stats[level] = 0
        
        group_total_rows.append(group_total_stats)
        group_stats_dict[group] = {
            'total': group_total_stats,
            'samples': group_sample_cols
        }
        
        # 计算分组平均值
        num_samples = len(group_sample_cols)
        if num_samples > 0:
            group_ave_stats = {
                'Type': f'{group}_ave',
                'Raw_read_count': round(group_read_count / num_samples),
                'ASV_count': round(group_asv_count / num_samples)
            }
            
            # 计算各个分类层级的平均值（分组总数 / 样本数）
            for level in taxonomy_levels:
                if level in level_dfs:
                    level_df = level_dfs[level]
                    level_group_cols = [col for col in group_sample_cols if col in level_df.columns]
                    if level_group_cols:
                        # 统计分组内所有样本中至少出现一次的分类单元数
                        level_count = (level_df[level_group_cols].sum(axis=1) > 0).sum()
                        group_ave_stats[level] = round(level_count / num_samples)
                    else:
                        group_ave_stats[level] = 0
                else:
                    group_ave_stats[level] = 0
            
            group_ave_rows.append(group_ave_stats)
            group_stats_dict[group]['ave'] = group_ave_stats
            
            # 计算分组覆盖率（组总和 / Total总和）
            group_coverage_stats = {
                'Type': f'{group}_coverage',
                'Raw_read_count': round(group_total_stats['Raw_read_count'] / total_stats['Raw_read_count'], 3) if total_stats['Raw_read_count'] > 0 else 0,
                'ASV_count': round(group_total_stats['ASV_count'] / total_stats['ASV_count'], 3) if total_stats['ASV_count'] > 0 else 0
            }
            
            for level in taxonomy_levels:
                if total_stats[level] > 0:
                    group_coverage_stats[level] = round(group_total_stats[level] / total_stats[level], 3)
                else:
                    group_coverage_stats[level] = 0
            
            group_coverage_rows.append(group_coverage_stats)
    
    # 先追加所有 total，然后所有 ave，最后所有 coverage
    results.extend(group_total_rows)
    results.extend(group_ave_rows)
    results.extend(group_coverage_rows)
    
    # 转换为 DataFrame
    summary_df = pd.DataFrame(results)
    
    # 确保列顺序正确
    columns_order = ['Type', 'Raw_read_count', 'ASV_count'] + taxonomy_levels
    summary_df = summary_df[columns_order]
    
    # 输出结果
    write_output_df(summary_df, output_summary_file, index=False)
    logger.success(f"总结表已保存到: {output_summary_file}")
    
    return summary_df


def main():
    args = parse_input()
    input_file, output_file = args.input, args.output
    output_dir = args.output_dir
    
    # 读取输入文件
    df = load_table(input_file)
    df['taxonomy'] = df['taxonomy'].str.replace(' ', '')
    os.makedirs(output_dir, exist_ok=True)
    levels_order = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    for tax_level in levels_order:
        tax_level_file = os.path.join(output_dir, f'{tax_level}.txt')
        tax_level_df = df.drop(columns=['ASV_ID']).copy()
        tax_level_df['index'] = tax_level_df["taxonomy"].str.extract(f'(^.*?{str.lower(tax_level[0])}__[^;]+)')
        tax_level_df = tax_level_df.drop(columns=["taxonomy"])
        # 只对数值列做 sum，分类列保留 first
        num_cols = tax_level_df.select_dtypes(include="number").columns
        tax_level_df = (
            tax_level_df
            .groupby('index', as_index=False)
            .agg({**{col: "sum" for col in num_cols},
                **{"index": "first"}})
        )
        tax_level_df.set_index('index', inplace=True)
        tax_level_df = convert_numeric_columns(tax_level_df, exclude_columns=['index'], convert_dtype='int')
        write_output_df(tax_level_df, tax_level_file)
    taxonomy_df = df[['ASV_ID', 'taxonomy']].copy()
    
    # 创建分类表
    logger.info("正在拆分分类信息...")
    taxonomy_df['taxonomy'] = taxonomy_df['taxonomy'].apply(normalize_taxonomy)

    # 展开 taxonomy 字典为多列
    expanded = pd.json_normalize(taxonomy_df['taxonomy'])
    expected_cols = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    for col in expected_cols:
        if col not in expanded.columns:
            expanded[col] = ''
    expanded = expanded[expected_cols]
    taxonomy_df = pd.concat([taxonomy_df[['ASV_ID']], expanded], axis=1)
    
    # 输出分类表
    # taxonomy_df.reset_index(inplace=True)
    write_output_df(taxonomy_df, output_file, index=False)
    
    # 如果提供了原始 taxonomy、样本描述文件和总结表输出路径，生成总结表
    if args.raw_taxonomy and args.samples_described and args.summary_output and output_dir:
        generate_summary_table(
            input_file=input_file,
            raw_taxonomy_file=args.raw_taxonomy,
            output_dir=output_dir,
            samples_described_file=args.samples_described,
            output_summary_file=args.summary_output
        )


if __name__ == '__main__':
    main()