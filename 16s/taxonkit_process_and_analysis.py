#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/05/14 14:00
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from distribution_table_plot import plot_sample_level_percentages

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser(description='Convert taxonomy table to 16S analysis format')
    p.add_argument('-e', '--expression', dest='expression', default='table_ASV.tsv', help='表达量文件路径')
    p.add_argument('-b', '--blast', dest='blast', default='nr2.blast', help='blast文件路径（去重的文件）')
    p.add_argument('-o', '--output-dir', dest='output_dir', default='.', help='输出目录路径')
    p.add_argument('--id-exclude', dest='id_exclude', default='', help='排除关键词(BLAST_ID)，多个关键词用逗号分隔')
    p.add_argument('--id-include', dest='id_include', default='', help='仅保留关键词(BLAST_ID)，多个关键词用逗号分隔')
    p.add_argument('--tax-exclude', dest='tax_exclude', default='', help='排除分类路径，格式如 d__Eukaryota;p__Basidiomycota，多条用逗号分隔')
    p.add_argument('--tax-include', dest='tax_include', default='', help='仅保留分类路径，格式如 d__Eukaryota;p__Basidiomycota，多条用逗号分隔')
    
    args = p.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, 'Prep_files'), exist_ok=True)
    
    return args


def taxonkit_output_to_16s_analysis_input(input_file):
    # 定义分类层级映射关系（前缀 + 列名）
    taxonomy_levels = [
        ('d__', 'domain'),
        ('p__', 'phylum'),
        ('c__', 'class'),
        ('o__', 'order'),
        ('f__', 'family'),
        ('g__', 'genus'),
        ('s__', 'species')
    ]
    
    try:
        df = load_table(input_file)

        # 处理每个分类层级
        formatted_lines = []
        for _, row in df.iterrows():
            parts = []
            # 添加分类信息
            for prefix, col in taxonomy_levels:
                value = row.get(col, '')
                if pd.isna(value) or str(value).strip() == '':
                    formatted = f"{prefix}__"
                else:
                    # 替换空格为下划线并添加前缀
                    formatted = f"{prefix}{str(value).replace(' ', '_')}"
                parts.append(formatted)
            
            formatted_lines.append(f"{row[df.columns[0]]}\t" + ";".join(parts))

        # 将字符串列表转换为两列数据
        data = [line.split('\t') for line in formatted_lines]
        result_df = pd.DataFrame(data, columns=['NCBI_taxon_ID', 'Taxon_ID'])
        
        return result_df
            
    except Exception as e:
        print(f"处理文件时发生错误: {str(e)}")
        raise


def run_taxonkit(taxid_file, output_file, data_dir='/home/train/.taxonkit'):
    """
    使用 taxonkit 处理分类信息
    
    Args:
        taxid_file (str): 包含 taxid 的输入文件路径
        output_file (str): 输出文件路径
        data_dir (str): taxonkit 数据目录路径
    """
    try:
        # 构建 taxonkit lineage 命令
        lineage_cmd = ['taxonkit', 'lineage', taxid_file, '--data-dir', data_dir]
        
        # 构建 taxonkit reformat 命令
        reformat_cmd = ['taxonkit', 'reformat', '-F', '-f', '{p}\t{c}\t{o}\t{f}\t{g}\t{s}']
        
        lineage_process = subprocess.Popen(lineage_cmd, stdout=subprocess.PIPE)
        reformat_process = subprocess.Popen(reformat_cmd, stdin=lineage_process.stdout, stdout=open(output_file, 'w'))
        
        # 等待命令执行完成
        lineage_process.stdout.close()
        reformat_process.communicate()
        
        if reformat_process.returncode != 0:
            raise Exception(f"taxonkit 命令执行失败，返回码: {reformat_process.returncode}")
            
    except Exception as e:
        print(f"执行 taxonkit 命令时发生错误: {str(e)}")
        raise


def filter_dataframe_by_blast_id(df, include_terms=None, exclude_terms=None, blast_id_column='Blast_ID'):
    """
    根据Blast_ID列的内容过滤DataFrame
    
    Args:
        df (pd.DataFrame): 要过滤的DataFrame
        include_terms (str, optional): 仅保留包含这些关键词的行，多个关键词用逗号分隔
        exclude_terms (str, optional): 排除包含这些关键词的行，多个关键词用逗号分隔
        blast_id_column (str): Blast_ID列的名称，默认为'Blast_ID'
    
    Returns:
        pd.DataFrame: 过滤后的DataFrame副本
    """
    if not include_terms and not exclude_terms:
        return df.copy()
    
    # 初始化过滤掩码
    if include_terms:
        # 仅保留模式：初始化为False，匹配到的设为True
        filter_mask = pd.Series(False, index=df.index)
        include_term_list = [term.strip() for term in include_terms.split(',') if term.strip()]
        for term in include_term_list:
            filter_mask |= df[blast_id_column].str.contains(term, case=False, na=False)
    if exclude_terms:
        # 排除模式：初始化为True，匹配到的设为False
        filter_mask = pd.Series(True, index=df.index)
        exclude_term_list = [term.strip() for term in exclude_terms.split(',') if term.strip()]
        for term in exclude_term_list:
            filter_mask &= ~df[blast_id_column].str.contains(term, case=False, na=False)
    
    return df[filter_mask].copy()


def _parse_taxon_path(path_str):
    """
    将类似于 d__Eukaryota;p__Basidiomycota 的路径解析为字典
    返回: {'domain': 'Eukaryota', 'phylum': 'Basidiomycota'}
    """
    if not isinstance(path_str, str) or path_str.strip() == '':
        return {}
    prefix_to_col = {
        'd__': 'domain',
        'p__': 'phylum',
        'c__': 'class',
        'o__': 'order',
        'f__': 'family',
        'g__': 'genus',
        's__': 'species'
    }
    result = {}
    for part in path_str.split(';'):
        part = part.strip()
        if not part:
            continue
        for prefix, col in prefix_to_col.items():
            if part.startswith(prefix):
                value = part[len(prefix):].strip('_')
                if value:
                    result[col] = value
                break
    return result


def filter_taxonomy_by_path(df, include_paths=None, exclude_paths=None):
    """
    基于分类路径字符串过滤 taxonomy DataFrame（包含列 domain/phylum/class/order/family/genus/species）。

    include_paths/exclude_paths: 逗号分隔的路径列表，如 'd__Eukaryota;p__Basidiomycota,g__Escherichia'
    规则:
      - include: 行需匹配任意一条包含路径（路径中每个层级都需相等，大小写不敏感，空值视为不匹配）
      - exclude: 行匹配任意一条排除路径则剔除
    两者同时给定时，先执行 include 再执行 exclude。
    """
    if not include_paths and not exclude_paths:
        return df.copy()

    work_df = df.copy()

    def row_matches_path(row, path_dict):
        for col, expected in path_dict.items():
            val = row.get(col, '')
            if pd.isna(val) or str(val).strip() == '':
                return False
            if str(val).strip().lower() != str(expected).strip().lower():
                return False
        return True

    mask = pd.Series(True, index=work_df.index)

    if include_paths:
        include_list = [s.strip() for s in include_paths.split(',') if s.strip()]
        include_dicts = [_parse_taxon_path(s) for s in include_list]
        inc_mask = pd.Series(False, index=work_df.index)
        for path_dict in include_dicts:
            if not path_dict:
                continue
            inc_mask |= work_df.apply(lambda r: row_matches_path(r, path_dict), axis=1)
        mask &= inc_mask

    if exclude_paths:
        exclude_list = [s.strip() for s in exclude_paths.split(',') if s.strip()]
        exclude_dicts = [_parse_taxon_path(s) for s in exclude_list]
        exc_mask = pd.Series(False, index=work_df.index)
        for path_dict in exclude_dicts:
            if not path_dict:
                continue
            exc_mask |= work_df.apply(lambda r: row_matches_path(r, path_dict), axis=1)
        mask &= ~exc_mask

    return work_df[mask].copy()


def generate_taxonomy_level_files(species_df, sample_col, output_dir):
    """
    生成各个分类级别的文件并统计非零值
    
    Args:
        species_df (pd.DataFrame): 包含物种级别数据的DataFrame
        sample_col (list): 样本列名列表
        output_dir (str): 输出目录路径
        
    Returns:
        pd.DataFrame: 包含各分类级别非零值统计的DataFrame
    """
    # 创建一个字典来存储每个分类级别的非零值统计
    non_zero_stats = {}
    
    # 生成各个分类级别的文件
    # 创建taxonomy_levels文件夹
    taxonomy_levels_dir = os.path.join(output_dir, 'taxonomy_levels')
    os.makedirs(taxonomy_levels_dir, exist_ok=True)
    
    # 创建species_df的副本以避免修改原始数据
    working_species_df = species_df.copy()
    working_species_df.drop(columns=['NCBI_taxon_ID', 'ASV_ID', 'Blast_ID'], inplace=True)
    write_output_df(working_species_df, os.path.join(taxonomy_levels_dir, 'Species.txt'), index=False)
    non_zero_stats['Species'] = (working_species_df[sample_col] > 0).sum()
    
    # Genus级别
    working_species_df = working_species_df.reset_index(drop=True)
    genus_df = working_species_df.copy()
    genus_df['Taxon_ID'] = genus_df['Taxon_ID'].str.replace(';s__.*', '', regex=True)
    genus_df = genus_df.groupby('Taxon_ID').agg({
        col: 'sum' if col in sample_col else 'first'
        for col in genus_df.columns
    })
    write_output_df(genus_df, os.path.join(taxonomy_levels_dir, 'Genus.txt'), index=False)
    non_zero_stats['Genus'] = (genus_df[sample_col] > 0).sum()
    
    # Family级别
    family_df = working_species_df.copy()
    family_df['Taxon_ID'] = family_df['Taxon_ID'].str.replace(';g__.*', '', regex=True)
    family_df = family_df.groupby('Taxon_ID').agg({
        col: 'sum' if col in sample_col else 'first'
        for col in family_df.columns
    })
    write_output_df(family_df, os.path.join(taxonomy_levels_dir, 'Family.txt'), index=False)
    non_zero_stats['Family'] = (family_df[sample_col] > 0).sum()
    
    # Order级别
    order_df = working_species_df.copy()
    order_df['Taxon_ID'] = order_df['Taxon_ID'].str.replace(';f__.*', '', regex=True)
    order_df = order_df.groupby('Taxon_ID').agg({
        col: 'sum' if col in sample_col else 'first'
        for col in order_df.columns
    })
    write_output_df(order_df, os.path.join(taxonomy_levels_dir, 'Order.txt'), index=False)
    non_zero_stats['Order'] = (order_df[sample_col] > 0).sum()
    
    # Class级别
    class_df = working_species_df.copy()
    class_df['Taxon_ID'] = class_df['Taxon_ID'].str.replace(';o__.*', '', regex=True)
    class_df = class_df.groupby('Taxon_ID').agg({
        col: 'sum' if col in sample_col else 'first'
        for col in class_df.columns
    })
    write_output_df(class_df, os.path.join(taxonomy_levels_dir, 'Class.txt'), index=False)
    non_zero_stats['Class'] = (class_df[sample_col] > 0).sum()
    
    # Phylum级别
    phylum_df = working_species_df.copy()
    phylum_df['Taxon_ID'] = phylum_df['Taxon_ID'].str.replace(';c__.*', '', regex=True)
    phylum_df = phylum_df.groupby('Taxon_ID').agg({
        col: 'sum' if col in sample_col else 'first'
        for col in phylum_df.columns
    })
    write_output_df(phylum_df, os.path.join(taxonomy_levels_dir, 'Phylum.txt'), index=False)
    non_zero_stats['Phylum'] = (phylum_df[sample_col] > 0).sum()
    
    # 创建统计结果DataFrame并绘图
    stats_df = pd.DataFrame(non_zero_stats).T
    plot_sample_level_percentages(stats_df, output_dir)
    
    return stats_df


def run_flattening_r(asv_file, output_file):
    """
    运行 R 脚本进行抽平分析
    
    Args:
        asv_file (str): ASV 文件路径
        output_file (str): 输出文件路径
    """
    try:
        # 构建 R 脚本命令
        r_script = '/home/colddata/qinqiang/script/16s/flattening.r'
        cmd = ['Rscript', r_script, '-i', os.path.abspath(asv_file), '-o', os.path.abspath(output_file)]
        
        # 执行 R 脚本
        process = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        if process.returncode != 0:
            print(f"R脚本错误输出：{process.stderr}")
            raise Exception(f"R 脚本执行失败，错误信息：{process.stderr}")
            
        print("抽平分析完成")
        print(process.stdout)
            
    except Exception as e:
        print(f"执行 R 脚本时发生错误: {str(e)}")
        raise


def main():
    args = parse_input()
    # ====== Load data ======
    expression_df = load_table(args.expression, dtype={"ASV_ID": str})
    # nr column qseqid -> ASV_ID; stitle -> Blast_ID; staxids -> NCBI_taxon_ID
    blast_df = load_table(args.blast, usecols=[0, 12, 13], header=None, names=['ASV_ID', 'Blast_ID', 'NCBI_taxon_ID'])
    
    # ====== 数据表预处理 ======
    expression_def_df = pd.merge(expression_df, blast_df, on='ASV_ID', how='inner')
    write_output_df(expression_def_df, os.path.join(args.output_dir, 'expression_def_df.txt'), index=False)
    
    # 获取所有列名，排除 'ASV_ID' 和 'Blast_ID'
    sample_col = [col for col in expression_df.columns if col not in ['ASV_ID', 'Blast_ID']]
    blast_df['NCBI_taxon_ID'] = blast_df['NCBI_taxon_ID'].str.split(';').str[0]

    # ====== 应用 BLAST_ID 过滤条件 ======
    expression_def_df = filter_dataframe_by_blast_id(
        expression_def_df, 
        include_terms=args.id_include if args.id_include else None,
        exclude_terms=args.id_exclude if args.id_exclude else None
    )
    
    # 按照 NCBI_taxon_ID 分组，样本列求和，Gene ID 保留每行 sum 中最高的 GeneID
    expression_def_df['row_sum'] = expression_def_df[sample_col].sum(axis=1)
    staxid_grouped_df = expression_def_df.groupby('NCBI_taxon_ID').agg({
        col: 'sum' if col in sample_col else lambda x: x.iloc[expression_def_df.loc[x.index, 'row_sum'].argmax()]
        for col in expression_def_df.columns if col not in ['row_sum', 'NCBI_taxon_ID']
    })
    # 重置索引，将 NCBI_taxon_ID 从索引转换为列
    staxid_grouped_df = staxid_grouped_df.reset_index()
    
    
    # ====== 运行 taxonkit ======
    taxid_file = os.path.join(args.output_dir, 'Prep_files', 'taxid_list.txt')
    taxonkit_raw_result_file = os.path.join(args.output_dir, 'Prep_files', 'taxonomy_raw.txt')
    taxonkit_formatted_file = os.path.join(args.output_dir, 'taxonomy_result.txt')
    taxonkit_filtered_file = os.path.join(args.output_dir, 'taxonomy_filtered.txt')
    staxid_grouped_df[['NCBI_taxon_ID']].to_csv(taxid_file, index=False, sep='\t', header=False)
    run_taxonkit(taxid_file, taxonkit_raw_result_file)
    taxonomy_df = load_table(
        taxonkit_raw_result_file,
        header=None, 
        names=['NCBI_taxon_ID', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    )
    taxonomy_df['domain'] = taxonomy_df['domain'].str.split(';').str[1]
    write_output_df(taxonomy_df, taxonkit_formatted_file, index=False)
    
    # ====== 基于分类路径进行过滤 ======
    if args.tax_include or args.tax_exclude:
        taxonomy_df_filtered = filter_taxonomy_by_path(
            taxonomy_df,
            include_paths=args.tax_include if args.tax_include else None,
            exclude_paths=args.tax_exclude if args.tax_exclude else None
        )
        allowed_taxids = set(taxonomy_df_filtered['NCBI_taxon_ID'].astype(str))
        staxid_grouped_df = staxid_grouped_df[staxid_grouped_df['NCBI_taxon_ID'].astype(str).isin(allowed_taxids)].copy()
        taxonomy_df = taxonomy_df[taxonomy_df['NCBI_taxon_ID'].astype(str).isin(allowed_taxids)].copy()
        write_output_df(taxonomy_df, taxonkit_filtered_file, index=False)
        
    # ====== 格式化分类信息
    taxonkit_formatted_df = taxonkit_output_to_16s_analysis_input(taxonkit_formatted_file)
    taxonkit_formatted_df['NCBI_taxon_ID'] = taxonkit_formatted_df['NCBI_taxon_ID'].astype(str)
    staxid_grouped_df['NCBI_taxon_ID'] = staxid_grouped_df['NCBI_taxon_ID'].astype(str)
    
    
    # ====== 合并数据
    species_df = pd.merge(taxonkit_formatted_df, staxid_grouped_df, on='NCBI_taxon_ID', how='inner')
    # taxonkit 运行之后 specie 会有重复，选择最高的去重复
    species_df['row_sum'] = species_df[sample_col].sum(axis=1)
    species_df = species_df.groupby('Taxon_ID').agg({
        col: 'sum' if col in sample_col else lambda x: x.iloc[species_df.loc[x.index, 'row_sum'].argmax()]
        for col in species_df.columns if col != 'row_sum'
    })
    # species_df.drop(columns=['row_sum'], inplace=True)
    write_output_df(species_df, os.path.join(args.output_dir, 'species_staxids_expression_result.txt'), index=False)
    # 生成各个分类级别的文件并统计非零值
    stats_df = generate_taxonomy_level_files(species_df, sample_col, args.output_dir)
    
    
    # 进行抽平
    flattening_dir = os.path.join(args.output_dir, 'flattening')
    os.makedirs(flattening_dir, exist_ok=True)
    asv_txt = os.path.join(flattening_dir, 'asv.txt')
    write_output_df(staxid_grouped_df[['ASV_ID'] + sample_col], asv_txt, index=False)
    # 运行 flattening 函数
    run_flattening_r(asv_txt, os.path.join(flattening_dir, 'asv_flattening.txt'))
    # 读取抽平后的结果
    flattened_df = load_table(os.path.join(flattening_dir, 'asv_flattening.txt'))
    # 获取非样本列
    non_sample_cols = [col for col in staxid_grouped_df.columns if col not in sample_col]
    # 合并抽平后的数据
    staxid_grouped_df = pd.merge(
        staxid_grouped_df[non_sample_cols],
        flattened_df,
        on='ASV_ID',
        how='inner'
    )
    write_output_df(staxid_grouped_df, os.path.join(args.output_dir, 'staxid_grouped_df.txt'), index=False)
       


if __name__ == '__main__':
    main()



