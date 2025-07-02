#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/05/13 17:04
# Author        : William GoGo
import os
import sys
import pandas as pd
import numpy as np
import argparse
import subprocess
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from data_check import df_drop_row_sum_eq_zero
from data_check import df_drop_element_side_space
from Rscript import draw_multigroup_heatmap
from load_input import load_table, write_output_df



def parse_input():
    parser = argparse.ArgumentParser()
    # 输入文件
    parser.add_argument('-i', '--targetgene', dest="target_gene_file", type=str, required=True,
                           help='输入文件，target gene 文件，至少包含两列，GeneID 和 Ontology')
    parser.add_argument('--fpkm', type=str, required=True, help='输入 fpkm_matrix.txt 文件')
    parser.add_argument('--kns', type=str, help='anova_p 文件添加定义输出')
    parser.add_argument('--deg-data-dir', dest='deg_data_dir', type=str, help='[必须]输入 DEG_data.txt 文件')
    parser.add_argument('-s', '--samplesinfo', type=str, required=True,
                        help='输入样品信息文件')
    # 运行参数
    parser.add_argument('--mean', default=False, type=bool, help='使用每组中的平均数画 heatmap')
    parser.add_argument('--specie-type', dest='specie_type',
                        choices=['homo_sapiens', 'mus_musculus', 'other'], default='other',
                        help='输入运行的物种类型，人，小鼠，或者其他，人小鼠使用 GeneSymbol 作为索引画图')
    # 输出
    parser.add_argument('-o', '--output', type=str, default=os.getcwd(), help='输出目录')
    
    args = parser.parse_args()

    return args


def group_vs_group_heatmap(group_target_gene_file, samples_file, specie_type, output_dir):
    df = load_table(group_target_gene_file, dtype={"GeneID": str})
    if specie_type in ['homo_sapiens', 'mus_musculus']:
        index_col = 'GeneSymbol'
    else:
        index_col = 'GeneID'
    df.dropna(subset=index_col, inplace=True)
    fpkm_df = df[[index_col] + [col for col in df.columns if col.endswith('FPKM')]]
    no_FPKM_str_columns = [col.replace('_FPKM', '') for col in fpkm_df.columns[1:]]
    fpkm_df.columns = [index_col] + no_FPKM_str_columns
    ontology_df = df[[index_col, 'SubOntology', 'Ontology']]
    
    os.makedirs(os.path.join(output_dir, 'Prep_files'), exist_ok=True)
    group_vs_group_heatmap_fname = os.path.join(
        output_dir,
        'Prep_files',
        os.path.basename(group_target_gene_file).replace('_data.txt', '_heatmap.xlsx')
    )
    group_vs_group_heatmap_pname = group_target_gene_file.replace('_data.txt', '_heatmap.jpeg')
    
    group1 = df.head(1)['sampleA'].values[0]
    group2 = df.head(1)['sampleB'].values[0]
    samples_df = pd.read_csv(samples_file, sep='\t', usecols=[0, 1])
    samples_df = samples_df[['sample', 'group']]
    samples_df = samples_df[samples_df['group'].isin([group1, group2])]
    
    with pd.ExcelWriter(group_vs_group_heatmap_fname, engine='openpyxl') as writer:
        fpkm_df.to_excel(writer, sheet_name='Sheet1', index=False)
        samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
        ontology_df.to_excel(writer, sheet_name='Sheet3', index=False)
    
    draw_multigroup_heatmap(group_vs_group_heatmap_fname, group_vs_group_heatmap_pname, other_args='--no-cluster-rows')


def deg_target_gene_summary(target_gene_data_list):
    """将每组目标基因相关的 DEG data 有表达的（Up 和 Down）合并

    Args:
        target_gene_data_list (list): target_gene_data list

    Returns:
        pd.DataFrame: 汇总文件
    """
    filtered_df_list = []
    for each_target_gene_data in target_gene_data_list:
        target_gene_data_df = load_table(each_target_gene_data, dtype={'GeneID': str})
        if 'regulation' not in target_gene_data_df.columns:
            logger.warning(f"文件 {each_target_gene_data} 中缺少 regulation 列")
            continue
        target_gene_data_df = target_gene_data_df[target_gene_data_df['regulation'].str.lower() != 'nosignificant']
        # 去掉列名包含 _FPKM 和 _reads 的列
        columns_to_keep = [col for col in target_gene_data_df.columns if '_FPKM' not in col and '_reads' not in col]
        target_gene_data_df = target_gene_data_df[columns_to_keep]
        filtered_df_list.append(target_gene_data_df)
    
    summary_df = pd.concat(filtered_df_list)
    
    return summary_df


def genesymbol_data_pre_process(df: pd.DataFrame) -> pd.DataFrame:
    # 如果 GeneSymbol 有 NA 则会去除一些
    df_source_lines_num = df.shape[0]
    df.dropna(subset='GeneSymbol', inplace=True)
    df['GeneSymbol'] = df['GeneSymbol'].replace(['NA', 'N/A', '-', '---', 'na', 'n/a', ''], np.nan)
    df_dropedna_lines_num = df.shape[0]
    df.drop_duplicates(subset='GeneSymbol', inplace=True)
    df_droped_duplicates_lines_num = df.shape[0]
    if df_dropedna_lines_num != df_droped_duplicates_lines_num:
        logger.info(f'对输入表预处理之前为 {df_source_lines_num}')
        logger.info(f'去 NA 之后为 {df_dropedna_lines_num}')
        logger.info(f'去重复之后为 {df_droped_duplicates_lines_num}')
        logger.warning(f'输入文件有 GeneSymbol 为空或重复，将会去重再分析画图')
    
    return df


def target_gene_heatmap(target_gene_df, fpkm_matrix_df, samples_df, specie_type, group_mean, output_pic):
    sample_columns = samples_df['sample'].tolist()
    
    # 目标基因添加 fpkm 值，准备画图文件
    gene_fpkm_df = pd.merge(target_gene_df, fpkm_matrix_df, on='GeneID', how='inner')
    gene_fpkm_df.sort_values(by=['Ontology', 'SubOntology'], inplace=True)
    
    # 先画全部基因的热图（原有功能）
    all_gene_heatmap_filename = output_pic.replace('.jpeg', '.xlsx')
    # multigroup_heatmap_data_df = gene_fpkm_df.loc[:, (gene_fpkm_df.columns != 'Ontology') & (gene_fpkm_df.columns != 'SubOntology')].copy()
    if specie_type in ['homo_sapiens', 'mus_musculus']:
        multigroup_heatmap_data_df = gene_fpkm_df[['GeneSymbol'] + sample_columns].copy()
        multigroup_heatmap_sheet3_df = gene_fpkm_df[['GeneSymbol', 'SubOntology', 'Ontology']].copy()
    else:
        multigroup_heatmap_data_df = gene_fpkm_df[['GeneID'] + sample_columns].copy()
        multigroup_heatmap_sheet3_df = gene_fpkm_df[['GeneID', 'SubOntology', 'Ontology']].copy()
    
    # 根据 samplesinfo 每组中的平均数画 heatmap
    if group_mean:
        for each_group in samples_df['group'].unique():
            group_samples = samples_df[samples_df['group'] == each_group]['sample']
            multigroup_heatmap_data_df[each_group] = multigroup_heatmap_data_df[group_samples].mean(axis=1)
        multigroup_heatmap_data_df.drop(columns=samples_df['sample'], inplace=True)
        samples_df = samples_df.drop(columns=['sample']).drop_duplicates(subset=['group'])
        samples_df['sample'] = samples_df['group']
        samples_df = samples_df[['sample', 'group']]

    with pd.ExcelWriter(all_gene_heatmap_filename, engine='openpyxl') as writer:
        multigroup_heatmap_data_df.to_excel(writer, sheet_name="Sheet1", index=False)
        samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
        multigroup_heatmap_sheet3_df.to_excel(writer, sheet_name='Sheet3', index=False)
    draw_multigroup_heatmap(all_gene_heatmap_filename, output_pic, other_args='--no-cluster-rows')

    # 每个 Ontology 单独画图
    for ontology, sub_df in gene_fpkm_df.groupby('Ontology'):
        if specie_type in ['homo_sapiens', 'mus_musculus']:
            sub_df = sub_df.dropna(subset=['GeneSymbol'])
            heatmap_data_df = sub_df[['GeneSymbol'] + sample_columns].copy()
            heatmap_sheet3_df = sub_df[['GeneSymbol', 'SubOntology', 'Ontology']].copy()
        else:
            heatmap_data_df = sub_df[['GeneID'] + sample_columns].copy()
            heatmap_sheet3_df = sub_df[['GeneID', 'SubOntology', 'Ontology']].copy()
        sub_samples_df = samples_df.copy()
        if group_mean:
            for each_group in sub_samples_df['group'].unique():
                group_samples = sub_samples_df[sub_samples_df['group'] == each_group]['sample']
                heatmap_data_df[each_group] = heatmap_data_df[group_samples].mean(axis=1)
            heatmap_data_df.drop(columns=sub_samples_df['sample'], inplace=True)
            sub_samples_df = sub_samples_df.drop(columns=['sample']).drop_duplicates(subset=['group'])
            sub_samples_df['sample'] = sub_samples_df['group']
            sub_samples_df = sub_samples_df[['sample', 'group']]
        # 处理文件名中的非法字符
        safe_ontology = str(ontology).replace('/', '_').replace(' ', '_')
        ontology_excel = all_gene_heatmap_filename.replace('.xlsx', f'_{safe_ontology}.xlsx')
        ontology_pic = output_pic.replace('.jpeg', f'_{safe_ontology}.jpeg')
        with pd.ExcelWriter(ontology_excel, engine='openpyxl') as writer:
            heatmap_data_df.to_excel(writer, sheet_name="Sheet1", index=False)
            sub_samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
            heatmap_sheet3_df.to_excel(writer, sheet_name='Sheet3', index=False)
        draw_multigroup_heatmap(ontology_excel, ontology_pic, other_args='--no-cluster-rows')


def main():
    # 数据加载预处理
    args = parse_input()
    output_dir = args.output
    target_gene_file, samples_file, fpkm_matrix_file = args.target_gene_file, args.samplesinfo, args.fpkm
    
    # 加载为 DataFrame
    samples_df = load_table(samples_file, usecols=[0, 1], dtype=str)
    target_gene_def_df = load_table(target_gene_file, dtype={'GeneID': str})
    
    # 对文件预处理
    samples_df = samples_df[['sample', 'group']]
    target_gene_def_df = df_drop_element_side_space(target_gene_def_df)
    if args.specie_type.lower() in ['homo_sapiens', 'mus_musculus']:
        target_gene_def_df = genesymbol_data_pre_process(target_gene_def_df)
        
    fpkm_matrix_df = load_table(fpkm_matrix_file, dtype={'GeneID': str})
    fpkm_matrix_df = df_drop_row_sum_eq_zero(fpkm_matrix_df)
    
    all_gene_heatmap_pic = os.path.join(output_dir, target_gene_file.split(os.sep)[-1].replace('.txt', '_heatmap.jpeg'))
    target_gene_heatmap(target_gene_def_df, fpkm_matrix_df, samples_df, args.specie_type, args.mean, all_gene_heatmap_pic)

    
    result_target_gene_data_list = []
    if args.deg_data_dir:
        deg_data_list = [x for x in os.listdir(args.deg_data_dir) if x.endswith('_DEG_data.txt')]
        for deg_data_file in deg_data_list:
            compare_name = os.path.basename(deg_data_file).replace('_DEG_data.txt', '')
            logger.info(f'正在找相关基因添加定义 {compare_name}')
            deg_data_file = os.path.join(args.deg_data_dir, deg_data_file)
            deg_data_df = load_table(deg_data_file, dtype={'GeneID': str})

            # result_df = pd.merge(deg_data_df, right=target_gene_def_df, on='GeneID', how='left', suffixes=('_df1', '_df2'))
            result_df = pd.merge(target_gene_def_df, right=deg_data_df, on='GeneID', how='inner', suffixes=('_df1', '_df2'))
            cols_to_drop = [col for col in result_df.columns if col.endswith('_df1')]
            result_df.drop(columns=cols_to_drop, inplace=True)
            result_df.columns = [col.replace('_df2', '') for col in result_df.columns]
            result_df.to_csv(os.path.join(args.output, f'{compare_name}_target_gene_data.txt'), sep='\t', index=False)
            result_target_gene_data_list.append(os.path.join(args.output, f'{compare_name}_target_gene_data.txt'))
            
            logger.info(f'正在画 {compare_name} heatmap')
            group_vs_group_heatmap(
                os.path.join(args.output, f'{compare_name}_target_gene_data.txt'),
                samples_file,
                args.specie_type,
                args.output
            )
    
    if len(result_target_gene_data_list) > 0:
        logger.info('正在对结果汇总')
        target_gene_summary_df = deg_target_gene_summary(result_target_gene_data_list)
        write_output_df(target_gene_summary_df, os.path.join(args.output, 'Target_gene_summary_data.txt'), index=False)
    
    logger.success("Done")


if __name__ == '__main__':
    main()