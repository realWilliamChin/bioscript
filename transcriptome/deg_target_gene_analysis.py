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
from data_check import df_replace_illegal_folder_chars
from Rscript import draw_multigroup_heatmap
from load_input import load_table, write_output_df



def parse_input():
    parser = argparse.ArgumentParser()
    # 输入文件
    parser.add_argument('-i', '--target-gene-file', required=True,
                        help='输入文件，target gene 文件，至少包含三列，GeneID 和 Ontology, SubOntology')
    parser.add_argument('--fpkm', required=True, help='输入 fpkm_matrix.txt 文件')
    parser.add_argument('-d', '--deg-data-dir', help='[必须]输入 DEG_data.txt 文件')
    parser.add_argument('-s', '--samplesinfo', required=True, help='输入样品信息文件')
    # 运行参数
    parser.add_argument('--mean', default=False, type=bool, help='使用每组中的平均数画 heatmap')
    parser.add_argument('--genesymbol', action='store_true', help='使用 GeneSymbol 作为索引，人和小鼠通常使用')
    # 输出
    parser.add_argument('-o', '--output', default=os.getcwd(), help='输出目录')
    
    args = parser.parse_args()

    return args


def genesymbol_data_pre_process(df: pd.DataFrame) -> pd.DataFrame:
    # 如果 GeneSymbol 有 NA 则会去除一些
    df_source_lines_num = df.shape[0]
    df['GeneSymbol'] = df['GeneSymbol'].replace(['NA', 'N/A', '-', '---', 'na', 'n/a', ''], np.nan)
    df.dropna(subset='GeneSymbol', inplace=True)
    df_dropedna_lines_num = df.shape[0]
    df.drop_duplicates(subset='GeneSymbol', inplace=True)
    df_droped_duplicates_lines_num = df.shape[0]
    if df_dropedna_lines_num != df_droped_duplicates_lines_num:
        logger.info(f'对输入表预处理之前为 {df_source_lines_num}')
        logger.info(f'去 NA 之后为 {df_dropedna_lines_num}')
        logger.info(f'去重复之后为 {df_droped_duplicates_lines_num}')
        logger.warning(f'输入文件有 GeneSymbol 为空或重复，将会去重再分析画图')
    
    return df


def apply_group_mean(df: pd.DataFrame, samples_df: pd.DataFrame):
    """
    对 df 按 samples_df 的 group 分组取均值，并调整 samples_df 结构。
    返回新的 df 和 samples_df。
    """
    for each_group in samples_df['group'].unique():
        group_samples = samples_df[samples_df['group'] == each_group]['sample']
        df[each_group] = df[group_samples].mean(axis=1)
    df = df.drop(columns=samples_df['sample'])
    samples_df = samples_df.drop(columns=['sample']).drop_duplicates(subset=['group'])
    samples_df['sample'] = samples_df['group']
    samples_df = samples_df[['sample', 'group']]
    return df, samples_df


def target_gene_heatmap(target_gene_df, fpkm_matrix_df, samples_df, index_col, group_mean, output_pic):
    sample_columns = samples_df['sample'].tolist()
    
    # 目标基因添加 fpkm 值，准备画图文件
    gene_fpkm_df = pd.merge(target_gene_df, fpkm_matrix_df, on='GeneID', how='inner')
    gene_fpkm_df.sort_values(by=['Ontology', 'SubOntology'], inplace=True)
    
    # 先画全部基因的热图（原有功能）
    all_gene_heatmap_filename = output_pic.replace('.jpg', '.xlsx')
    # multigroup_heatmap_data_df = gene_fpkm_df.loc[:, (gene_fpkm_df.columns != 'Ontology') & (gene_fpkm_df.columns != 'SubOntology')].copy()
    multigroup_heatmap_data_df = gene_fpkm_df[[index_col] + sample_columns].copy()
    multigroup_heatmap_sheet3_df = gene_fpkm_df[[index_col, 'SubOntology', 'Ontology']].copy()
    
    # if group mean 根据 samplesinfo 每组中的平均数画 heatmap
    if group_mean:
        multigroup_heatmap_data_df, samples_df = apply_group_mean(multigroup_heatmap_data_df, samples_df)

    with pd.ExcelWriter(all_gene_heatmap_filename, engine='openpyxl') as writer:
        multigroup_heatmap_data_df.to_excel(writer, sheet_name="Sheet1", index=False)
        samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
        multigroup_heatmap_sheet3_df.to_excel(writer, sheet_name='Sheet3', index=False)
    draw_multigroup_heatmap(all_gene_heatmap_filename, output_pic, other_args='--no-cluster-rows')

    # 每个 Ontology 单独画图
    for ontology, sub_df in gene_fpkm_df.groupby('Ontology'):
        logger.info(f'正在画 {ontology} heatmap')
        heatmap_data_df = sub_df[[index_col] + sample_columns].copy()
        heatmap_sheet3_df = sub_df[[index_col, 'SubOntology', 'Ontology']].copy()
        
        # if group mean 根据 samplesinfo 每组中的平均数画 heatmap
        if group_mean:
            heatmap_data_df, sub_samples_df = apply_group_mean(heatmap_data_df, sub_samples_df)
        
        ontology_excel = all_gene_heatmap_filename.replace('.xlsx', f'_{ontology}.xlsx')
        ontology_pic = output_pic.replace('.jpg', f'_{ontology}.jpg')
        with pd.ExcelWriter(ontology_excel, engine='openpyxl') as writer:
            heatmap_data_df.to_excel(writer, sheet_name="Sheet1", index=False)
            samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
            heatmap_sheet3_df.to_excel(writer, sheet_name='Sheet3', index=False)
        draw_multigroup_heatmap(ontology_excel, ontology_pic, other_args='--no-cluster-rows')


def deg_target_gene_summary(df_list, samples_info_df):
    """将每组目标基因相关的 DEG data 有表达的（Up 和 Down）合并

    Args:
        target_gene_data_list (list): target_gene_data list

    Returns:
        pd.DataFrame: 汇总文件
    """

    processed_df_list = []
    max_samples_number = samples_info_df.groupby('group').size().max()
    for df in df_list:
        df = df[df['regulation'].str.lower() != 'nosignificant']  # 只保留有表达的
        treat = df['sampleA'].values.tolist()[0]
        control = df['sampleB'].values.tolist()[0]
        treat_samples = samples_info_df[samples_info_df['group'] == treat]['sample'].values.tolist()
        control_samples = samples_info_df[samples_info_df['group'] == control]['sample'].values.tolist()
        
        # 获取原始的 FPKM 列和 reads 列
        df_treat_fpkm_column = [f'{x}_FPKM' for x in treat_samples]
        df_control_fpkm_column = [f'{x}_FPKM' for x in control_samples]
        df_treat_reads_column = [f'{x}_raw_reads' for x in treat_samples]
        df_control_reads_column = [f'{x}_raw_reads' for x in control_samples]
        
        # Treat FPKM
        for i in range(1, max_samples_number + 1):
            new_treat_fpkm = f"Treat_{i}_FPKM"
            if i <= len(df_treat_fpkm_column):
                df = df.rename(columns={df_treat_fpkm_column[i-1]: new_treat_fpkm})
            else:
                df[new_treat_fpkm] = 'N/A'

        # Control FPKM
        for i in range(1, max_samples_number + 1):
            new_control_fpkm = f"Control_{i}_FPKM"
            if i <= len(df_control_fpkm_column):
                df = df.rename(columns={df_control_fpkm_column[i-1]: new_control_fpkm})
            else:
                df[new_control_fpkm] = 'N/A'

        # Treat reads
        for i in range(1, max_samples_number + 1):
            new_treat_reads = f"Treat_{i}_reads"
            if i <= len(df_treat_reads_column):
                df = df.rename(columns={df_treat_reads_column[i-1]: new_treat_reads})
            else:
                df[new_treat_reads] = 'N/A'

        # Control reads
        for i in range(1, max_samples_number + 1):
            new_control_reads = f"Control_{i}_reads"
            if i <= len(df_control_reads_column):
                df = df.rename(columns={df_control_reads_column[i-1]: new_control_reads})
            else:
                df[new_control_reads] = 'N/A'
        
        processed_df_list.append(df)
    
    output_summary_df = pd.concat(processed_df_list)
    
    return output_summary_df



def deg_target_gene_heatmap(target_gene_def_df, samples_df, deg_data_dir, index_col, output_dir):
    result_target_gene_data_list = []
    samples_df = samples_df[['sample', 'group']]
    deg_data_list = [x for x in os.listdir(deg_data_dir) if x.endswith('_DEG_data.txt')]
    for deg_data_file in deg_data_list:
        comparison_name = os.path.basename(deg_data_file).replace('_DEG_data.txt', '')
        logger.info(f'正在画 {comparison_name} heatmap')
        
        group1 = comparison_name.split('-vs-')[0]
        group2 = comparison_name.split('-vs-')[1]
        comparison_samples_df = samples_df[samples_df['group'].isin([group1, group2])]
        comparison_samples_list = comparison_samples_df['sample'].values.tolist()
        comparison_dir = os.path.join(output_dir, comparison_name)
        os.makedirs(comparison_dir, exist_ok=True)
        
        logger.info(f'正在找相关基因添加定义 {comparison_name}')
        deg_data_file = os.path.join(deg_data_dir, deg_data_file)
        deg_data_df = load_table(deg_data_file, dtype={'GeneID': str})
        
        # 目标基因和 DEG_data 合并，文件预处理
        comparison_target_gene_df = pd.merge(target_gene_def_df, right=deg_data_df, on=index_col, how='inner', suffixes=('_df1', '_df2'))
        comparison_target_gene_df.dropna(subset=index_col, inplace=True)  # 去掉 NA 行
        cols_to_drop = [col for col in comparison_target_gene_df.columns if col.endswith('_df1')]
        comparison_target_gene_df.drop(columns=cols_to_drop, inplace=True)
        comparison_target_gene_df.columns = [col.replace('_df2', '') for col in comparison_target_gene_df.columns]
        # 添加到 list，输出汇总文件
        result_target_gene_data_list.append(comparison_target_gene_df.copy())
        
        comparison_target_gene_df.columns = [col[:-5] if col.endswith('_FPKM') else col for col in comparison_target_gene_df.columns]
        comparison_target_gene_file = os.path.join(comparison_dir, f'{comparison_name}_target_gene_data.xlsx')
        write_output_df(comparison_target_gene_df, comparison_target_gene_file, index=False)
        
        fpkm_df = comparison_target_gene_df[[index_col] + comparison_samples_list]
        ontology_df = comparison_target_gene_df[[index_col, 'SubOntology', 'Ontology']]
        
        group_vs_group_heatmap_fname = os.path.join(comparison_dir, f'{comparison_name}_target_gene_heatmap.xlsx')
        group_vs_group_heatmap_pname = os.path.join(comparison_dir, f'{comparison_name}_target_gene_heatmap.jpg')

        with pd.ExcelWriter(group_vs_group_heatmap_fname, engine='openpyxl') as writer:
            fpkm_df.to_excel(writer, sheet_name='Sheet1', index=False)
            comparison_samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
            ontology_df.to_excel(writer, sheet_name='Sheet3', index=False)
        
        draw_multigroup_heatmap(group_vs_group_heatmap_fname, group_vs_group_heatmap_pname, other_args='--no-cluster-rows')
        
        # 每个组间的每个 Ontology 画一张图
        for ontology, sub_df in comparison_target_gene_df.groupby('Ontology'):
            logger.info(f'正在画 {comparison_name} {ontology} heatmap')
            ontology_heatmap_data_df = sub_df[[index_col] + comparison_samples_list].copy()
            ontology_def_df = sub_df[[index_col, 'SubOntology', 'Ontology']].copy()

            ontology_excel_name = os.path.join(comparison_dir, f'{comparison_name}_{ontology}_heatmap.xlsx')
            ontology_pic_name = os.path.join(comparison_dir, f'{comparison_name}_{ontology}_heatmap.jpg')
            with pd.ExcelWriter(ontology_excel_name, engine='openpyxl') as writer:
                ontology_heatmap_data_df.to_excel(writer, sheet_name="Sheet1", index=False)
                comparison_samples_df.to_excel(writer, sheet_name="Sheet2", index=False)
                ontology_def_df.to_excel(writer, sheet_name="Sheet3", index=False)

            draw_multigroup_heatmap(ontology_excel_name, ontology_pic_name, other_args='--no-cluster-rows')
            
    # 超过两组比较, 汇总结果
    if len(result_target_gene_data_list) >= 2:
        logger.info('正在对结果汇总')
        target_gene_summary_df = deg_target_gene_summary(result_target_gene_data_list, samples_df)
        write_output_df(target_gene_summary_df, os.path.join(output_dir, 'Target_gene_summary_data.xlsx'), index=False)
    else:
        logger.info('跳过结果汇总，比较组少于 2 个')


def main():
    # 数据加载预处理
    args = parse_input()
    output_dir = args.output
    target_gene_file, samples_file, fpkm_matrix_file = args.target_gene_file, args.samplesinfo, args.fpkm
    
    # 加载为 DataFrame
    samples_df = load_table(samples_file, usecols=[0, 1], dtype=str)
    target_gene_def_df = load_table(target_gene_file, dtype={'GeneID': str})
    target_gene_def_df['Ontology'] = target_gene_def_df['Ontology'].str.replace('/', '_').str.replace(' ', '_')
    target_gene_def_df = df_drop_element_side_space(target_gene_def_df)
    target_gene_def_df = df_replace_illegal_folder_chars(target_gene_def_df, ['Ontology', 'SubOntology'])
    
    # 对文件预处理
    samples_df = samples_df[['sample', 'group']]
    if args.genesymbol:
        index_col = 'GeneSymbol'
        target_gene_def_df = genesymbol_data_pre_process(target_gene_def_df)
    else:
        index_col = 'GeneID'
        
    fpkm_matrix_df = load_table(fpkm_matrix_file, dtype={'GeneID': str})
    fpkm_matrix_df = df_drop_row_sum_eq_zero(fpkm_matrix_df)
    
    all_gene_heatmap_pic = os.path.join(output_dir, target_gene_file.split(os.sep)[-1].replace('.txt', '_heatmap.jpg'))
    target_gene_heatmap(target_gene_def_df, fpkm_matrix_df, samples_df, index_col, args.mean, all_gene_heatmap_pic)
    if args.deg_data_dir:
        deg_target_gene_heatmap(target_gene_def_df, samples_df, args.deg_data_dir, index_col, args.output)
    
    logger.success("Done")


if __name__ == '__main__':
    main()