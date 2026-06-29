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

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from data_check import df_drop_row_sum_eq_zero
from data_check import df_drop_element_side_space
from data_check import df_replace_illegal_folder_chars
from r_wrapper import smart_heatmap
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
    parser.add_argument('--genesymbol', help='使用 GeneSymbol 作为索引画热图，人和大鼠小鼠通常使用，输入包含 GeneID GeneSymbol 两列的文件')
    
    # 输出
    parser.add_argument('-o', '--output-dir', default=os.getcwd(), help='输出目录')
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)

    return args


def apply_group_mean(df: pd.DataFrame, samples_df: pd.DataFrame):
    """
    对 df 按 samples_df 的 group 分组取均值，并调整 samples_df 结构。
    返回新的 df 和 samples_df。
    """
    for each_group in samples_df['group'].unique():
        group_samples = samples_df[samples_df['group'] == each_group]['sample'].tolist()
        # 只使用在 df 中存在的样本列
        existing_samples = [s for s in group_samples if s in df.columns]
        if not existing_samples:
            logger.warning(f'组 {each_group} 的样本列在数据中不存在，跳过该组')
            continue
        df[each_group] = df[existing_samples].mean(axis=1)
    # 只删除在 df 中存在的样本列
    existing_sample_cols = [col for col in samples_df['sample'].tolist() if col in df.columns]
    if existing_sample_cols:
        df = df.drop(columns=existing_sample_cols)
    samples_df = samples_df.drop(columns=['sample']).drop_duplicates(subset=['group'])
    samples_df['sample'] = samples_df['group']
    samples_df = samples_df[['sample', 'group']]
    return df, samples_df


def target_gene_heatmap(target_gene_df, fpkm_matrix_df, samples_df, group_mean, output_dir, genesymbol_df=None):
    sample_columns = samples_df['sample'].tolist()
    # 保存原始的 samples_df，因为 apply_group_mean 会修改它
    original_samples_df = samples_df.copy()
    
    # 目标基因添加 fpkm 值，准备画图文件
    gene_fpkm_df = pd.merge(target_gene_df, fpkm_matrix_df, on='GeneID', how='inner')
    
    # 如果提供了 genesymbol_df，则替换 GeneID 为 GeneSymbol（仅用于画 heatmap）
    if genesymbol_df is not None:
        # 如果 gene_fpkm_df 中已经存在 'GeneSymbol' 列，先删除它以避免合并时产生后缀
        if 'GeneSymbol' in gene_fpkm_df.columns:
            gene_fpkm_df = gene_fpkm_df.drop(columns=['GeneSymbol'])
        # 合并 GeneSymbol 映射
        gene_fpkm_df = pd.merge(gene_fpkm_df, genesymbol_df, on='GeneID', how='inner')
        # 将 GeneID 列替换为 GeneSymbol
        gene_fpkm_df['GeneID'] = gene_fpkm_df['GeneSymbol']
        gene_fpkm_df.drop(columns=['GeneSymbol'], inplace=True)
    gene_fpkm_df.sort_values(by=['Ontology', 'SubOntology'], inplace=True)
    
    # 先画全部基因的热图（原有功能）
    all_gene_heatmap_filename = os.path.join(output_dir, '00_All_target_gene_heatmap.xlsx')
    multigroup_heatmap_data_df = gene_fpkm_df[['GeneID'] + sample_columns].copy()
    multigroup_heatmap_sheet3_df = gene_fpkm_df[['GeneID', 'Ontology']].copy()
    
    # if group mean 根据 samplesinfo 每组中的平均数画 heatmap
    if group_mean:
        multigroup_heatmap_data_df, samples_df = apply_group_mean(multigroup_heatmap_data_df, samples_df)

    with pd.ExcelWriter(all_gene_heatmap_filename, engine='openpyxl') as writer:
        multigroup_heatmap_data_df.to_excel(writer, sheet_name="Sheet1", index=False)
        samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
        multigroup_heatmap_sheet3_df.to_excel(writer, sheet_name='Sheet3', index=False)
    smart_heatmap(
        input_file=all_gene_heatmap_filename,
        output_file=os.path.join(output_dir, '00_All_target_gene_heatmap.jpg'),
        annotation_col=2,
        annotation_row=3,
        cluster_rows=False,
        cluster_cols=False,
        scale='row'
    )   

    # 每个 Ontology 单独画图
    for ontology, sub_df in gene_fpkm_df.groupby('Ontology'):
        logger.info(f'正在画 {ontology} heatmap')
        heatmap_data_df = sub_df[['GeneID'] + sample_columns].copy()
        heatmap_sheet3_df = sub_df[['GeneID', 'SubOntology']].copy()
        
        # if group mean 根据 samplesinfo 每组中的平均数画 heatmap
        # 使用原始的 samples_df，因为 heatmap_data_df 包含的是原始样本列名
        if group_mean:
            heatmap_data_df, sub_samples_df = apply_group_mean(heatmap_data_df, original_samples_df.copy())
        else:
            sub_samples_df = original_samples_df
        
        ontology_excel = os.path.join(output_dir, f'{ontology}_target_gene_heatmap.xlsx')
        ontology_pic = os.path.join(output_dir, f'{ontology}_target_gene_heatmap.jpg')
        with pd.ExcelWriter(ontology_excel, engine='openpyxl') as writer:
            heatmap_data_df.to_excel(writer, sheet_name="Sheet1", index=False)
            sub_samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
            heatmap_sheet3_df.to_excel(writer, sheet_name='Sheet3', index=False)
        smart_heatmap(
            input_file=ontology_excel,
            output_file=ontology_pic,
            annotation_col=2,
            annotation_row=3,
            cluster_rows=False,
            cluster_cols=False,
            scale='row',
            main=f'Ontology: {ontology}'
        )


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
        # 若过滤后无显著差异，跳过该比较组
        if df.empty:
            logger.warning('本比较组无显著差异目标基因，跳过汇总')
            continue
        if 'sampleA' not in df.columns or 'sampleB' not in df.columns:
            logger.warning('缺少 sampleA 或 sampleB 列，跳过该比较组汇总')
            continue
        treat = str(df['sampleA'].iloc[0])
        control = str(df['sampleB'].iloc[0])
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
    
    if len(processed_df_list) == 0:
        return pd.DataFrame()
    output_summary_df = pd.concat(processed_df_list)
    # 删掉一列全部是 'N/A' 的（只有两组且数量不一样会有这种情况）
    output_summary_df = output_summary_df.loc[:, ~(output_summary_df == 'N/A').all()]    
    return output_summary_df


def deg_target_gene_heatmap(target_gene_def_df, samples_df, deg_data_dir, output_dir, genesymbol_df=None):
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
        
        # 目标基因和 DEG_data 合并，文件预处理（使用 GeneID）
        comparison_target_gene_df = pd.merge(target_gene_def_df, right=deg_data_df, on='GeneID', how='inner', suffixes=('_df1', ''))
        comparison_target_gene_df.dropna(subset='GeneID', inplace=True)  # 去掉 NA 行
        cols_to_drop = [col for col in comparison_target_gene_df.columns if col.endswith('_df1')]
        comparison_target_gene_df.drop(columns=cols_to_drop, inplace=True)
        # 添加到 list，输出汇总文件
        result_target_gene_data_list.append(comparison_target_gene_df.copy())
        
        comparison_target_gene_df.columns = [col[:-5] if col.endswith('_FPKM') else col for col in comparison_target_gene_df.columns]
        comparison_target_gene_file = os.path.join(comparison_dir, f'00_{comparison_name}_target_gene_data.xlsx')
        write_output_df(comparison_target_gene_df, comparison_target_gene_file, index=False)
        
        # 如果提供了 genesymbol_df，则替换 GeneID 为 GeneSymbol（仅用于 heatmap）
        comparison_target_gene_df_for_heatmap = comparison_target_gene_df.copy()
        if genesymbol_df is not None:
            # 如果 comparison_target_gene_df_for_heatmap 中已经存在 'GeneSymbol' 列，先删除它以避免合并时产生后缀
            if 'GeneSymbol' in comparison_target_gene_df_for_heatmap.columns:
                comparison_target_gene_df_for_heatmap = comparison_target_gene_df_for_heatmap.drop(columns=['GeneSymbol'])
            # 合并 GeneSymbol 映射
            comparison_target_gene_df_for_heatmap = pd.merge(comparison_target_gene_df_for_heatmap, genesymbol_df, on='GeneID', how='inner')
            if comparison_target_gene_df_for_heatmap.shape[0] < 2:
                logger.warning(f'{comparison_name} 相关 genesymbol 数量小于 2，跳过画 Heatmap 图')
                continue
            # 将 GeneID 列替换为 GeneSymbol
            comparison_target_gene_df_for_heatmap['GeneID'] = comparison_target_gene_df_for_heatmap['GeneSymbol']
            comparison_target_gene_df_for_heatmap.drop(columns=['GeneSymbol'], inplace=True)
        
        # 使用 GeneID 作为索引列
        fpkm_df = comparison_target_gene_df_for_heatmap[['GeneID'] + comparison_samples_list]
        fpkm_df = df_drop_row_sum_eq_zero(fpkm_df)
        ontology_df = comparison_target_gene_df_for_heatmap[['GeneID', 'SubOntology', 'Ontology']]
        
        group_vs_group_heatmap_fname = os.path.join(comparison_dir, f'{comparison_name}_target_gene_heatmap.xlsx')
        group_vs_group_heatmap_pname = os.path.join(comparison_dir, f'{comparison_name}_target_gene_heatmap.jpg')

        with pd.ExcelWriter(group_vs_group_heatmap_fname, engine='openpyxl') as writer:
            fpkm_df.to_excel(writer, sheet_name='Sheet1', index=False)
            comparison_samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
            ontology_df.to_excel(writer, sheet_name='Sheet3', index=False)
        
        smart_heatmap(
            input_file=group_vs_group_heatmap_fname,
            output_file=group_vs_group_heatmap_pname,
            annotation_col=2,
            annotation_row=3,
            cluster_rows=False,
            cluster_cols=False,
            scale='row'
        )
        
        # 每个组间的每个 Ontology 画一张图
        for ontology, sub_df in comparison_target_gene_df_for_heatmap.groupby('Ontology'):
            logger.info(f'正在画 {comparison_name} {ontology} heatmap')
            ontology_heatmap_data_df = sub_df[['GeneID'] + comparison_samples_list].copy()
            ontology_heatmap_data_df = df_drop_row_sum_eq_zero(ontology_heatmap_data_df)
            ontology_def_df = sub_df[['GeneID', 'SubOntology']].copy()
            ontology_def_df = ontology_def_df[ontology_def_df['GeneID'].isin(ontology_heatmap_data_df['GeneID'])]
            if ontology_heatmap_data_df.shape[0] < 2:
                logger.warning(f'{comparison_name} {ontology} 中目标基因数量小于 2，不画 heatmap')
                continue
            ontology_excel_name = os.path.join(comparison_dir, f'{comparison_name}_{ontology}_heatmap.xlsx')
            ontology_pic_name = os.path.join(comparison_dir, f'{comparison_name}_{ontology}_heatmap.jpg')
            with pd.ExcelWriter(ontology_excel_name, engine='openpyxl') as writer:
                ontology_heatmap_data_df.to_excel(writer, sheet_name="Sheet1", index=False)
                comparison_samples_df.to_excel(writer, sheet_name="Sheet2", index=False)
                ontology_def_df.to_excel(writer, sheet_name="Sheet3", index=False)

            smart_heatmap(
                input_file=ontology_excel_name,
                output_file=ontology_pic_name,
                annotation_col=2,
                annotation_row=3,
                cluster_rows=False,
                cluster_cols=False,
                scale='row',
                main=f'Ontology: {ontology}'
            )
            
    logger.info('正在对结果汇总')
    target_gene_summary_df = deg_target_gene_summary(result_target_gene_data_list, samples_df)
    write_output_df(target_gene_summary_df, os.path.join(output_dir, '00_Target_gene_significant_DEG_summary.xlsx'), index=False)


def _each_ontology_deg_data_summary(df_list, samples_info_df):
    """合并所有组的 Ontology DEG 数据，统一列名格式
    
    Args:
        df_list: 包含所有比较组的 Ontology DEG 数据框列表
        samples_info_df: 样本信息数据框，包含 sample 和 group 列
        
    Returns:
        pd.DataFrame: 合并后的汇总数据框
    """
    processed_df_list = []
    max_samples_number = samples_info_df.groupby('group').size().max()
    for df in df_list:
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


def each_ontology_deg_data(target_gene_def_df, samples_df, deg_data_dir, output_dir):
    """
    功能:
        针对给定的目标基因定义表（target_gene_def_df），遍历差异表达基因（DEG）数据目录中的每一个 DEG_data.txt 文件，
        按 Ontology 对目标基因进行分组，并筛选出每个本体（Ontology）下在差异表达基因文件中出现的目标基因，
        合并所有组的数据，然后将对应的 DEG 数据分别输出为 txt 文件。

    输入参数:
        target_gene_def_df : pd.DataFrame
            目标基因定义表，至少包含 'GeneID' 列及 'Ontology' 列。
        samples_df : pd.DataFrame
            样本信息数据框，包含 sample 和 group 列。
        deg_data_dir : str
            DEG_data.txt 文件所在目录。
        output_dir : str
            输出目录，生成各组本体 DEG 文件。

    输出:
        每个 Ontology 生成一个 txt 文件，文件名格式为：
            {ontology}.txt
        每个文件包含所有比较组下属于该本体的所有目标基因 DEG 信息。
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 获取所有 DEG 数据文件
    deg_data_files = [f for f in os.listdir(deg_data_dir) if f.endswith("_DEG_data.txt") and not f.startswith("~")]
    
    # 按 Ontology 分组处理
    ontology_list = list(set(target_gene_def_df['Ontology'].tolist()))
    
    for ontology in ontology_list:
        logger.info(f'正在处理 Ontology: {ontology}')
        ontology_df = target_gene_def_df[target_gene_def_df['Ontology'] == ontology].copy()      
        
        output_summary_df_list = []
        # 遍历所有比较组
        for deg_file in deg_data_files:
            deg_file_path = os.path.join(deg_data_dir, deg_file)
            comparison_name = deg_file.replace("_DEG_data.txt", '')
            deg_data_df = load_table(deg_file_path, dtype={"GeneID": str})
            
            if 'sampleA' not in deg_data_df.columns or 'sampleB' not in deg_data_df.columns:
                logger.warning(f"{deg_file} 没有包含 sampleA 或 sampleB 列，跳过")
                continue
            
            # 筛选出属于当前 Ontology 的基因
            ontology_deg_data_df = pd.merge(ontology_df, deg_data_df, on='GeneID', how='inner')
            
            if ontology_deg_data_df.shape[0] == 0:
                logger.warning(f'{comparison_name} 的 {ontology} 未找到相关基因')
                continue
            
            output_summary_df_list.append(ontology_deg_data_df)
        
        # 输出汇总
        if len(output_summary_df_list) == 0:
            logger.warning(f'Ontology {ontology} 在组间中没有任何目标基因')
        else:
            output_summary_df = _each_ontology_deg_data_summary(output_summary_df_list, samples_df)
            
            # 生成文件名
            output_filename = os.path.join(output_dir, f"{ontology}.txt")
            write_output_df(output_summary_df, output_filename, index=False)
            logger.info(f'已输出 Ontology {ontology} 的数据到文件: {output_filename}')


def main():
    # 数据加载预处理
    args = parse_input()
    
    # 加载为 DataFrame
    samples_df = load_table(args.samplesinfo, usecols=[0, 1], dtype=str)
    target_gene_def_df = load_table(args.target_gene_file, dtype={'GeneID': str})
    target_gene_def_df['Ontology'] = target_gene_def_df['Ontology'].str.replace('/', '_').str.replace(' ', '_')
    target_gene_def_df = df_drop_element_side_space(target_gene_def_df)
    target_gene_def_df = df_replace_illegal_folder_chars(target_gene_def_df, ['Ontology', 'SubOntology'])
    
    # 对文件预处理
    samples_df = samples_df[['sample', 'group']]
    
    # 如果提供了 genesymbol 参数，读取 GeneID 和 GeneSymbol 映射文件（仅用于 heatmap）
    genesymbol_df = None
    if args.genesymbol:
        genesymbol_df = load_table(args.genesymbol)
        # 只保留 GeneID 和 GeneSymbol 两列
        genesymbol_df = genesymbol_df[['GeneID', 'GeneSymbol']].copy()
        # 删除空值
        genesymbol_df.dropna(subset=['GeneID', 'GeneSymbol'], inplace=True)
        logger.info(f'读取到 {genesymbol_df.shape[0]} 个 GeneID-GeneSymbol 映射')
        
    fpkm_matrix_df = load_table(args.fpkm, dtype={'GeneID': str})
    fpkm_matrix_df = df_drop_row_sum_eq_zero(fpkm_matrix_df)
    
    output_dir = os.path.join(args.output_dir, '00_Each_Target_gene_ontology_heatmap')
    os.makedirs(output_dir, exist_ok=True)
    target_gene_heatmap(target_gene_def_df, fpkm_matrix_df, samples_df, args.mean, os.path.join(args.output_dir, '00_Each_Target_gene_ontology_heatmap'), genesymbol_df)
    if args.deg_data_dir:
        deg_target_gene_heatmap(target_gene_def_df, samples_df, args.deg_data_dir, args.output_dir, genesymbol_df)
        # 如果 target_gene_def_df 中存在 Ontology 列则运行 each_ontology_deg_data 函数
        if 'Ontology' in target_gene_def_df.columns:
            each_ontology_deg_data_dir = os.path.join(args.output_dir, '00_Each_Target_gene_ontology_DEG_data')
            os.makedirs(each_ontology_deg_data_dir, exist_ok=True)
            each_ontology_deg_data(target_gene_def_df, samples_df, args.deg_data_dir, each_ontology_deg_data_dir)
        else:
            logger.warning('target_gene_def_df 中不存在 Ontology 列，跳过拆分每个本体 DEG 数据文件')
    else:
        logger.warning('deg_data_dir 为空，跳过画组间 heatmap 和拆分每个本体 DEG 数据文件')

    
    logger.success("Done")


if __name__ == '__main__':
    main()
