#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/29 17:19
# Author        : William GoGo
import os
import re
import sys
import pandas as pd
import argparse
import subprocess
import openpyxl
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/Rscript/')
from Rscript import smart_heatmap
sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df
from data_check import df_drop_row_sum_eq_zero, df_drop_element_side_space, df_replace_illegal_folder_chars


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', default='target_go_df.txt', help='GO_ID 的一个 list 文件, 列名包含 GO_ID')
    p.add_argument('-g', '--genego', help='swiss 注释解析出来的 gene_go 文件')
    p.add_argument('-f', '--fpkm', default='fpkm_matrix_filtered.txt', help='fpkm_matrix')
    p.add_argument('-s', '--samples', default='samples_described.txt', help='samples_described.txt 样本描述文件')
    p.add_argument('--genesymbol', help='使用 GeneSymbol 作为索引画热图，人和大鼠小鼠通常使用，输入包含 GeneID GeneSymbol 两列的文件')
    p.add_argument('-o', '--outputdir', default='00_Each_GO_pathway_heatmap', help='所有 heatmap 输出文件夹')
    p.add_argument('-e', '--expression-data', help='[暂时不用] reads_fpkm_matrix_def.txt')
    args = p.parse_args()
    return args


def each_go_gene_expression(target_go_df, gene_go_df, expression_data, output_dir='00_Each_GO_pathway_expression_data'):
    for _, row in target_go_df.iterrows():
        go_id = row['GO_ID']
        go_def = row['GO_def']
        each_go_id_df = gene_go_df[gene_go_df['GO_ID'] == go_id]
        if each_go_id_df.shape[0] < 1:
            logger.warning(f'没有 {go_id} 相关基因')
            continue
        each_go_id_gene_expression_df = pd.merge(each_go_id_df, expression_data, on='GeneID', how='inner')
        each_go_id_gene_expression_df.drop(columns=['GO_ID'], inplace=True)
        # 将任何不适合创建文件的字符（包括空格）变成 _
        go_replace_name = go_id.replace(':', '_') + '_' + re.sub(r'[\\/:*?"<>|\s]', '_', str(go_def))
        go_id_gene_expression_fn = os.path.join(output_dir, f'{go_replace_name}_gene_expression.xlsx')
        write_output_df(each_go_id_gene_expression_df, go_id_gene_expression_fn, index=False)


def each_go_gene_heatmap(target_go_df, gene_go_df, fpkm_matrix_df, samples_df, output_dir='00_Each_GO_pathway_heatmap', genesymbol_file=None):
    samples_df = samples_df[['sample', 'group']]
    samples_list = ['GeneID'] + samples_df['sample'].values.tolist()
    
    # 如果提供了 genesymbol 参数，读取 GeneID 和 GeneSymbol 映射文件
    genesymbol_df = None
    if genesymbol_file:
        genesymbol_df = load_table(genesymbol_file)
        # 只保留 GeneID 和 GeneSymbol 两列
        genesymbol_df = genesymbol_df[['GeneID', 'GeneSymbol']].copy()
        # 删除空值
        genesymbol_df.dropna(subset=['GeneID', 'GeneSymbol'], inplace=True)
        logger.info(f'读取到 {genesymbol_df.shape[0]} 个 GeneID-GeneSymbol 映射')
    
    for _, row in target_go_df.iterrows():
        go_pathway_id = row['GO_pathway_ID']
        go_pathway_def = row['GO_pathway_def']
        each_go_pathway_id_df = gene_go_df[gene_go_df['GO_pathway_ID'] == go_pathway_id]
        if each_go_pathway_id_df.shape[0] == 0:
            logger.warning(f'没有 {go_pathway_id} 相关基因')
            continue
        each_go_pathway_id_gene_fpkm_df = pd.merge(each_go_pathway_id_df, fpkm_matrix_df, on='GeneID', how='inner')
        each_go_pathway_id_gene_fpkm_df.dropna(how='any', axis=0, inplace=True)
        each_go_pathway_id_gene_fpkm_df = df_drop_row_sum_eq_zero(each_go_pathway_id_gene_fpkm_df)
        if each_go_pathway_id_gene_fpkm_df.shape[0] <= 2:
            logger.warning(f'{go_pathway_id} 相关基因表达量小于 2，跳过画 Heatmap 图')
            continue
        
        # 如果提供了 genesymbol_df，则替换 GeneID 为 GeneSymbol
        if genesymbol_df is not None:
            # 合并 GeneSymbol 映射
            each_go_pathway_id_gene_fpkm_df = pd.merge(each_go_pathway_id_gene_fpkm_df, genesymbol_df, on='GeneID', how='inner')
            if each_go_pathway_id_gene_fpkm_df.shape[0] < 2:
                logger.warning(f'{go_pathway_id} 相关 genesymbol 数量小于 2，跳过画 Heatmap 图')
                continue
            # 将 GeneID 列替换为 GeneSymbol
            each_go_pathway_id_gene_fpkm_df['GeneID'] = each_go_pathway_id_gene_fpkm_df['GeneSymbol']
            each_go_pathway_id_gene_fpkm_df.drop(columns=['GeneSymbol'], inplace=True)
        
        go_pathway_replace_name = go_pathway_id.replace(':', '_') + '_' + re.sub(r'[\\/:*?"<>|\s]', '_', str(go_pathway_def))
        # 准备 multigroup heatmap 输入文件
        go_pathway_id_gene_heatmap_filename = os.path.join(output_dir, f'{go_pathway_replace_name}_gene_heatmap.xlsx')
        with pd.ExcelWriter(go_pathway_id_gene_heatmap_filename, engine='openpyxl') as writer:
            each_go_pathway_id_gene_fpkm_df[samples_list].to_excel(writer, sheet_name='Sheet1', index=False)
            samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
            
            # # 如果 target_go_df 存在 Ontology 或 SubOntology 列，添加到第三个 Sheet
            # if 'Ontology' in target_go_df.columns or 'SubOntology' in target_go_df.columns:
            #     # 创建行注释 DataFrame
            #     row_annotation_data = {'GeneID': each_go_pathway_id_gene_fpkm_df['GeneID'].values}
                
            #     # 添加 Ontology 列
            #     if 'Ontology' in target_go_df.columns:
            #         ontology_value = row.get('Ontology', '')
            #         row_annotation_data['Ontology'] = [ontology_value] * len(each_go_pathway_id_gene_fpkm_df)
                
            #     # 添加 SubOntology 列
            #     if 'SubOntology' in target_go_df.columns:
            #         subontology_value = row.get('SubOntology', '')
            #         row_annotation_data['SubOntology'] = [subontology_value] * len(each_go_pathway_id_gene_fpkm_df)
                
            #     row_annotation_df = pd.DataFrame(row_annotation_data)
            #     row_annotation_df.to_excel(writer, sheet_name='Sheet3', index=False)
        
        # heatmap 画图
        go_pathway_id_heatmap_filename = os.path.join(output_dir, os.path.basename(go_pathway_id_gene_heatmap_filename).replace('.xlsx', '.jpeg'))
        heatmap_result = smart_heatmap(
            go_pathway_id_gene_heatmap_filename,
            go_pathway_id_heatmap_filename,
            annotation_col=2,
            cluster_rows=True,
            cluster_cols=False,
            scale="row"
        )
        if not heatmap_result:
            logger.error(f'{go_pathway_id} draw_multigroup_heatmap 程序失败')
    

def main():
    args = parse_input()
    genego_df = load_table(args.genego, header=None, names=['GeneID', 'GO_pathway_ID', 'GO_pathway_def'])
    # target_go_df 预处理
    gene_go_df = load_table(args.genego, header=None, names=['GeneID', 'GO_pathway_ID', 'GO_pathway_def'], dtype={"GeneID": str})
    go_id_def_df = gene_go_df[['GO_pathway_ID', 'GO_pathway_def']].drop_duplicates()
    target_go_df = load_table(args.input)
    target_go_df = df_drop_element_side_space(target_go_df)
    target_go_df['GO_pathway_ID'] = target_go_df['GO_pathway_ID'].str.split('_').str[0]
    
    # 合并 target_go_df 和 gene_go_df
    num_before_merge = target_go_df.shape[0]
    merged_target_go_df = pd.merge(left=go_id_def_df, right=target_go_df, on='GO_pathway_ID', how='inner')
    num_after_merge = merged_target_go_df.shape[0]
    logger.info(f"GO分析: 合并前 target_go_df 行数为: {num_before_merge}, 合并后为: {num_after_merge}")
    if num_after_merge < num_before_merge:
        lost_go_df = target_go_df[~target_go_df['GO_pathway_ID'].isin(merged_target_go_df['GO_pathway_ID'])]
        lost_go_log_dir = os.path.join(args.output, 'log')
        os.makedirs(lost_go_log_dir, exist_ok=True)
        lost_go_file = os.path.join(lost_go_log_dir, 'target_go_lost_in_merge_gene_go.txt')
        write_output_df(lost_go_df, lost_go_file, index=False)
        logger.warning(f"有 {num_before_merge-num_after_merge} 个 GO_pathway_ID 在 gene_go_df 中未找到, 详细列表见 {lost_go_file}")
    target_go_df = merged_target_go_df
    target_go_df = df_replace_illegal_folder_chars(target_go_df, ['GO_pathway_def', 'Ontology', 'SubOntology'])
    
    if args.fpkm:
        samples_df = load_table(args.samples)
        fpkm_df = load_table(args.fpkm)
        each_go_gene_heatmap(target_go_df, genego_df, fpkm_df, samples_df, args.outputdir, args.genesymbol)
    
    if args.expression_data:
        expression_df = load_table(args.expression_data)
        each_go_gene_expression(target_go_df, genego_df, expression_df, args.outputdir)
    
    logger.success('Done')


if __name__ == '__main__':
    main()