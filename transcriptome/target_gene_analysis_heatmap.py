#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/05/13 17:04
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse
import subprocess
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
# sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from Rscript import draw_multigroup_heatmap
from Rscript import anova_analysis
# from genedf_add_knsdef import add_kns_def


def parse_input():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--targetgene', dest="target_gene_file", type=str, required=True,
                           help='输入文件，target gene 文件，至少包含两列，GeneID 和 Ontology')
    parser.add_argument('--fpkm', type=str, required=True, help='输入 fpkm_matrix.txt 文件')
    parser.add_argument('--kns', type=str, help='anova_p 文件添加定义输出')
    parser.add_argument('--deg-data-dir', dest='deg_data_dir', type=str, help='[必须]输入 DEG_data.txt 文件')
    parser.add_argument('-s', '--samplesinfo', type=str, required=True,
                           help='输入样品信息文件')
    
    args = parser.parse_args()

    return args


def group_vs_group_heatmap(group_target_gene_file, samples_file):
    df = pd.read_csv(group_target_gene_file, sep='\t')
    fpkm_df = df[['GeneID'] + [col for col in df.columns if col.endswith('FPKM')]]
    fpkm_df.columns = ['GeneID'] + [col.replace('_FPKM', '') for col in fpkm_df.columns[1:]]
    ontology_df = df[['GeneID', 'Ontology']]
    group_vs_group_heatmap_fname = group_target_gene_file.replace('_def.txt', '_heatmap.xlsx')
    group_vs_group_heatmap_pname = group_target_gene_file.replace('_def.txt', '_heatmap.jpeg')
    
    group1 = df.iloc[0, 3]
    group2 = df.iloc[0, 4]
    samples_df = pd.read_csv(samples_file, sep='\t', usecols=[0, 1])
    samples_df = samples_df[['sample', 'group']]
    samples_df = samples_df[samples_df['group'].isin([group1, group2])]
    
    with pd.ExcelWriter(group_vs_group_heatmap_fname, engine='openpyxl') as writer:
        fpkm_df.to_excel(writer, sheet_name='Sheet1', index=False)
        samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
        ontology_df.to_excel(writer, sheet_name='Sheet3', index=False)
    
    draw_multigroup_heatmap(group_vs_group_heatmap_fname, group_vs_group_heatmap_pname, other_args='--no-cluster-rows')


def main():
    args = parse_input()
    target_gene_file, samples_file, fpkm_matrix_file = args.target_gene_file, args.samplesinfo, args.fpkm
    
    target_gene_def_df = pd.read_csv(target_gene_file, sep='\t', dtype={'GeneID': str})
    target_gene_def_df['Ontology'] = target_gene_def_df['Ontology'].str.strip()
    
    target_gene_df = target_gene_def_df[['GeneID', 'Ontology']]
    s_df_count = target_gene_df.shape[0]
    target_gene_df = target_gene_df.drop_duplicates(subset=['GeneID'])
    
    if target_gene_df.shape[0] != s_df_count:
        logger.warning(f"输入文件 ID 有重复，已进行去重，数量 {s_df_count - target_gene_df.shape[0]}")
        
    fpkm_matrix_df = pd.read_csv(fpkm_matrix_file, sep='\t', dtype={'GeneID': str})
    gene_fpkm_df = pd.merge(target_gene_df, fpkm_matrix_df, on='GeneID', how='left')
    gene_fpkm_df.drop(columns=['Ontology'], inplace=True)
    anova_file_name = target_gene_file.replace('.txt', '_anova_p.txt')
    gene_fpkm_df.to_csv(anova_file_name, sep='\t', index=False)
    anova_analysis(anova_file_name, samples_file, anova_file_name)
    anova_gene_fpkm_df = pd.read_csv(anova_file_name, sep='\t', dtype={'GeneID': str})
    # if args.kns:
    #     logger.info(f'正在对相关基因添加定义')
    #     anova_gene_fpkm_def_df = add_kns_def(anova_gene_fpkm_df, kns_file=args.kns)
    #     anova_gene_fpkm_def_df.to_csv(anova_file_name, sep='\t', index=False)
    
    # 添加定义改为使用原有定义(2024_06_05:张老师)
    anova_gene_fpkm_def_df = pd.merge(anova_gene_fpkm_df, target_gene_def_df, on='GeneID', how='left')
    anova_gene_fpkm_def_df.to_csv(anova_file_name, sep='\t', index=False)

    # 基因都是挑出来的基因，不需要对 p 值筛选画图 (2024_06_18:张老师)
    # anova_gene_fpkm_df = anova_gene_fpkm_df[anova_gene_fpkm_df['p_value'] <= 0.05]
    anova_gene_fpkm_df.drop(columns=['p_value', 'BH_p_value'], inplace=True)
    anova_gene_fpkm_df = pd.merge(anova_gene_fpkm_df, target_gene_df, on='GeneID', how='left')
    anova_gene_fpkm_df = anova_gene_fpkm_df.sort_values(by=['Ontology'])
    
    samples_df = pd.read_csv(samples_file, sep='\t', usecols=[0, 1])
    samples_df = samples_df[['sample', 'group']]

    # multigroup_heatmap 输入文件
    all_gene_ko_heatmap_filename = target_gene_file.replace('.txt', '_heatmap.xlsx')
    heatmap_filename = all_gene_ko_heatmap_filename.replace('.xlsx', '.jpeg')
    with pd.ExcelWriter(all_gene_ko_heatmap_filename, engine='openpyxl') as writer:
        anova_gene_fpkm_df.loc[:, anova_gene_fpkm_df.columns != 'Ontology'].to_excel(writer, sheet_name="Sheet1", index=False)
        samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
        anova_gene_fpkm_df[['GeneID', 'Ontology']].to_excel(writer, sheet_name='Sheet3', index=False)
    
    draw_multigroup_heatmap(all_gene_ko_heatmap_filename, heatmap_filename, other_args='--no-cluster-rows')
    
    
    deg_data_list = os.listdir(args.deg_data_dir)
    for deg_data_file in deg_data_list:
        
        if not deg_data_file.endswith('_DEG_data.txt'):
            continue

        compare_name = os.path.basename(deg_data_file).replace('_DEG_data.txt', '')
        logger.info(f'正在找相关基因添加定义 {compare_name}')
        deg_data_file = os.path.join(args.deg_data_dir, deg_data_file)
        deg_data_df = pd.read_csv(deg_data_file, sep='\t', dtype={'GeneID': str})

        # result_df = pd.merge(deg_data_df, right=target_gene_def_df, on='GeneID', how='left', suffixes=('_df1', '_df2'))
        result_df = pd.merge(target_gene_def_df, right=deg_data_df, on='GeneID', how='inner', suffixes=('_df1', '_df2'))
        cols_to_drop = [col for col in result_df.columns if col.endswith('_df1')]
        result_df.drop(columns=cols_to_drop, inplace=True)
        result_df.columns = [col.replace('_df2', '') for col in result_df.columns]
        result_df.set_index('GeneID')
        # result_df.dropna(inplace=True)
        result_df.to_csv(f'{compare_name}_target_gene_def.txt', sep='\t', index=False)
        
        logger.info(f'正在画 {compare_name} heatmap')
        group_vs_group_heatmap(f'{compare_name}_target_gene_def.txt', samples_file)
    
    logger.success("Done")


if __name__ == '__main__':
    main()