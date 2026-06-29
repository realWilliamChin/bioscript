#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/03/31 11:49
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from r_wrapper import draw_multigroup_heatmap
from r_wrapper import anova_analysis
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--targetgene', dest="target_gene_file", type=str, required=True,
                           help='输入文件，target gene 文件，至少包含两列，GeneID 和 Ontology')
    parser.add_argument('--fpkm', type=str, required=True, help='输入 fpkm_matrix.txt 文件')
    parser.add_argument('-d', '--data-dir', dest='data_dir', type=str, required=True,
                        help='输入，geneid 文件目录')
    parser.add_argument('--mean', default=False, type=bool, help='使用每组中的平均数画 heatmap')
    parser.add_argument('-s', '--samplesinfo', type=str, required=True,
                           help='输入样品信息文件')
    parser.add_argument('-o', '--output', type=str, default=os.getcwd(), help='输出目录')
    
    args = parser.parse_args()

    return args


def main():
    args = parse_input()
    output_dir = args.output
    target_gene_file, samples_file, fpkm_matrix_file = args.target_gene_file, args.samplesinfo, args.fpkm
    samples_df = load_table(samples_file, usecols=[0, 1])
    target_gene_def_df = load_table(target_gene_file, dtype={'GeneID': str})
    target_gene_def_df['Ontology'] = target_gene_def_df['Ontology'].str.strip()
    target_gene_def_df['SubOntology'] = target_gene_def_df['SubOntology'].str.strip()
    
    target_gene_df = target_gene_def_df[['GeneID', 'SubOntology', 'Ontology']]
    s_df_count = target_gene_df.shape[0]
    target_gene_df = target_gene_df.drop_duplicates(subset=['GeneID'])
    
    if target_gene_df.shape[0] != s_df_count:
        logger.warning(f"输入文件 ID 有重复，已进行去重，数量 {s_df_count - target_gene_df.shape[0]}")
        
    fpkm_matrix_df = load_table(fpkm_matrix_file, dtype={'GeneID': str})
    gene_fpkm_df = pd.merge(target_gene_df, fpkm_matrix_df, on='GeneID', how='left')
    gene_fpkm_df.drop(columns=['Ontology', 'SubOntology'], inplace=True)
    anova_file_name = os.path.join(output_dir, target_gene_file.split(os.sep)[-1].replace('.txt', '_anova_p.txt'))
    gene_fpkm_df.to_csv(anova_file_name, sep='\t', index=False)
    anova_analysis(anova_file_name, samples_file, anova_file_name)
    anova_gene_fpkm_df = pd.read_csv(anova_file_name, sep='\t', dtype={'GeneID': str})
    
    anova_gene_fpkm_def_df = pd.merge(anova_gene_fpkm_df, target_gene_def_df, on='GeneID', how='left')
    anova_gene_fpkm_def_df.to_csv(anova_file_name, sep='\t', index=False)

    anova_gene_fpkm_df.drop(columns=['p_value', 'BH_p_value'], inplace=True)
    anova_gene_fpkm_df = pd.merge(anova_gene_fpkm_df, target_gene_df, on='GeneID', how='left')
    anova_gene_fpkm_df = anova_gene_fpkm_df.sort_values(by=['Ontology', 'SubOntology'])
    
    # samples_df = pd.read_csv(samples_file, sep='\t', usecols=[0, 1])
    samples_df = samples_df[['sample', 'group']]

    # multigroup_heatmap 输入文件
    all_gene_ko_heatmap_filename = os.path.join(output_dir, target_gene_file.split(os.sep)[-1].replace('.txt', '_heatmap.xlsx'))
    heatmap_filename = all_gene_ko_heatmap_filename.replace('.xlsx', '.jpeg')
    multigroup_heatmap_data_df = anova_gene_fpkm_df.loc[:, anova_gene_fpkm_df.columns != 'Ontology']
    multigroup_heatmap_data_df = multigroup_heatmap_data_df.loc[:, multigroup_heatmap_data_df.columns != 'SubOntology']
    multigroup_heatmap_sheet3_df = anova_gene_fpkm_df[['GeneID', 'SubOntology', 'Ontology']].copy()
    
    # 根据 samplesinfo 每组中的平均数画 heatmap
    if args.mean:
        for each_group in samples_df['group'].unique():
            group_samples = samples_df[samples_df['group'] == each_group]['sample']
            multigroup_heatmap_data_df[each_group] = multigroup_heatmap_data_df[group_samples].mean(axis=1)
        multigroup_heatmap_data_df.drop(columns=samples_df['sample'], inplace=True)
        samples_df = samples_df.drop(columns=['sample']).drop_duplicates(subset=['group'])
        samples_df['sample'] = samples_df['group']
        samples_df = samples_df[['sample', 'group']]

    
    with pd.ExcelWriter(all_gene_ko_heatmap_filename, engine='openpyxl') as writer:
        multigroup_heatmap_data_df.to_excel(writer, sheet_name="Sheet1", index=False)
        samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
        multigroup_heatmap_sheet3_df.to_excel(writer, sheet_name='Sheet3', index=False)
    
    draw_multigroup_heatmap(all_gene_ko_heatmap_filename, heatmap_filename, other_args='--no-cluster-rows')
    

    data_list = os.listdir(args.data_dir)
    for each_data in data_list:
        item_name = each_data.replace('_gene.txt', '')
        logger.info(f'正在找相关基因添加定义 {item_name}')
        data_df = load_table(
            os.path.join(args.data_dir, each_data),
            header=None,
            usecols=[0,],
            names=['GeneID',],
            dtype={'GeneID': str}
        )

        data_def_df = pd.merge(data_df, target_gene_def_df, on='GeneID', how='inner')
        data_def_df = data_def_df.sort_values(by=['Ontology', 'SubOntology'])
        data_fpkm_df = pd.merge(data_def_df['GeneID'], fpkm_matrix_df, on='GeneID', how='inner')
        
        if data_fpkm_df.shape[0] == 0:
            logger.warning(f'目标基因在模块 {item_name} 中没有出现，跳过')
            continue
        
        # result_df.set_index('GeneID')
        write_output_df(
            data_def_df,
            os.path.join(args.output, f'{item_name}_target_gene_def.csv'),
            index=False
        )
        
        logger.info(f'正在画 {item_name} heatmap')
        
        with pd.ExcelWriter(os.path.join(args.output, f'{item_name}_target_gene_heatmap.xlsx'), engine='openpyxl') as writer:
            data_fpkm_df.to_excel(writer, sheet_name='Sheet1', index=False)
            samples_df[['sample', 'group']].to_excel(writer, sheet_name='Sheet2', index=False)
            data_def_df[['GeneID', 'SubOntology', 'Ontology']].to_excel(writer, sheet_name='Sheet3', index=False)
        
        if data_fpkm_df.shape[0] < 2:
            logger.warning(f'目标基因在模块 {item_name} 中出现 1 个，跳过画 heatmap 图')
            continue
        
        draw_multigroup_heatmap(
            os.path.join(args.output, f'{item_name}_target_gene_heatmap.xlsx'),
            os.path.join(args.output, f'{item_name}_target_gene_heatmap.jpeg'),
            other_args='--no-cluster-rows'
        )
        
    
    logger.success("Done")


if __name__ == '__main__':
    main()
