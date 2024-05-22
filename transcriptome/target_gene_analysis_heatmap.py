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
from Rscript import draw_multigroup_heatmap
from Rscript import anova_analysis


def parse_input():
    argparser = argparse.ArgumentParser(description='')
    argparser.add_argument('-i', '--targetgene', dest="target_gene_file", type=str, required=True,
                           help='输入文件，target gene 文件，至少包含两列，GeneID 和 Ontology')
    argparser.add_argument('--fpkm', type=str, required=True, help='输入 fpkm_matrix.txt 文件')
    argparser.add_argument('-s', '--samplesinfo', type=str, required=True,
                           help='输入样品信息文件')
    
    args = argparser.parse_args()

    return args


def target_gene_draw_multigroup_heatmap(target_gene_file, samples_file, fpkm_matrix_file):
    target_gene_df = pd.read_csv(target_gene_file, sep='\t', dtype={'GeneID': str})
    target_gene_df = target_gene_df[['GeneID', 'Ontology']]
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
    anova_gene_fpkm_df = anova_gene_fpkm_df[anova_gene_fpkm_df['p_value'] <= 0.05]
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


def main():
    args = parse_input()
    target_gene_draw_multigroup_heatmap(args.target_gene_file, args.samplesinfo, args.fpkm)
    
    logger.success("Done")


if __name__ == '__main__':
    main()