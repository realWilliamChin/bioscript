#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/05 12:26
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd


def parse_input():
    p = argparse.ArgumentParser(description="""对每个单独的组的样本进行筛选，每组中至少有一个大于 min default=50 计为有效
计算每组中 kegg_pathway map 上的基因数量和总共（所有组中至少有一 1 个样本大于 min）
""")
    # p.add_argument(dest='')
    p.add_argument('-k', '--keggclean', help='输入注释出来的 KEGG_clean 文件')
    p.add_argument('-s', '--samplesfile', help='输入 samples_described.txt 文件')
    p.add_argument('-i', '--inputreadstable', dest='input_file',
                   help='输入 gene_count_matrix.txt')
    p.add_argument('-o', '--output', default='KEGG_Summary.txt',
                   help='输出文件名，默认 KEGG_mapped_count_summary_table.txt')
    p.add_argument('-m', '--min', default=50, type=int,
                   help='输入每组中至少有一个样本的 reads count 大于这个数才算有效')
    
    args = p.parse_args()
    
    return args


if __name__ == '__main__':
    args = parse_input()
    
    kegg_ref_df = pd.read_csv(args.keggclean, sep='\t', usecols=[0, 1], header=None, names=['GeneID', 'KEGG_Pathway'])
    
    gene_df = pd.read_csv(args.input_file, sep='\t')
    samples_df = pd.read_csv(args.samplesfile, sep='\t', usecols=[0, 1])
    group_list = samples_df['group'].drop_duplicates().tolist()
    
    result_df = pd.DataFrame()
    result_df['GeneID'] = gene_df['GeneID'].copy()
    
    for group in group_list:
        cur_group_samples = samples_df[samples_df['group'] == group]['sample'].tolist()
        result_df[group] = gene_df[cur_group_samples].apply(max, axis=1)

    result_df['Total'] = result_df.max(axis=1, numeric_only=True)
    for sample in (['Total'] + group_list):
        sample_df = result_df[['GeneID', sample]].copy()
        sample_df = sample_df[sample_df[sample] > args.min]
        kegg_ref_df = pd.merge(kegg_ref_df, sample_df, on='GeneID', how='left')

    kegg_ref_df.fillna(0).to_csv(f'{args.output}.gene.txt', sep='\t', index=False)
    kegg_ref_df = kegg_ref_df.drop(columns='GeneID')

    kegg_ref_df = kegg_ref_df.groupby('KEGG_Pathway').count()
    kegg_ref_df.to_csv(args.output, sep='\t')
