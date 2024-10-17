#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/06 17:23
# Author        : William GoGo
import pandas as pd
import os, sys
import argparse
from loguru import logger
from get_passed_path import passed_path


def parse_input():
    p = argparse.ArgumentParser(description="")
    p.add_argument('--kid', help='输入的 KEGG_ID 文件')
    p.add_argument('-f', '--files', nargs='+', help='每个样本的 reads count 文件', required=True)
    p.add_argument('--keggclean', help='KEGG_clean 文件')
    
    args = p.parse_args()
    
    return args


if __name__ == '__main__':
    args = parse_input()
    kegg_df = pd.read_csv(args.keggclean, sep='\t', header=None, usecols=[0, 1, 4], names=['GeneID', 'KO', 'KEGG_ID'])
    kegg_df['KO'] = kegg_df['KO'].str.split(':').str[0]
    path_file = "/home/colddata/qinqiang/script/Rscript/pathview/passed_path.txt"
    path_df = pd.read_csv(path_file, sep='\t', header=None, names=['KO', 'Def'])
    # all_sample_df = pd.read_csv("gene_count_matrix.txt",sep='\t')
    # all_sample_list = all_sample_df.columns[1:]
    all_sample_list = args.files

    for sample_file in all_sample_list:
        
        # 这块儿筛选出每个组的 reads_count 
        # sample_df = all_sample_df[['GeneID', sample]].copy()
        # sample_df[sample].astype(float)
        # source_rows_number = sample_df.shape[0]
        # sample_df = sample_df[sample_df[sample] > 50]
        # sample_df = sample_df['GeneID']
        # filter_rows_number = sample_df.shape[0]
        # print(f'{sample} 原始行数 {source_rows_number}，过滤完 {args.min} 之后的行数 {filter_rows_number}')
        
        sample_name = sample_file.replace('.txt', '')
        os.mkdir(sample_name)
        sample_df = pd.read_csv(sample_file, sep='\t')

        result_df = pd.merge(sample_df, kegg_df, how='left', on='GeneID')
        result_df['regulation'] =  1
        result_df = result_df.drop_duplicates(subset=['KEGG_ID'])
        result_df = result_df.drop(columns=['GeneID'])
        result_df = result_df.dropna(subset=['KEGG_ID'])
        if args.kid:
            passed_path(args.kid, args.kid.replace('.txt', '_passed_path.txt'))
        else:
            passed_path_df = pd.merge(result_df['KO'], path_df, on='KO', how='left')
            passed_path_df.dropna(subset=['Def']).to_csv(f'{sample_name}/passed_path.txt', sep='\t', index=False, header=False)
        result_df = result_df.drop(columns=['KO'])
        
        result_df[['KEGG_ID', 'regulation']].to_csv(f'{sample_name}/regulation.txt', sep='\t', index=False)