#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/10/09 16:24
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', dest='gene_file', help='输入第一列为 GeneID 的文件，需要有 GeneID 列名')
    p.add_argument('-o', dest='output_file', help='输出文件')
    p.add_argument('--header', dest='header', action='store_true', help='输入文件是否有 GeneID header')
    p.add_argument('-k', '--keggclean', help='KEGG_clean 文件')
    
    args = p.parse_args()
    
    return args


def main():
    args = parse_input()
    if args.header:
        gene_df = load_table(args.gene_file)
    else:
        gene_df = load_table(args.gene_file, header=None, names=['GeneID'])
    kegg_df = load_table(args.keggclean, header=None, usecols=[0, 1, 4], names=['GeneID', 'KO_ID', 'KEGG_ID'])
    kegg_df['KO_ID'] = kegg_df['KO_ID'].str.split(':').str[0]
    
    ko_def_df = pd.merge(gene_df, kegg_df, on='GeneID', how='left')
    regulation_df = ko_def_df.copy()
    regulation_df['regulation'] = '1'
    write_output_df(regulation_df[['KEGG_ID', 'regulation']].drop_duplicates().dropna(), 'regulation.txt', index=False)
    ko_def_df.drop(columns=['GeneID', 'KEGG_ID'], inplace=True)
    ko_def_df.drop_duplicates(subset=['KO_ID'], inplace=True)
    ko_def_df.dropna(inplace=True)
    # ko_def_df = pd.read_csv(ko_file, sep='\t', names=['KEGG_ID'], dtype=str, usecols=[0], header=0)
    passed_path_file = '/home/colddata/qinqiang/script/Analysis/pathview/passed_path.txt'
    passed_path_df = pd.read_csv(passed_path_file, sep='\t', names=['KO_ID', 'Ko_Def'], dtype=str)
    
    ko_def_df = passed_path_df[passed_path_df['KO_ID'].isin(ko_def_df['KO_ID'])]
    ko_def_df = ko_def_df.drop_duplicates(subset=['KO_ID'])
    write_output_df(ko_def_df, args.output_file, index=False, header=False)


if __name__ == '__main__':
    main()