#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/10/09 16:24
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument(dest='gene_file', help='输入第一列为 GeneID 的文件')
    p.add_argument(dest='output_file', help='输出文件')
    p.add_argument('-k', '--keggclean', help='KEGG_clean 文件')
    
    args = p.parse_args()
    
    return args


def passed_path(ko_file, output_file):
    ko_def_df = pd.read_csv(ko_file, sep='\t', names=['KEGG_ID'], dtype=str, usecols=[0], header=0)
    passed_path_file = '/home/colddata/qinqiang/script/Rscript/pathview/passed_path.txt'
    passed_path_df = pd.read_csv(passed_path_file, sep='\t', names=['KEGG_ID', 'Ko_Def'], dtype=str)
    ko_def_df = passed_path_df[passed_path_df['KEGG_ID'].isin(ko_def_df['KEGG_ID'])]
    ko_def_df = ko_def_df.drop_duplicates(subset=['KEGG_ID'])
    ko_def_df.to_csv(output_file, sep='\t', index=False, header=False)


if __name__ == '__main__':
    args = parse_input()
    
    passed_path(args.ko_file, args.output_file)