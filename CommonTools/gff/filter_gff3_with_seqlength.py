#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/08/19 11:57
# Author        : William GoGo
import os, sys
import re
import argparse
from loguru import logger
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='input', help='gff3 输入文件')
    parser.add_argument(dest='output', help='输出文件')
    parser.add_argument('-m', type=int, help='最低长度')
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_input()
    input_file = args.input
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute']
    df = pd.read_csv(input_file, sep='\t', header=None, names=gff3_columns)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    
    df['ID'] = df['attribute'].str.split(';').str[0].str.replace('ID=', '').str.split('.').str[0]
    long_genes_list = df[(df['type'] == 'gene') & ((df['end'] - df['start']) > args.m)]['ID'].to_list()
    print(len(long_genes_list))
    print(df.head(10))
    df = df[df['ID'].isin(long_genes_list)]
    df = df.drop(columns=['ID'])
    df.to_csv(args.output, sep='\t', index=False)
    