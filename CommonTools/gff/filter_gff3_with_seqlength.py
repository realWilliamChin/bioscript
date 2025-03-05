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
    parser.add_argument('-i', dest='input_file', help='gff3 输入文件')
    parser.add_argument('-o', dest='output_file', help='输出文件')
    parser.add_argument('--iheader', action='store_true', help='输入文件是否包含 header')
    
    subparsers = parser.add_subparsers(dest='command', help='子命令帮助')
    
    # 按照序列长度过滤
    seq_length_parser = subparsers.add_parser('seqlength', help='按照 sequnce length 过滤 gff3')
    seq_length_parser.add_argument('-m', type=int, help='最低长度')
    
    # 按照得分过滤
    score_parser = subparsers.add_parser('score', help='按照 score 得分过滤 gff3')
    score_parser.add_argument('--min', type=float, default=0, help='最低得分')
    score_parser.add_argument('--max', type=float, default=0, help='最高得分')

    args = parser.parse_args()
    
    # # 如果没有子命令，打印帮助信息
    if args.command is None:
        parser.print_help()
        return None
    
    # # 调用子命令对应的处理函数
    # args.func(args)

    return args


def seq_length_filter(gff3_df):

    gff3_df['start'] = gff3_df['start'].astype(int)
    gff3_df['end'] = gff3_df['end'].astype(int)
    
    gff3_df['ID'] = gff3_df['attribute'].str.split(';').str[0].str.replace('ID=', '').str.split('.').str[0]
    long_genes_list = gff3_df[(gff3_df['type'] == 'gene') & ((gff3_df['end'] - gff3_df['start']) > args.m)]['ID'].to_list()
    print(len(long_genes_list))
    print(gff3_df.head(10))
    gff3_df = gff3_df[gff3_df['ID'].isin(long_genes_list)]
    gff3_df = gff3_df.drop(columns=['ID'])
    
    return gff3_df


def score_filter(gff3_df, socre_min, score_max):
    source_row_number = gff3_df.shape[0]
    
    if socre_min != 0:
        gff3_df = gff3_df[gff3_df['score'] >= socre_min]
    
    if score_max != 0:
        gff3_df = gff3_df[gff3_df['score'] <= score_max]
        
    filtered_row_number = gff3_df.shape[0]
    
    logger.info(f'根据 socre 过滤前 {source_row_number}, 过滤后 {filtered_row_number}')
    
    return gff3_df
    

if __name__ == '__main__':
    args = parse_input()
    gff3_file = args.input_file
    output_file = args.output_file
    iheader = args.iheader
    
    
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute']
    
    if iheader:
        gff3_df = pd.read_csv(gff3_file, sep='\t')
    else:
        gff3_df = pd.read_csv(gff3_file, sep='\t', skiprows=1, names=gff3_columns)
    
    if args.command == 'seqlength':
        logger.info('通过序列长度处理')
        gff3_df['start'] = gff3_df['start'].astype(int)
        gff3_df['end'] = gff3_df['end'].astype(int)
        gff3_result_df = seq_length_filter(gff3_df)
    elif args.command == 'score':
        logger.info('通过 socre 列进行处理')
        gff3_df['socre'] = gff3_df['score'].astype(float)
        gff3_result_df = score_filter(gff3_df, args.min, args.max)
    
    
    gff3_result_df.to_csv(output_file, sep='\t', index=False)
    