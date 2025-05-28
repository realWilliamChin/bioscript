#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/10/09 16:24
# Author        : William GoGo
import os, sys
import pandas as pd
import argparse
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input-dir', dest='input_dir', help='Target_KEGG_analysis 目标分析目录')
    p.add_argument('-o', '--output', default='ko03000_data_summary.csv',
                   help='输出汇总文件，默认为 ko03000_data_summary.csv')
    
    args = p.parse_args()
    
    return args


def main():
    args = parse_input()
    ko03000_df_list = []
    for each_dir in os.listdir(args.input_dir):
        if '-vs-' not in each_dir:
            continue
        ko03000_data_df = load_table(
            os.path.join(args.input_dir, each_dir, 'DEG_expression_data', f'{each_dir}_ko03000_DEG_data.csv'),
            header = 0,
            usecols = ['GeneID', 'sampleA', 'sampleB', 'pvalue', 'padj', 'regulation', 'FC']
        )
        ko03000_data_df = ko03000_data_df[ko03000_data_df['regulation'].str.lower() != 'nosignificant']
        ko03000_df_list.append(ko03000_data_df)
    ko03000_df = pd.concat(ko03000_df_list)
    write_output_df(ko03000_df, args.output, index=False)
    
    logger.success(f'Success {args.output} output!')


if __name__ == '__main__':
    main()