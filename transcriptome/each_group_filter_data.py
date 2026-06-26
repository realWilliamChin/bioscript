#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/23 16:02
# Author        : William GoGo
import os, sys
import pandas as pd
import argparse
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', help='输入文件，通常是 fpkm_matrix 或者是 reads_matrix 文件')
    parser.add_argument('-o', '--output-dir', help='输出目录，输出文件按照组名称命名')
    parser.add_argument('-s', '--samples', help='samples_described.txt')
    parser.add_argument('-f', '--filter-value', type=float, default=50, help='过滤值，默认大于等于 50')
    
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    args = parse_input()
    df = load_table(args.input_file)
    df[df.columns[0]].astype(str)
    samples_df = load_table(args.samples, usecols=[0, 1])
    
    all_group_list = samples_df['group'].drop_duplicates().tolist()
    grouped_samples_df = samples_df.groupby(by='group')
    for each_group in all_group_list:
        each_group_sample_list = grouped_samples_df.get_group(each_group)['sample'].values.tolist()
        each_group_sample_df = df[[df.columns[0]] + each_group_sample_list].copy()
        each_group_sample_df['max'] = each_group_sample_df.max(axis=1, numeric_only=True)
        each_group_sample_df = each_group_sample_df[each_group_sample_df['max'] >= args.filter_value].drop(columns=['max'])
        logger.info(f'{each_group}, {each_group_sample_df.shape[0]}')
        each_group_sample_df_filename = os.path.join(args.output_dir, f'{each_group}.txt')
        write_output_df(each_group_sample_df, each_group_sample_df_filename, index=False)
