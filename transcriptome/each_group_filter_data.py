#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/23 16:02
# Author        : William GoGo
import os, sys
import pandas as pd
import argparse
import subprocess
from loguru import logger


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='input_file', help='输入文件，通常是 fpkm_matrix 或者是 reads_matrix 文件')
    parser.add_argument(dest='output_dir', help='输出目录，输出文件按照组名称命名')
    parser.add_argument('-s', '--samples', help='samples_described.txt')
    parser.add_argument('--filter-value', dest='filter_value', type=float, default=50, help='过滤值，默认大于等于 50')
    parser.add_argument('--venn', action='store_true', help='加上这个参数，出 Venn 图')
    
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    args = parse_input()
    output_dir = args.output_dir
    df = pd.read_csv(args.input_file, sep='\t')
    df['GeneID'].astype(str)
    samples_df = pd.read_csv(args.samples, sep='\t', usecols=[0, 1])

    all_group_list = samples_df.drop_duplicates(subset='group')['group'].tolist()
    grouped_samples_df = samples_df.groupby(by='group')

    dataframe_list = []
    for each_group in all_group_list:
        each_group_sample_list = grouped_samples_df.get_group(each_group)['sample'].tolist()
        each_group_sample_df = df[['GeneID'] + each_group_sample_list].copy()
        each_group_sample_df['max'] = each_group_sample_df.max(axis=1, numeric_only=True)
        each_group_sample_df = each_group_sample_df[each_group_sample_df['max'] >= args.filter_value]
        each_group_sample_df = each_group_sample_df.drop(columns=['max'])
        logger.info(f'{each_group}, {each_group_sample_df.shape[0]}')
        each_group_sample_df_filename = f'{each_group}.txt'
        each_group_sample_df.to_csv(os.path.join(output_dir, each_group_sample_df_filename), sep='\t', index=False)
        each_group_sample_df = each_group_sample_df.drop(columns=each_group_sample_list)
        each_group_sample_df = each_group_sample_df.rename(columns={"GeneID": each_group})
        dataframe_list.append(each_group_sample_df)
        

    names_list = [x.columns[0] for x in dataframe_list]
    venn_df = pd.concat(dataframe_list, axis=1, names=names_list)
    print(venn_df)
    venn_df.to_csv('venn.txt', sep='\t', index=False)
    logger.success('Done')

