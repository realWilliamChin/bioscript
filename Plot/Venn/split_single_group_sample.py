#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/08/04 11:58
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df
from data_check import df_drop_row_sum_eq_zero


def parse_input():
    p = argparse.ArgumentParser(description='按照每组拆分成单个文件，并且筛选大于 0 的')
    p.add_argument('-i', help='输入多样本表文件')
    p.add_argument('-s', help='输入 samples_described.txt 文件')
    p.add_argument('-o', help='输出文件目录', default=os.getcwd())

    args = p.parse_args()
    
    return args


def main():
    args = parse_input()
    df = load_table(args.i)
    samples_df = load_table(args.s, usecols=[0, 1])

    group_samples = samples_df.groupby('group')['sample'].apply(list)
    for group_name, samples in group_samples.items():
        logger.info(f'Processing group: {group_name}, samples: {samples}')

        group_df = df[[df.columns[0]] + samples].copy()
        group_df = df_drop_row_sum_eq_zero(group_df)
        
        write_output_df(group_df[df.columns[0]], os.path.join(args.o, f'{group_name}.txt'), index=False)      


if __name__ == '__main__':
    main()