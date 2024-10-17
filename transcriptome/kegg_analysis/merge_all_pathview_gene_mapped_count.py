#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/25 15:34
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-f', '--files', nargs='+', help='每个样本的 mapped_count 文件', required=True)
    p.add_argument('-n', '--totalnames', nargs='+', help='每个样本的 Total 列的名字', required=True)
    p.add_argument('-o', '--output', help='输出文件')
    
    return p.parse_args()


def merge_mappped_file(file_list, total_name_list, output_file):
    df = pd.read_csv(file_list[0], sep='\t')
    df = df.rename(columns={'Total': total_name_list[0]})
    for f, total_name in zip(file_list[1:], total_name_list[1:]):
        f_df = pd.read_csv(f, sep='\t')
        f_df = f_df.rename(columns={'Total': total_name})
        df = pd.merge(df, f_df, on='KEGG_Pathway', how='outer')
    df.to_csv(output_file, sep='\t', index=False)


def main():
    args = parse_input()
    merge_mappped_file(args.files, args.totalnames, args.output)


if __name__ == '__main__':
    main()