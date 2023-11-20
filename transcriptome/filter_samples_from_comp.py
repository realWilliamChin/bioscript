#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/10/31 14:56
# Author        : William GoGo
"""
从定义好的 compare_info.txt 中，去掉 samples_described.txt 中多余的组
用来多组 compare_info.txt 的 multideseq 的自动化
"""
import os
import pandas as pd
import argparse


def parse_input():
    args = argparse.ArgumentParser()
    args.add_argument('-s', '--samples', dest="samples",
                      help='samples_described.txt, 默认是当前目录的 samples_described.txt，执行程序会进行替换')
    args.add_argument('-c', '--compare', help='compare_info.txt，默认是当前目录的 compare*.txt')
    args.add_argument('-o', '--output', help='输出文件')
    
    args = args.parse_args()
    
    if not args.samples:
        args.samples = [x for x in os.listdir() if x.startswith('samples_described')][0]
    if not args.compare:
        args.compare = [x for x in os.listdir() if x.startswith('compare_info')][0]
    if not args.output:
        args.output = args.samples
    
    return args


def filter_samples_from_comp(samples, compares, output_file):
    samples_df = pd.read_csv(samples, sep='\t', usecols=[0, 1])
    compares_df_1 = pd.read_csv(compares, sep='\t', usecols=[0], skiprows=1, names=['group'])
    compares_df_2 = pd.read_csv(compares, sep='\t', usecols=[1], skiprows=1, names=['group'])
    compares_df = pd.concat([compares_df_1, compares_df_2])
    compares_df = compares_df.drop_duplicates(subset='group')
    
    samples_df = samples_df[samples_df['group'].isin(compares_df['group'])]
    
    if set(samples_df['group']) != set(compares_df['group']):
        samples_not_in_compare = set(compares_df['group']) - set(samples_df['group'])
        compare_not_in_samples = set(samples_df['group']) - set(compares_df['group'])
        print(f'samples 不在 compare {samples_not_in_compare}')
        print(f'compare 不在 samples {compare_not_in_samples}')
        exit(1)
    samples_df.to_csv(output_file, sep='\t', index=False)


def main():
    args = parse_input()
    filter_samples_from_comp(args.samples, args.compare, args.output)


if __name__ == '__main__':
    main()
