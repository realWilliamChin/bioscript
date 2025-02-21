#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/05/20 21:41
# Author        : William GoGo
import os
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    argparser = argparse.ArgumentParser(description='生成 Trinity 拼接的样本文件')
    argparser.add_argument('-s', '--samples_described', default='samples_described.txt', required=True,
                           help='samples described file')
    argparser.add_argument('-d', '--datadir', help='input dir', default='00_Cleandata')
    argparser.add_argument('-o', '--output', help='output file', default='samples_trinity.txt')
    argparser.add_argument('-t', '--type', choices=['all', 'max', 'custom'], default='max',
                           help="all: all samples；max: 每组中选文件最长的(指定-s)；custom: custom samples")
    args = argparser.parse_args()
    if args.type == 'max':
        if not args.samples_described:
            logger.warning('Please input samples described file.')
    return args


def main():
    args = parse_input()
    # change args.input_dir to real path
    data_dir = os.path.abspath(args.datadir)
    samples_file = args.samples_described
    # samples_file 读取为字典
    samples_df = pd.read_csv(samples_file, sep='\t', skiprows=1, names=['group', 'sample', 'R1', 'R2'])
    samples_df['R1'] = samples_df['R1'].apply(lambda x: os.path.join(data_dir, x))
    samples_df['R2'] = samples_df['R2'].apply(lambda x: os.path.join(data_dir, x))
    if args.type == 'all':
        samples_df.to_csv(args.output, sep='\t', index=False, header=False)
    elif args.type == 'max':
        samples_df['length'] = samples_df.apply(lambda x: max(os.path.getsize(x['f1']), os.path.getsize(x['f2'])), axis=1)
        samples_df = samples_df.sort_values(by=['group', 'length'], ascending=False)
        samples_df = samples_df.drop_duplicates(subset=['group'], keep='first')
        samples_df.drop(columns=['length'], inplace=True)
        samples_df.to_csv(args.output, sep='\t', index=False, header=False)
    logger.info(samples_df)
    logger.success('Done.')
    

if __name__ == '__main__':
    main()
