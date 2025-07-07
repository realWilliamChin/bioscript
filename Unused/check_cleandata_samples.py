#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/12/12 16:43
# Author        : William GoGo
import os, sys
import pandas as pd
import argparse
from loguru import logger


def parse_input():
    args = argparse.ArgumentParser()
    args.add_argument('-s', '--samples', help='samples_described.txt', default='samples_described.txt')
    args.add_argument('--cd', help='cleandata file dir')
    
    return args.parse_args()


def check_data(samples_file, cleandata_dir):
    """检查 cleandata 文件是否存在

    Args:
        samples_file (str): samples_described.txt
        cleandata_dir (str): cleandata 文件路径
    """
    samples_df = pd.read_csv(samples_file, sep='\t', usecols=['sample', 'R1', 'R2'])
    for index, row in samples_df.iterrows():
        samples_name = row['sample']
        file1 = os.path.join(cleandata_dir, row['R1'])
        file2 = os.path.join(cleandata_dir, row['R2'])
        if not os.access(file1, os.R_OK) or not os.access(file2, os.R_OK):
            logger.error(f'文件 {file1} 或 {file2} 不存在')
            continue
    logger.info('cleandata 文件检查完毕')


def main():
    args = parse_input()
    check_data(args.samples, args.cd)


if __name__ == '__main__':
    main()