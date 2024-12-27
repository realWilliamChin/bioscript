#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/20 22:51
# Author        : William GoGo
import os
import sys
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    argparser = argparse.ArgumentParser(description='校验 compare_info.txt 和 samples_described.txt 是否一致')
    argparser.add_argument('-c', '--compare', help='compare_info.txt 文件', default='compare_info.txt')
    argparser.add_argument('-s', '--sample', help='samples_described.txt 文件', default='samples_described.txt')
    argparser.add_argument('-d', '--cleandata', help='cleandata 文件夹')

    args = argparser.parse_args()
    return args


def check_sample(samples_file, cleandata_dir):
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
            logger.error(f'文件 {file1} 或 {file2} 不存在，将跳过此样本 {samples_name} 的比对')
            continue
    logger.info('cleandata 文件检查完毕')


def check_sample_comp(samples_file, comp_file):
    samples = pd.read_csv(samples_file, sep='\t')
    comp = pd.read_csv(comp_file, sep='\t')
    samples_group = samples['group'].tolist()
    comp_treat = comp['Treat'].tolist()
    comp_control = comp['Control'].tolist()
    treat_err = [(i+1, j) for i, j in enumerate(comp_treat) if j not in samples_group]
    control_err = [(i+1, j) for i, j in enumerate(comp_control) if j not in samples_group]
    if treat_err or control_err:
        logger.error(f'treat err in group lst: {str(treat_err)}')
        logger.error(f'control err in group lst: {control_err}')
        return 1
    elif treat_err == [] and control_err == []:
        logger.success('samples_described.txt and Compare_info.txt is Match!')
        return 0


def main():
    args = parse_input()
    if args.compare:
        rep = check_sample_comp(args.sample, args.compare)
        if rep == 1:
            sys.exit(1)
    if args.cleandata:
        check_sample(args.sample, args.cleandata)


if __name__ == '__main__':
    main()
