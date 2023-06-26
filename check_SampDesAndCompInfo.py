#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/20 22:51
# Author        : William GoGo
import os
import argparse
import pandas as pd


def parse_input():
    argparser = argparse.ArgumentParser(description='校验 compare_info.txt 和 samples_described.txt 是否一致')
    argparser.add_argument('-c', '--compare', help='compare_info.txt 文件', default='compare_info.txt')
    argparser.add_argument('-s', '--sample', help='samples_described.txt 文件', default='samples_described.txt')

    args = argparser.parse_args()
    return args


def main():
    args = parse_input()
    group_lst = pd.read_csv(args.sample, sep='\t')['group'].tolist()
    treat_lst = pd.read_csv(args.compare, sep='\t')['Treat'].tolist()
    control_lst = pd.read_csv(args.compare, sep='\t')['Control'].tolist()
    
    treat_err_in_group_lst = [(i, j) for i, j in enumerate(treat_lst) if j not in group_lst]
    control_err_in_group_lst = [(i, j) for i, j in enumerate(control_lst) if j not in group_lst]
    
    print(f'treat err in group lst: {treat_err_in_group_lst}')
    print(f'control err in group lst: {control_err_in_group_lst}')


if __name__ == '__main__':
    main()
