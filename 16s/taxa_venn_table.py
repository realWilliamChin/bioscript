#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/29 13:58
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    p = argparse.ArgumentParser(description='08/taxa_bar_output 下面的表处理成 venn 输入数据')
    p.add_argument(dest='input_file', help='输入文件')
    p.add_argument(dest='output_file', help='输出文件')
    p.add_argument('-m', '--min-valid-value', default=2, type=int, help='最小有效值')
    
    return p.parse_args()


def taxa_venn_table(input_file, output_file):
    df = pd.read_csv(input_file)
    
    gene_ids = df.iloc[:, 0]
    for i in range(1, df.shape[1]):
        df.iloc[:, i] = df.iloc[:, i].apply(lambda x: gene_ids[df.index[df.iloc[:, i] == x][0]] if x >= 2 else None)


if __name__ == '__main__':
    args = parse_input()
    taxa_venn_table(args.input_file, args.output_file)
    
    