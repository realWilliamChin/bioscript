#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/29 13:58
# Author        : William GoGo
import os, sys
import pandas as pd
import argparse
from loguru import logger


def parse_input():
    p = argparse.ArgumentParser(description='每列的和的值为 100%，其他数值是除以总值百分比，输出新的结果')
    p.add_argument(dest='input_file', help='输入文件，有 Header，第一列是 row names')
    p.add_argument(dest='output_file', help='输出文件')
    
    return p.parse_args()


def data2fpkm(input_file, output_file):
    df = pd.read_csv(input_file)
    
    if not df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').notnull().all().all():
        logger.critical("输入文件中包含非数值数据，请检查。")
        sys.exit(1)
    
    row_names = df.iloc[:, 0]
    column_sums = df.iloc[:, 1:].sum(axis=0)
    normalized_df = df.iloc[:, 1:].div(column_sums, axis=1) * 100
    normalized_df.insert(0, df.columns[0], row_names)
    
    # 输出到文件
    normalized_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    args = parse_input()
    data2fpkm(args.input_file, args.output_file)
    