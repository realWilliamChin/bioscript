#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/29 13:58
# Author        : William GoGo
import os, sys
import numpy as np
import pandas as pd
import argparse
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser(description='每列的和的值为 100%，其他数值是除以总值百分比，输出新的结果')
    p.add_argument(dest='input_file', help='输入文件，有 Header，第一列是 row names')
    p.add_argument(dest='output_file', help='输出文件')
    
    return p.parse_args()


def data2fpkm(df):
    if not df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').notnull().all().all():
        logger.critical("输入文件中包含非数值数据，请检查。")
        sys.exit(1)
    
    row_names = df.iloc[:, 0]
    
    column_sums = df.iloc[:, 1:].sum(axis=0)
    column_sums = np.where(column_sums == 0, 1, column_sums)
    normalized_df = df.iloc[:, 1:].div(column_sums, axis=1) * 100
    normalized_df.insert(0, df.columns[0], row_names)
    
    return normalized_df


if __name__ == "__main__":
    args = parse_input()
    df = load_table(args.input_file)
    out_df = data2fpkm(df)
    
    write_output_df(out_df, args.output_file, index=False)