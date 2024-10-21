#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/29 14:02
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger
import pandas as pd


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument("input_file", help="输入文件的路径")
    p.add_argument("output_file", help="输出文件的路径")
    
    p.add_argument("--read-options", nargs="*", help="传递给 pd.read_* 的其他选项参数, 格式: key=value")
    p.add_argument("--write-options", nargs="*", help="传递给 pd.write_* 的其他选项参数, 格式: key=value")

    return p.parse_args()
    

def load_table(table_file, *args, **kwargs):
    """根据末尾格式读取输入的 table 文件

    末尾格式包括, txt, csv, tsv, xls, xlsx
    可自定义更多选项参数, 解析 header=*, usecols, names 等等

    Args:
        table_file (str): 输入文件的路径

    Returns:
        pd.DataFrame: 读取后的 DataFrame 对象
    """

    # 根据文件扩展名选择读取函数
    ext = table_file.split('.')[-1]
    if ext in ['csv', 'tsv']:
        reader = pd.read_csv
    elif ext in ['xls', 'xlsx']:
        reader = pd.read_excel
    elif ext == 'txt':
        reader = lambda file, **kwargs: pd.read_csv(file, sep='\t', **kwargs)
    else:
        raise ValueError(f"不支持的文件格式: {ext}")

    # 读取数据, 并将 *args, **kwargs 传递给读取函数
    df = reader(table_file, *args, **kwargs)
    return df


def write_output_df(df, output_file, *args, **kwargs):
    """根据末尾格式输出 table 文件

    末尾格式包括, txt, csv, tsv, xls, xlsx
    可自定义更多选项参数, 解析 header=*, usecols, names 等等

    Args:
        df (_type_): _description_
        output_file (_type_): _description_
    """
    ext = output_file.split('.')[-1]
    if ext in ['csv', 'tsv']:
        df.to_csv(output_file, *args, **kwargs)
    elif ext in ['xls', 'xlsx']:
        df.to_excel(output_file, *args, **kwargs)
    elif ext == 'txt':
        df.to_csv(output_file, sep='\t', *args, **kwargs)
    else:
        raise ValueError(f"不支持的文件格式: {ext}")


def main():
    args = parse_input()
    df = load_table(args.input_file, *args.read_options)
    write_output_df(df, args.output_file, *args.write_options)


if __name__ == '__main__':
    main()
