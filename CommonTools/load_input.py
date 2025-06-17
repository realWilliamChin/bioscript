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
    if ext in ['xls', 'xlsx']:
        reader = lambda file, **kwargs: pd.read_excel(file, **kwargs, engine='openpyxl')
    elif ext in ['txt', 'tsv', 'blast', 'gff', 'gff3', 'gtf']:
        reader = lambda file, **kwargs: pd.read_csv(file, sep='\t', **kwargs)
    else:
        reader = lambda file, **kwargs: pd.read_csv(file, **kwargs)

    # 读取数据, 并将 *args, **kwargs 传递给读取函数
    df = reader(table_file, *args, **kwargs)
    
    # 检查列名中是否包含 geneid/GeneID/gene_id, 如果有则将其转换为字符串类型（经常忘，写在这里）
    gene_id_cols = [col for col in df.columns if 'geneid' in str(col).lower() or 'gene_id' in str(col).lower()]
    if gene_id_cols:
        df[gene_id_cols] = df[gene_id_cols].astype(str)
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
    if ext in ['xls', 'xlsx']:
        df.to_excel(output_file, engine='openpyxl', *args, **kwargs)
    elif ext in ['txt', 'tsv', 'blast', 'gff', 'gff3', 'gtf']:
        df.to_csv(output_file, sep='\t', *args, **kwargs)
    else:
        df.to_csv(output_file, *args, **kwargs)


def main():
    args = parse_input()
    df = load_table(args.input_file, *args.read_options)
    write_output_df(df, args.output_file, *args.write_options)


if __name__ == '__main__':
    main()
