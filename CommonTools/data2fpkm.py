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
    p = argparse.ArgumentParser(description='将 reads 转换为 fpkm 值')
    p.add_argument(dest='input_file', help='输入文件，有 Header，第一列是 row names')
    p.add_argument(dest='output_file', help='输出文件')
    
    p.add_argument('--geneinfo', help='kns 文件或 basicinfo 文件')
    
    return p.parse_args()


def replace_negative_values(df, columns=None):
    """
    将数据框中的负数值替换为0.001
    
    参数:
    df: pandas DataFrame
    columns: 需要处理的列名列表，如果为 None 则处理所有数值列
    
    返回:
    处理后的 DataFrame
    """
    df_copy = df.copy()
    if columns is None:
        # 获取所有数值类型的列
        columns = df_copy.select_dtypes(include=['float64', 'int64']).columns
    
    for col in columns:
        # df_copy[col] = df_copy[col].apply(lambda x: 0.000001 if x < 0 else int(x))
        df_copy[col] = df_copy[col].apply(lambda x: 0.000001 if x < 0 else x)
    
    return df_copy


def rawreads_convert_fpkm(reads_df, geneinfo_df):
    """
    将raw reads count数据转换为FPKM值
    fpkm <- (raw_count * 10^9) / (gene_number * gene_length)
    
    参数:
    reads_df: pandas DataFrame, 包含GeneID列和样本的raw count数据
    kns_df: pandas DataFrame, 包含GeneID和 Start 列 和 End 的基因长度信息
    
    返回:
    pandas DataFrame: FPKM标准化后的表达量数据
    """
    
    geneinfo_df['Gene_Length'] = (geneinfo_df['End'] - geneinfo_df['Start']).abs() + 1
    
    index_column = reads_df.columns.tolist()[0]
    reads_df = reads_df[reads_df[index_column].isin(geneinfo_df[index_column].values.tolist())]
    sample_cols = [col for col in reads_df.columns if col != index_column]
    merged_df = pd.merge(reads_df, geneinfo_df[[index_column, 'Gene_Length']], on=index_column, how='inner')
    fpkm_df = merged_df[[index_column]].copy()

    for col in sample_cols:
        # X FPKM = (reads_count * 10^9) / (total_reads * gene_length)
        # √ FPKM = (reads_count * 10^9) / (gene_number * gene_length)
        fpkm_df[col] = (merged_df[col] * 1e9) / (merged_df.shape[0] * merged_df['Gene_Length'])

    # 将负值或0值替换为一个很小的正数
    fpkm_df = replace_negative_values(fpkm_df, columns=sample_cols)
    new_reads_df = reads_df[reads_df[index_column].isin(fpkm_df[index_column].values)]
    
    return fpkm_df, new_reads_df


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
    if args.geneinfo:
        geneinfo_df = load_table(args.geneinfo)
        fpkm_df, new_reads_df = rawreads_convert_fpkm(df, geneinfo_df)
        write_output_df(fpkm_df, args.output_file, index=False)
        write_output_df(new_reads_df, f'{args.output_file}.reads.txt', index=False)
    else:
        out_df = data2fpkm(df)
        write_output_df(out_df, args.output_file, index=False)