#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/11/27 14:06
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='输入表格格式文件，包含 chrom region_start 和 region_end，其他列自动保留')
    p.add_argument('-b', '--basicinfo', help='输入 basicinfo 文件')
    p.add_argument('-o', '--output', help='输出文件')
    
    return p.parse_args()


def main():
    args = parse_input()
    
    result_df_lst = []
    
    df = pd.read_csv(args.input)
    df_columns = df.columns.tolist()
    # GeneID Chr_Number Start End Strand Gene_Def
    basic_info_df = pd.read_csv(args.basicinfo, sep='\t')
    basic_info_df['Chr_Number'] = basic_info_df['Chr_Number'].astype(str)
    df['CHROM'] = df['CHROM'].astype(str)
    
    # print(df.dtypes)
    # print(basic_info_df.dtypes)
    
    for each_row in df.itertuples():
        each_row_df = basic_info_df[
            (each_row.region_start <= basic_info_df['Start']) &
            (each_row.region_end >= basic_info_df['End']) &
            (each_row.CHROM == basic_info_df['Chr_Number'])
        ].copy()
        # 追加其他所有信息
        for column in df_columns:
            each_row_df[column] = getattr(each_row, column)
        result_df_lst.append(each_row_df)
    
    result_df = pd.concat(result_df_lst)
    
    if result_df.empty:
        logger.error('返回为空')
        sys.exit(1)
    
    result_df.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()