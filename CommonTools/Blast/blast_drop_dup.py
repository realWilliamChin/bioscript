#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/10/17 17:47
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_arguments():
    parser = argparse.ArgumentParser(description='BLAST结果去重')
    
    # 输入参数组
    input_group = parser.add_argument_group('输入参数')
    input_group.add_argument('-i', '--infile', required=True,
                              help='输入BLAST文件路径')
    input_group.add_argument('--i-has-header', dest='i_has_header', action='store_true',
                              help='输入文件是否包含表头')
    input_group.add_argument('--i-columns', dest='i_columns', required=True,
                               help='输入文件的列名，用逗号分隔，没有 header 也要输入')
    
    # 输出参数组
    output_group = parser.add_argument_group('输出参数')
    output_group.add_argument('-o', '--outfile',
                               help='输出文件路径，默认在原文件名后添加_uniq后缀')
    output_group.add_argument('--o-columns', dest='o_columns',
                                help='输出文件的列名，使用逗号分隔，默认为输入的 columns，通常对输出文件重命名 header 时使用')
    output_group.add_argument('--o-no-header', dest='o_no_header', action='store_true',
                               help='不输出表头，默认输出表头')
    
    # 去重参数组
    dedup_group = parser.add_argument_group('去重参数')
    dedup_group.add_argument('-c', '--column', type=str, default='qacc', 
                              help='去重的列名称，通常是 qacc 或者 GeneID 列')
    dedup_group.add_argument('--sort-by', dest='sort_by', type=str,
                              help='用于排序的列名称')
    dedup_group.add_argument('--keep', choices=['first', 'last', 'max', 'min'], default='first', required=True,
                              help='去重时保留的策略：first(第一个), last(最后一个), max(最大值), min(最小值)')
    
    args = parser.parse_args()
    
    # 处理输入输出列名
    args.i_columns = [x.strip() for x in args.i_columns.split(',')]
    if args.o_columns:
        args.o_columns = [x.strip() for x in args.o_columns.split(',')]
    else:
        args.o_columns = args.i_columns
    
    # 设置输出文件名
    if not args.outfile:
        args.outfile = str(args.infile).replace('.blast', '_uniq.blast')
    
    return args


def main():
    args = parse_arguments()
    
    # 读取输入文件
    df = load_table(args.infile, header=0 if args.i_has_header else None, names=args.i_columns)
    
    # 根据keep策略选择排序方式
    if args.keep in ['max', 'min']:
        ascending = args.keep == 'min'
    
    # 排序和去重
    if args.sort_by:
        df.sort_values(by=[args.column, args.sort_by], 
                      ascending=[True, ascending], inplace=True)
    
    if args.keep == 'first':
        df.drop_duplicates(subset=[args.column], keep='first', inplace=True)
    elif args.keep == 'last':
        df.drop_duplicates(subset=[args.column], keep='last', inplace=True)

    elif args.keep in ['max', 'min']:
        df = df.groupby(args.column).apply(
            lambda x: x.loc[x[args.sort_by].idxmax() if args.keep == 'max' 
                          else x[args.sort_by].idxmin()]
        ).reset_index(drop=True)
    
    # 输出文件
    write_output_df(df, args.outfile, index=False, header=False if args.o_no_header else True)


if __name__ == '__main__':
    main()
    