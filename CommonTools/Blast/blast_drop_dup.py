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


def validate_arguments(args):
    """
    验证命令行参数的有效性
    
    Args:
        args: 解析后的参数对象
    
    Raises:
        argparse.ArgumentError: 当参数验证失败时抛出
    """
    dedup_params = [args.column, args.sort_by, args.keep]
    filter_params = [args.filter_by, args.filter_condition, args.filter_value]
    
    # 验证参数组的互斥性
    if any(dedup_params) and any(filter_params):
        raise argparse.ArgumentError(None, 
            "去重参数和过滤参数不能同时使用，请选择其中一种功能")
    
    if not any(dedup_params) and not any(filter_params):
        raise argparse.ArgumentError(None, 
            "必须使用去重参数或过滤参数中的一种功能")
    
    # 验证过滤参数完整性
    if any(filter_params):
        if not all(filter_params):
            raise argparse.ArgumentError(None, 
                "过滤参数不完整，需要同时提供 --filter-by、--filter-condition 和 --filter-value")
    
    # 验证排序参数
    if args.keep in ['max', 'min'] and not args.sort_by:
        raise argparse.ArgumentError(None, 
            "使用 max/min 保留策略时必须指定 --sort-by 参数")


def parse_arguments():
    parser = argparse.ArgumentParser(description='BLAST结果去重或过滤')
    
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
                               help='输出文件路径，默认在原文件名后添加_uniq或_filtered后缀')
    output_group.add_argument('--o-columns', dest='o_columns',
                                help='输出文件的列名，使用逗号分隔，默认为输入的 columns，通常对输出文件重命名 header 时使用')
    output_group.add_argument('--o-no-header', dest='o_no_header', action='store_true',
                               help='不输出表头，默认输出表头')
    
    # 去重参数组
    dedup_group = parser.add_argument_group('去重参数（与过滤参数互斥）')
    dedup_group.add_argument('-c', '--column', default='qacc', 
                              help='去重的列名称，通常是 qacc 或者 GeneID 列')
    dedup_group.add_argument('--sort-by', dest='sort_by',
                              help='用于排序的列名称')
    dedup_group.add_argument('--keep', choices=['first', 'last', 'max', 'min'], default='first',
                              help='去重时保留的策略：first(第一个), last(最后一个), max(最大值), min(最小值)')
    
    # 过滤参数组
    filter_group = parser.add_argument_group('过滤参数（与去重参数互斥）')
    filter_group.add_argument('--filter-by', dest='filter_by',
                              help='用于过滤的列名称')
    filter_group.add_argument('--filter-condition', dest='filter_condition', choices=['>', '<', '>=', '<=', '=='],
                              help='过滤条件：>, <, >=, <=, ==')
    filter_group.add_argument('--filter-value', dest='filter_value', type=float,
                              help='过滤的阈值')
    
    args = parser.parse_args()
    
    # 检查参数组的互斥性
    dedup_params = [args.column, args.sort_by, args.keep]
    filter_params = [args.filter_by, args.filter_condition, args.filter_value]
    
    if any(dedup_params) and any(filter_params):
        parser.error("去重参数和过滤参数不能同时使用，请选择其中一种功能")
    
    if not any(dedup_params) and not any(filter_params):
        parser.error("必须使用去重参数或过滤参数中的一种功能")
    
    # 处理输入输出列名
    args.i_columns = [x.strip() for x in args.i_columns.split(',')]
    if args.o_columns:
        args.o_columns = [x.strip() for x in args.o_columns.split(',')]
    else:
        args.o_columns = args.i_columns
    
    # 设置输出文件名
    if not args.outfile:
        if any(filter_params):
            args.outfile = str(args.infile).replace('.blast', '_filtered.blast')
        else:
            args.outfile = str(args.infile).replace('.blast', '_uniq.blast')
    
    return args


def filter_blast_results(df, filter_by, filter_condition, filter_value):
    """
    过滤BLAST结果
    
    Args:
        df (pd.DataFrame): 输入的DataFrame
        filter_by (str): 过滤的列名
        filter_condition (str): 过滤条件
        filter_value (float): 过滤阈值
    
    Returns:
        pd.DataFrame: 过滤后的DataFrame
    """
    condition = f"{filter_by} {filter_condition} {filter_value}"
    filtered_df = df.query(condition)
    logger.info(f"应用过滤条件: {condition}")
    logger.info(f"过滤后剩余行数: {len(filtered_df)}")
    return filtered_df


def deduplicate_blast_results(df, column, sort_by, keep):
    """
    对BLAST结果进行去重
    
    Args:
        df (pd.DataFrame): 输入的DataFrame
        column (str): 去重的列名
        sort_by (str): 排序的列名
        keep (str): 保留策略
    
    Returns:
        pd.DataFrame: 去重后的DataFrame
    """
    if keep in ['max', 'min']:
        ascending = keep == 'min'
        df = df.sort_values(by=[column, sort_by], 
                           ascending=[True, ascending])
        df = df.groupby(column).first().reset_index()
    else:
        # 排序和去重
        if sort_by:
            df.sort_values(by=[column, sort_by], 
                          ascending=[True, ascending], inplace=True)
        
        if keep == 'first':
            df.drop_duplicates(subset=[column], keep='first', inplace=True)
        elif keep == 'last':
            df.drop_duplicates(subset=[column], keep='last', inplace=True)
    
    logger.info(f"去重后剩余行数: {len(df)}")
    return df


def main():
    args = parse_arguments()
    
    # 读取输入文件
    df = load_table(args.infile, header=0 if args.i_has_header else None, names=args.i_columns)
    
    # 检查使用的是哪个功能组
    dedup_params = [args.column, args.sort_by, args.keep]
    filter_params = [args.filter_by, args.filter_condition, args.filter_value]
    
    if any(filter_params):
        # 执行过滤功能
        df = filter_blast_results(df, args.filter_by, args.filter_condition, args.filter_value)
    else:
        # 执行去重功能
        df = deduplicate_blast_results(df, args.column, args.sort_by, args.keep)
    
    # 输出文件
    write_output_df(df, args.outfile, index=False, header=False if args.o_no_header else True)


if __name__ == '__main__':
    main()
    