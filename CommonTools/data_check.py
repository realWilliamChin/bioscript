#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/04/01 11:06
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger
import re


def parse_input():
    """解析命令行参数，使用子命令结构"""
    p = argparse.ArgumentParser(description='数据检查和清理工具')
    subparsers = p.add_subparsers(dest='command', help='可用命令', required=True)
    
    # 通用参数函数
    def add_common_args(parser):
        parser.add_argument('--input', '-i', required=True, help='输入 table，任意格式')
        parser.add_argument('--output', '-o', required=True, help='输出 table，任意格式')
        parser.add_argument('--no-input-header', action='store_true',
                           help='输入文件不包含 header，则输出文件也不包含 header')
        parser.add_argument('--sep', default=None, help='输入文件分隔符，默认自动检测')
        parser.add_argument('--output-sep', default=None, help='输出文件分隔符，默认与输入相同')
    
    # 子命令: drop_element_side_space
    parser_drop_space = subparsers.add_parser('drop_element_side_space',
                                               help='去除 DataFrame 中每个元素两侧的空格')
    add_common_args(parser_drop_space)
    
    # 子命令: drop_row_sum_eq_zero
    parser_drop_zero = subparsers.add_parser('drop_row_sum_eq_zero',
                                             help='删除所有数值型列和为0的行')
    add_common_args(parser_drop_zero)
    
    # 子命令: replace_illegal_folder_chars
    parser_replace_chars = subparsers.add_parser('replace_illegal_folder_chars',
                                                 help='替换指定列中不能作为文件夹名的非法字符')
    add_common_args(parser_replace_chars)
    parser_replace_chars.add_argument('--columns', '-c', required=True, nargs='+',
                                      help='需要处理的列名列表')
    parser_replace_chars.add_argument('--replace-with', default='_',
                                      help='替换成的字符，默认为下划线')
    
    # 子命令: convert_numeric_columns
    parser_convert_numeric = subparsers.add_parser('convert_numeric_columns',
                                                   help='转换数字但类型是字符串的列为数值类型')
    add_common_args(parser_convert_numeric)
    parser_convert_numeric.add_argument('--exclude-columns', nargs='+', default=None,
                                        help='指定这些列不做转换')
    
    args = p.parse_args()
    return args


def read_dataframe(input_file, no_header=False, sep=None):
    """读取 DataFrame，支持自动检测分隔符"""
    if sep is None:
        # 尝试自动检测分隔符
        with open(input_file, 'r', encoding='utf-8') as f:
            first_line = f.readline()
            if '\t' in first_line:
                sep = '\t'
            elif ',' in first_line:
                sep = ','
            else:
                sep = None  # 让 pandas 自动检测
    
    if no_header:
        df = pd.read_csv(input_file, sep=sep, header=None)
    else:
        df = pd.read_csv(input_file, sep=sep)
    
    return df, sep


def write_dataframe(df, output_file, no_header=False, sep=None, output_sep=None):
    """写入 DataFrame"""
    if output_sep is None:
        output_sep = sep if sep is not None else '\t'
    
    if no_header:
        df.to_csv(output_file, sep=output_sep, index=False, header=False)
    else:
        df.to_csv(output_file, sep=output_sep, index=False)
    
    logger.info(f'结果已保存到: {output_file}')


def file_check():
    pass


def dir_check():
    pass


def dataframe_check():
    pass


def df_drop_element_side_space(df):
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df


def df_drop_row_sum_eq_zero(df: pd.DataFrame) -> pd.DataFrame:
    """
    建议先运行 convert_numeric_columns 强制转换所有可能为数值列的类型
    删除所有数值型列和为0的行
    输入: pd.DataFrame
    返回: pd.DataFrame
    """
    numeric_cols = df.select_dtypes(include='number').columns  # 选择数值列
    row_sums = df[numeric_cols].sum(axis=1)                   # 计算每行和
    before = df.shape[0]
    filtered_df = df[row_sums != 0]
    after = filtered_df.shape[0]
    if before != after:
        logger.info(f'去除数值型行和为0的行: 原始{before}行, 过滤后{after}行, 去除{before-after}行')
    return filtered_df


def df_replace_illegal_folder_chars(df, columns, replace_with="_"):
    """
    替换指定列中不能作为文件夹名的非法字符
    :param df: pandas.DataFrame
    :param columns: 需要处理的列名列表
    :param replace_with: 替换成的字符，默认为下划线
    :return: 处理后的 DataFrame
    """
    # Windows 和 Linux 下常见非法字符
    illegal_chars = r'[ \/\\:\*\?\'"<>\|\(\)（）【】]'
    for col in columns:
        if col in df.columns:
            # 替换非法字符为下划线，并合并连续下划线为一个下划线
            # 将 _-_ 替换为 -
            df[col] = df[col].astype(str).apply(
                lambda x: re.sub(r'_-_', '-', re.sub(r'_+', '_', re.sub(illegal_chars, replace_with, x)))
            )
        else:
            logger.warning(f'输入的 {col} 没有在输入文件中，不会对当前列处理非法字符问题')
    return df


def str_replace_illegal_folder_chars(str, replace_with="_"):
    """
    替换指定字符串中不能作为文件夹名的非法字符
    :param str: 需要处理的字符串
    :param replace_with: 替换成的字符，默认为下划线
    :return: 处理后的字符串
    """
    illegal_chars = r'[ \/\\:\*\?\'"<>\|\(\)（）【】]'
    return re.sub(illegal_chars, replace_with, str)


def convert_numeric_columns(df: pd.DataFrame, exclude_columns=None, convert_dtype='float') -> pd.DataFrame:
    """
    转换数字但类型是字符串的列为数值类型。
    空字符串会被视为缺失值允许保留为 NaN。
    :param df: 待处理的DataFrame
    :param exclude_columns: 可迭代的列名，指定这些列不做转换（可为None）
    :param dtype: 转换后的数值类型，'float' 或 'int'，默认为 'float'
    :return: 转换后的DataFrame
    """
    if exclude_columns is None:
        exclude_columns = []
    df_converted = df.copy()
    for col in df_converted.columns:
        if col in exclude_columns:
            logger.debug(f'列 {col} 在排除列表中，跳过数值转换')
            continue
        ser = df_converted[col].astype(str).str.strip()
        ser = ser.replace({"": None})
        # 尝试数值转换（无法转换的字符串会变为 NaN）
        converted = pd.to_numeric(ser, errors='coerce')
        # 获取非空值的掩码（排除原本就是空值的位置）
        non_empty_mask = ser.notna()
        # 只有当所有非空值都能成功转换为数值时，才转换该列
        # 如果存在无法转换的字符串（如 "abc"），converted 中对应位置会是 NaN
        # 此时 converted[non_empty_mask].notna().all() 会返回 False，该列不会被转换
        if non_empty_mask.any() and converted[non_empty_mask].notna().all():
            if convert_dtype == 'int':
                # 只在所有非空都成功转换后且无缺失值时可转为 int
                if converted.isna().any():
                    # 存在缺失值不能强转为 int，否则会报错，只能保持 float
                    df_converted[col] = converted
                else:
                    # 能完整转为 int
                    df_converted[col] = converted.astype(int)
            else:
                df_converted[col] = converted
        elif not non_empty_mask.any():
            # 如果列全部为空，也可以转换为数值类型（全部为 NaN 的数值列），但不能是 int（全为 NaN 的 int 是不合法的）
            df_converted[col] = converted  # 保持 float 全 NaN
        else:
            # 存在无法转换为数值的非空字符串，记录详细日志方便排查
            bad_mask = non_empty_mask & converted.isna()
            bad_values = ser[bad_mask].unique()
            # 只展示部分示例，避免日志过长
            example_values = bad_values[:5]
            logger.warning(
                f"列 {col} 中存在无法转换为数值的非空值，共 {bad_mask.sum()} 个，例如: {list(example_values)}；"
                "该列将保持原始类型，不进行数值转换"
            )
    return df_converted


def main():
    args = parse_input()
    
    # 读取输入文件
    df, input_sep = read_dataframe(args.input, args.no_input_header, args.sep)
    logger.info(f'读取文件: {args.input}, 形状: {df.shape}')
    
    # 根据子命令执行相应操作
    if args.command == 'drop_element_side_space':
        df = df_drop_element_side_space(df)
    
    elif args.command == 'drop_row_sum_eq_zero':
        # 建议先转换数值列
        df = convert_numeric_columns(df)
        df = df_drop_row_sum_eq_zero(df)
    
    elif args.command == 'replace_illegal_folder_chars':
        df = df_replace_illegal_folder_chars(df, args.columns, args.replace_with)
    
    elif args.command == 'convert_numeric_columns':
        df = convert_numeric_columns(df, exclude_columns=args.exclude_columns)
        logger.info(f'数值列转换完成，转换类型: {args.convert_dtype}')
    
    # 写入输出文件
    output_sep = args.output_sep if hasattr(args, 'output_sep') else None
    write_dataframe(df, args.output, args.no_input_header, input_sep, output_sep)


if __name__ == '__main__':
    main()
