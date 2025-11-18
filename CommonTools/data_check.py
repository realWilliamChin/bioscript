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
    p = argparse.ArgumentParser()
    p.add_argument('-i', help='输入 table，任意格式，必须包含 Header')
    p.add_argument('--no-input-header', action='store_true',
                   help='输入文件不包含 header，则输出文件也不包含 header')
    p.add_argument('-o', help='输出 table，任意格式')

    p.add_argument('-r', '--run', choices=['drop_element_side_space', 'drop_row_sum_eq_zero', 'replace_illegal_folder_chars'])
    
    args = p.parse_args()
    
    return args


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
    删除所有数值型列和为0的行。
    输入: pd.DataFrame
    返回: pd.DataFrame
    """
    df = convert_numeric_columns(df)
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
    illegal_chars = r'[ \/\\:\*\?"<>\|]'
    for col in columns:
        if col in df.columns:
            df[col] = df[col].astype(str).apply(lambda x: re.sub(illegal_chars, replace_with, x))
        else:
            logger.warning(f'输入的 {col} 没有在输入文件中，不会对当前列处理非法字符问题')
    return df


def convert_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    遍历所有列：若某列中所有非空值均可成功转为数值，则将该列转换为数值类型。
    空字符串会被视为缺失值允许保留为 NaN。
    """
    df_converted = df.copy()
    for col in df_converted.columns:
        # 先将值标准化为字符串，去除首尾空白，并将空字符串视为缺失
        ser = df_converted[col].astype(str).str.strip()
        ser = ser.replace({"": None})
        # 尝试数值转换（不可转的变为 NaN）
        converted = pd.to_numeric(ser, errors='coerce')
        non_empty_mask = ser.notna()
        # 所有非空项均可成功转为数值则接受转换
        if converted[non_empty_mask].notna().all():
            df_converted[col] = converted
    return df_converted


def main():
    args = parse_input()
    


if __name__ == '__main__':
    main()
