#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/04/01 11:06
# Author        : William GoGo
import os, sys
from loguru import logger
import re


def file_check():
    pass


def dir_check():
    pass


def dataframe_check():
    pass


def df_drop_element_side_space(df):
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df


def df_drop_row_sum_eq_zero(df):
    numeric_cols = df.select_dtypes(include='number').columns  # 选择数值列
    row_sums = df[numeric_cols].sum(axis=1)                   # 计算每行和
    filtered_df = df[row_sums != 0]
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


