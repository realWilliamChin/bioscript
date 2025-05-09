#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/04/01 11:06
# Author        : William GoGo
import os, sys
from loguru import logger


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