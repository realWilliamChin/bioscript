#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from typing import Any

import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df

def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='输入文件')
    parser.add_argument('-o', '--output', help='输出文件')
    return parser.parse_args()


def proportion(df: pd.DataFrame) -> pd.DataFrame:
    """
    将原始丰度矩阵转换为百分比相对丰度矩阵。

    行: 样本 / 特征 (与输入保持一致, 默认行为: 按行求和后标准化为 100%)

    Args:
        df (pd.DataFrame): 原始 count 矩阵

    Returns:
        pd.DataFrame: 相对丰度(%)矩阵
    """
    # 仅对数值型样本列做运算, 避免把基因ID等字符串列参与除法
    numeric_df = df.select_dtypes(include="number")
    meta_df = df.drop(columns=numeric_df.columns, errors="ignore")

    # 行方向求和并标准化为 100%
    numeric_t = numeric_df.T
    numeric_t["count"] = numeric_t.sum(axis=1)

    # 防止 0 求和导致除以 0
    zero_mask = numeric_t["count"] == 0
    if zero_mask.any():
        logger.warning(f"检测到 {zero_mask.sum()} 行总和为 0, 将其相对丰度设为 0")
        numeric_t.loc[zero_mask, "count"] = 1

    numeric_t = numeric_t.div(numeric_t["count"], axis=0).mul(100)
    numeric_t = numeric_t.drop(columns=["count"])

    result = numeric_t.T

    # 若存在非数值信息列则拼回去
    if not meta_df.empty:
        result = pd.concat([meta_df, result], axis=1)

    return result


def main():
    """
    使用 CommonTools.load_input 读取输入表格,
    计算相对丰度后再写出到输出文件。
    """
    args = parse_input()

    logger.info(f"读取输入文件: {args.input}")
    df = load_table(args.input)

    logger.info("开始计算相对丰度(%)矩阵")
    df_prop = proportion(df)

    logger.info(f"写出结果到: {args.output}")
    write_output_df(df_prop, args.output, index=False)


if __name__ == "__main__":
    main()
