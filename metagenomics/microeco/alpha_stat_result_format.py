#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/03/17 11:42
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df

# 定义常见度量指标名称映射
MEASURE_MAPPING = {
    'Observed': 'observed_features',
    'Chao1': 'chao1',
    'ACE': 'ace',
    'Shannon': 'shannon',
    'Simpson': 'simpson',
    'InvSimpson': 'invsimpson',
    'Fisher': 'fisher',
    'Pielou': 'pielou_e',
    'Coverage': 'goods_coverage',
    'Faith_pd': 'faith_pd',
    'FaithPD': 'faith_pd'
}


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Convert microeco alpha diversity stat result from long format to wide format'
    )
    parser.add_argument('-i', '--input', required=True, help='Alpha_stats_single.xlsx')
    parser.add_argument('-o', '--output', help='Output file path')
    parser.add_argument('--no-rename', action='store_true', help='Do not rename measure names, keep original names')
    args = parser.parse_args()
    return args


def format_conversion(input_file, output_file, no_rename=False):
    """执行格式转换"""
    logger.info(f"Reading input file: {input_file}")
    # 读取输入文件，支持tab分隔
    df = load_table(input_file, usecols=['sample', 'Measure', 'Mean'])

    logger.info(f"Found {df['sample'].nunique()} samples and {df['Measure'].nunique()} measures")

    # 透视表转换：行=sample，列=Measure，值=Mean
    logger.info("Pivoting table to wide format")
    wide_df = df.pivot(index='sample', columns='Measure', values='Mean')

    # 重置索引，将sample变为一列
    wide_df = wide_df.reset_index()
    # 将sample列重命名为sample-iddd
    wide_df = wide_df.rename(columns={'sample': 'sample-id'})

    # 如果需要，重命名度量指标列
    if not no_rename:
        rename_dict = {}
        for col in wide_df.columns:
            if col in MEASURE_MAPPING:
                rename_dict[col] = MEASURE_MAPPING[col]
        if rename_dict:
            logger.info(f"Renaming {len(rename_dict)} measure columns")
            wide_df = wide_df.rename(columns=rename_dict)

    # 保存结果
    logger.info(f"Saving result to: {output_file}")
    write_output_df(wide_df, output_file, index=False, float_format='%.2f')
    logger.info(f"Conversion completed. Output shape: {wide_df.shape}")
    logger.info(f"Output columns: {', '.join(wide_df.columns)}")

    return wide_df


def main():
    args = parse_args()
    format_conversion(args.input, args.output, args.no_rename)


if __name__ == '__main__':
    main()