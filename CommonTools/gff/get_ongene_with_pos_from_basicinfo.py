#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/11/28 16:24
# Author        : William GoGo
import argparse
import os
import sys

import pandas as pd
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(description='根据 POS 与 basicinfo start/end 判断是否 On_Gene')
    parser.add_argument('-i', '--input', required=True, help='输入表格文件，至少包含 CHROM、POS 列')
    parser.add_argument('-r', '--reference', required=True, help='basicinfo 或包含 Chromosome/Start/End/Strand 列的参考文件')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径')
    return parser.parse_args()


def ensure_columns(df, required_cols, label):
    missing = required_cols - set(df.columns)
    if missing:
        logger.error(f'{label} 缺少必需列: {missing}')
        sys.exit(1)


def prepare_reference(reference_df):
    ref_columns = reference_df.columns.tolist()
    ref_df = reference_df.copy()
    ref_df['Start'] = pd.to_numeric(ref_df['Start'], errors='coerce')
    ref_df['End'] = pd.to_numeric(ref_df['End'], errors='coerce')
    invalid = ref_df['Start'].isna() | ref_df['End'].isna()
    if invalid.any():
        logger.warning(f'参考文件有 {invalid.sum()} 行缺失 start/end，将被跳过')
        ref_df = ref_df[~invalid].copy()
    ref_df['Interval_Low'] = ref_df[['Start', 'End']].min(axis=1)
    ref_df['Interval_High'] = ref_df[['Start', 'End']].max(axis=1)
    grouped = {str(chrom): sub_df for chrom, sub_df in ref_df.groupby('Chromosome')}
    return grouped, ref_columns


def merge_on_gene(input_df, ref_groups, ref_columns):
    result_rows = []
    columns = input_df.columns.tolist()
    for row in input_df.itertuples(index=False, name=None):
        row_dict = dict(zip(columns, row))
        pos = row_dict['POS']
        chrom_val = str(row_dict['CHROM'])

        def append_off_gene():
            combined = {**row_dict}
            for col in ref_columns:
                combined.setdefault(col, pd.NA)
            combined['On_Gene_Status'] = 'Off_Gene'
            result_rows.append(combined)

        if pd.isna(pos):
            append_off_gene()
            continue

        chrom_ref = ref_groups.get(chrom_val)
        if chrom_ref is None:
            append_off_gene()
            continue

        matches = chrom_ref[(chrom_ref['Interval_Low'] <= pos) & (chrom_ref['Interval_High'] >= pos)]
        if matches.empty:
            append_off_gene()
            continue

        for _, match in matches.iterrows():
            combined = {**row_dict}
            for col in ref_columns:
                combined[col] = match[col]
            combined['On_Gene_Status'] = 'On_Gene'
            result_rows.append(combined)

    return pd.DataFrame(result_rows)


def main():
    args = parse_input()
    input_df = load_table(args.input, low_memory=False)
    reference_df = load_table(args.reference, low_memory=False)

    ensure_columns(input_df, {'CHROM', 'POS'}, '输入文件')
    ensure_columns(reference_df, {'Chromosome', 'Start', 'End', 'Strand'}, '参考文件')

    input_df['POS'] = pd.to_numeric(input_df['POS'], errors='coerce')
    if input_df['POS'].isna().any():
        nan_rows = input_df['POS'].isna().sum()
        logger.warning(f'有 {nan_rows} 行 POS 无法转换为数字，将被视为 Off_Gene')

    ref_groups, ref_columns = prepare_reference(reference_df)
    merged_df = merge_on_gene(input_df, ref_groups, ref_columns)
    write_output_df(merged_df, args.output, index=False)
    logger.success(f'输出完成: {args.output}')


if __name__ == '__main__':
    main()