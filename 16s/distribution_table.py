#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/10/16 11:10
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument(dest='input_table_file', help='输入文件 csv xlsx txt 格式')
    p.add_argument(dest='output_table_file', help='输出文件')
    p.add_argument('-t', '--taxa-type', dest='taxa_type', default='genus', help='taxon type')
    
    return p.parse_args()


def process_distribution_table(input_df, taxa_type):
    taxa_type = taxa_type.lower()
    taxa_map = {
        'genus': 'g__',
        'species': 's__',
        'phylum': 'p__',
        'class': 'c__',
        'order': 'o__',
        'family': 'f__'
    }
    if taxa_type not in taxa_map:
        raise ValueError('taxa type must be one of genus, species, phylum, class, order, family')
    str_contain = taxa_map[taxa_type]
    index_col = input_df.columns[0]

    # 1. 筛选包含 taxa_type 并去除 unclassified，排序
    has_strcontain_df = input_df[
        input_df[index_col].str.contains(str_contain) &
        ~input_df[index_col].str.contains(f'{str_contain}unclassified__{taxa_type}')
    ].copy()
    has_strcontain_df[index_col] = has_strcontain_df[index_col].str.split(str_contain).str[-1]
    has_strcontain_df['row_sum'] = has_strcontain_df.sum(axis=1, numeric_only=True)
    has_strcontain_df = has_strcontain_df.sort_values(by='row_sum', ascending=False)
    has_strcontain_df.drop(columns=['row_sum'], inplace=True)

    # 2. 筛选不包含 taxa_type 的行
    no_strcontain_df = input_df[~input_df[index_col].str.contains(str_contain)].copy()

    # 3. 合并
    all_df = pd.concat([has_strcontain_df, no_strcontain_df], ignore_index=True)

    # 4. 处理 top15 和 Others
    if all_df.shape[0] > 15:
        top15_df = all_df.iloc[:15]
        others_sum = all_df.iloc[15:].sum(numeric_only=True)
        others_row = pd.DataFrame([["Others"] + list(others_sum)], columns=all_df.columns)
        output_df = pd.concat([top15_df, others_row], ignore_index=True)
    else:
        output_df = all_df.copy()

    # 5. 如果有 unclassified_{taxa_type}，合并到 Others
    unclassified_mask = input_df[index_col].str.contains(f'{str_contain}unclassified__{taxa_type}')
    if unclassified_mask.any():
        unclassified_sum = input_df[unclassified_mask].sum(numeric_only=True)
        if (output_df[index_col] == 'Others').any():
            output_df.loc[output_df[index_col] == 'Others', output_df.columns[1]:] += unclassified_sum.values
        else:
            others_row = pd.DataFrame([["Others"] + list(unclassified_sum)], columns=output_df.columns)
            output_df = pd.concat([output_df, others_row], ignore_index=True)

    return output_df


if __name__ == '__main__':
    args = parse_input()
    input_df = load_table(args.input_table_file)
    output_df = process_distribution_table(input_df, taxa_type=args.taxa_type)
    write_output_df(output_df, args.output_table_file, index=False)