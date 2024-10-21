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
    p.add_argument(dest='input_table_file', help='输入文件')
    p.add_argument(dest='output_table_file', help='输出文件')
    p.add_argument('-t', '--taxa-type', dest='taxa_type', default='genus', help='taxon type')
    
    return p.parse_args()


def process_distribution_table(input_df, taxa_type):
    if taxa_type.lower() == 'genus':
        str_contain = 'g__'
    elif taxa_type.lower() == 'species':
        str_contain = 's__'
    elif taxa_type.lower() == 'phylum':
        str_contain = 'p__'
    elif taxa_type.lower() == 'class':
        str_contain = 'c__'
    elif taxa_type.lower() == 'order':
        str_contain = 'o__'
    elif taxa_type.lower() == 'family':
        str_contain = 'f__'
    else:
        raise ValueError('taxa type must be one of genus, species, phylum, class, order, family')
    # 选取索引字符串包含 str_contain 的
    has_strcontain_df = input_df[input_df['index'].str.contains(str_contain)].copy()
    has_strcontain_df['index'] = has_strcontain_df['index'].str.split(';').str[-1].str.replace(str_contain, '')
    has_strcontain_df['row_sum'] = has_strcontain_df.sum(axis=1, numeric_only=True)
    has_strcontain_df = has_strcontain_df.sort_values(by='row_sum', ascending=False)
    has_strcontain_df.drop(columns=['row_sum'], inplace=True)
    no_genus_df = input_df[~input_df['index'].str.contains(str_contain)].copy()
    all_df = pd.concat([has_strcontain_df, no_genus_df])
    all_df.set_index('index', inplace=True)
    # 选出前 15 行，其他按照列，计算总和 row_names 为 Others
    top15_df = all_df.iloc[:15, :]
    all_df.loc['Others', :] = all_df.iloc[15:, :].sum(axis=0, numeric_only=True)
    all_df = all_df.tail(1)
    output_df = pd.concat([top15_df, all_df])
    return output_df


if __name__ == '__main__':
    args = parse_input()
    input_df = load_table(args.input_table_file)
    output_df = process_distribution_table(input_df, taxa_type=args.taxa_type)
    write_output_df(output_df, args.output_table_file, index=True)