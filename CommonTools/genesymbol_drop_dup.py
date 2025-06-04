#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/05/09 16:02
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(description='Drop duplicate genesymbols')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
    return parser.parse_args()


def drop_dup_genesymbol(df):
    df = df.sort_values(by=['GeneID', 'GeneSymbol'], ascending=True)
    filter_before_number = df.shape[0]
    # 对每个重复的 GeneID 组进行处理
    def handle_duplicates(group):
        # 如果只有一行，直接返回
        if len(group) == 1:
            return group.iloc[0]
            
        # 检查 genesymbol 中是否包含 'LOC' 或 'AS1'
        loc_mask = group['GeneSymbol'].str.contains('LOC', na=False)
        as1_mask = group['GeneSymbol'].str.contains('AS1', na=False)
        
        # 如果有些含LOC/AS1有些不含，保留不含的第一个
        if (loc_mask.any() or as1_mask.any()) and not (loc_mask.all() or as1_mask.all()):
            filtered_group = group[~(loc_mask | as1_mask)]
            if len(filtered_group) > 0:
                return filtered_group.iloc[0]
        
        # 如果都含LOC/AS1或都不含，保留第一个
        return group.iloc[0]

    # 按GeneID分组并应用处理函数，添加 include_groups=False 参数
    df = df.groupby('GeneID', as_index=False).apply(handle_duplicates, include_groups=False)
    filter_after_number = df.shape[0]
    logger.info(f'Filter {filter_before_number - filter_after_number} duplicate genesymbols')
    return df


def main():
    args = parse_input()
    df = load_table(args.input)
    df = drop_dup_genesymbol(df)
    write_output_df(df, args.output, index=False)


if __name__ == '__main__':
    main()

