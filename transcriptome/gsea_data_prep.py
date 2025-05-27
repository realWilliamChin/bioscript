#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/10/24 09:30
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-e', '--expression', type=str, help='fpkm_matrix_filtered.txt')
    p.add_argument('-s', '--samples', type=str, help='samples_described.txt')
    p.add_argument('-c', '--compare', type=str, help='compare_info.txt')
    p.add_argument('-g', '--genesymbol', dest='genesymbol',
                   help='genesymbol.txt, 如果没有 genesymbol 文件，不输出 _genesymbol.txt 文件')
    
    args = p.parse_args()
    return args


def process_genesymbol(genesymbol_file):
    if not genesymbol_file:
        return None
    
    genesymbol_df = load_table(genesymbol_file, dtype=str, header=0, names=['GeneID', 'GeneSymbol'])
    return genesymbol_df.dropna(subset=['GeneSymbol'])


def write_expression_data(compare_df, compare, genesymbol_df=None):
    """写入表达数据文件"""
    # 写入 GeneID 版本
    fpkm_data_name = f'{compare}_GeneID.txt'
    write_output_df(compare_df, fpkm_data_name, index=False)
    
    # 如果存在 genesymbol，写入 GeneSymbol 版本
    if genesymbol_df is not None:
        genesymbol_data_df = pd.merge(genesymbol_df, compare_df, on='GeneID', how='inner')
        genesymbol_data_df.drop(columns=['GeneID'], inplace=True)
        genesymbol_data_name = f'{compare}_GeneSymbol.txt'
        write_output_df(genesymbol_data_df, genesymbol_data_name, index=False)


def get_expression_data(samples, compares, expression, genesymbol=None):
    # 加载所有输入数据
    genesymbol_df = process_genesymbol(genesymbol)
    samples_df = load_table(samples, dtype=str)
    compare_df = load_table(compares, dtype=str)
    expression_df = load_table(expression, dtype={'GeneID': str})
    
    # 处理每个比较组
    compare_df['compare'] = compare_df['Treat'] + '_vs_' + compare_df['Control']
    for compare in compare_df['compare']:
        treat_group, control_group = compare.split('_vs_')
        compares_sample_list = (
            samples_df[samples_df['group'] == treat_group]['sample'].tolist() +
            samples_df[samples_df['group'] == control_group]['sample'].tolist()
        )
        compare_data = expression_df[['GeneID'] + compares_sample_list].copy()
        write_expression_data(compare_data, compare, genesymbol_df)


def main():
    args = parse_input()
    get_expression_data(args.samples, args.compare, args.expression, args.genesymbol)
    logger.success('Done!')


if __name__ == '__main__':
    main()