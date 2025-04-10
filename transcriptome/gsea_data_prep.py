#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/10/24 09:30
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser(help='有 GeneSymbol 使用 GeneSymbol 文件')
    p.add_argument('-e', '--expression', type=str, help='fpkm_matrix_filtered.txt')
    p.add_argument('-s', '--samples', type=str, help='samples_described.txt')
    p.add_argument('-c', '--compare', type=str, help='compare_info.txt')
    p.add_argument('-k', '--kegg-shortname', dest='kegg_shortname', type=str, help='kegg_shortname.txt')
    
    return p.parse_args()


def get_expression_data(samples, compares, kegg_shortname, expression):
    kegg_df = load_table(kegg_shortname, skiprows=1, dtype=str, header=None, names=['GeneID', 'Gene_shortname'])
    kegg_df.dropna(subset=['Gene_shortname'], inplace=True)
    samples_df = load_table(samples)
    compare_df = load_table(compares)
    expression_df = load_table(expression)
    
    compare_df['compare'] = compare_df['Treat'] + '_vs_' + compare_df['Control']
    compare_list = list(compare_df['compare'])
    for compare in compare_list:
        compares_sample_treat_list = samples_df[samples_df['group'] == compare.split('_vs_')[0]]['sample'].tolist()
        compares_sample_control_list = samples_df[samples_df['group'] == compare.split('_vs_')[1]]['sample'].tolist()
        compares_sample_list = compares_sample_treat_list + compares_sample_control_list
        compare_df = expression_df[['GeneID'] + compares_sample_list].copy()
        fpkm_data_name = f'{compare}_GeneID.txt'
        write_output_df(compare_df, fpkm_data_name, index=False)
        
        genesymbol_df = pd.merge(kegg_df, compare_df, on='GeneID', how='inner')
        genesymbol_df.drop(columns=['GeneID'], inplace=True)
        genesymbol_data_name = f'{compare}_genesymbol.txt'
        write_output_df(genesymbol_df, genesymbol_data_name, index=False)


def main():
    args = parse_input()
    get_expression_data(args.samples, args.compare, args.kegg_shortname, args.expression)


if __name__ == '__main__':
    main()