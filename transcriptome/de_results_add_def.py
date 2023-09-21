#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/4/28 14:46
# Author        : WilliamGoGo
import os
import pandas as pd
from genedf_add_def import add_kns_def


def parse_input():
    """
    解析输入参数，输入 kegg, nr, swiss file 的路径
    """
    import argparse
    parser = argparse.ArgumentParser(description='输入 kegg, nr, swiss file 的路径')
    parser.add_argument('-k', '--kegg', type=str, help='kegg file')
    parser.add_argument('-n', '--nr', type=str, help='nr file')
    parser.add_argument('-s', '--swiss', type=str, help='swiss file')
    parser.add_argument('--kns', type=str, help='输入 kns_def.txt，则不用输入上面的三个')

    return parser.parse_args()


def process_deresults(de_results_file, kegg_file, nr_file, swiss_file, kns_file):
    de_df = pd.read_csv(de_results_file, sep='\t', dtype={"GeneID": str})
    de_matrix_df = pd.read_csv(de_results_file + '_readCounts.matrix', sep='\t', dtype={"GeneID": str})
    de_df = pd.merge(left=de_df, right=de_matrix_df, on='GeneID', how='left')
    # 排序 down，up，NOsig。down 的 FC 值从小到大，up 的 FC 值从大到小
    # 先分三份，再合并
    de_df_down = de_df[de_df['regulation'] == 'Down'].copy()
    de_df_down.sort_values(by='FC', ascending=True, inplace=True)
    de_df_up = de_df[de_df['regulation'] == 'Up'].copy()
    de_df_up.sort_values(by='FC', ascending=False, inplace=True)
    de_df_nosig = de_df[de_df['regulation'] == 'NoSignificant'].copy()
    # 合并
    de_df = pd.concat([de_df_down, de_df_up, de_df_nosig])
    # 添加注释
    de_df = add_kns_def(de_df, kegg_file, nr_file, swiss_file, kns_file)
    de_df.to_csv(de_results_file.replace('DE_results', 'DEG_data.txt'), sep='\t', index=False)


def main():
    args = parse_input()

    for de_results_file in os.listdir():
        if de_results_file.endswith('DE_results'):
            print(f'processing---{de_results_file}')
            process_deresults(de_results_file, args.kegg, args.nr, args.swiss, args.kns)


if __name__ == '__main__':
    main()


