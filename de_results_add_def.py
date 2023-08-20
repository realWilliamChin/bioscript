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
    args = parser.parse_args()
    return args.kegg, args.nr, args.swiss


def process_deresults(de_results_file, kegg_file, nr_file, swiss_file):
    de_df = pd.read_csv(de_results_file, sep='\t')
    de_matrix_df = pd.read_csv(de_results_file + '_readCounts.matrix', sep='\t')
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
    de_df = add_kns_def(de_df, kegg_file, nr_file, swiss_file)
    de_df.to_csv(de_results_file.replace('DE_results', 'DEG_data.txt'), sep='\t', index=False)


def main():
    file = parse_input()
    kegg_file = file[0]
    nr_file = file[1]
    swiss_file = file[2]
    # kegg_file = '../def/Medicago_sativa_KEGG_gene_def.txt'
    # nr_file = '../def/Medicago_sativa_ZhongmuNo1_cds_rename_nr_diamond_nr_gene_def.txt'
    # swiss_file = '../def/Medicago_sativa_ZhongmuNo1_cds_rename_swiss_gene_def.txt'

    for de_results_file in os.listdir():
        if de_results_file.endswith('DE_results'):
            process_deresults(de_results_file, kegg_file, nr_file, swiss_file)


if __name__ == '__main__':
    main()


