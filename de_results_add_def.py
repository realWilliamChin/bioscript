# -*- coding: UTF-8 -*-
# Created Time  : 2023/4/28 14:46
# Author        : WilliamGoGo
import os
import pandas as pd


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


def add_kns_def(df, kegg_file, nr_file, swiss_file):
    """
    添加 kegg nr swiss def
    """
    kegg_df = pd.read_csv(kegg_file, sep='\t', skiprows=1, names=['GeneID', 'KO_Number', 'KEGG_Def'])
    kegg_df['KEGG_Def'] = kegg_df['KO_Number'] + '::' + kegg_df['KEGG_Def']
    kegg_df.drop(columns=['KO_Number'], inplace=True)
    nr_df = pd.read_csv(nr_file, sep='\t', skiprows=1, usecols=[0, 2], names=['GeneID', 'NR_Def'])
    swiss_df = pd.read_csv(swiss_file, sep='\t', skiprows=1, usecols=[0, 2], names=['GeneID', 'Swiss_Def'])
    result_df = pd.merge(left=df, right=kegg_df, on='GeneID', how='left')
    result_df = pd.merge(left=result_df, right=nr_df, on='GeneID', how='left')
    result_df = pd.merge(left=result_df, right=swiss_df, on='GeneID', how='left')
    return result_df


def process_deresults():
    file = parse_input()
    kegg_file = file[0]
    nr_file = file[1]
    swiss_file = file[2]
    # kegg_file = '../def/Medicago_sativa_KEGG_gene_def.txt'
    # nr_file = '../def/Medicago_sativa_ZhongmuNo1_cds_rename_nr_diamond_nr_gene_def.txt'
    # swiss_file = '../def/Medicago_sativa_ZhongmuNo1_cds_rename_swiss_gene_def.txt'

    for de_results_file in os.listdir():
        if de_results_file.endswith('DE_results'):
            de_df = pd.read_csv(de_results_file, sep='\t')
            de_matrix_df = pd.read_csv(de_results_file + '_readCounts.matrix', sep='\t')
            de_df = pd.merge(left=de_df, right=de_matrix_df, on='GeneID', how='left')
            de_df = add_kns_def(de_df, kegg_file, nr_file, swiss_file)
            de_df.to_csv(de_results_file + '_def.txt', sep='\t', index=False)


if __name__ == '__main__':
    process_deresults()
