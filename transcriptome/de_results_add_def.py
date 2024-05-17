#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/4/28 14:46
# Author        : WilliamGoGo
import os, sys
import pandas as pd
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/transcriptome/'))
from genedf_add_knsdef import add_kns_def


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
    de_reads_df = pd.read_csv(de_results_file + '_readCounts.matrix', sep='\t', dtype={"GeneID": str})
    de_reads_df.columns = ['GeneID'] + [x + '_raw_reads' for x in de_reads_df.columns.tolist() if x.lower() != 'geneid']
    de_df = pd.merge(left=de_df, right=de_reads_df, on='GeneID', how='left')
    # 排序 down，up，NOsig。down 的 FC 值从小到大，up 的 FC 值从大到小
    # 先分三份，再合并
    de_df_down = de_df[de_df['regulation'] == 'Down'].copy()
    de_df_down.sort_values(by='FC', ascending=True, inplace=True)
    de_df_up = de_df[de_df['regulation'] == 'Up'].copy()
    de_df_up.sort_values(by='FC', ascending=False, inplace=True)
    de_df_nosig = de_df[de_df['regulation'] == 'NoSignificant'].copy()
    # 合并
    de_df = pd.concat([de_df_down, de_df_up, de_df_nosig])
    if kns_file or kegg_file or nr_file or swiss_file:
        # 添加注释
        de_df = add_kns_def(de_df, kegg_file, nr_file, swiss_file, kns_file)
    de_df.to_csv(de_results_file.replace('DE_results', 'DEG_data.txt'), sep='\t', index=False)


def main():
    args = parse_input()

    for de_results_file in os.listdir():
        if de_results_file.endswith('DE_results'):
            logger.info(f'processing---{de_results_file}')
            process_deresults(de_results_file, args.kegg, args.nr, args.swiss, args.kns)
    if args.kns or args.swiss or args.kegg or args.nr:
        for up_down_id_file in os.listdir("DEG_analysis_results"):
            up_down_id_file = os.path.join("DEG_analysis_results", up_down_id_file)
            if up_down_id_file.endswith('Down_ID.txt') or up_down_id_file.endswith('Up_ID.txt'):
                up_down_id_df = pd.read_csv(up_down_id_file, sep='\t', names=['GeneID'], dtype={"GeneID": str})
                result_df = add_kns_def(up_down_id_df, args.kegg, args.nr, args.swiss, args.kns)
                result_df.to_csv(up_down_id_file.replace('.txt', '_def.txt'), sep='\t', index=False)

    logger.success('Done!')


if __name__ == '__main__':
    main()


