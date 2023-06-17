# -*- coding: UTF-8 -*-
# Created Time  : 2023/4/26 9:44
# Author        : WilliamGoGo
import os
import pandas as pd


def format_tbl(tbl_file):
    """
    格式化 hmmscan 的 tbl 文件
    """
    tbl_names = ['target_name', 'accession', 'query_name', 'accession', 'fs_E-value', 'fs_score', 'fs_bias', 'b1d_E-value', 'b1d_score', 'b1d_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target']
    tbl_lst = []
    with open(tbl_file, 'r') as f:
        for line in f.readlines():
            if '#' in line:
                continue
            lst = line.strip().split()
            lst = lst[:18] + [' '.join(lst[18:])]
            tbl_lst.append(lst)
    tbl_df = pd.DataFrame(tbl_lst, columns=tbl_names)

    return tbl_df


def parse_input():
    """
    解析输入参数
    """
    import argparse
    parser = argparse.ArgumentParser(description='Format hmmscan tbl file, and output a txt file. '
                                                 'Not input file name, default all tbl file in current dir'
                                                 'File must with tbl suffix')
    parser.add_argument('-i', '--input', type=str, help='hmmscan tbl file')
    args = parser.parse_args()
    if args.input:
        file_name = args.input
    else:
        file_name = 'all'
    return file_name


if __name__ == '__main__':
    tbl_file = parse_input()
    if tbl_file == 'all':
        tbl_lst = [x for x in os.listdir() if 'tbl' in x]
        for tbl_file in tbl_lst:
            tbl = format_tbl(tbl_file)
            tbl.to_csv(tbl_file.replace('.tbl', '_tbl.txt'), sep='\t', index=False)
    else:
        tbl = format_tbl(tbl_file)
        tbl.to_csv(tbl_file.replace('.tbl', '_tbl.txt'), sep='\t', index=False)

