# -*- coding: UTF-8 -*-
# Created Time  : 3/6/2023 3:25 PM
# Author        : WilliamGoGo
# 根据 samples 对一些文件列排序，主要针对 fpkm_matrix_filterd.txt 和 reads_matrix_filterd.txt 优化
import argparse
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser(description='根据 samples 对一些文件列排序，主要针对 fpkm_matrix_filterd.txt 和 reads_matrix_filterd.txt 优化')
    parser.add_argument('-s', '--sample', help='指定 samples_described.txt 文件', required=True)
    parser.add_argument('-f', '--file', required=True, help='输入需要重新对列排序的文件')
    parser.add_argument('-r', '--replacefile', action='store_true', help='是否对源文件替换，不替换则会重新生成 _realign.txt')
    args = parser.parse_args()
    return args.sample, args.file, args.replacefile


def check_columns_name(lst, df_columns):
    for i in lst:
        if i == df_columns[0]:
            continue
        if not i in df_columns:
            print(f'{i} not in {df_columns}')
            exit(1)
    
    
def reindex(lst, file, replace_or_not):
    df = pd.read_csv(file, sep='\t')
    lst.insert(0, df.columns[0])
    if replace_or_not:
        df.reindex(columns=lst).to_csv(file, sep='\t', index=False)
    else:
        df.reindex(columns=lst).to_csv(file.replace('.txt', '_realign.txt'), sep='\t', index=False)


def main():
    parse_input = parse_input()
    samples_file = parse_input[0]
    realign_file = parse_input[1]
    replace_ornot = parse_input[2]
    
    sample_arr = list(pd.read_csv(samples_file, sep='\t')['sample'])
    sample_arr
    check_columns_name(sample_arr, pd.read_csv(realign_file, sep='\t').columns)
    reindex(sample_arr, realign_file, replace_ornot)


if __name__ == '__main__':
    main()
