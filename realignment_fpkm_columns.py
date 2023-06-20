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
    

def reindex(lst, file, replace_or_not):
    df = pd.read_csv(file, sep='\t')
    lst.insert(0, df.columns[0])
    if replace_or_not:
        df.reindex(columns=lst).to_csv(file, sep='\t', index=False)
    else:
        df.reindex(columns=lst).to_csv(file.replace('.txt', '_realign.txt'), sep='\t', index=False)


if __name__ == '__main__':
    parse_input = parse_input()
    samples_file = parse_input[0]
    realign_file = parse_input[1]
    replace_ornot = parse_input[2]
    
    sample_arr = list(pd.read_csv(samples_file, sep='\t')['sample'])
    sample_arr
    reindex(sample_arr, realign_file, replace_ornot)
