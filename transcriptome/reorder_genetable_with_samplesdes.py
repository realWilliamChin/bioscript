#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Created Time  : 3/6/2023 3:25 PM
# Author        : WilliamGoGo
import argparse
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser(description='根据 samples_described.txt 对表文件列重新排序')
    parser.add_argument('-s', '--sample', default='samples_described.txt',
                        help='指定 samples_described.txt 文件, 默认当前文件夹下的 samples_described.txt')
    parser.add_argument('-f', '--file', required=True, help='输入需要重新对列排序的文件')
    parser.add_argument('-r', '--isreplace', action='store_true', help='是否对源文件替换，不替换则会重新生成 _realign.txt')
    args = parser.parse_args()
    return args


def check_columns_name(lst, df_columns):
    error_lst = []
    for i in df_columns:
        if i == df_columns[0]:
            continue
        if i not in lst:
            error_lst.append(i)
    if error_lst:
        message = f"{', '.join(error_lst)} 未包含，将忽略"
        print(message)
    
    
def reindex(lst, file, replace_or_not):
    df = pd.read_csv(file, sep='\t', dtype=str)
    lst.insert(0, df.columns[0])
    if replace_or_not:
        df.reindex(columns=lst).to_csv(file, sep='\t', index=False)
    else:
        df.reindex(columns=lst).to_csv(file.replace('.txt', '_realign.txt'), sep='\t', index=False)


def main():
    args = parse_input()
    
    sample_arr = list(pd.read_csv(args.sample, sep='\t')['sample'])
    check_columns_name(sample_arr, pd.read_csv(args.file, sep='\t').columns.values)
    reindex(sample_arr, args.file, args.isreplace)

    print('\nDone!\n')


if __name__ == '__main__':
    main()
