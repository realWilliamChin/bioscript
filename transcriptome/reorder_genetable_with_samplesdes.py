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
    parser.add_argument('-f', '--file', required=True,
                        help='输入需要重新对列排序的文件, 可以针对 csv 和 txt 和 xlsx 格式文件进行处理')
    parser.add_argument('-o', '--output', default='output.txt',
                        help='生成的文件名，默认在当前文件夹下生成 output.txt, 也可以指定其他格式, csv 或者 xlsx')
    args = parser.parse_args()
    return args
    
    
def reindex(lst, input_file, output_file):
    if input_file.endswith('.csv'):
        df = pd.read_csv(input_file, dtype=str)
    elif input_file.endswith('.xlsx'):
        df = pd.read_excel(input_file, dtype=str, engine='openpyxl')
    else:
        df = pd.read_csv(input_file, sep='\t', dtype=str)
    df_columns = df.columns.values[1:]
    error_lst = []
    
    # 合并检查和退出逻辑
    if any(i not in df.columns[1:] for i in lst): 
        print(f"数据表不包含此列：{', '.join(i for i in lst if i not in df.columns[1:])}")
        exit(1)

    # 合并打印错误列表逻辑
    error_lst = [i for i in df_columns if i not in lst]
    if error_lst: print(f"{', '.join(error_lst)} 未包含，将忽略")
        
    
    df = df.reindex(columns=[df.columns[0],] + lst)
    if len(set(df.columns) - set(lst)) - 1 == 0:
        print('输入的 title list 和输出的表 title 数量匹配正确')
    else:
        print('输入的 title list 和输出的表 title 数量不匹配，请检查')
        exit(1)

    if output_file.endswith('.csv'):
        df.to_csv(output_file, index=False)
    elif output_file.endswith('.xlsx'):
        df.to_excel(output_file, index=False, engine='openpyxl')
    else:
        df.to_csv(output_file, sep='\t', index=False)


def main():
    args = parse_input()
    
    sample_arr = list(pd.read_csv(args.sample, sep='\t')['sample'])
    reindex(sample_arr, args.file, args.output)

    print('\nDone!\n')


if __name__ == '__main__':
    main()
