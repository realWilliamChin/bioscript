#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Created Time  : 3/6/2023 3:25 PM
# Author        : WilliamGoGo
import os, sys
import argparse
import pandas as pd

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


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
    
    
def convert_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    遍历所有列：若某列中所有非空值均可成功转为数值，则将该列转换为数值类型。
    空字符串会被视为缺失值允许保留为 NaN。
    """
    df_converted = df.copy()
    for col in df_converted.columns:
        # 先将值标准化为字符串，去除首尾空白，并将空字符串视为缺失
        ser = df_converted[col].astype(str).str.strip()
        ser = ser.replace({"": None})
        # 尝试数值转换（不可转的变为 NaN）
        converted = pd.to_numeric(ser, errors='coerce')
        non_empty_mask = ser.notna()
        # 所有非空项均可成功转为数值则接受转换
        if converted[non_empty_mask].notna().all():
            df_converted[col] = converted
    return df_converted


def reindex(lst, input_file, output_file):
    df = load_table(input_file, dtype=str)
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

    # 在输出前将可转换为数值的列转换为数值类型
    df = convert_numeric_columns(df)
    write_output_df(df, output_file, index=False)


def main():
    args = parse_input()
    
    sample_arr = list(pd.read_csv(args.sample, sep='\t')['sample'])
    reindex(sample_arr, args.file, args.output)

    print('\nDone!\n')


if __name__ == '__main__':
    main()
