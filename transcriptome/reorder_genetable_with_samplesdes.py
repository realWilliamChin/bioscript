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
    parser = argparse.ArgumentParser(description='根据 samples_described.txt 对表文件列重新排序，支持通过 old_sample 列重命名矩阵列')
    parser.add_argument('-s', '--sample', default='samples_described.txt',
                        help='指定 samples_described.txt 文件, 默认当前文件夹下的 samples_described.txt')
    parser.add_argument('-f', '--file', required=True,
                        help='输入需要重新对列排序的文件, 可以针对 csv 和 txt 和 xlsx 格式文件进行处理')
    parser.add_argument('-o', '--output', default='output.txt',
                        help='生成的文件名，默认在当前文件夹下生成 output.txt, 也可以指定其他格式, csv 或者 xlsx')
    parser.add_argument('--use-old-sample', action='store_true', default=False,
                        help='是否使用 old_sample 列重命名矩阵列，重命名后再按 sample 列排序')
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
    
    # 读取样本描述表
    sample_df = pd.read_csv(args.sample, sep='\t', dtype=str)
    
    # 检查 sample 列是否存在
    if 'sample' not in sample_df.columns:
        print("样本描述表缺少必要的 'sample' 列")
        exit(1)
    
    sample_arr = list(sample_df['sample'])
    
    # 如果使用 old_sample 重命名列
    if args.use_old_sample:
        # 检查 old_sample 列是否存在
        if 'old_sample' not in sample_df.columns:
            print("样本描述表缺少 'old_sample' 列，无法使用 --use-old-sample 参数")
            exit(1)
        
        # 建立 old_sample 到 sample 的映射
        sample_map = dict(zip(sample_df['old_sample'], sample_df['sample']))
        print(f"加载 {len(sample_map)} 个样本映射关系")
        
        # 检查重复的 old_sample
        duplicated_old = sample_df[sample_df['old_sample'].duplicated(keep=False)]
        if not duplicated_old.empty:
            print(f"警告：样本描述表中存在重复的old_sample: {', '.join(duplicated_old['old_sample'].unique())}，将使用第一个出现的映射关系")
        
        # 加载矩阵文件
        df = load_table(args.file, dtype=str)
        print(f"原始矩阵包含 {df.shape[0]} 行, {df.shape[1]} 列")
        
        # 找出匹配的 old_sample 列
        existing_old_samples = [col for col in sample_map.keys() if col in df.columns]
        missing_old_samples = [col for col in sample_map.keys() if col not in df.columns]
        
        if missing_old_samples:
            print(f"警告：矩阵中不存在以下old_sample列: {', '.join(missing_old_samples)}")
        
        if not existing_old_samples:
            print("错误：矩阵中没有找到任何匹配的old_sample列，无法继续")
            exit(1)
        
        print(f"在矩阵中找到 {len(existing_old_samples)} 个匹配的old_sample列，将进行重命名")
        
        # 重命名列
        df = df.rename(columns=sample_map)
        
        # 将重命名后的 df 保存到临时文件，供 reindex 函数读取
        import tempfile
        tmp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
        tmp_path = tmp_file.name
        write_output_df(df, tmp_path, index=False)
        tmp_file.close()
        
        # 使用临时文件进行排序
        reindex(sample_arr, tmp_path, args.output)
        
        # 删除临时文件
        os.unlink(tmp_path)
    else:
        # 原有逻辑，直接排序
        reindex(sample_arr, args.file, args.output)

    print('\nDone!\n')


if __name__ == '__main__':
    main()
