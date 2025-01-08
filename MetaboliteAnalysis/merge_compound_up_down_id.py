#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/12/30 15:50
# Author        : William GoGo
"""
代码用来对，代谢目录 05_Enrich 目录，所有的 Up 和
Donw compound ID 合并成 pathview 代码运行的格式
"""
import os, sys
import argparse
import pandas as pd


def parse_input():
    args = argparse.ArgumentParser()
    args.add_argument('-i', '--input', help='输入文件夹')
    args.add_argument('-o', '--output', help='输出合并后文件的目录')

    return args.parse_args()


def main():
    args = parse_input()
    input_dir, output_dir = args.input, args.output
    
    vs_group_list = [x for x in os.listdir(input_dir)]
    
    for vs in vs_group_list:
        up_id_file = os.path.join(input_dir, vs, f'{vs}_Up_Compound_ID.txt')
        down_id_file = os.path.join(input_dir, vs, f'{vs}_Down_Compound_ID.txt')
        
        up_id_df = pd.read_csv(up_id_file, sep='\t', header=None, names=['Compound_ID'])
        down_id_df = pd.read_csv(down_id_file, sep='\t', header=None, names=['Compound_ID'])
        
        up_id_df['regulation'] = 1
        down_id_df['regulation'] = -1
        
        os.makedirs(os.path.join(output_dir, vs), exist_ok=True)
        result_file = os.path.join(output_dir, vs, f'{vs}_compound_id_regulation.txt')
        result_df = pd.concat([up_id_df, down_id_df])
        result_df.to_csv(result_file, sep='\t', index=False, header=False)
        

if __name__ == '__main__':
    main()