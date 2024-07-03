#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/07/03 10:49
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import subprocess
from loguru import logger


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='input_file', help='组间分析中的 Differential_metabolite_count_summary.xlsx 文件')
    parser.add_argument(dest='comp_file', help='两组之间画 venn 图的比较文件')
    
    return parser.parse_args()


def draw_venn_graph(input_file, output_file):
    cmd = f'Rscript /home/colddata/qinqiang/script/Plot/Venn/venn2.r {input_file} {output_file}'
    rep = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if rep != 0:
        logger.error(f'运行 {cmd} 错误')
        logger.error(f'标准错误：{rep.stderr.decode()}')
        logger.error(f'标准输出：{rep.stdout.decode()}')
    else:
        logger.success(f'运行 {cmd} 成功')
    

def main():
    args = parse_input()
    
    df = pd.read_excel(args.input_file, engine='openpyxl')
    df = df.set_index('group')
    comp_df = pd.read_csv(args.comp_file, sep='\t', header=None, names=['Group1', 'Group2'])
    for row in comp_df.itertuples():
        group1 = row.Group1
        group2 = row.Group2
        group1_up_list = df.loc[group1]['Up_list'].split(',')
        group2_up_list = df.loc[group2]['Up_list'].split(',')
        
        group1_down_list = df.loc[group1]['Down_list'].split(',')
        group2_down_list = df.loc[group2]['Down_list'].split(',')
        
        Up_venn_df = pd.DataFrame.from_dict({group1: group1_up_list, group2: group2_up_list}, orient='index').transpose()
        Down_venn_df = pd.DataFrame.from_dict({group1: group1_down_list, group2: group2_down_list}, orient='index').transpose()
        
        Up_venn_df.to_csv(f'{group1}_and_{group2}_Up_venn.txt', sep='\t', index=False)
        Down_venn_df.to_csv(f'{group1}_and_{group2}_Down_venn.txt', sep='\t', index=False)
        
        draw_venn_graph(f'{group1}_and_{group2}_Up_venn.txt', f'{group1}_and_{group2}_Up_venn.jpeg')
        draw_venn_graph(f'{group1}_and_{group2}_Down_venn.txt', f'{group1}_and_{group2}_Down_venn.jpeg')
        


if __name__ == '__main__':
    main()