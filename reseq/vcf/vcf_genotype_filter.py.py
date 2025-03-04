#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/12/25 15:34
# Author        : William GoGo
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    argparser = argparse.ArgumentParser()
    argparser.add_argument(dest='input', help='输入文件')
    argparser.add_argument(dest='output', help='输出文件')
    argparser.add_argument('-m', '--male', help='父本列名')
    argparser.add_argument('-f', '--female', help='母本列名')
    return argparser.parse_args()


def skip_rows(input_file, skip_str):
    # vcf 文件读取前处理
    skip_rows = 0
    with open(input_file, "r") as file:
        for line in file:
            if line.startswith(skip_str):
                skip_rows += 1
            else:
                break
    logger.info(f'跳过 {input_file} 的前 {skip_rows} 行')
    return skip_rows



def filter_parent(vcf_df, male_parent, female_parent):
    mask = (
        # 排除父本和母本同为纯合的情况（包括 0/0, 0|0, 1/1, 1|1 组合）
        ~(
            vcf_df[male_parent].str.split(':').str[0].str.match(r'^0[/|]0$') & 
            vcf_df[female_parent].str.split(':').str[0].str.match(r'^0[/|]0$')
        ) 
        & 
        ~(
            vcf_df[male_parent].str.split(':').str[0].str.match(r'^1[/|]1$') & 
            vcf_df[female_parent].str.split(':').str[0].str.match(r'^1[/|]1$')
        )
    )
    
    filtered_df = vcf_df[mask]
    droped_df = vcf_df[~mask]
    
    return filtered_df, droped_df


def main():
    args = parse_input()
    in_file = args.input
    skiprows = skip_rows(in_file, '##')
                
    df = pd.read_csv(in_file, sep='\t', skiprows=skiprows, low_memory=False)
    
    filter_before = df.shape[0]
    filtered_df, droped_df = filter_parent(df, args.male, args.female)
    filter_after = filtered_df.shape[0]
    logger.info(f'过滤之前 {filter_before} 行，过滤之后 {filter_after} 行')
    
    # 去掉父本母本
    filtered_df.drop(columns=[args.male, args.female], inplace=True)
    
    # 读取原始文件的注释行
    with open(in_file, 'r') as f:
        header_lines = [line for line in f if line.startswith('##')]

    # 保存过滤后的结果（包含注释行）
    with open(args.output, 'w') as f:
        f.writelines(header_lines)
        filtered_df.to_csv(f, sep='\t', index=False, header=True)   
    
    droped_df.to_csv('droped.txt', sep='\t', index=False)


if __name__ == '__main__':
    main()