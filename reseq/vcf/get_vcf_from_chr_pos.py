#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/29 15:34
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger
import pysam

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-i', '--input', help='输入文件，第一列 CHROM, 第二列 POS')
    argparser.add_argument('-v', '--vcf', help='参考文件，被提取的 vcf 文件')
    argparser.add_argument('-o', '--output', help='输出文件')
    
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


def main():
    args = parse_input()
    skiprows = skip_rows(args.vcf, '##')   
    vcf_df = pd.read_csv(args.vcf, sep='\t', skiprows=skiprows, low_memory=False)
    df = load_table(args.input, usecols=['Marker'])
    vcf_df['Marker'] = vcf_df['#CHROM'].astype(str) + '_' + vcf_df['POS'].astype(str)
    
    filtered_df = pd.merge(df, vcf_df, how='left', on='Marker')
    filtered_df = filtered_df.drop(columns=['Marker'])
    # 读取原始文件的注释行
    with open(args.vcf, 'r') as f:
        header_lines = [line for line in f if line.startswith('##')]

    # 保存过滤后的结果（包含注释行）
    with open(args.output, 'w') as f:
        f.writelines(header_lines)
        filtered_df.to_csv(f, sep='\t', index=False, header=True)   


if __name__ == '__main__':
    main()