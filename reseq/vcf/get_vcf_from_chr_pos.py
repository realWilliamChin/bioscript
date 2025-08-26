#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/29 15:34
# Author        : William GoGo
import argparse
import pandas as pd
from loguru import logger
import pysam


def parse_input():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-i', '--input', help='输入文件')
    argparser.add_argument('-r', '--ref', help='参考文件，没有表头，第一列chrom，第二列pos')
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
    in_file = args.input
    skiprows = skip_rows(in_file, '##')   
    vcf_df = pd.read_csv(in_file, sep='\t', skiprows=skiprows, low_memory=False)
    ref_df = pd.read_csv(args.ref, sep='\t', header=None, usecols=[0, 1], names=['#CHROM', 'POS'])
    vcf_df['chr_pos'] = vcf_df['#CHROM'].astype(str) + '_' + vcf_df['POS'].astype(str)
    ref_df['chr_pos'] = ref_df['#CHROM'].astype(str) + '_' + ref_df['POS'].astype(str)
    ref_df = ref_df.drop(columns=['#CHROM', 'POS'])
    
    filtered_df = pd.merge(ref_df, vcf_df, how='left', on='chr_pos')
    filtered_df = filtered_df.drop(columns=['chr_pos'])

    # 读取原始文件的注释行
    with open(in_file, 'r') as f:
        header_lines = [line for line in f if line.startswith('##')]

    # 保存过滤后的结果（包含注释行）
    with open(args.output, 'w') as f:
        f.writelines(header_lines)
        filtered_df.to_csv(f, sep='\t', index=False, header=True)   


if __name__ == '__main__':
    main()