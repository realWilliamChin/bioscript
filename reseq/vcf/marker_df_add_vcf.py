#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/29 15:34
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    description = """
    该脚本用于将 VCF 文件中的信息合并到包含 Marker 列的输入文件中。
    
    功能说明：
    1. 读取包含 Marker 列的输入文件（数据框格式）
    2. 读取 VCF 文件，自动跳过以 '##' 开头的注释行
    3. 从 VCF 文件的 #CHROM 和 POS 列生成 Marker（格式：CHROM_POS）
    4. 基于 Marker 列将 VCF 文件的信息左连接（left join）到输入文件中
    5. 自动处理重复列，保留原始输入文件的列，添加 VCF 文件中的新列
    6. 输出合并后的结果文件
    """
    argparser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    argparser.add_argument('-i', '--input', help='输入文件，包含 Marker')
    argparser.add_argument('-v', '--vcf', help='vcf 文件')
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
    vcf_file = args.vcf
    skiprows = skip_rows(vcf_file, '##')

    df = load_table(in_file)
    vcf_df = pd.read_csv(vcf_file, sep='\t', skiprows=skiprows, low_memory=False)
    vcf_df['Marker'] = vcf_df['#CHROM'].astype(str) + '_' + vcf_df['POS'].astype(str)
    
    result_df = pd.merge(df, vcf_df, how='left', on='Marker', suffixes=('', '_dup'))
    # 去掉被 merge 出来的重复列，仅保留前面的
    dup_cols = [col for col in result_df.columns if col.endswith('_dup')]
    result_df = result_df.drop(columns=dup_cols)
    write_output_df(result_df, args.output, index=False)


if __name__ == '__main__':
    main()