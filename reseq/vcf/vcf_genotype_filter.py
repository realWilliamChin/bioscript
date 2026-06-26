#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/12/25 15:34
# Author        : William GoGo
#
# 【项目模板脚本】本脚本中的样本名/亲本组合等过滤逻辑是按具体项目需求硬编码的，
# 每次新项目使用时需要按需修改内部的样本列表与基因型组合规则，不是一个通用 CLI 工具。
import argparse
import pandas as pd
from loguru import logger
import pysam


def parse_input():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-i', '--input', help='输入文件')
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



def filter_vcf(vcf_df, fertile_20008, fertile_20515, fertile_20549, fertile_20594, fertile_20705, fertile_20879, sterile_285):
    fertile_20008_gt = vcf_df[fertile_20008].str.split(':').str[0]
    fertile_20515_gt = vcf_df[fertile_20515].str.split(':').str[0]
    fertile_20549_gt = vcf_df[fertile_20549].str.split(':').str[0]
    fertile_20594_gt = vcf_df[fertile_20594].str.split(':').str[0]
    fertile_20705_gt = vcf_df[fertile_20705].str.split(':').str[0]
    fertile_20879_gt = vcf_df[fertile_20879].str.split(':').str[0]
    sterile_285_gt = vcf_df[sterile_285].str.split(':').str[0]
    
    hom_ref = r'^0[/|]0$'      # 纯合参考型 (0/0 或 0|0)
    hom_alt = r'^1[/|]1$'      # 纯合变异型 (1/1 或 1|1)
    hom_het = r'^0[/|]1$'      # 杂合型 (0/1 或 0|1)
    hom_het2 = r'^1[/|]0$'      # 杂合型 (1/0 或 1|0)
    missing = r'^\.[/|]\.$'    # 缺失数据 (./. 或 .|.)
    
    # 通用纯合变异 (支持任意 ALT 索引)
    hom_alt_general = r'^[1-9]\d*[/|][1-9]\d*$'  # 如 3/3, 5|5
    # 通用杂合 (支持任意 ALT 组合)
    het_general = r'^[0-9]\d*[/|][0-9]\d*$'      # 如 0/3, 2|4  
    
    # 找出 fertile* 共有的（必须两个一模一样），sterile_285 没有的
    mask = (
        (fertile_20008_gt.str.match(hom_ref) & fertile_20515_gt.str.match(hom_ref) & fertile_20549_gt.str.match(hom_ref) & fertile_20594_gt.str.match(hom_ref) & fertile_20705_gt.str.match(hom_ref) & fertile_20879_gt.str.match(hom_ref)) |
        (fertile_20008_gt.str.match(hom_het) & fertile_20515_gt.str.match(hom_het) & fertile_20549_gt.str.match(hom_het) & fertile_20594_gt.str.match(hom_het) & fertile_20705_gt.str.match(hom_het) & fertile_20879_gt.str.match(hom_het)) |
        (fertile_20008_gt.str.match(hom_het2) & fertile_20515_gt.str.match(hom_het2) & fertile_20549_gt.str.match(hom_het2) & fertile_20594_gt.str.match(hom_het2) & fertile_20705_gt.str.match(hom_het2) & fertile_20879_gt.str.match(hom_het2)) |
        (fertile_20008_gt.str.match(hom_alt) & fertile_20515_gt.str.match(hom_alt) & fertile_20549_gt.str.match(hom_alt) & fertile_20594_gt.str.match(hom_alt) & fertile_20705_gt.str.match(hom_alt) & fertile_20879_gt.str.match(hom_alt))
    )
    # # sterile_285 和其他都不一样的
    # mask_sterile = (
    #     ~mask & ~sterile_285_gt.str.match(hom_ref) & ~sterile_285_gt.str.match(hom_het) & ~sterile_285_gt.str.match(hom_het2) & ~sterile_285_gt.str.match(hom_alt)
    # )
    
    return vcf_df[mask], vcf_df[~mask]


def main():
    args = parse_input()
    in_file = args.input
    skiprows = skip_rows(in_file, '##')
                
    df = pd.read_csv(in_file, sep='\t', skiprows=skiprows, low_memory=False)
    
    filter_before = df.shape[0]
    filtered_df, droped_df = filter_vcf(df, 'fertile_20008', 'fertile_20515', 'fertile_20549', 'fertile_20594', 'fertile_20705', 'fertile_20879', 'sterile_285')
    filter_after = filtered_df.shape[0]
    logger.info(f'过滤之前 {filter_before} 行，过滤之后 {filter_after} 行')
    
    # 去掉父本母本
    # filtered_df.drop(columns=[args.male, args.female], inplace=True)
    
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