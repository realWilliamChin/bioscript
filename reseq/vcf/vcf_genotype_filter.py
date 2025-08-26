#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/12/25 15:34
# Author        : William GoGo
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



def filter_vcf(vcf_df, son, father, mother):
    son_gt = vcf_df[son].str.split(':').str[0]
    father_gt = vcf_df[father].str.split(':').str[0]
    mother_gt = vcf_df[mother].str.split(':').str[0]
    
    hom_ref = r'^0[/|]0$'      # 纯合参考型 (0/0 或 0|0)
    hom_alt = r'^1[/|]1$'      # 纯合变异型 (1/1 或 1|1)
    missing = r'^\.[/|]\.$'    # 缺失数据 (./. 或 .|.)
    
    # 通用纯合变异 (支持任意 ALT 索引)
    hom_alt_general = r'^[1-9]\d*[/|][1-9]\d*$'  # 如 3/3, 5|5
    # 通用杂合 (支持任意 ALT 组合)
    het_general = r'^[0-9]\d*[/|][0-9]\d*$'      # 如 0/3, 2|4  
    
    # son_hom = son_gt.str.match(hom_ref) | son_gt.str.match(hom_alt)
    # father_hom = father_gt.str.match(hom_ref) | father_gt.str.match(hom_alt)
    # mother_hom = mother_gt.str.match(hom_ref) | mother_gt.str.match(hom_alt)
    # father_missing = father_gt.str.match(missing)
    # mother_missing = mother_gt.str.match(missing)
    
    # mask = (
    #     # 儿子是纯合
    #     son_hom &
    #     # 父母都不缺失
    #     ~father_missing & 
    #     ~mother_missing &
    #     # 父亲或母亲中至少有一个不是纯合
    #     (~father_hom | ~mother_hom)
    # )
    
    son_only_mask = (
        ~(son_gt.str.match(hom_ref) & father_gt.str.match(hom_ref))
        & ~(son_gt.str.match(hom_ref) & mother_gt.str.match(hom_ref))
        & ~(son_gt.str.match(hom_alt) & father_gt.str.match(hom_alt))
        & ~(son_gt.str.match(hom_alt) & mother_gt.str.match(hom_alt))
        & ~(son_gt.str.match(missing) & father_gt.str.match(missing))
        & ~(son_gt.str.match(missing) & mother_gt.str.match(missing))
        & ~(son_gt.str.match(r'^0[/|]1$') & father_gt.str.match(r'^0[/|]1$'))
        & ~(son_gt.str.match(r'^0[/|]1$') & mother_gt.str.match(r'^0[/|]1$'))
        & ~(son_gt.str.match(r'^1[/|]0$') & father_gt.str.match(r'^1[/|]0$'))
        & ~(son_gt.str.match(r'^1[/|]0$') & mother_gt.str.match(r'^1[/|]0$'))
    )

    son_mother_mask = (
        (son_gt.str.match(r'^1[/|]0$') & mother_gt.str.match(r'^1[/|]0$')) |
        (son_gt.str.match(r'^1[/|]0$') & mother_gt.str.match(r'^0[/|]1$')) |
        (son_gt.str.match(r'^0[/|]1$') & mother_gt.str.match(r'^0[/|]1$')) |
        (son_gt.str.match(r'^0[/|]1$') & mother_gt.str.match(r'^1[/|]0$'))
    )
    son_father_mask = (
        (son_gt.str.match(r'^1[/|]0$') & father_gt.str.match(r'^1[/|]0$')) |
        (son_gt.str.match(r'^1[/|]0$') & father_gt.str.match(r'^0[/|]1$')) |
        (son_gt.str.match(r'^0[/|]1$') & father_gt.str.match(r'^0[/|]1$')) |
        (son_gt.str.match(r'^0[/|]1$') & father_gt.str.match(r'^1[/|]0$'))
    )
    
    
    filtered_df = vcf_df[son_only_mask]
    droped_df = vcf_df[~son_only_mask]
    
    return filtered_df, droped_df


def main():
    args = parse_input()
    in_file = args.input
    skiprows = skip_rows(in_file, '##')
                
    df = pd.read_csv(in_file, sep='\t', skiprows=skiprows, low_memory=False)
    
    filter_before = df.shape[0]
    filtered_df, droped_df = filter_vcf(df, 'Child', 'Father', 'Mother')
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
    
    droped_df.to_csv('son_only_droped.txt', sep='\t', index=False)


if __name__ == '__main__':
    main()