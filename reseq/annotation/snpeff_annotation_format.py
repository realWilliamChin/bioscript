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
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-i', '--input', help='snpeff 的结果 vcf 文件')
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


def parse_info(vcf_df):
    # 解析 INFO 字段，找到 END 字段和 SVTYPE 字段，生成新列
    if 'INFO' not in vcf_df.columns:
        logger.warning('VCF文件缺少INFO列，无法解析 END/SVTYPE')
        vcf_df['END'] = None
        vcf_df['SVTYPE'] = None
        return vcf_df

    # 提取 END 字段
    vcf_df['END'] = vcf_df['INFO'].str.extract(r'(?:^|;| )END=([^;]+)')[0]
    # 提取 SVTYPE 字段
    vcf_df['SVTYPE'] = vcf_df['INFO'].str.extract(r'(?:^|;| )SVTYPE=([^;]+)')[0]

    return vcf_df


def parse_snpeff(vcf_df):
    # 取出来 EFF= 后面的 第一个 ( 前面的词作为突变类型
    vcf_df['Mutation_type'] = vcf_df['INFO'].str.extract(r';EFF=([^(]+)\(')[0]
    # 取第一个 ( 后第一个 | 前的内容作为 Significance
    # 例: ';EFF=missense_variant(MODERATE|MISSENSE|...' -> 'MODERATE'
    vcf_df['Significance'] = vcf_df['INFO'].str.extract(r';EFF=[^()]*\(([^|]*)\|')[0]
    return vcf_df


def parse_format(vcf_df):
    # 解析 VCF 的 FORMAT 字段，并将格式及每个样本列展开为新列
    if 'FORMAT' not in vcf_df.columns:
        logger.warning('VCF文件缺少FORMAT列，无法解析')
        return vcf_df

    # 假定标准 VCF 字段
    std_vcf_cols = ['#CHROM', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    # 找出样本名字列（即 FORMAT 后面所有列）
    std_col_idx = [i for i, c in enumerate(vcf_df.columns) if c in std_vcf_cols]
    max_idx = max(std_col_idx) if std_col_idx else (vcf_df.columns.get_loc('FORMAT') if 'FORMAT' in vcf_df.columns else 8)
    sample_cols = vcf_df.columns[(max_idx+1):]
    if len(sample_cols) == 0:
        logger.warning('VCF文件无样本列，无需解析FORMAT成新列')
        return vcf_df

    # 针对每一行，解析 FORMAT 及样本内容，并展开
    format_fields = vcf_df['FORMAT'].str.split(':')
    for sample in sample_cols:
        sample_values = vcf_df[sample].str.split(':')
        # 预存字典用于组装 (行号, 字段) -> 值
        out_dict = {}
        for idx, (fields, values) in enumerate(zip(format_fields, sample_values)):
            if isinstance(fields, list) and isinstance(values, list):
                for f, v in zip(fields, values):
                    out_dict.setdefault(f, {})[idx] = v
        # 新增列，列名: sample 名 + "_" + format 字段名
        for f in set().union(*format_fields.dropna()):
            v = pd.Series(out_dict.get(f, {}))
            v.index = vcf_df.index  # 对齐索引
            vcf_df[f"{sample}_{f}"] = v

    return vcf_df


def main():
    args = parse_input()
    skiprows = skip_rows(args.input, '##')
    vcf_df = pd.read_csv(args.input, sep='\t', skiprows=skiprows)
    vcf_df = parse_info(vcf_df)
    vcf_df = parse_snpeff(vcf_df)
    vcf_df = parse_format(vcf_df)

    write_output_df(vcf_df, args.output, index=False, sep='\t')


if __name__ == '__main__':
    main()