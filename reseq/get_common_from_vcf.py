#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/07 14:19
# Author        : William GoGo
import argparse
import pandas as pd
import os
from common import skip_rows


def parse_input():
    argparser = argparse.ArgumentParser(description='')
    argparser.add_argument('-i', '--input', help='输入文件')
    argparser.add_argument('-o', '--output', help='输出文件')
    return argparser.parse_args()


def drop_nonvarin_line(df):
    for samples in df.columns.tolist()[9:]:
        df = df[df[samples].str.split(':').str[0] != '0|0']
        df = df[df[samples].str.split(':').str[0] != '0/0']
        df = df[df[samples].str.split(':').str[0] != './.']
    return df
    

def merge_from_bcftools():
    file_lst = [x for x in os.listdir() if x.endswith('.vcf')]
    skiprows = skip_rows(file_lst[0], '##')
    df = pd.read_csv(file_lst[0], sep='\t', skiprows=skiprows, low_memory=False)
    for file in file_lst[1:]:
        skiprows = skip_rows(file, '##')
        df_2 = pd.read_csv(file, sep='\t', skiprows=skiprows, low_memory=False)
        df = pd.merge(left=df, right=df_2, on=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], how='outer')

    return df


def diffset():
    left_file = 'BSA_Sg-Z_WT_INDELs.vcf'
    right_file = '0000.vcf'
    skiprows = skip_rows(left_file, '##')
    left_df = pd.read_csv(left_file, sep='\t', skiprows=skiprows, low_memory=False)
    skiprows = skip_rows(right_file, '##')
    right_df = pd.read_csv(right_file, sep='\t', skiprows=skiprows, low_memory=False, usecols=[0, 1, 9])
    df = pd.merge(left=left_df, right=right_df, on=['#CHROM', 'POS'], how='inner')
    df.to_csv('diffset.vcf', sep='\t', index=False)
    

def main():
    args = parse_input()
    file = args.input
    skiprows = skip_rows(file, '##')
    df = pd.read_csv(file, sep='\t', skiprows=skiprows, low_memory=False)
    df = drop_nonvarin_line(df)
    df.to_csv(args.output, sep='\t', index=False)
    # merge_from_bcftools().to_csv('merge.vcf', sep='\t', index=False)
    # diffset()


if __name__ == '__main__':
    main()