#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/07/14 10:50
# Author        : William GoGo
import argparse
import pandas as pd


def parse_input():
    argparser = argparse.ArgumentParser(description='Get specific sample from vcf file')
    argparser.add_argument('-v', '--vcf', help='input vcf file')
    argparser.add_argument('-l', '--sample_list', help='sample list file')
    argparser.add_argument('-o', '--output', default='specific_sample.vcf', help='output vcf file')
    return argparser.parse_args()


def reorder(vcf, sample_lst, output_file):
    # read headers startwith '#' and not startwith '##' as a list
    with open(vcf, 'r') as f, open(output_file, 'a') as out_f:
        column_name_lst = []
        for line in f:
            if line[:2] == '##':
                out_f.write(line)
            elif line[:1] == '#':
                column_name_lst = line.strip().split('\t')
                break
    vcf_df = pd.read_csv(vcf, sep='\t', comment='#', header=None, names=column_name_lst, low_memory=False)
    vcf_sample_lst = vcf_df.columns.tolist()[:9] + sample_lst
    result_df = vcf_df[vcf_sample_lst]
    result_df.to_csv(output_file, sep='\t', index=False, mode='a')


def main():
    args = parse_input()
    with open(args.sample_list, 'r') as f:
        sample_lst = f.read().strip().split('\n')
    reorder(args.vcf, sample_lst, args.output)


if __name__ == '__main__':
    main()
    