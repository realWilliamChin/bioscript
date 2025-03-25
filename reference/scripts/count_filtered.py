#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/03/11 16:17
# Author        : William GoGo
import os, sys
import pandas as pd
import numpy as np
import argparse
from loguru import logger
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools'))
from load_input import load_table

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/transcriptome'))
from reorder_genetable_with_samplesdes import reindex


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default='gene_count_matrix.csv',
                        help='输入 gene_count_matrix.csv 文件路径')
    parser.add_argument('-o', '--output', default='reads_matrix_filtered.txt',
                        help='输出过滤后的 count 的文件路径，默认 reads_matrix_filtered.txt')
    parser.add_argument('-n', type=int, default=50,
                        help='输入 counts 每行最低出现的数量阈值，把低于这个值的基因去掉')
    parser.add_argument('-s', default='samples_described.txt',
                        help='输入 sampels_described.txt 样本描述文件，用来对输出文件进行排序')
    parser.add_argument('--id', type=str, help='all_gene_id.txt 文件，输出只在 all_gene_id.txt 中存在的基因')

    return parser.parse_args()


def parse_raw_file(gene_count_matrix_file, parsed_raw_file):
    # os.system("sed 's/,/\t/g' gene_count_matrix.csv > gene_count_matrix.txt")
    gene_count_matrix_df = pd.read_csv(gene_count_matrix_file)
    gene_count_matrix_df.rename(columns={"gene_id": "GeneID"}, inplace=True)
    if gene_count_matrix_df['GeneID'].str.contains('|').mean() > 0.6:
        if gene_count_matrix_df['GeneID'].str.startswith('gene-').mean() > 0.6:
            gene_count_matrix_df['GeneID'] = np.where(
                gene_count_matrix_df['GeneID'].str.startswith('gene-'),
                gene_count_matrix_df['GeneID'].str.extract('gene-(\S+)|')[0],
                gene_count_matrix_df['GeneID']
            )
            gene_count_matrix_df['GeneID'] = np.where(
                gene_count_matrix_df['GeneID'].str.startswith('LOC'),
                gene_count_matrix_df['GeneID'].str.replace('LOC', ''),
                gene_count_matrix_df['GeneID']
            )
        gene_count_matrix_df['GeneID'] = np.where(
            gene_count_matrix_df['GeneID'].str.contains('|'),
            gene_count_matrix_df['GeneID'].str.split('|').str[0],
            gene_count_matrix_df['GeneID']
        )
    if gene_count_matrix_df['GeneID'].str.contains('gene:').mean() > 0.6:
        gene_count_matrix_df['GeneID'] = np.where(
            gene_count_matrix_df['GeneID'].str.contains('gene:'),
            gene_count_matrix_df['GeneID'].str.split('gene:').str[1],
            gene_count_matrix_df['GeneID']
        )
    gene_count_matrix_df.to_csv(parsed_raw_file, sep='\t', index=False)


def matrix_filter(gene_count_matrix_file, output_file, filter_number: int):
    f2=open(output_file,'w')
    num=0
    num1=0

    with open(gene_count_matrix_file,'r') as f1:
        for line in f1.readlines():
            line=line.strip()
            if(num==0):
                each_part=line.split('\t')
                f2.write('GeneID'+ '\t')
                #for each in each_part[1:]:
                f2.write(str('\t'.join(each_part[1:])))
                print(str('\t'.join(each_part[1:])))
                f2.write('\n')
            else:
                each_part = line.split('\t')
                for each in each_part[1:]:
                    if(int(float(each))>filter_number):
                        line1=line.replace('\t','\t')
                        f2.write(line1)
                        f2.write('\n')
                        num1+=1
                        break
            num+=1
    f2.close()
    logger.info('total_filtered(not include head line):' + str(num1))
    logger.info('total(include head line):'+str(num))


def main():
    args = parse_input()
    parsed_raw_file = args.input.replace('.csv', '.txt')
    parse_raw_file(args.input, parsed_raw_file)
    
    if args.id:
        all_geneid_df = pd.read_csv(args.id, sep='\t', header=None, dtype=str)
        all_geneid_df.rename(columns={0: 'GeneID'}, inplace=True)
        before_count = gene_count_matrix_df.shape[0]
        gene_count_matrix_df = pd.merge(left=all_geneid_df, right=gene_count_matrix_df, on='GeneID', how='left')
        gene_count_matrix_df.fillna(0, inplace=True)
        gene_count_matrix_df = gene_count_matrix_df.astype({'GeneID': str})
        gene_count_matrix_df.to_csv(parsed_raw_file, sep='\t', index=False)
        after_count = gene_count_matrix_df.shape[0]
        logger.info(f'before: {before_count}, after: {after_count}，未包含在 all_gene_id.txt 中的基因数：{before_count - after_count}')        

    matrix_filter(parsed_raw_file, args.output, args.n)
    
    sample_lst = load_table(args.s)['sample'].tolist()
    reindex(sample_lst, parsed_raw_file, parsed_raw_file)
    reindex(sample_lst, args.output, args.output)


if __name__ == '__main__':
    main()