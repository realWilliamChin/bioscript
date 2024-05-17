#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/5/29 19:43
# Author        : WilliamGoGo
# 合并 fpkm 和 reads 然后添加 def
import os, sys
import argparse
import pandas as pd

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/transcriptome/'))
from genedf_add_knsdef import add_kns_def


def parse_input():
    parser = argparse.ArgumentParser(description='合并 fpkm 和 reads matrix, 然后加上三个注释的定义')
    parser.add_argument('-f', '--fpkm', type=str, required=True, help='指定 fpkm 文件')
    parser.add_argument('-r', '--reads', type=str, required=True, help='指定 reads 文件')
    parser.add_argument('--kns', help='输入 kns_def 文件，添加 KEGG_ID, KEGG_GeneID, NR_Def, Swiss_protein_ID')
    parser.add_argument('-k', '--kegg', type=str, help='指定 kegg_gene_def file')
    parser.add_argument('-n', '--nr', type=str, help='指定 nr_gene_def file')
    parser.add_argument('-s', '--swiss', type=str, help='指定 swiss_gene_def file')
    parser.add_argument('-o', '--output', type=str, help='指定 output 文件名，不指定默认 fpkm_and_reads_matrix_filtered_data_def.txt')

    return parser.parse_args()


def merge_fpkm_reads(fpkm_file, reads_file):
    fpkm_df = pd.read_csv(fpkm_file, sep='\t', dtype={"GeneID": str})
    reads_df = pd.read_csv(reads_file, sep='\t', dtype={"GeneID": str})
    # 对 GeneID 列进行重命名，如果是用其他方式写的 gene_id geneid 等等
    if 'gene' in fpkm_df.columns[0].lower() and 'id' in fpkm_df.columns[0].lower():
        fpkm_df.rename(columns={fpkm_df.columns[0]: 'GeneID'}, inplace=True)
    if 'gene' in reads_df.columns[0].lower() and 'id' in reads_df.columns[0].lower():
        reads_df.rename(columns={reads_df.columns[0]: 'GeneID'}, inplace=True)

    # 列名添加 _fpkm _reads，先把列名提取出来加上，然后给表换上新的列名
    fpkm_column_lst = list(each + '_fpkm' for each in list(fpkm_df.columns)[1:])
    reads_column_lst = list(each + '_reads' for each in list(reads_df.columns)[1:])
    fpkm_column_lst.insert(0, fpkm_df.columns[0])
    reads_column_lst.insert(0, fpkm_column_lst[0])
    fpkm_df.columns = fpkm_column_lst
    reads_df.columns = reads_column_lst
    
    # 合并
    reads_fpkm_df = pd.merge(left=fpkm_df, right=reads_df, on=fpkm_df.columns[0], how='left')
    return reads_fpkm_df


def main():
    args = parse_input()
    merged_df = merge_fpkm_reads(args.fpkm, args.reads)
    result_df = add_kns_def(merged_df, args.kegg, args.nr, args.swiss, args.kns)
    if args.output:
        result_df.to_csv(args.output, sep='\t', index=False)
    else:
        result_df.to_csv('fpkm_and_reads_matrix_filtered_data_def.txt', sep='\t', index=False)
        
    print('\nDone!\n')


if __name__ == '__main__':
    main()
