#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/14 17:05
# Author        : William GoGo
"""
nr 注释程序
生成 nr.blast 文件，nr_gene_def.txt 文件 和 nr_TF_def.txt 文件
"""
import argparse
import os
from unittest import result
import pandas as pd
from sympy import use


def parse_input():
    args = argparse.ArgumentParser(description='输入 nr.blast 文件的路径')
    args.add_argument('-n', '--nr-blast-file', type=str, dest='input_file',
                      help='nr.blast file')
    args.add_argument('-b', '--basicinfo', type=str, dest='basicinfo',
                           help='gff 类型, embl or ncbi，默认自动检测，检测失败手动输入')

    return args.parse_args()


def nr(nr_blast_file, gene_basicinfo_file):
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend','sstart','send','evalue','bitscore','stitle']

    data_frame = pd.read_csv(nr_blast_file, sep='\t', names=columns, dtype=str)

    data_frame = data_frame.sort_values(by=['qseqid', 'pident'], ascending=[True, False])
    data_frame = data_frame.drop_duplicates(subset=['qseqid'], keep='first')
    data_frame.to_csv(nr_blast_file.replace('.blast', '_uniq.blast'), sep='\t', index=False)
    nr_gene_def_df = data_frame[['qseqid', 'sseqid', 'stitle']].copy()
    nr_gene_def_df.columns = ['GeneID', 'NCBI_ID', 'NR_Def']
    nr_gene_def_df['NR_Def'] = nr_gene_def_df['NR_Def'].str.split(n=1).str[1]
    print(f'注释上的基因数量是 {nr_gene_def_df.shape[0]} 个')
    
    # 有参基因中没有注释上的 gene 的 biotype 类型，如果是无参不需要此步骤
    if gene_basicinfo_file and os.path.exists(gene_basicinfo_file):
        nr_gene_def_df = nr_def_add_not_protein_coding(nr_gene_def_df, gene_basicinfo_file)
    else:
        print('没有检测到 gene_basicinfo 文件，或文件输入有误，如没有输入 -b 参数，忽略此条消息')
        
    nr_gene_def_df.to_csv(nr_blast_file.replace('.blast', '_gene_def.txt'), sep='\t', index=False)
    
    nr_TF_def_df = nr_gene_def_df[nr_gene_def_df['NR_Def'].str.contains('transcription')]
    nr_TF_def_df = nr_TF_def_df.sort_values(by='NR_Def', key=lambda x: x.str.lower())
    nr_TF_def_df.to_csv(nr_blast_file.replace('.blast', '_TF_def.txt'), sep='\t', index=False)


def nr_def_add_not_protein_coding(nr_gene_def_df, gene_basicinfo_file):
    gff_basicinfo_df = pd.read_csv(gene_basicinfo_file, sep='\t', usecols=['GeneID', 'Gene_Def'])
    print(f'总基因数量是 {gff_basicinfo_df.shape[0]} 个')
    
    # 没有注释上的
    gene_protein_coding_df = gff_basicinfo_df[gff_basicinfo_df['Gene_Def'].str.contains('protein_coding')].copy()
    non_annotationed_df = gene_protein_coding_df[~gene_protein_coding_df['GeneID'].isin(nr_gene_def_df['GeneID'])]
    non_annotationed_df = non_annotationed_df.rename(columns={'Gene_Def': 'NR_Def'})
    print(f'protein coding 没有注释上的有 {non_annotationed_df.shape[0]} 个')
    
    # 不是 protein_coding 的
    gene_non_protein_coding_df = gff_basicinfo_df[~gff_basicinfo_df['Gene_Def'].str.contains('protein_coding')].copy()
    gene_non_protein_coding_df = gene_non_protein_coding_df.rename(columns={'Gene_Def': 'NR_Def'})
    print(f'不是 protein_coding 没有注释上的有 {gene_non_protein_coding_df.shape[0]} 个')
    
    result = pd.concat([nr_gene_def_df, non_annotationed_df, gene_non_protein_coding_df], axis=0)
    result = result.drop_duplicates(subset='GeneID', keep='first')
    result = result.sort_values(by='GeneID')
    print(f'加上没有注释上的和不是 protein_coding 的有 {result.shape[0]} 个')
    if gff_basicinfo_df.shape[0] == result.shape[0]:
        print('\n注释结果正确！')
    
    return result


def main():
    args = parse_input()
    if not args.input_file:
        args.input_file = [x for x in os.listdir() if x.endswith('nr.blast')][0]
    nr(args.input_file, args.basicinfo)
    print('\nDone!\n')


if __name__ == '__main__':
    main()