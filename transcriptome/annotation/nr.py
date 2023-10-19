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
import pandas as pd


def parse_input():
    args = argparse.ArgumentParser(description='输入 nr.blast 文件的路径')
    args.add_argument('-i', '--input', type=str, help='nr.blast file')
    
    return args.parse_args()


def nr(nr_blast_file):
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend','sstart','send','evalue','bitscore','stitle']

    data_frame = pd.read_csv(nr_blast_file, sep='\t', names=columns)

    data_frame = data_frame.sort_values(by=['qseqid', 'pident'], ascending=[True, False])
    data_frame = data_frame.drop_duplicates(subset=['qseqid'], keep='first')
    nr_gene_def_df = data_frame[['qseqid', 'sseqid', 'stitle']].copy()
    nr_gene_def_df.columns = ['GeneID', 'NCBI_ID', 'NR_Def']
    nr_gene_def_df['NR_Def'] = nr_gene_def_df['NR_Def'].str.split(n=1).str[1]
    nr_gene_def_df.to_csv(nr_blast_file.replace('.blast', '_gene_def.txt'), sep='\t', index=False)
    nr_TF_def_df = nr_gene_def_df[nr_gene_def_df['NR_Def'].str.contains('transcription')]
    nr_TF_def_df = nr_TF_def_df.sort_values(by='NR_Def', key=lambda x: x.str.lower())
    nr_TF_def_df.to_csv(nr_blast_file.replace('.blast', '_TF_def.txt'), sep='\t', index=False)


def main():
    args = parse_input()
    if not args.input:
        args.input = [x for x in os.listdir() if x.endswith('nr.blast')][0]
    nr(args.input)
    print('\nDone!\n')


if __name__ == '__main__':
    main()

