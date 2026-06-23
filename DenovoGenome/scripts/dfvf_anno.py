#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/07 14:41
# Author        : William GoGo
import os, sys
import argparse
import subprocess
import pandas as pd

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


# blastp -outfmt 6 的列名(与下方 dfvf_blast 的 -outfmt 保持一致)
BLAST_COLUMNS = ['qacc', 'sacc', 'qcovhsp', 'ppos', 'length', 'mismatch', 'gapopen',
                 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='输入需要注释的 fasta 文件')
    p.add_argument('-o', '--output', help='输出文件')

    args = p.parse_args()

    return args


def dfvf_blast(fasta_file, output_file):
    cmd = [
        "/home/data/opt/biosoft/ncbi-blast-2.9.0+/bin/blastp",
        "-db", "/home/colddata/qinqiang/script/DenovoGenome/libs/DFVF/DFVF",
        "-query", fasta_file,
        "-out", output_file,
        "-evalue", "1e-5",
        "-num_threads", "18",
        "-outfmt", "6 qacc sacc qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    ]
    subprocess.run(cmd, check=True)


def add_def(blast_file, output_file):
    def_df = load_table('/home/colddata/qinqiang/script/DenovoGenome/libs/DFVF/All_genes_def.txt')
    blast_df = load_table(blast_file, header=None, names=BLAST_COLUMNS)
    result_df = pd.merge(blast_df, def_df, left_on='sacc', right_on='GeneID', how='left')
    result_df = result_df.drop(columns=['GeneID'])
    write_output_df(result_df, output_file, index=False)


def main():
    args = parse_input()
    blast_file = f'{args.input}.dfvf.blast'
    dfvf_blast(args.input, blast_file)
    add_def(blast_file, args.output)


if __name__ == '__main__':
    main()