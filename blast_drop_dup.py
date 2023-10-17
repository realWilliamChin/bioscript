#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/10/17 17:47
# Author        : William GoGo

import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description='biogrid，')
    parser.add_argument('-i', '--infile', required=True, help='输入 blast 文件')
    parser.add_argument('-o', '--outfile', help='输出 uniq blast 的文件名')
    return parser.parse_args()


def drop_dup(blast_file, out_file):
    # blast_names = ['qacc','sacc','qcovhsp','ppos','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','stitle']
    # df = pd.read_csv(blast_file, sep='\t', names=blast_names)
    df = pd.read_csv(blast_file, sep='\t', header=None)
    print(df.columns.values)
    df.sort_values(by=[0, 2], ascending=[True, False], inplace=True)
    df.drop_duplicates(subset=0, keep='first', inplace=True)
    df.to_csv(out_file, sep='\t', index=False, header=False)
    

def main():
    args = parse_arguments()
    if args.outfile:
        outfile = args.outfile
    else:
        outfile = str(args.infile).replace('.blast', '_uniq.blast')
    drop_dup(args.infile, outfile)
    

if __name__ == '__main__':
    main()
    