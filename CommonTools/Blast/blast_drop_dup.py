#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/10/17 17:47
# Author        : William GoGo
import sys
import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description='按照 qacc 第一列和倒数第二列 bitscore 的得分进行去重')
    parser.add_argument('-i', '--infile', required=True, help='输入 blast 文件')
    parser.add_argument('-o', '--outfile', help='输出 uniq blast 的文件名')
    parser.add_argument('-n', '--column-names', help='输出 blast 列名 格式这么写，qacc,sacc,pident', dest='column_names',
                        default="qacc sacc pident qcovs qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle")
    
    args = parser.parse_args()
    args.column_names = args.column_names.split()

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.outfile:
        outfile = args.outfile
    else:
        outfile = str(args.infile).replace('.blast', '_uniq.blast')
    
    df = pd.read_csv(args.infile, sep='\t', header=None)
    print(df.columns.values)
    # 计算倒数第二列的正整数索引
    index_last_second = df.shape[1] - 2  # 列总数减2得到倒数第二列的索引
    # 按第一列升序和倒数第二列降序排序
    df.sort_values(by=[df.columns[0], df.columns[index_last_second]], ascending=[True, False], inplace=True)
    df.drop_duplicates(subset=[0], keep='first', inplace=True)
    if args.column_names:
        column_names = [x.strip() for x in args.column_names.split(',')]
        df.to_csv(outfile, sep='\t', index=False, header=column_names)
    else:
        df.to_csv(outfile, sep='\t', index=False, header=None)

if __name__ == '__main__':
    main()
    