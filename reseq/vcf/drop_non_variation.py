#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/10/13 16:16
# Author        : William GoGo

import argparse
import pandas as pd


def parse_input():
    arg = argparse.ArgumentParser(description="")
    arg.add_argument('-i', '--input', required=True, help='输入需要过滤的vcf文件')
    arg.add_argument('-o', '--output', required=True, help='输入过滤后的文件名')
    return arg.parse_args()


def main():
    args = parse_input()
    with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
        for line in infile:
            # 保留头部信息
            if line.startswith("#"):
                outfile.write(line)
                continue
            
            # 分割行以获取基因型信息
            split_line = line.strip().split("\t")
            genotype_info = split_line[9].split(":")[0]  # 假设基因型在第10列
            
            # 过滤基因型
            if genotype_info not in {"0|0", "0/0", ".|.", "./."}:
                outfile.write(line)
                
    print('\nDone!\n')


if __name__ == '__main__':
    main()