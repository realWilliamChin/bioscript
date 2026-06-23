#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/06/28 10:21
# Author        : William GoGo
# Description   : Convert cmsearch output rfam_out.tab to gff format
import os, sys
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    parser = argparse.ArgumentParser(description="Convert cmsearch output rfam_out.tab to gff")
    parser.add_argument('-i', '--input', type=str, dest='input_file', help="Input file")
    parser.add_argument('-o', '--output', type=str, help="Output file")
    # parser.add_argument('-f', '--filter', type=int, help="score 小于多少的扔掉")
    
    return parser.parse_args()


def process_original_gff3(input_file, output_file):
    title_list = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute']
    with open(input_file, 'r') as f, open(output_file, 'a') as w:
        w.write('\t'.join(title_list) + '\n')
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith('#'):
                continue
            else:
                line_list = line.split()
                w.write('\t'.join(line_list) + '\n')


def main():
    args = parse_input()
    input_file = args.input_file
    output_file = args.output
    process_original_gff3(input_file, output_file)
    

if __name__ == '__main__':
    main()
