#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/06/28 16:35
# Author        : William GoGo
# Version       : 1.0
# Description   : 对 gff3 文件找出有序列的一行，提取出 fasta 文件
import os, sys
import re
import argparse
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser(description="Extract fasta sequence from gff3 file")
    parser.add_argument('-g', dest="gff3", help="gff3 file", default='tRNA.gff3')
    return parser.parse_args()


def main():
    args = parse_input()
    gff3 = args.gff3
    output_file = gff3.replace('.gff3', '.fasta')
    with open(gff3, "r") as f:
        lines = f.readlines()
    with open(output_file, "w") as f1, open(gff3.replace('.gff3', '_rmseq.gff3'), "w") as f2:
        for line in lines:
            if 'Sequence=' in line:
                gene_id = re.findall(r'ID=(.*?);', line)[0]
                seq = re.findall(r'Sequence=([A-Za-z]+)', line)[0]
                f1.write(f">{gene_id}\n{seq}\n")
                f2.write(line.replace(f"Sequence={seq}", '').replace(';;', ';'))
            else:
                f2.write(line)
                
    print("Done!")


if __name__ == "__main__":
    main()