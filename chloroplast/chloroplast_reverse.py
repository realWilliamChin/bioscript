#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/05 19:00
# Author        : William GoGo
import argparse
from Bio.Seq import reverse_complement


def parse_input():
    argparser = argparse.ArgumentParser(description='')
    argparser.add_argument('-i', '--input', help='输入文件',
                           default='embplant_pt.K105.scaffolds.graph1.1.path_sequence.fasta')
    argparser.add_argument('-o', '--output', help='输出文件')
    return argparser.parse_args()
    

def main():
    args = parse_input()

    f2=open(args.output,'w')
    with open(args.input,'r') as f, open(args.output, 'w') as out_f:

        for each_line in f:
            if '>' not in each_line:
                line=each_line
                break
            else:
                out_f.write(each_line)
        my_seq = reverse_complement(line)
        out_f.write(my_seq)
    print('Done!')


if __name__ == '__main__':
    main()
