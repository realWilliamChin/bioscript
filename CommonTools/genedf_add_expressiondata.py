#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/04/30 09:24
# Author        : William GoGo
import argparse
import pandas as pd

def parse_input():
    parser = argparse.ArgumentParser(description="根据输入文件 GeneID 列合并 expression 表")
    parser.add_argument('-i', '--input', type=str, help='输入文件')
    parser.add_argument('-o', '--output', type=str, help='输出文件')
    parser.add_argument('-e', '--expression', type=str, help='表达量文件，合并的那个文件')
    parser.add_argument('--col-on', type='str', help='根据哪个列合并')
    parser.add_argument('--only-fpkm', dest="only_fpkm", action="store_true", help='只合并 fpkm')
    parser.add_argument('--only-reads', dest="only_reads", action="store_true", help='只合并 reads')
    parser.add_argument('--def', action="store_true", help='合并加上定义')

    return parser.parse_args()


def main():
    # TODO: 继续写，没写完，写了一半，不知道用什么方法写好了。
    args = parse_input()
    input_df = pd.read_csv(args.input, args.expression)
    if args.only_fpkm:
        result_df = pd.merge(input_df, args.expression, on="GeneID", how='left')
    result_df.to_csv(args.output, sep='\t', index=False)
    

if __name__ == '__main__':
    main()