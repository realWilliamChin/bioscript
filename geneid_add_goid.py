#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/27 11:00
# Author        : William GoGo
import os
import pandas as pd
import argparse
from get_wego_pic import get_wego_pic

def parse_input():
    argparser = argparse.ArgumentParser(description='')
    argparser.add_argument('-a', '--all', default='all_geneid_Go.txt',
                           help='指定所有基因的 GO 文件, 默认当前文件夹下的 _all_geneid_Go 文件')
    argparser.add_argument('-p', '--picture', action='store_true', default=False,
                           help='是否生成图片，默认不生成')
    args = argparser.parse_args()
    return args


def geneid_goid(all_geneid_goid_file, updownid_file):
    """
    生成每个样本的上下调基因的 GO 文件
    """
    updown_go_filename = updownid_file.replace('_ID.txt', '_Go.txt')
    # 读取 updownid_file 为 list
    with open(updownid_file, 'r') as file:
        id_lst = list(file.read().splitlines())
    # 循环读取 all_geneid_goid_file，如果 GeneID 在 id_lst 中，就写入 updown_go_filename
    with open(all_geneid_goid_file, 'r') as file:
        for line in file:
            geneid = line.split('\t')[0]
            if geneid in id_lst:
                with open(updown_go_filename, 'a') as updown_go_file:
                    updown_go_file.write(line)


def main():
    args = parse_input()
    for each_file in os.listdir():
        if 'ID.txt' not in each_file:
            continue
        if 'Up' or 'Down' in each_file:
            geneid_goid(args.all, each_file)
    if args.picture:
        for each_file in os.listdir():
            if each_file.endswith('Go.txt'):
                get_wego_pic(each_file)
    

if __name__ == '__main__':
    main()
    