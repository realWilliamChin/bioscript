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
    argparser.add_argument('-a', '--all', default='all_geneid_GO.txt',
                           help='指定所有基因的 GO 文件, 默认当前文件夹下的 _all_geneid_Go 文件')
    argparser.add_argument('-f', '--file',
                           help='指定基因的 ID 文件(结尾必须是_ID.txt), 默认当前文件夹下的 Up Down_ID.txt 文件')
    argparser.add_argument('-t', '--generate_type', default='wego_genomics', choices=['wego_genomics', 'tbtools'],
                           help='指定生成类型, wego_genomics(tab 分隔) 或者 tbtools(逗号分隔), 默认 wego_web)')
    argparser.add_argument('-p', '--picture', action='store_true', default=False,
                           help='是否生成图片，默认不生成')
    argparser.add_argument('-l', '--level', type=str, default='3',
                           help='指定 GO 的 level, 默认 3, 与 -p 一块儿使用')
    args = argparser.parse_args()
    return args


def geneid_goid(all_geneid_goid_file, updownid_file, generate_type):
    """
    生成每个样本的上下调基因的 GO 文件
    """
    updown_go_filename = updownid_file.replace('_ID.txt', '_GO.txt')
    # 读取 updownid_file 为 list
    print(f'reading {updownid_file}, generating {updown_go_filename} ...')
    with open(updownid_file, 'r') as id_lst_file:
        id_lst = list(id_lst_file.read().splitlines())
    # 循环读取 all_geneid_goid_file，如果 GeneID 在 id_lst 中，就写入 updown_go_filename
    with open(all_geneid_goid_file, 'r') as input_file, open(updown_go_filename, 'w') as updown_go_file:
        for line in input_file:
            geneid = line.split('\t')[0]
            if geneid in id_lst:
                # 如果是 wego_genomics，就直接写入，如果是 tbtools，就把除第一个 tab 换成逗号
                if generate_type == 'tbtools':
                    line = line.replace('\t', ',')
                    line = line.replace(',', '\t', 1)
                updown_go_file.write(line)
    return updown_go_filename


def main():
    args = parse_input()
    if args.file:
        geneid_goid(args.all, args.file, args.generate_type)
    else:
        for each_file in os.listdir():
            if 'ID.txt' not in each_file:
                continue
            if 'Up' or 'Down' in each_file:
                Go_filename = geneid_goid(args.all, each_file, args.generate_type)
            if args.picture:
                get_wego_pic(Go_filename, args.level)
    

if __name__ == '__main__':
    main()
    