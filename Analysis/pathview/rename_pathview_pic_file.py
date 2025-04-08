#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/01 23:06
# Author        : William GoGo

import os, sys
import argparse
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--passed-path-file', dest='input_file', default='/home/colddata/qinqiang/script/Analysis/pathview/passed_path.txt',
                        help='KEGG_ID 和 KEGG_DEF 文件，默认 passed_path.txt')
    parser.add_argument('-i', '--input-dir', dest='input_dir', default='./', help='所有图片目录，默认当前目录')
    parser.add_argument('-o', '--output-dir', dest='output_dir', help='输出目录，默认当前目录', default='./')
    
    args = parser.parse_args()
    
    return args


def rename(input_dir, ref_file):
    ref_df = pd.read_csv(ref_file, sep='\t', header=None, names=['KO', 'Def'])
    ref_df = ref_df.drop_duplicates(subset='KO')
    ko_png_list = [x for x in os.listdir(input_dir) if x.endswith('ko.data.png') and '_' not in x]
    for png in ko_png_list:
        png_ko_name = os.path.basename(png).split('.')[0]
        replace_name = ref_df[ref_df['KO'] == png_ko_name]['Def'].to_list()
        if len(replace_name) == 0:
            continue
        # print(type(replace_name), replace_name)
        # exit(0)
        replace_name = f'{png_ko_name}_{replace_name[0]}.ko.data.png'
        print(png, png_ko_name, replace_name)
        # exit(0)
        os.rename(os.path.join(input_dir, png), os.path.join(input_dir, replace_name))


if __name__ == '__main__':
    args = parse_input()
    rename(args.input_dir, args.input_file)