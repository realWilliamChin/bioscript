#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/01/22 16:07
# Author        : William GoGo
"""
输入文件是一个 GeneID list，包含 t1 t2，包含 header
结果生成一个包含替换新的 GeneID 的文件 TargetGeneID
"""
import os, sys
import re
import argparse
from loguru import logger
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='input', help='gff3 输入文件')
    parser.add_argument(dest='output', help='输出文件')
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_input()
    
    result_list = []
    with open(args.input, 'r') as i, open(args.output, 'a') as o:
        o.write(f'GeneID\tTargetGeneID\n')
        line_list = i.readlines()
        for each_line in line_list:
            if each_line == line_list[0]:
                result_list.append([each_line.strip()])
            elif each_line.split('.t')[0] == result_list[-1][0].split('.t')[0]:
                result_list[-1] += [each_line.strip()]
            else:
                result_list.append([each_line.strip()])
        
        for each_id_list in result_list:
            num = 1
            for each_id in each_id_list:
                each_new_id = each_id.split('.t')[0] + '.t' + str(num)
                num += 1
                o.write(f'{each_id}\t{each_new_id}\n')