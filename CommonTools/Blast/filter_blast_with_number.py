#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/25 10:06
# Author        : William GoGo
import argparse
import pandas as pd
import re

def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('input_file', help='输入 .blast 结尾的文件')
    p.add_argument('output_file', help='输出过滤之后的 .blast 文件')
    p.add_argument('--filter', help='过滤选项（列名）', required=True)
    p.add_argument('--filter-value', help='过滤值，格式为 >20', required=True)
    p.add_argument('--header', help='如果没有列名，请输入列名，以 "," 分隔', default=None)
    
    return p.parse_args()


def parse_filter_value(filter_value):
    match = re.match(r'([><=]=?)(\d+\.?\d*)', filter_value)
    if match:
        operator, value = match.groups()
        return operator, float(value)
    else:
        raise ValueError(f"无效的过滤值格式: {filter_value}")


def main(input_file, output_file, filter_column, filter_value, header=None):
    if header:
        header_list = header.split(',')
        df = pd.read_csv(input_file, sep='\t', header=None, names=header_list)
    else:
        df = pd.read_csv(input_file, sep='\t')

    operator, value = parse_filter_value(filter_value)
    operations = {
        '>': lambda x, y: x > y,
        '>=': lambda x, y: x >= y,
        '=': lambda x, y: x == y,
        '<': lambda x, y: x < y,
        '<=': lambda x, y: x <= y
    }

    if filter_column not in df.columns:
        raise ValueError(f"列 {filter_column} 不存在于数据中")

    df_filtered = df[operations[operator](df[filter_column], value)]
    df_filtered.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    args = parse_input()
    main(args.input_file, args.output_file, args.filter, args.filter_value, args.header)
