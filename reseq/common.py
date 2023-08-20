#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/07 11:35
# Author        : William GoGo


def skip_rows(file_name, skip_str):
    # 找出 vcf 文件中以 '##' 开头的行数，读取时跳过这些行
    skip_rows = 0
    with open(file_name, "r") as file:
        for line in file:
            if line.startswith(skip_str):
                skip_rows += 1
            else:
                break
    return skip_rows