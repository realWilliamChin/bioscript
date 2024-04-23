#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/11/14 22:31
# Author        : William GoGo
import os
import argparse
import pandas as pd


def sort_excel_sheets(input_file, output_file):
    xls = pd.ExcelFile(input_file, engine='openpyxl')
    writer = pd.ExcelWriter(output_file, engine='openpyxl')

    for sheet_name in xls.sheet_names:
        data = pd.read_excel(xls, sheet_name=sheet_name, header=None, engine='openpyxl')
        comments = data.iloc[:9, :]
        content = data.iloc[9:, :]
        sorted_content = content.sort_values(by=6)
        final_data = pd.concat([comments, sorted_content])
        final_data.to_excel(writer, sheet_name=sheet_name, index=False, header=False)
    writer.close()


def main():
    os.mkdir('output')
    for f in os.listdir():
        if f.endswith('.xlsx'):
            print(f"Sorting {f} ...")
            sort_excel_sheets(f, "output/"+f)


if __name__ == '__main__':
    main()
