#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/02 01:25
# Author        : William GoGo

import os, sys
import argparse
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='input_file', help='input file')
    parser.add_argument(dest='output_file', help='output_file')
    
    args = parser.parse_args()
    
    return args


def fake_data(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    fake_columns = df.columns[1:]
    
    for col in fake_columns:
        df[col].astype(float)
        df[f'{col}_1'] = df[col]
        df[f'{col}_2'] = df[col] * 0.95
        df[f'{col}_3'] = df[col] * 1.05
        df[f'{col}_1'] = df[f'{col}_1'].astype(int)
        df[f'{col}_2'] = df[f'{col}_2'].astype(int)
        df[f'{col}_3'] = df[f'{col}_3'].astype(int)
        df = df.drop(columns=[col])
    
    df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    args = parse_input()
    fake_data(args.input_file, args.output_file)
    