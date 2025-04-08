#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/29 17:19
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument(dest='ko_file', help='输入包含 KEGG_Pathway_ID 列文件')
    p.add_argument(dest='output_file', help='输出文件 passed_path.txt')
    
    args = p.parse_args()
    
    return args


def passed_path(ko_file, output_file):
    ko_def_df = pd.read_csv(ko_file, sep='\t', names=['KEGG_Pathway_ID'], dtype=str, usecols=[0], header=0)
    passed_path_file = '/home/colddata/qinqiang/script/Analysis/pathview/passed_path.txt'
    passed_path_df = pd.read_csv(passed_path_file, sep='\t', names=['KEGG_Pathway_ID', 'Ko_Def'], dtype=str)
    ko_def_df = passed_path_df[passed_path_df['KEGG_Pathway_ID'].isin(ko_def_df['KEGG_Pathway_ID'])]
    ko_def_df = ko_def_df.drop_duplicates(subset=['KEGG_Pathway_ID'])
    ko_def_df.to_csv(output_file, sep='\t', index=False, header=False)


if __name__ == '__main__':
    args = parse_input()
    
    passed_path(args.ko_file, args.output_file)