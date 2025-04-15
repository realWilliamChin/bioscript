#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/01/10 11:03
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='输入 table_depth_ASV_w_tax.txt 文件')
    p.add_argument('-o', '--output', help='输出文件夹，包含 2 个文件，Species.txt Species_taxon.txt')
    
    return p.parse_args()



def main():
    args = parse_input()
    input_file, output_dir = args.input, args.output
    
    df = pd.read_csv(input_file, sep='\t', skiprows=1)

    # df.rename(columns={'#OTU ID', 'OTU'}, inplace=True)
    df['OTU'] = 'OTU_' + (pd.RangeIndex(start=0, stop=len(df)).astype(str).str.zfill(3))
    
    taxon_df = df[['OTU', 'taxonomy']].copy()
    
    df.drop(columns=['taxonomy', '#OTU ID'], inplace=True)
    print(df)
    df.set_index('OTU', inplace=True)
    
    
    df.to_csv(os.path.join(output_dir, 'Species.txt'), sep='\t')
    


if __name__ == '__main__':
    main()
