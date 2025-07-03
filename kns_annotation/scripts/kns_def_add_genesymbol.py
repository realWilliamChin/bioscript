#!/usr/bin/env python
"""

"""
import os
import argparse
import pandas as pd


def parse_input():
    args = argparse.ArgumentParser()
    args.add_argument('-i', '--input', help='输入 kns_def 文件, 默认是当前文件夹下第一个 kns_def.txt')
    args.add_argument('-g', '--genesymbol', help='输入 gene_symbol 文件, 默认是当前文件夹下第一个 gene_symbol.txt')

    return args.parse_args()


def shortname2genesymbol(kns_def_file, genesymbol_file):
    columns = ['GeneID', 'NR_ID', 'NR_ID_Des', 'Swiss_ID', 'Swiss_ID_Des', 'KEGG_ID', 'GeneSymbol', 'EC_Number', 'KEGG_Description']
    kns_df = pd.read_csv(kns_def_file, sep='\t', skiprows=1, names=columns)
    kns_df = kns_df.drop(columns=['GeneSymbol'])
    genesymbol_df = pd.read_csv(genesymbol_file, sep='\t')
    
    # kns_df 的 KEGG_shortname 列替换为 genesymbol_df 的 GeneSymbol 列
    kns_df = pd.merge(left=kns_df, right=genesymbol_df, on='GeneID', how='left')
    # 对 kns_df 的列排序，按照 columns 的顺序
    kns_df = kns_df[columns]
    kns_df = kns_df.fillna('NA')
    kns_df = kns_df.drop_duplicates(subset='GeneID', keep='first')
    kns_df.to_csv(kns_def_file, sep='\t', index=False)
    

def main():
    args = parse_input()
    if not args.input:
        try:
            args.input = [x for x in os.listdir() if x.endswith('kns_def.txt')][0]
        except IndexError:
            print('没有找到 kns_def.txt 文件')
            exit(1)
    if not args.genesymbol:
        try:
            args.genesymbol = [x for x in os.listdir() if x.endswith('GeneSymbol.txt')][0]
        except IndexError:
            print('没有找到 GeneSymbol.txt 文件')
            exit(1)
    
    shortname2genesymbol(args.input, args.genesymbol)
    
    print('\nDone\n')


if __name__ == '__main__':
    main()