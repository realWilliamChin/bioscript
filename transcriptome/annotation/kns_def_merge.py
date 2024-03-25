#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/18 09:22
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/transcriptome/'))
from genedf_add_def import add_kns_def


def parse_input():
    argparser = argparse.ArgumentParser(description='输入 k，n，s，def 文件，添加 KEGG_ID, KEGG_GeneID, NR_Def, Swiss_protein_ID')
    argparser.add_argument('-k', '--kegg', help='KEGG_gene_def 文件，会添加 KEGG_ID 和 KEGG_Shortname 列')
    argparser.add_argument('-n', '--nr', help='NR_gene_def 文件')
    argparser.add_argument('-s', '--swiss', help='Swiss_gene_def 文件')
    argparser.add_argument('-b', '--basicinfo', help="basic_info.txt 文件，如果指定则添加到 Kns 中")
    argparser.add_argument('-i', '--input', required=True, help='输入 all_geneid 文件')
    argparser.add_argument('-o', '--output', help='输出文件，支持 xlsx, csv, txt 格式输出')
    
    args = argparser.parse_args()
    
    if not args.output:
        args.output = args.input.split(os.sep)[-1].split('_all')[0] + '_kns_def.txt'
    
    return args


def main():
    args = parse_input()
        
    df = pd.read_csv(args.input, sep='\t', names=['GeneID'], dtype={"GeneID": str})
    
    if args.basicinfo:
        basic_df = pd.read_csv(args.basicinfo, sep='\t', dtype={"GeneID": str})
        result_df = pd.merge(df, basic_df, on="GeneID", how="left")
    
    result_df = add_kns_def(result_df, args.kegg, args.nr, args.swiss)
    
    

    if args.output.endswith('.csv'):
        result_df.to_csv(args.output, index=False)
    elif args.output.endswith('.xlsx'):
        result_df.to_excel(args.output, index=False, engine='openpyxl')
    else:
        result_df.to_csv(args.output, sep='\t', index=False)

    print('\nDone!\n')


if __name__ == '__main__':
    main()



