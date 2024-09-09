#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/4/28 14:46
# Author        : WilliamGoGo
import os, sys
import argparse
from loguru import logger
import pandas as pd
import openpyxl


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='输入 group.txt 文件')
    p.add_argument('-d', '--datadir', help='Enrichment 文件目录')
    p.add_argument('-o', '--outputdir', help='输出文件目录')
    
    args = p.parse_args()
    
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
        
    return args


if __name__ == '__main__':
    args = parse_input()
    df = pd.read_csv(args.input, sep='\t')
    
    up_go_suffix = '_Up_EnrichmentGO.xlsx'
    down_go_suffix = '_Down_EnrichmentGO.xlsx'
    up_kegg_suffix = '_Down_EnrichmentKEGG.xlsx'
    down_kegg_suffix = '_Down_EnrichmentKEGG.xlsx'

    grouped = df.groupby('group')
    
    for group_name, group in grouped:
        compare_list = group['compare'].tolist()
        logger.info(f'{group_name} ---- {compare_list}')
        go_df = pd.DataFrame()
        kegg_df = pd.DataFrame()
        for compare in compare_list:
            compare_up_go_file = os.path.join(args.datadir, f'{compare}{up_go_suffix}')
            compare_down_go_file = os.path.join(args.datadir, f'{compare}{down_go_suffix}')
            compare_up_kegg_file = os.path.join(args.datadir, f'{compare}{up_kegg_suffix}')
            compare_down_kegg_file = os.path.join(args.datadir, f'{compare}{down_kegg_suffix}')
            
            compare_up_go_df = pd.read_excel(compare_up_go_file)
            compare_down_go_df = pd.read_excel(compare_down_go_file)
            compare_up_kegg_df = pd.read_excel(compare_up_kegg_file)
            compare_down_kegg_df = pd.read_excel(compare_down_kegg_file)
            
            compare_up_go_df['Regulation'] = 'Up'
            compare_down_go_df['Regulation'] = 'Down'
            compare_go_df = pd.concat([compare_up_go_df, compare_down_go_df], ignore_index=True)
            
            compare_up_kegg_df['Regulation'] = 'Up'
            compare_down_kegg_df['Regulation'] = 'Down'
            compare_kegg_df = pd.concat([compare_up_kegg_df, compare_down_kegg_df], ignore_index=True)

            if go_df.empty:
                go_df = compare_go_df[(compare_go_df['Count'] > 1) & (compare_go_df['pvalue'] < 0.05)]
                kegg_df = compare_kegg_df[(compare_kegg_df['Count'] > 1) & (compare_kegg_df['pvalue'] < 0.05)]
            else:
                compare_go_df = compare_go_df[(compare_go_df['Count'] > 1) & (compare_go_df['pvalue'] < 0.05)]
                compare_kegg_df = compare_kegg_df[(compare_kegg_df['Count'] > 1) & (compare_kegg_df['pvalue'] < 0.05)]
                go_df = pd.concat([go_df, compare_go_df], ignore_index=True)
                kegg_df = pd.concat([kegg_df, compare_kegg_df], ignore_index=True)
        go_df_filename = os.path.join(args.outputdir, f'{group_name}_GO_Enrichment.csv')
        kegg_df_filename = os.path.join(args.outputdir, f'{group_name}_KEGG_Enrichment.csv')
        go_df.to_csv(go_df_filename, index=False)
        kegg_df.to_csv(kegg_df_filename, index=False)
    
    