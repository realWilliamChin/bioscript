#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2025/4/21 17:38
# Author        : WilliamGoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input-dir', dest='input_dir', default=os.getcwd(), type=str,
                   help='input data dir')
    p.add_argument('-o', '--output-file', dest='output_file', type=str,
                   default='00_Enrichment_significant_pathway_summary.xlsx',
                   help='output file endswith .xlsx, GO in sheet1, KEGG in sheet2')
    p.add_argument('--kegg-pvalue-cutoff', dest='kegg_pvalue_cutoff', type=float,
                   default=0.05, help='< kegg pvalue cutoff')
    p.add_argument('--kegg-count-cutoff', dest='kegg_count_cutoff', type=int,
                   default=2, help='> kegg count cutoff')
    p.add_argument('--go-padjust-cutoff', dest='go_padjust_cutoff', type=float,
                   default=0.05, help='< go padjust cutoff')
    p.add_argument('--go-count-cutoff', dest='go_count_cutoff', type=int,
                   default=2, help='> go count cutoff')
    
    return p.parse_args()


def main():
    args = parse_input()
    go_enrich_df_list = []
    kegg_enrich_df_list = []
    
    for enrich_file in os.listdir(args.input_dir):
        sample_name = enrich_file.split('_Enrichment')[0]
        if enrich_file.endswith('_EnrichmentGO.xlsx'):
            go_enrich_df = load_table(os.path.join(args.input_dir, enrich_file))
            go_enrich_df['Sample'] = sample_name
            go_enrich_df_list.append(go_enrich_df)
        elif enrich_file.endswith('_EnrichmentKEGG.xlsx'):
            kegg_enrich_df = load_table(os.path.join(args.input_dir, enrich_file))
            kegg_enrich_df['Sample'] = sample_name
            kegg_enrich_df_list.append(kegg_enrich_df)
            
    go_summary_df = pd.concat(go_enrich_df_list)
    kegg_summary_df = pd.concat(kegg_enrich_df_list)
    
    go_summary_df = go_summary_df[go_summary_df['p.adjust'] < args.go_padjust_cutoff]
    go_summary_df = go_summary_df[go_summary_df['Count'] > args.go_count_cutoff]
    kegg_summary_df = kegg_summary_df[kegg_summary_df['pvalue'] < args.kegg_pvalue_cutoff]
    kegg_summary_df = kegg_summary_df[kegg_summary_df['Count'] > args.kegg_count_cutoff]
        
    with pd.ExcelWriter(args.output_file, engine='openpyxl') as w:
        go_summary_df.to_excel(w, sheet_name='GO_analysis', index=False)
        kegg_summary_df.to_excel(w, sheet_name='KEGG_analysis', index=False)
    

if __name__ == '__main__':
    main()
