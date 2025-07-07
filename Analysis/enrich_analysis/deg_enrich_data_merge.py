#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/4/28 14:46
# Author        : WilliamGoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input-dir', default=os.getcwd(), help='input data dir')
    p.add_argument('-c', '--compare-file',  help='input compare_info.txt file')
    p.add_argument('-o', '--output-file', default='DEG_enrichment_significant_pathway_summary.xlsx',
                   help='output file endswith .xlsx, GO in sheet1, KEGG in sheet2')
    
    return p.parse_args()


def deg_enrich_data_merge(input_dir, comp_file, output_file):
    go_enrich_df_down_list = []
    go_enrich_df_up_list = []
    kegg_enrich_df_down_list = []
    kegg_enrich_df_up_list = []
    
    compare_df = load_table(comp_file)
    for i, compare in compare_df.iterrows():
        compare_info = compare['Treat'] + "-vs-" + compare['Control']

        for enrich_file in os.listdir(input_dir):
            if enrich_file.startswith(compare_info) and enrich_file.endswith('Down_EnrichmentGO.xlsx'):
                go_enrich_df_down = load_table(os.path.join(input_dir, enrich_file))
                go_enrich_df_down['Comparison_ID'] = compare_info
                go_enrich_df_down['Regulation'] = 'Down'
                go_enrich_df_down_list.append(go_enrich_df_down)
            elif enrich_file.startswith(compare_info) and enrich_file.endswith('Up_EnrichmentGO.xlsx'):
                go_enrich_df_up = load_table(os.path.join(input_dir, enrich_file))
                go_enrich_df_up['Comparison_ID'] = compare_info
                go_enrich_df_up['Regulation'] = 'Up'
                go_enrich_df_up_list.append(go_enrich_df_up)
            elif enrich_file.startswith(compare_info) and enrich_file.endswith('Down_EnrichmentKEGG.xlsx'):
                kegg_enrich_df_down = load_table(os.path.join(input_dir, enrich_file))
                kegg_enrich_df_down['Comparison_ID'] = compare_info
                kegg_enrich_df_down['Regulation'] = 'Down'
                kegg_enrich_df_down_list.append(kegg_enrich_df_down)
            elif enrich_file.startswith(compare_info) and enrich_file.endswith('Up_EnrichmentKEGG.xlsx'):
                kegg_enrich_df_up = load_table(os.path.join(input_dir, enrich_file))
                kegg_enrich_df_up['Comparison_ID'] = compare_info
                kegg_enrich_df_up['Regulation'] = 'Up'
                kegg_enrich_df_up_list.append(kegg_enrich_df_up)
                
                
        go_down_summary = pd.concat(go_enrich_df_down_list)
        go_down_summary['p.adjust'] = go_down_summary['p.adjust'].astype(float)
        go_down_summary = go_down_summary[go_down_summary['Count'] >= 2]
        go_down_summary = go_down_summary[go_down_summary['p.adjust'] < 0.05]
        go_down_summary.sort_values(by='Description', key=lambda x: x.str.lower(), inplace=True)
        
        go_up_summary = pd.concat(go_enrich_df_up_list)
        go_up_summary['p.adjust'] = go_up_summary['p.adjust'].astype(float)
        go_up_summary = go_up_summary[go_up_summary['Count'] >= 2]
        go_up_summary = go_up_summary[go_up_summary['p.adjust'] < 0.05]
        go_up_summary.sort_values(by='Description', key=lambda x: x.str.lower(), inplace=True)
        
        kegg_down_summary = pd.concat(kegg_enrich_df_down_list)
        kegg_down_summary['pvalue'] = kegg_down_summary['pvalue'].astype(float)
        kegg_down_summary = kegg_down_summary[kegg_down_summary['Count'] >= 2]
        kegg_down_summary = kegg_down_summary[kegg_down_summary['pvalue'] < 0.05]
        kegg_down_summary.sort_values(by='Description', key=lambda x: x.str.lower(), inplace=True)
        
        kegg_up_summary = pd.concat(kegg_enrich_df_up_list)
        kegg_up_summary['pvalue'] = kegg_up_summary['pvalue'].astype(float)
        kegg_up_summary = kegg_up_summary[kegg_up_summary['Count'] >= 2]
        kegg_up_summary = kegg_up_summary[kegg_up_summary['pvalue'] < 0.05]
        kegg_up_summary.sort_values(by='Description', key=lambda x: x.str.lower(), inplace=True)
        
        go_summary = pd.concat([go_down_summary, go_up_summary])
        kegg_summary = pd.concat([kegg_down_summary, kegg_up_summary])
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as w:
            go_summary.to_excel(w, sheet_name='GO_analysis', index=False)
            kegg_summary.to_excel(w, sheet_name='KEGG_analysis', index=False)
    


def main():
    args = parse_input()
    deg_enrich_data_merge(args.input_dir, args.compare_file, args.output_file)


if __name__ == '__main__':
    main()
    
    