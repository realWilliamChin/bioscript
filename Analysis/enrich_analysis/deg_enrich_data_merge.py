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
    p.add_argument('--filter-col', choices=['p.adjust', 'pvalue', 'none'], default='pvalue',
                   help='选择用于显著性筛选的列，GO分析默认p.adjust, KEGG分析默认pvalue')
    p.add_argument('--filter-value', type=float, default=0.05, help='显著性筛选的阈值，默认 0.05')
    p.add_argument('--count-threshold', type=int, default=2, help='Count筛选的最小值，默认 2')
    
    return p.parse_args()


def deg_enrich_data_merge(input_dir, comp_file, output_file, filter_col='pvalue', filter_value=0.05, count_threshold=2):
    go_enrich_df_down_list = []
    go_enrich_df_up_list = []
    kegg_enrich_df_down_list = []
    kegg_enrich_df_up_list = []
    
    compare_df = load_table(comp_file)
    logger.info(f'已读取 Compare_info，总共有 {compare_df.dtypes[0]} 行')
    for i, compare in compare_df.iterrows():
        compare_info = compare['Treat'] + "-vs-" + compare['Control']
        logger.info(f'读取 {compare_info} EnrichmentGO/KEGG 等文件中')
        for enrich_file in os.listdir(input_dir):
            if enrich_file.startswith(f'{compare_info}_Down_') and 'EnrichmentGO' in enrich_file:
                go_enrich_df_down = load_table(os.path.join(input_dir, enrich_file))
                go_enrich_df_down['Comparison_ID'] = compare_info
                go_enrich_df_down['Regulation'] = 'Down'
                go_enrich_df_down_list.append(go_enrich_df_down)
            elif enrich_file.startswith(f'{compare_info}_Up_') and 'EnrichmentGO' in enrich_file:
                go_enrich_df_up = load_table(os.path.join(input_dir, enrich_file))
                go_enrich_df_up['Comparison_ID'] = compare_info
                go_enrich_df_up['Regulation'] = 'Up'
                go_enrich_df_up_list.append(go_enrich_df_up)
            elif enrich_file.startswith(f'{compare_info}_Down_') and 'EnrichmentKEGG' in enrich_file:
                kegg_enrich_df_down = load_table(os.path.join(input_dir, enrich_file))
                kegg_enrich_df_down['Comparison_ID'] = compare_info
                kegg_enrich_df_down['Regulation'] = 'Down'
                kegg_enrich_df_down_list.append(kegg_enrich_df_down)
            elif enrich_file.startswith(f'{compare_info}_Up_') and 'EnrichmentKEGG' in enrich_file:
                kegg_enrich_df_up = load_table(os.path.join(input_dir, enrich_file))
                kegg_enrich_df_up['Comparison_ID'] = compare_info
                kegg_enrich_df_up['Regulation'] = 'Up'
                kegg_enrich_df_up_list.append(kegg_enrich_df_up)
    
    if len(go_enrich_df_down_list) == 0:
        logger.error('EnrichmentGO Down 没有找到任何文件')
    if len(go_enrich_df_up_list) == 0:
        logger.error('EnrichmentGO Up 没有找到任何文件')
    if len(kegg_enrich_df_down_list) == 0:
        logger.error('EnrichmentKEGG Down 没有找到任何文件')
    if len(kegg_enrich_df_up_list) == 0:
        logger.error('EnrichmentKEGG Up 没有找到任何文件')
    
                
    # go筛选
    logger.info('正在对 GO 进行过滤')
    go_down_summary = pd.concat(go_enrich_df_down_list)
    go_down_summary["p.adjust"] = go_down_summary["p.adjust"].astype(float)
    go_down_summary = go_down_summary[go_down_summary['Count'] >= count_threshold]
    go_down_summary = go_down_summary[go_down_summary[filter_col] < filter_value] if filter_col in go_down_summary.columns else go_down_summary
    go_down_summary.sort_values(by='Description', key=lambda x: x.str.lower(), inplace=True)
    
    go_up_summary = pd.concat(go_enrich_df_up_list)
    go_up_summary["p.adjust"] = go_up_summary["p.adjust"].astype(float)
    go_up_summary = go_up_summary[go_up_summary['Count'] >= count_threshold]
    go_up_summary = go_up_summary[go_up_summary[filter_col] < filter_value] if filter_col in go_up_summary.columns else go_up_summary
    go_up_summary.sort_values(by='Description', key=lambda x: x.str.lower(), inplace=True)
    
    # kegg筛选
    logger.info('正在对 KEGG 进行过滤')
    kegg_down_summary = pd.concat(kegg_enrich_df_down_list)
    kegg_down_summary["pvalue"] = kegg_down_summary["pvalue"].astype(float)
    kegg_down_summary = kegg_down_summary[kegg_down_summary['Count'] >= count_threshold]
    kegg_down_summary = kegg_down_summary[kegg_down_summary[filter_col] < filter_value] if filter_col in kegg_down_summary.columns else kegg_down_summary
    kegg_down_summary.sort_values(by='Description', key=lambda x: x.str.lower(), inplace=True)
    
    kegg_up_summary = pd.concat(kegg_enrich_df_up_list)
    kegg_up_summary["pvalue"] = kegg_up_summary["pvalue"].astype(float)
    kegg_up_summary = kegg_up_summary[kegg_up_summary['Count'] >= count_threshold]
    kegg_up_summary = kegg_up_summary[kegg_up_summary[filter_col] < filter_value] if filter_col in kegg_up_summary.columns else kegg_up_summary
    kegg_up_summary.sort_values(by='Description', key=lambda x: x.str.lower(), inplace=True)
    
    go_summary = pd.concat([go_down_summary, go_up_summary])
    kegg_summary = pd.concat([kegg_down_summary, kegg_up_summary])
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as w:
        go_summary.to_excel(w, sheet_name='GO_analysis', index=False)
        kegg_summary.to_excel(w, sheet_name='KEGG_analysis', index=False)
    logger.info(f'已汇总到 {output_file}')
    


def main():
    args = parse_input()
    deg_enrich_data_merge(args.input_dir, args.compare_file, args.output_file, 
                          args.filter_col, args.filter_value, args.count_threshold)


if __name__ == '__main__':
    main()
    
    