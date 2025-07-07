#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/11/19 15:41
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/Rscript/')
from Rscript import enrich_analysis
sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--compare', default='compare_info.txt',
                        help='compare_info.txt 文件，主要用来读取组间名称')
    parser.add_argument('--degdata-dir', dest='degdata_dir', default='DEG_analysis_results',
                        help='提供 compareinfo up 和 down gene id 的目录，通常是./DEG_analysis_results')
    parser.add_argument('-k', '--keggclean', help='kegg_clean.txt 文件')
    parser.add_argument('-g', '--genego', help='kegg_clean.txt 文件')
    parser.add_argument('-o', '--outputdir', default='Pathway_enrichment_analysis')
    
    args = parser.parse_args()
    
    os.makedirs(args.outputdir, exist_ok=True)
    
    return args


def deg_enrich(compare, degdata_dir, genego_file, keggclean_file, outputdir):
    compare_df = load_table(compare)
    compare_df['compare_name'] = compare_df['Treat'] + '-vs-' + compare_df['Control']
    for comp in compare_df['compare_name'].tolist():
        deg_Upid_file = os.path.join(degdata_dir, f'{comp}_Up_ID.txt')
        deg_Downid_file = os.path.join(degdata_dir, f'{comp}_Down_ID.txt')

        if os.path.exists(deg_Upid_file):
            logger.info(f'进行 {deg_Upid_file} 的 enrich 分析')
            enrich_analysis(deg_Upid_file, genego_file, keggclean_file, outputdir)
        if os.path.exists(deg_Downid_file):
            logger.info(f'进行 {deg_Downid_file} 的 enrich 分析')
            enrich_analysis(deg_Downid_file, genego_file, keggclean_file, outputdir)


def main():
    args = parse_input()
    deg_enrich(args.compare, args.degdata_dir, args.genego, args.keggclean, args.outputdir)


if __name__ == '__main__':
    main()