#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/11/19 15:41
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
from Rscript import enrich_analysis

def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--compare', default='compare_info.txt', type=str,
                        help='compare_info.txt 文件，主要用来读取组间名称')
    parser.add_argument('--degdata-dir', default=os.curdir(), dest='degdata_dir',
                        help='提供 compareinfo up 和 down gene id 的目录，通常是./DEG_analysis_results')
    parser.add_argument('-k', '--keggclean', help='kegg_clean.txt 文件')
    parser.add_argument('-g', '--genego', help='kegg_clean.txt 文件')
    parser.add_argument('-o', '--outputdir')


def main():
    args = parse_input()
    compare_df, degdata_dir, genego_file, keggclean_file, outputdir = args.compare, args.degdata_dir, args.genego, args.keggclean, args.outputdir
    compare_df['compare_name'] = compare_df['Treat'] + '_vs_' + compare_df['Control']
    for comp in compare_df['compare_name'].tolist():
        deg_Upid_file = os.path.join(degdata_dir, f'{comp}_Up_ID.txt')
        deg_Downid_file = os.path.join(degdata_dir, f'{comp}_Down_ID.txt')

        if os.path.exists(deg_Upid_file):
            enrich_analysis(deg_Upid_file, genego_file, keggclean_file, outputdir)
        if os.path.exists(deg_Downid_file):
            enrich_analysis(deg_Downid_file, genego_file, keggclean_file, outputdir)


if __name__ == '__main__':
    main()