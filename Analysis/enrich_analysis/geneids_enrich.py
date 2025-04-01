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
    parser.add_argument('-i', '--clusterid', help='提供 compareinfo up 和 down gene id 的目录，通常是./DEG_analysis_results')
    parser.add_argument('-k', '--keggclean', help='kegg_clean.txt 文件')
    parser.add_argument('-g', '--genego', help='kegg_clean.txt 文件')
    parser.add_argument('-o', '--outputdir', help='文件输出目录')

    return parser.parse_args()


def main():
    args = parse_input()
    cluster_id_file = args.clusterid
    cluster_id_df = pd.read_csv(cluster_id_file, sep='\t', usecols=[0])
    
    tmp_file = cluster_id_file.replace('.txt', '_ID.txt')
    cluster_id_df.to_csv(tmp_file, sep='\t', index=False, header=False)
    enrich_analysis(tmp_file, args.genego, args.keggclean, args.outputdir)


if __name__ == '__main__':
    main()
    
    