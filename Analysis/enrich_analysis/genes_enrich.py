#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/11/19 15:41
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
from Rscript import enrich_analysis

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df

def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--ids-dir', dest='ids_dir', required=True,
                   help='包含单 GeneID 列所有文件的文件夹')
    p.add_argument('--ids-header', dest='ids_header', action='store_true',
                        help='所有 id 文件是否有 header(GeneID)')
    p.add_argument('-k', '--kegg-clean', dest="kegg_clean", required=True,
                   help='kegg_clean.txt 文件')
    p.add_argument('-g', '--gene-go', dest="gene_go", required=True,
                   help='swiss 注释出来的的 gene_go.txt 文件')
    p.add_argument('-o', '--output-dir', dest='output_dir', default=os.getcwd(),
                   help='文件输出目录')
    
    args = p.parse_args()
    
    if args.ids_header:
        args.ids_header = 0
    else:
        args.ids_header = None

    return args


def main():
    args = parse_input()
    ids_dir, ids_header = args.ids_dir, args.ids_header
    for module_file in os.listdir(ids_dir):
        if not os.path.isfile(os.path.join(ids_dir, module_file)):
            continue
        logger.info(f'正在处理 {module_file}')
        module_df = load_table(
            os.path.join(ids_dir, module_file),
            header=ids_header,
            usecols=[0],
            names=['GeneID'],
            dtype={"GeneID": str}
        )
        
        module_name = module_file.split('.')[0]
        tmp_file = os.path.join(args.output_dir, f'{module_name}_ID.txt')
        write_output_df(module_df, tmp_file, index=False, header=False)
        enrich_analysis(tmp_file, args.gene_go, args.kegg_clean, args.output_dir)

    logger.success('Done!')

if __name__ == '__main__':
    main()
    
    