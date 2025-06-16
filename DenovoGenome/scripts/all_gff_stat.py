#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/07/01 14:23
# Author        : William GoGo
# Description   : 合并 braker.gff3 rrna.gff3 trna.gff3 还有 others.gff3 并计算统计信息
import os, sys
import argparse
import pandas as pd
import openpyxl
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(description="合并 braker.gff3 rrna.gff3 trna.gff3 还有 others.gff3 并计算统计信息")
    parser.add_argument("-b", "--braker", type=str, required=True, help="braker.gff3 文件")
    parser.add_argument("-r", "--rrna", type=str, required=True, help="rrna.gff3 文件")
    parser.add_argument("-t", "--trna", type=str, required=True, help="trna.gff3 文件")
    parser.add_argument("-o", "--others", type=str, required=True, help="others.gff3 文件")
    parser.add_argument("-p", "--prefix", type=str, required=True, help="输出文件前缀，可以加上目录")
    return parser.parse_args()


def main():
    args = parse_input()
    braker_file, rrna_file, trna_file, others_file = args.braker, args.rrna, args.trna, args.others
    
    braker_df = load_table(braker_file, header=None)
    rrna_df = load_table(rrna_file, header=None)
    trna_df = load_table(trna_file, header=None)
    others_df = load_table(others_file, header=None)
    
    all_df = pd.concat([braker_df, rrna_df, trna_df, others_df], ignore_index=True)
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    all_df.columns = gff3_columns
    write_output_df(all_df, f"{args.prefix}_all.gff3", index=False, header=True)
    
    # 统计每种类型的 gene id 个数，出一个统计表
    braker_id_counts = braker_df[braker_df[2] == 'gene'].shape[0]
    rrna_id_counts = rrna_df.shape[0]
    trna_id_counts = trna_df[trna_df[2] == 'gene'].shape[0]
    others_id_counts = others_df.shape[0]
    
    
    stat_results = pd.DataFrame({
        'type': ['braker', 'rrna', 'trna', 'others'],
        'id_counts': [braker_id_counts, rrna_id_counts, trna_id_counts, others_id_counts]
    })
    
    write_output_df(stat_results, f"{args.prefix}_gff_gene_type_counts.txt", index=False)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()