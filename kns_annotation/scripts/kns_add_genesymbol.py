#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2025/04/11 11:43
# Author        : WilliamGoGo
import os, sys
import pandas as pd
import argparse
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df


def parse_arguments():
    p = argparse.ArgumentParser(description='kns 添加 gene symbol，去掉 KEGG_shortname 列')
    p.add_argument('--kns', help='输入文件')
    p.add_argument('--symbol', help='输入文件')
    p.add_argument('-o', '--output', help='输出文件')
    args = p.parse_args()
    return args


def main():
    args = parse_arguments()
    kns_df = load_table(args.kns, dtype=str)
    symbol_df = load_table(args.symbol, dtype=str)
    
    # 记录合并前的行数
    before_count = len(kns_df)
    logger.info(f"合并前数据行数: {before_count}")
    
    kns_df = kns_df.merge(symbol_df, on='GeneID', how='left')
    
    # 记录去重前的行数
    before_dedup = len(kns_df)
    logger.info(f"去重前数据行数: {before_dedup}")
    
    kns_df.drop_duplicates(subset=['GeneID'], inplace=True)
    
    # 记录去重后的行数
    after_dedup = len(kns_df)
    logger.info(f"去重后数据行数: {after_dedup}")
    logger.info(f"去重过滤掉的行数: {before_dedup - after_dedup}")
    
    kns_df.fillna('NA', inplace=True)
    write_output_df(kns_df, args.output, index=False)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()

