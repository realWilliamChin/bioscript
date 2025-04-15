#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/04/15 11:03
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger
sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df

def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='输入 Species.txt 文件')
    p.add_argument('-o', '--output', default='Species_taxon.txt',
                   help='输出文件默认 Species_taxon.txt')
    
    return p.parse_args()


def split_taxonomy(taxonomy_str):
    """将分类学字符串拆分成各个分类级别，识别 k__, o__, s__ 等前缀"""
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    prefixes = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    parts = taxonomy_str.split('|')
    result = {}
    
    # 初始化所有级别为默认值
    for level, prefix in zip(levels, prefixes):
        result[level] = prefix
    
    # 处理每个部分
    for part in parts:
        for level, prefix in zip(levels, prefixes):
            if part.startswith(prefix):
                result[level] = part
                break
    
    return result


def main():
    args = parse_input()
    input_file, output_file = args.input, args.output
    df = load_table(input_file)
    # 生成一列 OTU，OTU_ID 为 OTU_1, OTU_2, OTU_3, ...
    df['OTU'] = df.index.map(lambda x: f'OTU_{x}')
    species_name_df = df[['OTU', 'Species']]
    df.set_index('OTU', inplace=True)
    species_value_df = df.drop(columns=['Species'])
    write_output_df(species_value_df, 'Species_aboun.txt', index_label='ID')
    
    taxonomy_df = pd.DataFrame(species_name_df['Species'].apply(split_taxonomy).tolist())
    taxonomy_df['ID'] = species_name_df['OTU']
    taxonomy_df.set_index('ID', inplace=True)
    write_output_df(taxonomy_df, args.output, index_label='ID')
    logger.info(f"处理完成，结果已保存至 {output_file}")


if __name__ == '__main__':
    main()