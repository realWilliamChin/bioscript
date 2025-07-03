#!/usr/bin/env python
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    args = argparse.ArgumentParser()
    args.add_argument('-i', '--input', help='输入 kns_def 文件, 默认是当前文件夹下第一个 kns_def.txt')
    args.add_argument('-g', '--genesymbol', help='输入 gene_symbol 文件, 默认是当前文件夹下第一个 gene_symbol.txt')
    args.add_argument('-o', '--output', help='输出文件名称，默认是 args.input')
    
    args = args.parse_args()
    
    if not args.output:
        args.output = args.input
    
    return args


# def shortname2genesymbol(kns_def_file, genesymbol_file):
#     columns = ['GeneID', 'NR_ID', 'NR_ID_Des', 'Swiss_ID', 'Swiss_ID_Des', 'KEGG_ID', 'GeneSymbol', 'EC_Number', 'KEGG_Description']
#     kns_df = pd.read_csv(kns_def_file, sep='\t', skiprows=1, names=columns)
#     kns_df = kns_df.drop(columns=['GeneSymbol'])
#     genesymbol_df = pd.read_csv(genesymbol_file, sep='\t')
    
#     # kns_df 的 KEGG_shortname 列替换为 genesymbol_df 的 GeneSymbol 列
#     kns_df = pd.merge(left=kns_df, right=genesymbol_df, on='GeneID', how='left')
#     # 对 kns_df 的列排序，按照 columns 的顺序
#     kns_df = kns_df[columns]
#     kns_df = kns_df.fillna('NA')
#     kns_df = kns_df.drop_duplicates(subset='GeneID', keep='first')
#     kns_df.to_csv(kns_def_file, sep='\t', index=False)

def kns_add_symbol(kns_def_file, genesymbol_file):
    kns_df = load_table(kns_def_file, dtype={"GeneID": str})
    genesymbol_df = load_table(genesymbol_file, dtype={"GeneID": str})

    result_df = pd.merge(kns_df, genesymbol_df, on='GeneID', how='left')
    source_shape = result_df.shape[0]
    result_df.drop_duplicates(subset='GeneID', inplace=True)
    filtered_shape = result_df.shape[0]
    
    if source_shape != filtered_shape:
        logger.warning(f'数据表合并之后有重复，已进行去重处理，请检查 genesymbol 或 kns 文件是否有重复')
        logger.warning(f'重复数量是 {source_shape - filtered_shape}')
    result_df.fillna(value='N/A', inplace=True)

    return result_df


def main():
    args = parse_input()
    # shortname2genesymbol(args.input, args.genesymbol)
    result_df = kns_add_symbol(args.input, args.genesymbol)
    write_output_df(result_df, args.output, index=False)
    
    logger.success(f'Done, 已输出到文件 {args.output}')


if __name__ == '__main__':
    main()