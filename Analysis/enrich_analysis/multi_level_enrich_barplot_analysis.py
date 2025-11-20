#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2025/11/03 11:10
# Author        : WilliamGoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df
sys.path.append('/home/colddata/qinqiang/script/transcriptome')
from genedf_add_expression_and_def import add_def


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', required=True, help='输入需要分析的 GO_ID 或者 KEGG_ID，两列 ID 和 Description')
    p.add_argument('-e', '--enrich', nargs='+', required=True,
                   help='输入所有需要分析的组 enrch 文件（可输入多个文件，空格分隔）')
    p.add_argument('-o', '--output', help='输出文件', default='summary.txt')
    
    args = p.parse_args()
    
    return args


def main():
    args = parse_input()
    
    input_df = load_table(args.input)
    input_df = input_df[['ID', 'Description']]
    
    df_list = []
    for f in args.enrich:
        # enrich_df = load_table(f)
        base_name = os.path.basename(f)
        comparison_name = base_name.split('_Up')[0] if '_Up' in base_name else base_name.split('_Down')[0]
        comparison_regulation = base_name.rsplit('_Enrichment', 1)[0].split('_')[-1]
        
        logger.info(f'Processing {f} - {comparison_name} - {comparison_regulation}')
        specific_enrich_df = add_def(
            file_df = input_df, 
            index_column = 'ID',
            def_file = f,
            merge_how = 'left',
            remove_duplicates = True,
            dup_col = 'left'
        )
        
        specific_enrich_df['Regulation'] = comparison_regulation
        specific_enrich_df['Comparison'] = comparison_name
        df_list.append(specific_enrich_df)
    
    summary_df = pd.concat(df_list, axis=0)
    # 避免 replace 的 downcasting 警告：使用 pd.to_numeric 处理数值转换
    summary_df['Count'] = summary_df['Count'].astype(str).replace('N/A', '0').fillna('0')
    summary_df['Count'] = pd.to_numeric(summary_df['Count'], errors='coerce').fillna(0).astype(int)
    summary_df['p.adjust'] = summary_df['p.adjust'].astype(str).replace('N/A', '1').fillna('1')
    summary_df['p.adjust'] = pd.to_numeric(summary_df['p.adjust'], errors='coerce').fillna(1).astype(float)
    write_output_df(summary_df, args.output, index=False)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()