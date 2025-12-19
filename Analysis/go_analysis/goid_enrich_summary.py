#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/12/18 10:43
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger
import datetime

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df
from data_check import df_drop_element_side_space
from data_check import df_replace_illegal_folder_chars
from logger_config import get_logger

if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    sys.exit(1)


def parse_input():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-i', '--input', help="输入文件，文件中至少包含三列，GO_pathway_ID、Ontology、SubOntology")
    argparser.add_argument('-c', '--compare', help='输入 compare_info.txt 文件')
    argparser.add_argument('-g', '--genego', help='swiss 注释出来的的 gene_go.txt')
    argparser.add_argument('-s', '--samples-described', help='输入 samples_described.txt 文件')
    argparser.add_argument('-e', '--enrich-data-dir', help='转录组 输入 enrich.r 运行出来的结果文件夹')
    argparser.add_argument('-o', '--output', default='00_Target_GO_Enrich_data_summary.xlsx', help='输出目录，默认为当前目录，输出目录如果不存在则会尝试创建')
    
    args = argparser.parse_args()
    
    
    # Data check
    for f in [args.input, args.genego]:
        if not os.access(f, os.R_OK):
            logger.critical(f"输入文件 {f} 未找到或不可读")
            sys.exit(1)
    
    return args


def goid_enrich_summary(target_go_df: pd.DataFrame, comparisons_df: pd.DataFrame, enrich_data_dir: str, output: str) -> None:
    target_go_list = target_go_df['GO_pathway_ID'].values.tolist()
    target_go_enrich_df_list = []
    comparisons_df['comparisons'] = comparisons_df['Treat'] + '-vs-' + comparisons_df['Control']
    comparison_list = comparisons_df['comparisons'].values.tolist()
    for comparison in comparison_list:
        up_enrich_df = load_table(os.path.join(enrich_data_dir, f'{comparison}_Up_EnrichmentGO.xlsx'))
        down_enrich_df = load_table(os.path.join(enrich_data_dir, f'{comparison}_Down_EnrichmentGO.xlsx'))
        target_go_enrich_up_df = up_enrich_df[up_enrich_df['ID'].isin(target_go_list)]
        target_go_enrich_down_df = down_enrich_df[down_enrich_df['ID'].isin(target_go_list)]
        
        target_go_enrich_up_df.insert(0, 'Group', comparison)
        target_go_enrich_up_df.insert(1, 'Regulation', 'Up')
        target_go_enrich_down_df.insert(0, 'Group', comparison)
        target_go_enrich_down_df.insert(1, 'Regulation', 'Down')
        
        target_go_enrich_df_list.append(target_go_enrich_up_df)
        target_go_enrich_df_list.append(target_go_enrich_down_df)
    
    target_go_enrich_summary_df = pd.concat(target_go_enrich_df_list)
    target_go_enrich_summary_df.drop(columns=['Ontology'], inplace=True)
    target_go_enrich_summary_df.rename(columns={'ID': 'GO_pathway_ID'}, inplace=True)
    if 'Ontology' in target_go_df.columns and 'SubOntology' in target_go_df.columns:
        target_go_enrich_summary_df = pd.merge(
            target_go_df[['GO_pathway_ID', 'Ontology', 'SubOntology']],
            target_go_enrich_summary_df,
            how='right',
            on='GO_pathway_ID'
        )
    target_go_enrich_summary_df.sort_values(
        by=['Ontology', 'SubOntology', 'Group', 'Regulation'],
        ascending=[True, True, True, False],
        inplace=True
    )
    write_output_df(target_go_enrich_summary_df, output, index=False)


if __name__ == '__main__':
    args = parse_input()
    logger.info(f'正在对输入数据进行预处理，去除两边空格，替换 Ontology 非法字符')
    # target_go_df 预处理
    target_go_df = load_table(args.input)
    if 'GO_def' not in target_go_df.columns:
        target_go_df['GO_def'] = target_go_df['GO_ID'].str.split('_').str[1]
    target_go_df['GO_ID'] = target_go_df['GO_ID'].str.split('_').str[0]
    target_go_df = df_drop_element_side_space(target_go_df)
    target_go_df = df_replace_illegal_folder_chars(target_go_df, ['GO_def', 'Ontology', 'SubOntology'])
    target_go_df['GO_ID'] = target_go_df['GO_ID'].str.split("_").str[0]  # 为了和 enrich.r 出来的文件的 ID 对应上，添加一列只包含 ID 的列， GO:0010111
    
    comparison_df = load_table(args.compare)
    ontology_list = list(set(target_go_df['Ontology'].tolist()))
    go_id_list = list(set(target_go_df['GO_ID'].tolist()))
    gene_go_df = load_table(args.genego, header=None, names=['GeneID', 'GO_ID'], dtype={"GeneID": str})
    
    logger.info(f'执行 go 分析')
    goid_enrich_summary(target_go_df, comparison_df, args.enrich_data_dir, args.output)
    
    logger.success('Done!')