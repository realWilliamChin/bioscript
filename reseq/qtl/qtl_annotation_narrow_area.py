#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/03/03 10:48
# Author        : William GoGo
"""
BSA R 脚本流程后的分析脚本
"""
import os, sys
import argparse
import pandas as pd
from loguru import logger

import target_gene_from_basicinfo as target_gene

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table
from load_input import write_output_df
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/transcriptome/'))
from genedf_add_expression_and_def import add_kns_def


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('--chr-region', dest='chr_region_file', help='chr_region file')
    p.add_argument('--qtl', dest='qtl_file', help='qtl file')

    p.add_argument('-t', '--threads', default=1, type=int, dest='threads',
                           help='输入进程，默认单个进程处理')
    p.add_argument('-r', '--pos-range', dest='pos_range', type=int, default=10000,
                           help='上下 10000(默认10k) 区间找 target gene')
    p.add_argument('-b', '--basicinfo', type=str, dest='basicinfo', required=True,
                           help='process_gff.py 处理 gff 文件出来的 gene_basic 文件')
    p.add_argument('--kns', type=str, dest='kns',
                           help='输入 kns_def.txt(all_gene_def.txt)，添加定义')
    
    p.add_argument('-o', '--output', type=str, dest='output',
                        help='输出文件名，target_gene_def.txt')

    
    return p.parse_args()


def prep_target_gene_input(chr_region_file, qtl_file):
    chr_region_df = load_table(chr_region_file, usecols=['CHROM', 'qtl', 'start', 'end', 'maxGprime', 'avgDeltaSNP'])
    qtl_df = load_table(qtl_file, usecols=['CHROM', 'POS', 'Gprime'])
    
    # 找到 chr_region 每个区域中(从 qtl_df 中找) Gprime 值最大的位点
    
    each_region_big_point_df_list = []
    
    for i, each_row in chr_region_df.iterrows():
        each_region_df = qtl_df[
            (each_row['start'] <= qtl_df['POS'])
            & (each_row['end'] >= qtl_df['POS']) 
            & (each_row['CHROM'] == qtl_df['CHROM'])
        ]
        biggest_in_region = each_region_df.loc[each_region_df['Gprime'].idxmax()]
        biggest_in_region = pd.DataFrame(biggest_in_region).T
        biggest_in_region['qtl'] = each_row['qtl']
        biggest_in_region['maxGprime'] = each_row['maxGprime']
        biggest_in_region['avgDeltaSNP'] = each_row['avgDeltaSNP']
        
        each_region_big_point_df_list.append(biggest_in_region)
    
    each_region_big_point_df = pd.concat(each_region_big_point_df_list, axis=0)
    
    each_region_big_point_df.drop(columns=['Gprime'], inplace=True)
    
    return each_region_big_point_df
    
    # write_output_df(each_region_big_point_df, output_file, index=False)


def main():
    args = parse_input()
    prep_target_gene_df = prep_target_gene_input(args.chr_region_file, args.qtl_file)
    # prep_target_gene_df.rename(columns={'CHROM': 'SNP_CHROM', 'POS': 'SNP_POS'}, inplace=True)
    
    # 找到 target_gene
    prep_target_gene_df['SNP_MARKER'] = prep_target_gene_df['CHROM'].astype(str) + '_' + prep_target_gene_df['POS'].astype(str)
    target_gene_result_df = target_gene.find_target_gene_multithreads(prep_target_gene_df, args.basicinfo, args.pos_range, args.threads)

    if target_gene_result_df.empty:
        logger.error(f'没有找到 target gene')
        sys.exit(0)
    if args.kns:
        # 添加 kns 定义
        target_gene_result_df.rename(columns={"Target_GeneID": "GeneID"}, inplace=True)
        no_kns_rows = target_gene_result_df.shape[0]
        target_gene_result_df = add_kns_def(target_gene_result_df, kns_file=args.kns)
        add_kns_rows = target_gene_result_df.shape[0]
        logger.debug(f'添加定义之前的行数 {no_kns_rows}, 添加定义之后的行数 {add_kns_rows}')
        target_gene_result_df.rename(columns={"GeneID": "Target_GeneID"}, inplace=True)

    target_gene_result_df = pd.concat([
        target_gene_result_df['SNP_MARKER'],
        target_gene_result_df.iloc[:, target_gene_result_df.columns != 'SNP_MARKER']
    ], axis=1)
    target_gene_result_df.rename(
        columns={
            'CHROM':'SNP_CHROM',
            'POS': 'SNP_POS',
            'avgDeltaSNP':'qtl_avgDeltaSNP',
            'maxGprime': 'qtl_maxGprime'
        },
        inplace=True
    )
    target_gene_result_df['CHROM'] = target_gene_result_df['SNP_CHROM']
    
    # 添加一列 CHROM 从 SNP_CHROM 复制，到 Target_GeneID 后面
    target_idx = target_gene_result_df.columns.get_loc('Target_GeneID')
    new_columns = (
        target_gene_result_df.columns[:target_idx + 1].tolist() +
        ['CHROM'] +
        target_gene_result_df.columns[target_idx + 1:].drop('CHROM').tolist()
    )
    target_gene_result_df = target_gene_result_df[new_columns]
    target_gene_result_df.sort_values(by=['qtl_maxGprime'], ascending=False, inplace=True)
    write_output_df(target_gene_result_df, args.output, index=False)
    
    logger.success(f'处理完成，结果文件为 {args.output}, 结果行数为 {target_gene_result_df.shape[0]}')
    

if __name__ == '__main__':
    main()

