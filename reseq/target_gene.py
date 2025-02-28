#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/07/04 18:09
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/transcriptome/'))
from genedf_add_expression_and_def import add_kns_def

def parse_input():
    argparser = argparse.ArgumentParser(description="target gene")
    argparser.add_argument('-i', '--input', required=True, type=str, dest='input_file',
                           help='输入要进行处理的表，vcf, csv 格式文件（至少需要 CHROM POS 两列）')
    argparser.add_argument('-t', '--threads', default=1, type=int, dest='threads',
                           help='输入进程，默认单个进程处理')
    argparser.add_argument('-r', '--pos-range', dest='pos_range', type=int, default=10000,
                           help='上下 10000(默认10k) 区间找 target gene')
    argparser.add_argument('-b', '--basicinfo', type=str, dest='basicinfo', required=True,
                           help='process_gff.py 处理 gff 文件出来的 gene_basic 文件')
    argparser.add_argument('--kns', type=str, dest='kns',
                           help='输入 kns_def.txt，添加定义')
    argparser.add_argument('-o', '--output', type=str, dest='output',
                           help='输出文件名，target_gene_def.txt')
    args = argparser.parse_args()
    
    if not os.path.isfile(args.basicinfo) or not os.path.isfile(args.input_file):
        logger.error('输入的文件不存在，请检查')
        exit(1)
    
    return args


# 是否在基因内部
def gene_check_on_off_gene(row):
    if row['Target_Start'] <= row['POS'] <= row['Target_End']:
        return 'On_Gene'
    else:
        return 'Off_Gene'


def utr_check(group):
    # 如果所有组内 On_Gene_Status 都是 Off_Gene
    if all(group['On_Gene_Status'] == 'Off_Gene'):
        conditions = [
            ((group['Target_Start'] - 1500 <= group['POS']) & (group['POS'] <= group['Target_Start'] - 1) & (group['Target_Gene_Strand'] == '+')),
            ((group['Target_End'] + 1 <= group['POS']) & (group['POS'] <= group['Target_End'] + 1500) & (group['Target_Gene_Strand'] == '-')),
            ((group['Target_Start'] - 1500 <= group['POS']) & (group['POS'] <= group['Target_Start'] - 1) & (group['Target_Gene_Strand'] == '-')),
            ((group['Target_End'] + 1 <= group['POS']) & (group['POS'] <= group['Target_End'] + 1500) & (group['Target_Gene_Strand'] == '+'))
        ]
        utr_labels = ['UTR_5_Prime', 'UTR_5_Prime', 'UTR_3_Prime', 'UTR_3_Prime']

        for condition, label in zip(conditions, utr_labels):
            group.loc[condition, 'On_Gene_Status'] = label

    return group


def process_sub_dataframe(sub_df, basicinfo_df):
    # 2. 如果 input_df 的 start 在 basicinfo 的 start 和 end 之间
    result_df_lst = []
    sub_df_columns = sub_df.columns.tolist()
    # eachrow_info_df = pd.DataFrame()
    for each_row in sub_df.itertuples():
        # merge input_df pos 上下 10k 在 basicinfo 的 start 和 end 之间的数据
        eachrow_df = basicinfo_df[
            ((each_row.start <= basicinfo_df['Target_Start']) & (each_row.end >= basicinfo_df['Target_Start']))
            | ((each_row.start <= basicinfo_df['Target_End']) & (each_row.end >= basicinfo_df['Target_End']))
            | ((basicinfo_df['Target_Start'] <= each_row.POS) & (basicinfo_df['Target_End'] >= each_row.POS))
        ].copy()
        # eachrow_df = eachrow_df.drop_duplicates(subset='Target_GeneID', keep='first')
        # 追加其他所有信息
        for column in sub_df_columns:
            if column in ['start', 'end']:
                continue
            eachrow_df[column] = getattr(each_row, column)
        # eachrow_df = eachrow_df.assign(POS=each_row.POS, REF=each_row.REF, ALT=each_row.ALT, CHROM=each_row.CHROM)
        result_df_lst.append(eachrow_df)
    
    # 3. 判断 pos 是否在基因内部，和 utr check
    result_df = pd.concat(result_df_lst)
    
    # 如何找不到 target gene 则返回空的 df
    if result_df.empty:
        return result_df
    
    result_df['On_Gene_Status'] = result_df.apply(gene_check_on_off_gene, axis=1)
    result_df = result_df.groupby('POS').apply(utr_check)
    return result_df


def find_target_gene_multithreads(input_df, gene_basic_info, pos_range, num_threads):
    """
    多进程处理
    input_df: 输入的表, 至少需要包含 POS 列
    """
    basicinfo_df = pd.read_csv(gene_basic_info, sep='\t', usecols=[0, 2, 3, 4], skiprows=1,
                         names=["Target_GeneID", "Target_Start", "Target_End", "Target_Gene_Strand"],
                         dtype={"Target_GeneID": str, "Target_Start": int, "Target_End": int, "Target_Gene_Strand": str})
    # 标出 pos 位置的上下 pos_range 默认10k，再去比对 gff start 和 end 判断是否包含在内，包含在内则把 basicinfo_df 的信息添加到 input_df 中
    # 1. 标出 pos 位置的上下 pos_range，命名为 start 和 end
    input_df['start'] = input_df['POS'].apply(lambda x: max(x - pos_range, 0))
    input_df['end'] = input_df['POS'] + pos_range
    if num_threads == 1:
        result_df = process_sub_dataframe(input_df, basicinfo_df)
    else:
        # 使用多进程进行处理
        sub_dfs = np.array_split(input_df, num_threads)

        with ProcessPoolExecutor(max_workers=num_threads) as executor:
            futures = []
            for sub_df in sub_dfs:
                sub_executor = executor.submit(process_sub_dataframe, sub_df, basicinfo_df)
                futures.append(sub_executor)
            
            # futures = [executor.submit(process_sub_dataframe, sub_df, basicinfo_df) for sub_df in sub_dfs]
            results = [future.result() for future in futures if future.result().shape[0] > 0]
            result_df = pd.concat(results)
        if len(results) == 0:
            return pd.DataFrame()
    
    sources_columns = input_df.columns.tolist()
    sources_columns.remove('start')
    sources_columns.remove('end')
    re_columns = sources_columns + ['Target_GeneID', 'Target_Start', 'Target_End', 'Target_Gene_Strand', 'On_Gene_Status']
    result_df = result_df.reindex(columns=re_columns)
    result_df = result_df.sort_values(by=['POS'])
    
    return result_df


def main():
    args = parse_input()
    
    if args.input_file.endswith('.vcf'):
        # vcf 文件读取前处理
        skip_rows = 0
        with open(args.input_file, "r") as file:
            for line in file:
                if line.startswith('#'):
                    skip_rows += 1
                else:
                    break
        logger.info(f'跳过 {args.input_file} 的前 {skip_rows} 行')
        
        vcf_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        input_df = pd.read_csv(args.input_file, sep='\t', skiprows=skip_rows, usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], 
                               low_memory=False, names=vcf_columns,
                               dtype={"CHROM": str, "POS": int, "ID": str, "REF": str, "ALT": str, "QUAL": str, "FILTER": str, "INFO": str, "FORMAT": str})
    
        input_df['Marker'] = input_df['CHROM'].astype(str) + '_' + input_df['POS'].astype(str)
        df = find_target_gene_multithreads(input_df, args.basicinfo, args.pos_range, args.threads)
        
        if df.empty:
            logger.error(f'没有找到 target gene')
            exit(0)
        
        # 添加 kns 定义
        df.rename(columns={"Target_GeneID": "GeneID"}, inplace=True)
        no_kns_rows = df.shape[0]
        df = add_kns_def(df, kns_file=args.kns)
        add_kns_rows = df.shape[0]
        logger.debug(f'添加定义之前的行数 {no_kns_rows}, 添加定义之后的行数 {add_kns_rows}')
        df.rename(columns={"GeneID": "Target_GeneID"}, inplace=True)
        
        df = pd.concat([df['Marker'], df.iloc[:, df.columns != 'Marker']], axis=1)
        # df_columns = df.columns.tolist() # 修改为上面那一句代码
        # df_columns = [x for x in df_columns if x not in vcf_columns + ['Marker']]
        # df = df.reindex(columns=['Marker', 'ID', 'CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + df_columns)
        
        df.to_csv(args.output, sep='\t', index=False)
        
    elif args.input_file.endswith('.csv'):
        input_df = pd.read_csv(args.input_file)
        input_df['Marker'] = input_df['CHROM'].astype(str) + '_' + input_df['POS'].astype(str)
        df = find_target_gene_multithreads(input_df, args.basicinfo, args.pos_range, args.threads)
        
        df.rename(columns={"Target_GeneID": "GeneID"}, inplace=True)
        no_kns_rows = df.shape[0]
        df = add_kns_def(df, kns_file=args.kns)
        add_kns_rows = df.shape[0]
        logger.debug(f'添加定义之前的行数 {no_kns_rows}, 添加定义之后的行数 {add_kns_rows}')
        df.rename(columns={"GeneID": "Target_GeneID"}, inplace=True)
        
        df = pd.concat([df['Marker'], df.iloc[:, df.columns != 'Marker']], axis=1)
        
        df.to_csv(args.output, sep='\t', index=False)
        if df.empty:
            logger.error(f'没有找到 target gene')
            exit(0)
    else:
        # 其他格式文件未做优化，可能会出现 bug
        logger.error('不支持其他文件格式')
        exit(1)

    logger.success(f'处理完成，结果文件为 {args.output}, 结果行数为 {df.shape[0]}')
    
if __name__ == '__main__':
    main()