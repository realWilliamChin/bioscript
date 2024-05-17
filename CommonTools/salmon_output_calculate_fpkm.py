#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/05/08 14:27
# Author        : William GoGo
import os
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--salmon-output', dest='salmon_output', type=str, help='salmon 的输出文件夹')
    parser.add_argument('--samples', type=str, help='samples_described.txt，如果不指定则不按照 samples 列的顺序整理表')
    parser.add_argument('--filter-num', dest='filter_num', type=int, default=50, help='reads 每行小于此数值过滤')
    
    return parser.parse_args()


def main():
    args = parse_input()
    
    samples_file, salmon_output, filter_num = args.samples, args.salmon_output, args.filter_num
    merged_fpkm = pd.DataFrame()
    merged_reads = pd.DataFrame()
    if args.samples:
        samples_df = pd.read_csv(samples_file, sep='\t')
        samples = samples_df['sample'].values.tolist()
    else:
        samples = [x for x in os.listdir(salmon_output)]
    for sample in samples:
        logger.info(f"正在处理 {sample}")
        quant_genes_sf = os.sep.join([salmon_output, sample, 'quant.genes.sf'])
        if not os.access(quant_genes_sf, os.R_OK):
            logger.error(f'{quant_genes_sf} 无法访问此文件')
        df = pd.read_csv(quant_genes_sf, sep='\t', header=0, names=['GeneID', 'Length', 'EffectiveLength', 'TPM', 'NumReads'])
        df['GeneID'] = df['GeneID'].astype(str)
        df['EffectiveLength'] = df['EffectiveLength'].astype(float)
        df['Length'] = df['Length'].astype(float)
        df['NumReads'] = df['NumReads'].astype(int)
        # 公式参考网址 https://zhuanlan.zhihu.com/p/688419066
        df['FPKM'] = (df['NumReads'] * (10**9)) / (df['NumReads'].sum() * df['EffectiveLength'])
        df.to_csv(quant_genes_sf.replace('.sf', '_fpkm.txt'), sep='\t', index=False)
        if merged_fpkm.empty:
            merged_fpkm['GeneID'] = df['GeneID']
            merged_fpkm[sample] = df['FPKM']
            merged_reads['GeneID'] = df['GeneID']
            merged_reads[sample] = df['NumReads']
        else:
            merged_fpkm = pd.merge(merged_fpkm, df[['GeneID', "FPKM"]], how='left', on='GeneID')
            merged_fpkm.rename(columns={"FPKM": sample}, inplace=True)

            merged_reads = pd.merge(merged_reads, df[['GeneID', "NumReads"]], how='left', on='GeneID')
            merged_reads.rename(columns={"NumReads": sample}, inplace=True)

    merged_fpkm.to_csv("fpkm_matrix.txt", sep='\t', index=False)
    merged_reads.to_csv("reads_matrix.txt", sep='\t', index=False)
    
    # filter
    num_cols = merged_reads.select_dtypes(include='number').columns
    merged_reads_filter = merged_reads[(merged_reads[num_cols] > filter_num).any(axis=1)]
    merged_fpkm_filter = pd.merge(merged_reads_filter[['GeneID']], merged_fpkm, on='GeneID', how='left')
    merged_reads_filter.to_csv("reads_matrix_filtered.txt", sep='\t', index=False)
    merged_fpkm_filter.to_csv("fpkm_matrix_filtered.txt", sep='\t', index=False)

    logger.info(f'根据 reads 每行至少有一个大于 {filter_num} 过滤前 {merged_fpkm.shape[0]}, 过滤后 {merged_fpkm_filter.shape[0]}')
    logger.info("Done!")

if __name__ == '__main__':
    main()