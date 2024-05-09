#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/05/08 14:27
# Author        : William GoGo
import os
import argparse
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--salmon-output', dest='salmon_output', type=str, help='salmon 的输出文件夹')
    
    return parser.parse_args()


def main():
    args = parse_input()
    merged_fpkm = pd.DataFrame()
    merged_reads = pd.DataFrame()
    salmon_output_list = [x for x in os.listdir(args.salmon_output)]
    for sample in salmon_output_list:
        quant_genes_sf = os.sep.join([sample,'quant.genes.sf'])
        # df = pd.read_csv(quant_genes_sf, sep='\t', header=1, names=['GeneID', 'Length', 'EffetiveLength', 'TPM', 'NumReads'])
        df = pd.read_csv(quant_genes_sf, sep='\t', header=1, names=['GeneID', 'Length', 'EffectiveLength', 'TPM', 'NumReads'])
        df['GeneID'] = df['GeneID'].astype(str)
        df['EffectiveLength'] = df['EffectiveLength'].astype(float)
        df['Length'] = df['Length'].astype(float)
        df['NumReads'] = df['NumReads'].astype(int)
        ExonMappedFragments = df['NumReads']
        TotalMappedFragments = df['NumReads'].sum()
        ExonLength = df['EffectiveLength']
        df['FPKM'] = (df['NumReads'] * 10^9) / (TotalMappedFragments * df['EffectiveLength'])
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
    merged_reads_filter = merged_reads[(merged_reads[num_cols] > 50).any(axis=1)]
    merged_fpkm_filter = pd.merge(merged_reads_filter[['GeneID']], merged_fpkm, on='GeneID', how='left')
    merged_reads_filter.to_csv("reads_matrix_filtered.txt", sep='\t', index=False)
    merged_fpkm_filter.to_csv("fpkm_matrix_filtered.txt", sep='\t', index=False)


if __name__ == '__main__':
    main()