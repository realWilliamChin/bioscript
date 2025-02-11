#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/02/11 15:49
# Author        : William GoGo
import os
import argparse
import pandas as pd
import numpy as np


def step1_commongene(group_name, common_gene_file, output_dir):
    common_df = pd.read_csv(common_gene_file, sep='\t')
    
    gene_df = pd.read_csv(f"{group_name}_DEG_gene.txt", sep='\t', usecols=['GeneID', 'log2FoldChange', 'pvalue'])
    gene_df.columns = ['GeneID', 'log2FC_RNA', 'pvalue_RNA']
    
    result_df = pd.merge(common_df, gene_df, how='left', on='GeneID')
    protein_df = pd.read_csv(f"{group_name}_DEG_protein.txt", sep='\t', usecols=['GeneID', 'log2FoldChange', 'pvalue'])
    protein_df.columns = ['GeneID', 'log2FC_protein', 'pvalue_protein']
    
    result_df = pd.merge(result_df, protein_df, how='left', on='GeneID')
    result_df.fillna(0, inplace=True)
    result_df.to_csv(f"{output_dir}/{group_name}_common_result.txt", sep='\t', index=False)


def step2_group(group_name, output_dir, genesymbol_file=None):
    rna_value = 0.26
    protein_value = 0.26
    
    df = pd.read_csv(f'{output_dir}/{group_name}_common_result.txt', sep='\t')
    df['log2FC_protein'] = df['log2FC_protein'].apply(lambda x: x if x <= 9 else 9)
    df['log2FC_protein'] = df['log2FC_protein'].apply(lambda x: x if x >= -9 else -9)
    df['log2FC_RNA'] = df['log2FC_RNA'].apply(lambda x: x if x <= 9 else 9)
    df['log2FC_RNA'] = df['log2FC_RNA'].apply(lambda x: x if x >= -9 else -9)

    group1 = df[(df['log2FC_protein'] > protein_value) & (df['log2FC_RNA'] < -rna_value)].copy()
    group2 = df[(df['log2FC_protein'] > protein_value) & (df['log2FC_RNA'] > -rna_value) & (df['log2FC_RNA'] < rna_value)].copy()
    group3 = df[(df['log2FC_protein'] > protein_value) & (df['log2FC_RNA'] > rna_value)].copy()
    group4 = df[(df['log2FC_protein'] < protein_value) & (df['log2FC_protein'] > -protein_value) & (df['log2FC_RNA'] < -rna_value)].copy()
    group5 = df[(df['log2FC_protein'] < protein_value) & (df['log2FC_protein'] > -protein_value) & (df['log2FC_RNA'] > -rna_value) & (df['log2FC_RNA'] < rna_value)].copy()
    group6 = df[(df['log2FC_protein'] < protein_value) & (df['log2FC_protein'] > -protein_value) & (df['log2FC_RNA'] > rna_value)].copy()
    group7 = df[(df['log2FC_protein'] < -protein_value) & (df['log2FC_RNA'] < -rna_value)].copy()
    group8 = df[(df['log2FC_protein'] < -protein_value) & (df['log2FC_RNA'] > -rna_value) & (df['log2FC_RNA'] < rna_value)].copy()
    group9 = df[(df['log2FC_protein'] < -protein_value) & (df['log2FC_RNA'] > rna_value)].copy()
    
    group1['Group'] = 'group1'
    group2['Group'] = 'group2'
    group3['Group'] = 'group3'
    group4['Group'] = 'group4'
    group5['Group'] = 'group5'
    group6['Group'] = 'group6'
    group7['Group'] = 'group7'
    group8['Group'] = 'group8'
    group9['Group'] = 'group9'
    
    group3['regulation'] = 'Up'
    group7['regulation'] = 'Down'
    
    if genesymbol_file:
        group1 = add_genesymbol(group1, genesymbol_file)
        group3 = add_genesymbol(group3, genesymbol_file)
        group7 = add_genesymbol(group7, genesymbol_file)
        group9 = add_genesymbol(group9, genesymbol_file)
        
        group1.loc[group1['log2FC_RNA'] > -2, 'GeneSymbol'] = np.nan
        group1.loc[group1['log2FC_protein'] < 2, 'GeneSymbol'] = np.nan
        group3.loc[group3['log2FC_RNA'] < 2, 'GeneSymbol'] = np.nan
        group3.loc[group3['log2FC_protein'] < 2, 'GeneSymbol'] = np.nan
        group7.loc[group7['log2FC_RNA'] > -2, 'GeneSymbol'] = np.nan
        group7.loc[group7['log2FC_protein'] > -2, 'GeneSymbol'] = np.nan
        group9.loc[group9['log2FC_RNA'] < 2, 'GeneSymbol'] = np.nan
        group9.loc[group9['log2FC_protein'] > -2, 'GeneSymbol'] = np.nan
        
    
    result_df = pd.concat([group1, group2, group3, group4, group5, group6, group7, group8, group9])
    # result_df.fillna(0, inplace=True)
    result_df.to_csv(f"{output_dir}/{group_name}_result_group.txt", sep='\t', index=False)
    
    group1.drop(columns=['Group']).to_csv(f"{output_dir}/{group_name}_result_group1.txt", sep='\t', index=False)
    group2.drop(columns=['Group']).to_csv(f"{output_dir}/{group_name}_result_group2.txt", sep='\t', index=False)
    group3.drop(columns=['Group']).to_csv(f"{output_dir}/{group_name}_result_group3.txt", sep='\t', index=False)
    group4.drop(columns=['Group']).to_csv(f"{output_dir}/{group_name}_result_group4.txt", sep='\t', index=False)
    group5.drop(columns=['Group']).to_csv(f"{output_dir}/{group_name}_result_group5.txt", sep='\t', index=False)
    group6.drop(columns=['Group']).to_csv(f"{output_dir}/{group_name}_result_group6.txt", sep='\t', index=False)
    group7.drop(columns=['Group']).to_csv(f"{output_dir}/{group_name}_result_group7.txt", sep='\t', index=False)
    group8.drop(columns=['Group']).to_csv(f"{output_dir}/{group_name}_result_group8.txt", sep='\t', index=False)
    group9.drop(columns=['Group']).to_csv(f"{output_dir}/{group_name}_result_group9.txt", sep='\t', index=False)


def add_genesymbol(df, genesymbol_file):
    genesymbol_df = pd.read_csv(genesymbol_file, sep='\t')
    
    df = pd.merge(df, genesymbol_df, how='left', on='GeneID')
    
    return df
    

def parse_input():
    p = argparse.ArgumentParser()
    # p.add_argument('-i', '--input', help='包含 DEG_gene 和 DEG_protein 的目录')
    p.add_argument('-o', '--output', default=os.getcwd(), help='输入输出目录')
    p.add_argument('--common', '--common-gene-list', help='Common_gene_list.txt 文件，输入文件')
    
    p.add_argument('--genesymbel', default=None, help='是否添加 genesymbol 列，画图的时候会')
    
    return p.parse_args()


def main():
    args = parse_input()
    # input_dir = args.input
    
    
    for each_file in os.listdir():
        if each_file.endswith('_DEG_gene.txt'):
            group_name = each_file.split('_DEG_gene.txt')[0]
            output_dir = os.path.join(args.output, f"{group_name}_result")
            os.makedirs(output_dir, exist_ok=True)
            step1_commongene(group_name, args.common, output_dir)
            step2_group(group_name, output_dir, args.genesymbel)


if __name__ == '__main__':
    main()
