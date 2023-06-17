# -*- coding: UTF-8 -*-
# Created Time  : 2023/5/29 19:43
# Author        : WilliamGoGo
# 合并 fpkm 和 reads 然后添加 def
import argparse
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser(description='合并 fpkm 和 reads matrix, 然后加上三个注释的定义')
    parser.add_argument('-f', '--fpkm', type=str, required=True, help='指定 fpkm 文件')
    parser.add_argument('-r', '--reads', type=str, required=True, help='指定 reads 文件')
    parser.add_argument('-k', '--kegg', type=str, required=True, help='指定 kegg_gene_def file')
    parser.add_argument('-n', '--nr', type=str, required=True, help='指定 nr_gene_def file')
    parser.add_argument('-s', '--swiss', type=str, required=True, help='指定 swiss_gene_def file')
    parser.add_argument('-o', '--output', type=str, help='指定 output 文件名，不指定默认 fpkm_and_reads_matrix_filtered_data_def.txt')
    args = parser.parse_args()
    return args.fpkm, args.reads, args.kegg, args.nr, args.swiss, args.output


def merge_fpkm_reads(fpkm_file, reads_file):
    fpkm_df = pd.read_csv(fpkm_file, sep='\t')
    reads_df = pd.read_csv(reads_file, sep='\t')
    # 对 GeneID 列进行重命名，如果是用其他方式写的 gene_id geneid 等等
    if 'gene' in fpkm_df.columns[0].lower() and 'id' in fpkm_df.columns[0].lower():
        fpkm_df.rename(columns={fpkm_df.columns[0]: 'GeneID'}, inplace=True)
    if 'gene' in reads_df.columns[0].lower() and 'id' in reads_df.columns[0].lower():
        reads_df.rename(columns={reads_df.columns[0]: 'GeneID'}, inplace=True)

    # 列名添加 _fpkm _reads，先把列名提取出来加上，然后给表换上新的列名
    fpkm_column_lst = list(each + '_fpkm' for each in list(fpkm_df.columns)[1:])
    reads_column_lst = list(each + '_reads' for each in list(reads_df.columns)[1:])
    fpkm_column_lst.insert(0, fpkm_df.columns[0])
    reads_column_lst.insert(0, fpkm_column_lst[0])
    fpkm_df.columns = fpkm_column_lst
    reads_df.columns = reads_column_lst
    
    # 合并
    reads_fpkm_df = pd.merge(left=fpkm_df, right=reads_df, on=fpkm_df.columns[0], how='left')
    return reads_fpkm_df


def add_kns_def(df, kegg_file, nr_file, swiss_file):
    """
    添加 kegg nr swiss def
    """
    kegg_df = pd.read_csv(kegg_file, sep='\t', skiprows=1, usecols=[0, 1], names=['GeneID', 'KEGG_ID'])
    nr_df = pd.read_csv(nr_file, sep='\t', skiprows=1, usecols=[0, 1, 2], names=['GeneID', 'NR_ID', 'NR_Def1'])
    # 没有 NCBI ID 直接加会丢失先 fillna，再去掉 NA::
    nr_df.fillna(value='NA', inplace=True)
    nr_df['NR_Def'] = nr_df['NR_ID'] + '::' + nr_df['NR_Def1']
    nr_df['NR_Def'] = nr_df['NR_Def'].str.replace('NA::', '', regex=False)
    nr_df.drop(columns=['NR_ID', 'NR_Def1'], inplace=True)
    swiss_df = pd.read_csv(swiss_file, sep='\t', skiprows=1, usecols=[0, 1], names=['GeneID', 'Swiss_protein_ID'])
    result_df = pd.merge(left=df, right=kegg_df, on='GeneID', how='left')
    result_df = pd.merge(left=result_df, right=nr_df, on='GeneID', how='left')
    result_df = pd.merge(left=result_df, right=swiss_df, on='GeneID', how='left')
    result_df.fillna(value='NA', inplace=True)
    return result_df


if __name__ == '__main__':
    parse_input = parse_input()
    fpkm, reads, kegg, nr, swiss, output_file = parse_input[0], parse_input[1], parse_input[2], parse_input[3], parse_input[4], parse_input[5]
    merged_df = merge_fpkm_reads(fpkm, reads)
    result_df = add_kns_def(merged_df, kegg, nr, swiss)
    if output_file:
        result_df.to_csv(output_file, sep='\t', index=False)
    else:
        result_df.to_csv('fpkm_and_reads_matrix_filtered_data_def.txt', sep='\t', index=False)
