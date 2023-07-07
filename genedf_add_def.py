#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/07/03 16:50
# Author        : William GoGo
import argparse
import pandas as pd


def parse_input():
    argparser = argparse.ArgumentParser(description='输入 k，n，s，def 文件，添加 KEGG_ID, KEGG_GeneID, NR_Def, Swiss_protein_ID')
    argparser.add_argument('-k', '--kegg', help='KEGG_gene_def 文件，会添加 KEGG_ID 和 KEGG_Shortname 列')
    argparser.add_argument('--kegg_id', action='store_true', help='只保留 KEGG_ID 列')
    argparser.add_argument('--kegg_short', action='store_true', help='只保留 KEGG_Shortname 列')
    argparser.add_argument('-n', '--nr', help='NR_gene_def 文件')
    argparser.add_argument('-s', '--swiss', help='Swiss_gene_def 文件')
    argparser.add_argument('-i', '--input', help='输入文件') 
    argparser.add_argument('--input_col', default=0, type=int, help='输入第几列作为 GeneID 列，默认第一列，从 0 开始数')
    argparser.add_argument('--input_sep', help='输入文件分隔符，默认制表符')
    argparser.add_argument('--input_format', default='txt', choices=['txt', 'csv', 'xlsx', 'other'],
                           help='输入文件格式，默认 csv(, 分隔符), xlsx, txt(tab 分隔符), other 指定 --input_sep 分隔符')
    argparser.add_argument('-o', '--output', default='output.txt', help='输出文件')
    return argparser.parse_args()


def add_kns_def(df, kegg_file=None, nr_file=None, swiss_file=None, kegg_id=False, kegg_short=False):
    """
    添加 KEGG_Pathway_ID, KEGG_Shortname, NR_Def, Swiss_protein_ID
    """
    if kegg_file:
        kegg_df = pd.read_csv(kegg_file, sep='\t', skiprows=1, usecols=[0, 1, 2], names=['GeneID', 'KEGG_ID', 'Gene_Shortname'])
        result_df = pd.merge(left=df, right=kegg_df, on='GeneID', how='left')
        if kegg_id:
            result_df.drop(columns=['Gene_Shortname'], inplace=True)
        elif kegg_short:
            result_df.drop(columns=['KEGG_ID'], inplace=True)
    if nr_file:
        nr_df = pd.read_csv(nr_file, sep='\t', skiprows=1, usecols=[0, 1, 2], names=['GeneID', 'NR_ID', 'NR_Def1'])
        # 没有 NCBI ID 直接加会丢失先 fillna，再去掉 NA::
        nr_df.fillna(value='NA', inplace=True)
        nr_df['NR_Def'] = nr_df['NR_ID'] + '::' + nr_df['NR_Def1']
        nr_df['NR_Def'] = nr_df['NR_Def'].str.replace('NA::', '', regex=False)
        nr_df.drop(columns=['NR_ID', 'NR_Def1'], inplace=True)
        result_df = pd.merge(left=result_df, right=nr_df, on='GeneID', how='left')
    if swiss_file:
        swiss_df = pd.read_csv(swiss_file, sep='\t', skiprows=1, usecols=[0, 1], names=['GeneID', 'Swiss_Protein_ID'])
        result_df = pd.merge(left=result_df, right=swiss_df, on='GeneID', how='left')
    result_df.fillna(value='NA', inplace=True)
    return result_df


def main():
    args = parse_input()
    if args.input_format == 'csv':
        df = pd.read_csv(args.input, sep=',')
    elif args.input_format == 'xlsx':
        df = pd.read_excel(args.input, engine='openpyxl')
    elif args.input_format == 'txt':
        df = pd.read_csv(args.input, sep='\t')
        
    df = pd.read_csv(args.input, sep='\t')
    source_key_col_name = df.columns[args.input_col]
    df.rename(columns={source_key_col_name: 'GeneID'}, inplace=True)
    result_df = add_kns_def(df, args.kegg, args.nr, args.swiss, args.kegg_id, args.kegg_short)
    result_df.rename(columns={'GeneID': source_key_col_name}, inplace=True)
    result_df.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
