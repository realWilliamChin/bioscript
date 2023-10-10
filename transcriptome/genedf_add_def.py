#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/07/03 16:50
# Author        : William GoGo
import argparse
import pandas as pd


def parse_input():
    argparser = argparse.ArgumentParser(description='输入 k，n，s，def 文件，添加 KEGG_ID, KEGG_GeneID, NR_Def, Swiss_protein_ID')
    argparser.add_argument('--kns', help='输入 kns_def 文件，添加 KEGG_ID, KEGG_GeneID, NR_Def, Swiss_protein_ID')
    argparser.add_argument('-k', '--kegg', help='KEGG_gene_def 文件，会添加 KEGG_ID 和 KEGG_Shortname 列')
    argparser.add_argument('-n', '--nr', help='NR_gene_def 文件')
    argparser.add_argument('-s', '--swiss', help='Swiss_gene_def 文件')
    argparser.add_argument('-i', '--input', help='输入文件')
    argparser.add_argument('--input_col', default="GeneID",
                           help='输入作为 GeneID 的列，可以是其他列，如果没有列名，则输入第几列，从 0 开始数')
    argparser.add_argument('--input_sep', help='输入文件分隔符，默认制表符')
    argparser.add_argument('--input_format', default='txt', choices=['txt', 'csv', 'xlsx', 'other'],
                           help='输入文件格式，默认 csv(, 分隔符), xlsx, txt(tab 分隔符), other 指定 --input_sep 分隔符')
    argparser.add_argument('--input_header', default='y', choices=['y', 'n'],
                           help="输入文件是否有列名，默认包含列名，如果没有列名，输入参数 n")
    argparser.add_argument('-o', '--output', default='output.txt', help='输出文件')
    return argparser.parse_args()


def add_kns_def(file_df, kegg_file=None, nr_file=None, swiss_file=None, kns_file=None):
    """
    添加 NR_Def, Swiss_protein_ID, KEGG_Pathway_ID, KEGG_Shortname, EC_number, KEGG_description 列
    """
    result_df = file_df
    source_shape = file_df.shape[0]

    if nr_file:
        nr_df = pd.read_csv(nr_file, sep='\t', skiprows=1, usecols=[0, 1, 2], names=['GeneID', 'NR_ID', 'NR_Def1'], dtype={'GeneID': str})
        # 没有 NCBI ID 直接加会丢失。先 fillna，再去掉 NA::
        nr_df.fillna(value='NA', inplace=True)
        nr_df['NR_ID_Des'] = nr_df['NR_ID'] + '::' + nr_df['NR_Def1']
        nr_df['NR_ID_Des'] = nr_df['NR_ID_Des'].str.replace('NA::', '', regex=False)
        # ID 列
        nr_id_df = nr_df[nr_df['NR_ID_Des'].str.contains('::')].copy()
        nr_id_df.drop(columns=['NR_Def1', 'NR_ID_Des'], inplace=True)
        result_df = pd.merge(left=result_df, right=nr_id_df, on='GeneID', how='left')
        # ID_Des 列
        nr_df.drop(columns=['NR_ID', 'NR_Def1'], inplace=True)
        result_df = pd.merge(left=result_df, right=nr_df, on='GeneID', how='left')

    if swiss_file:
        swiss_df = pd.read_csv(swiss_file, sep='\t', skiprows=1, usecols=[0, 1, 2], names=['GeneID', 'Swiss_ID', 'Swiss_Def'], dtype={'GeneID': str})
        swiss_df.fillna(value='NA', inplace=True)
        swiss_df['Swiss_ID_Des'] = swiss_df['Swiss_ID'] + '::' + swiss_df['Swiss_Def']
        swiss_df['Swiss_ID_Des'] = swiss_df['Swiss_ID_Des'].str.replace('NA::', '', regex=False)
        # ID 列
        swiss_id_df = swiss_df[swiss_df['Swiss_ID_Des'].str.contains('::')].copy()
        swiss_id_df.drop(columns=['Swiss_Def', 'Swiss_ID_Des'], inplace=True)
        result_df = pd.merge(left=result_df, right=swiss_id_df, on='GeneID', how='left')
        # ID_Des 列
        swiss_df.drop(columns=['Swiss_ID', 'Swiss_Def'], inplace=True)
        result_df = pd.merge(left=result_df, right=swiss_df, on='GeneID', how='left')
        
    if kegg_file:
        kegg_df = pd.read_csv(kegg_file, sep='\t', skiprows=1, names=['GeneID', 'KEGG_ID', 'KEGG_Shortname', 'EC_Number', 'KEGG_Description'], dtype={'GeneID': str})
        result_df = pd.merge(left=result_df, right=kegg_df, on='GeneID', how='left')

    if kns_file:
        kns_df = pd.read_csv(kns_file, sep='\t', dtype={'GeneID': 'str'})
        result_df = pd.merge(left=file_df, right=kns_df, on='GeneID', how='left')
        result_df['GeneID'].to_csv('test.list', sep='\t', index=False)

    result_df.fillna(value='NA', inplace=True)
    result_shape = result_df.shape[0]
    
    if source_shape != result_shape:
        print(f"原表行数{source_shape}, 结果表行数{result_shape}, 可能输入文件基因 ID 有重复，或注释文件有重复")
        result_df.drop_duplicates(subset='GeneID', inplace=True)
        print('已自动去重处理')

    return result_df


def main():
    args = parse_input()
    if args.input_format == 'csv':
        df = pd.read_csv(args.input, sep=',')
    elif args.input_format == 'xlsx':
        df = pd.read_excel(args.input, engine='openpyxl')
    elif args.input_format == 'txt':
        df = pd.read_csv(args.input, sep='\t')
    
    if args.input_header == 'y':
        df = pd.read_csv(args.input, sep='\t')
    elif args.input_header == 'n':
        df = pd.read_csv(args.input, sep='\t', header=None)
        args.input_col = int(args.input_col)
    # source_key_col_name = str(df.columns[args.input_col])
    df.rename(columns={args.input_col: 'GeneID'}, inplace=True)
    df['GeneID'] = df["GeneID"].astype(str)
    result_df = add_kns_def(df, args.kegg, args.nr, args.swiss, args.kns)
    
    if args.input_header == 'y':
        result_df.rename(columns={'GeneID': args.input_col}, inplace=True)
    
    result_df.to_csv(args.output, sep='\t', index=False)
    print('\nDone!\n')


if __name__ == '__main__':
    main()
