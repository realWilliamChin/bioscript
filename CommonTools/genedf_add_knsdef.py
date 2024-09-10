#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/07/03 16:50
# Author        : William GoGo
import argparse
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser(description='对各种表添加基因定义文件，指定 --kns 则无需指定 -k,-n,-s')
    parser.add_argument('--kns', help='输入 kns_def 文件，添加 KEGG_ID, KEGG_GeneID, NR_Def, Swiss_protein_ID')
    parser.add_argument('-k', '--kegg', help='KEGG_gene_def 文件，会添加 KEGG_ID 和 KEGG_Shortname 列')
    parser.add_argument('-n', '--nr', help='NR_gene_def 文件')
    parser.add_argument('-s', '--swiss', help='Swiss_gene_def 文件')
    parser.add_argument('-i', '--input', help='输入文件，可根据文件后缀格式类型读取，txt是tab分隔')
    parser.add_argument('--input-sep', dest='input_sep', help='自定义文件分隔符，默认制表符')
    parser.add_argument('--input-header', default='GeneID', dest='input_header',
                        help="默认 GeneID 添加定义，如果有其他列名，请写出列名，如果没有列名，输入列的位置，从 0 开始数，列名不可以是数字")
    parser.add_argument('-o', '--output', default='output.txt', help='输出文件')
    
    add_expression = parser.add_argument_group(title="指定添加 expression 的参数")
    add_expression.add_argument('-e', '--expression', type=str, help='表达量文件，通常是 fpkm_reads_matrix_data_def.txt')
    
    args = parser.parse_args()
    
    if args.input_sep:
        pass
    elif args.input.endswith('.txt'):
        args.input_sep = '\t'
    elif args.input.endswith('.csv'):
        args.input_sep = ','
    elif args.input.endswith('.tsv'):
        args.input_sep = '\t'
    
    if args.kns:
        args.kegg, args.nr, args.swiss = None, None, None
    
    return args


def add_kns_def(file_df, kegg_file=None, nr_file=None, swiss_file=None, kns_file=None):
    """
    传入表中必须包含 GeneID 列
    添加 NR_Def, Swiss_protein_ID, KEGG_ID, KEGG_Shortname, EC_number, KEGG_Description 列
    """
    result_df = file_df.copy()
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
        result_df = pd.merge(left=result_df, right=kns_df, on='GeneID', how='left')

    diff_col = list(set(result_df.columns) - set(file_df.columns))
    result_df[diff_col] = result_df[diff_col].fillna(value='NA')
    result_shape = result_df.shape[0]
    
    if source_shape != result_shape:
        print(f"原表行数{source_shape}, 结果表行数{result_shape}, 可能输入文件基因 ID 有重复，或注释文件有重复")
        result_df.drop_duplicates(subset='GeneID', inplace=True)
        print('已自动去重处理')

    return result_df


def add_expression_data(input_file_df, expression_data):
    expression_df = pd.read_csv(expression_data, sep='\t', dtype={'GeneID':str})
    expression_df_numeric_cols = expression_df.select_dtypes(include=['number']).columns
    expression_df_string_cols = expression_df.select_dtypes(include=['object']).columns
    expression_df[expression_df_numeric_cols] = expression_df[expression_df_numeric_cols].fillna(0)
    expression_df[expression_df_string_cols] = expression_df[expression_df_string_cols].fillna('NA')

    result_df = pd.merge(input_file_df, expression_df, on="GeneID", how='left')
    result_df.fillna(0, inplace=True)
    # result_df.to_csv(args.output, sep='\t', index=False)
    
    return result_df


def main():
    args = parse_input()
    
    try:
        args.input_header = int(args.input_header)
    except ValueError:
        pass
    
    if args.input.endswith('.xlsx'):
        df = pd.read_excel(args.input, engine='openpyxl')
    elif type(args.input_header) == str:
        df = pd.read_csv(args.input, sep=args.input_sep)
    elif type(args.input_header) == int:
        df = pd.read_csv(args.input, sep=args.input_sep, header=None)
    else:
        print('输入的 input_header 有误，请检查')
        exit(1)
    
    df.rename(columns={args.input_header: 'GeneID'}, inplace=True)
    df['GeneID'] = df["GeneID"].astype(str)
    
    if args.expression:
        result_df = add_expression_data(df, args.expression)
    else:
        result_df = add_kns_def(df, args.kegg, args.nr, args.swiss, args.kns)
    
    if args.output.endswith('.xlsx'):
        result_df = result_df.rename(columns={'GeneID': args.input_header})
        result_df.to_excel(args.output, index=False, engine='openpyxl')
    if type(args.input_header) == str:
        result_df = result_df.rename(columns={'GeneID': args.input_header})
        result_df.to_csv(args.output, sep='\t', index=False)
    elif type(args.input_header) == int:
        result_df.to_csv(args.output, sep='\t', index=False, header=True)
        
    print('\nDone!\n')


if __name__ == '__main__':
    main()
