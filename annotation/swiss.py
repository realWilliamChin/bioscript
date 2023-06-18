#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2/24/2023 5:30 PM
# Author        : WilliamGoGo
import os, sys
import time
import argparse
import pandas as pd
import numpy as np


def parse_input():
    argparser = argparse.ArgumentParser(description='')
    argparser.add_argument('-b', '--blast', help='指定 swiss.blast 文件, 默认当前文件夹下的 _swiss.blast 文件')
    argparser.add_argument('-p', '--prefix', required=True, help='生成文件的前缀')
    args = argparser.parse_args()
    return args


def _keep_goid(s):
    """
    处理表里的元素，只留下 GO number，格式 GO:xxxxxx
    """
    if not pd.isna(s):
        if 'O:' in s and ']' in s:
            return 'G' + s.split(']')[0]
        else:
            return ''


def process_go(swiss_df, ref_df):
    df_lst = []
    for each_go in ['GO_BP', 'GO_CC', 'GO_MF']:
        # 分成 bp，cc，mf 的
        df = pd.merge(left=swiss_df[['GeneID', 'GOID']], right=ref_df[['GOID', each_go]], on='GOID', how='left')
        df = df.dropna().drop(columns=['GOID'])

        # 去除掉每个基因少于 3 个的 go
        df = df[df[each_go].str.count('GO') >= 3]

        # 把后面的 GO 功能都加到第二列
        df_expand = df[each_go].str.split(';', expand=True)
        df = pd.concat([df['GeneID'], df_expand], axis=1).melt(id_vars=['GeneID'], value_name='GOID')

        df = df.drop(columns=['variable']).dropna()

        # 对后面的 GOID 格式改成 ID 在前，使用 _ 连接
        df_expand = df['GOID'].str.split('\[GO', expand=True)
        df['Def'] = df_expand[1].str.replace(':', 'GO:', regex=True).str.replace(']', '_', regex=True) + df_expand[0].str.strip()
        df = df.drop(columns=['GOID']).sort_values('GeneID')
        df_lst.append(df)
    return df_lst


def main():
    args = parse_input()
    swiss_file = args.blast if args.blast else [x for x in os.listdir() if '_swiss.blast' in x][0]
    # 读取 swiss 参考文件和 blast 文件，并初始化
    swiss_df = pd.read_csv(swiss_file, sep='\t', usecols=[0, 1, 15], names=['GeneID', 'GOID', 'Swiss_Def'])
    swiss_df = swiss_df.drop_duplicates(subset='GeneID', keep='first', inplace=False)
    siwss_df_expand = swiss_df['Swiss_def'].str.split(';', expand=True)
    swiss_df['Swiss_def'] = siwss_df_expand.iloc[:, 0]
    # 在 swiss_gene_def 中 Swiss_Def 删掉 RecName: Full=
    swiss_df['Swiss_def'] = swiss_df['Swiss_def'].str.replace('RecName: Full=', '')
    # 生成 _unigene_swiss_gene_def.txt
    swiss_df.to_csv(args.prefix + '_swiss_gene_def.txt', sep='\t', index=False, header=['GeneID', 'Swissprot_ID', 'Swiss_Def'])

    ref_file = '/home/data/ref_data/db/swiss_go_txt/Swiss_protein_go.txt'
    ref_df = pd.read_csv(ref_file, sep='\t', skiprows=1, names=['GOID', 'GO_BP', 'GO_CC', 'GO_MF'])

    # 生成 idNO_def 文件
    idNO_def_filename = args.prefix + '_swiss_idNo_def.txt'
    ref_df['merge_go'] = ref_df['GO_BP'] + '_' + ref_df['GO_CC'] + '_' + ref_df['GO_MF']
    ref_df['merge_go'].replace('', np.nan, regex=True, inplace=True)
    idNo_def = pd.merge(left=swiss_df.iloc[:, [0, 1]], right=ref_df.iloc[:, [0, 4]], on='GOID', how='left')
    idNo_def = idNo_def.dropna().drop(columns=['GOID'])
    idNo_def_expand = idNo_def['merge_go'].str.split('\[G', expand=True)
    idNo_def_expand = idNo_def_expand.applymap(_keep_goid)
    idNo_def = pd.concat([idNo_def['GeneID'], idNo_def_expand], axis=1).fillna('')
    idNo_def.to_csv(idNO_def_filename, sep='\t', index=False, header=False)
    with open(idNO_def_filename, 'r') as idNo_def_file:
        os.remove(idNO_def_filename)
        for line in idNo_def_file.readlines():
            line = line.replace('\t\t', '\t')
            line = line.strip() + '\n'
            with open(idNO_def_filename, 'a') as f:
                f.write(line)

    # 生成 wego 注释所需要的格式
    # go_bp_df_wego = go_bp_df_expand.applymap(_keep_goid)
    # go_cc_df_wego = go_cc_df_expand.applymap(_keep_goid)
    # go_mf_df_wego = go_mf_df_expand.applymap(_keep_goid)
    # go_bp_df_wego = pd.concat([go_bp_df['GeneID'], go_bp_df_wego], axis=1).fillna('')
    # go_cc_df_wego = pd.concat([go_cc_df['GeneID'], go_cc_df_wego], axis=1).fillna('')
    # go_mf_df_wego = pd.concat([go_mf_df['GeneID'], go_mf_df_wego], axis=1).fillna('')

    result_df = process_go(swiss_df, ref_df)
    # 保存 GO_BP GO_CC GO_MF 文件
    result_df[0].to_csv(swiss_file.replace('.blast', '_GO_BP_ID.txt'), sep='\t', index=False, header=False)
    result_df[1].to_csv(swiss_file.replace('.blast', '_GO_CC_ID.txt'), sep='\t', index=False, header=False)
    result_df[2].to_csv(swiss_file.replace('.blast', '_GO_MF_ID.txt'), sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()

