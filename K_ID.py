#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/29 17:19
# Author        : William GoGo
# ---------------------------------------
# 读取 de_results_add_def.py 生成的 DEG_data.txt 文件
# 把 regulation 变为 1,0,-1 对应 up,nosig,down
# 生成的文件只包括 K_ID 和 regulation
# 例如：K12345    0
# 一个 K_ID 对应多个 GeneID，有可能 up,down,nosig 都有
# 对 K_ID 去重，哪个 regulation 多留下哪个，一样多的变为 0
# 并生成一个中间文件，方便检查
# ---------------------------------------
import os
import pandas as pd
import argparse


def parse_input():
    argparser = argparse.ArgumentParser(description='读取 de_results_add_def.py 生成的文件')
    argparser.add_argument('--kegg-clean', type=str, required=True, dest='kegg_clean',
                           help="输入 kegg 注释出来的 'KEGG_clean 文件'")
    argparser.add_argument('--kegg-path', type=str, default='/home/colddata/qinqiang/script/lib/passed_path.txt', dest='kegg_pathway',
                           help="输入 passed_path.txt, 默认在 /home/colddata/qinqiang/script/lib/passed_path.txt 已经存在")
    return argparser.parse_args()


def k_id(kegg_file, kegg_passed_path):
    for each_file in os.listdir():
        if 'DEG_data.txt' in each_file:
            # 数据前处理
            deg_data_df = pd.read_csv(each_file, sep='\t', usecols=['GeneID', 'regulation', 'KEGG_ID'])
            deg_data_df = deg_data_df[deg_data_df['KEGG_ID'] != 'NA']
            deg_data_df['regulation'] = deg_data_df['regulation'].str.replace('Down', '-1', regex=False)
            deg_data_df['regulation'] = deg_data_df['regulation'].str.replace('Up', '1', regex=False)
            deg_data_df['regulation'] = deg_data_df['regulation'].str.replace('NoSignificant', '0', regex=False)
            
            # 生成中间文件 K_ID，然后 GeneID + regulation
            temp_file_df = deg_data_df[['KEGG_ID', 'GeneID', 'regulation']].copy()
            temp_file_df['GeneID_regulation'] = temp_file_df['GeneID'] + '_' + temp_file_df['regulation']
            temp_file_df = temp_file_df.groupby('KEGG_ID').agg(lambda x: ';'.join(x))
            temp_file_df[['GeneID_regulation']].to_csv(each_file.replace('_DEG_data.txt', '_K_GeneID.txt'), sep='\t')
            
            # 生成 KEGG_ID, regulation 文件
            deg_data_df['regulation'].astype(int)
            grouped = deg_data_df.groupby('KEGG_ID')['regulation'].value_counts().unstack(fill_value=0)
            up_reg, non_reg, down_reg = None, None, None
            for i, columns in enumerate(grouped.columns):
                if columns == '1':
                    up_reg = i
                elif columns == '0':
                    non_reg = i
                elif columns == '-1':
                    down_reg = i
            result_data = {'KEGG_ID': [], 'regulation': []}
            for kegg_id, row in grouped.iterrows():
                if row[up_reg] > row[down_reg]:
                    result_data['KEGG_ID'].append(kegg_id)
                    result_data['regulation'].append(1)
                elif row[up_reg] < row[down_reg]:
                    result_data['KEGG_ID'].append(kegg_id)
                    result_data['regulation'].append(-1)
                else:
                    result_data['KEGG_ID'].append(kegg_id)
                    result_data['regulation'].append(0)
            result_df = pd.DataFrame(result_data)
            result_df.sort_values(by='regulation', inplace=True, ascending=False)
            # 检查是否有空行
            test = result_df.dropna(axis=0)
            if result_df.shape == test.shape:
                print('OK!')
            else:
                print("Warning: 可能有空行", result_df.shape, test.shape)
            result_df.to_csv(each_file.replace('_DEG_data.txt', '_K_ID.txt'), sep='\t', index=False)
            
            
            # 输出 Kegg_valid_pathway
            kegg_clean_df = pd.read_csv(kegg_file, sep='\t', usecols=[1, 4], names=['ko_number', 'KEGG_ID'])
            kegg_clean_df.drop_duplicates(subset=['ko_number', 'KEGG_ID'], inplace=True)
            kegg_clean_df['ko_number'] = kegg_clean_df['ko_number'].str.split(':').str[0]
            kegg_clean_df = kegg_clean_df[kegg_clean_df['KEGG_ID'].isin(result_df['KEGG_ID'])]
            kegg_pathway_df = pd.read_csv(kegg_passed_path, sep='\t', names=['ko_number', 'Des'])
            kegg_pathway_df = kegg_pathway_df[kegg_pathway_df['ko_number'].isin(kegg_clean_df['ko_number'])]
            kegg_pathway_df.to_csv(each_file.replace('DEG_data.txt', 'passed_path.txt'), sep='\t', index=False, header=False)


def main():
    args = parse_input()
    k_id(args.kegg_clean, args.kegg_pathway)
    

if __name__ == "__main__":
    main()