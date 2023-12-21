#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import pandas as pd
import numpy as np
import argparse
from shutil import copyfile

def parse_input():
    parser = argparse.ArgumentParser(description='合并 fpkm 和 reads matrix, 然后加上三个注释的定义')
    parser.add_argument('-i', '--id', type=str, help='all_gene_id.txt 文件，输出只在 all_gene_id.txt 中存在的基因')

    return parser.parse_args()


def process(all_geneid_file=None):
    gene_fpkm_matrix_df = pd.read_csv('gene_fpkm_matrix.csv')
    gene_fpkm_matrix_df.rename(columns={"gene_id": "GeneID"}, inplace=True)

    # 如果 GeneID 中大部分包含 |，则取 | 前面的部分
    if gene_fpkm_matrix_df['GeneID'].str.contains('|').mean() > 0.6:
        if gene_fpkm_matrix_df['GeneID'].str.startswith('gene-').mean() > 0.6:
            gene_fpkm_matrix_df['GeneID'] = np.where(
                gene_fpkm_matrix_df['GeneID'].str.startswith('gene-'),
                gene_fpkm_matrix_df['GeneID'].str.extract('gene-(\S+)|')[0],
                gene_fpkm_matrix_df['GeneID']
            )
            gene_fpkm_matrix_df['GeneID'] = np.where(
                gene_fpkm_matrix_df['GeneID'].str.startswith('LOC'),
                gene_fpkm_matrix_df['GeneID'].str.replace('LOC', ''),
                gene_fpkm_matrix_df['GeneID']
            )
        gene_fpkm_matrix_df['GeneID'] = np.where(
            gene_fpkm_matrix_df['GeneID'].str.contains('|'),
            gene_fpkm_matrix_df['GeneID'].str.split('|').str[0],
            gene_fpkm_matrix_df['GeneID']
        )
    if gene_fpkm_matrix_df['GeneID'].str.contains('gene:').mean() > 0.6:
        gene_fpkm_matrix_df['GeneID'] = np.where(
            gene_fpkm_matrix_df['GeneID'].str.contains('gene:'),
            gene_fpkm_matrix_df['GeneID'].str.split('gene:').str[1],
            gene_fpkm_matrix_df['GeneID']
        )
    
    if all_geneid_file:
        all_geneid_df = pd.read_csv(all_geneid_file, sep='\t', header=None, dtype=str)
        all_geneid_df.rename(columns={0: 'GeneID'}, inplace=True)
        before_count = gene_fpkm_matrix_df.shape[0]
        gene_fpkm_matrix_df = pd.merge(left=all_geneid_df, right=gene_fpkm_matrix_df, on='GeneID', how='left')
        gene_fpkm_matrix_df.fillna(0, inplace=True)
        gene_fpkm_matrix_df = gene_fpkm_matrix_df.astype({'GeneID': str})
        gene_fpkm_matrix_df.to_csv('gene_fpkm_matrix.txt', sep='\t', index=False)
        after_count = gene_fpkm_matrix_df.shape[0]
        print(f'before: {before_count}, after: {after_count}，未包含在 all_gene_id.txt 中的基因数：{before_count - after_count}')
    else:
        gene_fpkm_matrix_df.to_csv('gene_fpkm_matrix.txt', sep='\t', index=False)

    num=0
    num1=0
    num3=0
    id_list=[]
    f3=open('reads_matrix_filtered.txt','r')
    for line in f3.readlines():
            id=line.split('\t')[0]
            id_list.extend([id])
            num1+=1
    print('num of count_matri_filter(include head line)：'+str(num1))
    f2=open('fpkm_matrix_filtered.txt','w')
    num=0
    with open('gene_fpkm_matrix.txt','r') as f1:
        for line in f1.readlines():
            if (num == 0):
                each_part = line.split('\t')
                f2.write('GeneID' + '\t')

                #for each in each_part[1:]:
                f2.write('\t'.join(each_part[1:]))


            else:
                line=line.strip()
                id_part=line.split('|')[0]
                if ('|' not in line):
                    id_part=line.split('\t',1)[0]

                if(id_part in id_list):
                    new_line=line.split('\t',1)[-1]
                    f2.write(id_part+'\t'+new_line+'\n')
                    num3+=1
            num+=1
    f2.close()
    f3.close()
    print('mapped fpkm list(not include head line):'+str(num3))
    print('num of fpkm_matrix(include head line)'+str(num))


def main():
    args = parse_input()
    process(args.id)


if __name__ == '__main__':
    main()