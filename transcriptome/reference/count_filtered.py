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
    # os.system("sed 's/,/\t/g' gene_count_matrix.csv > gene_count_matrix.txt")
    gene_count_matrix_df = pd.read_csv('gene_count_matrix.csv')
    gene_count_matrix_df.rename(columns={"gene_id": "GeneID"}, inplace=True)
    if gene_count_matrix_df['GeneID'].str.contains('|').mean() > 0.6:
        if gene_count_matrix_df['GeneID'].str.startswith('gene-').mean() > 0.6:
            gene_count_matrix_df['GeneID'] = np.where(
                gene_count_matrix_df['GeneID'].str.startswith('gene-'),
                gene_count_matrix_df['GeneID'].str.extract('gene-(\S+)|')[0],
                gene_count_matrix_df['GeneID']
            )
            gene_count_matrix_df['GeneID'] = np.where(
                gene_count_matrix_df['GeneID'].str.startswith('LOC'),
                gene_count_matrix_df['GeneID'].str.replace('LOC', ''),
                gene_count_matrix_df['GeneID']
            )
        gene_count_matrix_df['GeneID'] = np.where(
            gene_count_matrix_df['GeneID'].str.contains('|'),
            gene_count_matrix_df['GeneID'].str.split('|').str[0],
            gene_count_matrix_df['GeneID']
        )
    if gene_count_matrix_df['GeneID'].str.contains('gene:').mean() > 0.6:
        gene_count_matrix_df['GeneID'] = np.where(
            gene_count_matrix_df['GeneID'].str.contains('gene:'),
            gene_count_matrix_df['GeneID'].str.split('gene:').str[1],
            gene_count_matrix_df['GeneID']
        )

    if all_geneid_file:
        all_geneid_df = pd.read_csv(all_geneid_file, sep='\t', header=None, dtype=str)
        all_geneid_df.rename(columns={0: 'GeneID'}, inplace=True)
        before_count = gene_count_matrix_df.shape[0]
        gene_count_matrix_df = pd.merge(left=all_geneid_df, right=gene_count_matrix_df, on='GeneID', how='left')
        gene_count_matrix_df.fillna(0, inplace=True)
        gene_count_matrix_df = gene_count_matrix_df.astype({'GeneID': str})
        gene_count_matrix_df.to_csv('gene_count_matrix.txt', sep='\t', index=False)
        after_count = gene_count_matrix_df.shape[0]
        print(f'before: {before_count}, after: {after_count}，未包含在 all_gene_id.txt 中的基因数：{before_count - after_count}')
    else:
        gene_count_matrix_df.to_csv('gene_count_matrix.txt', sep='\t', index=False)


    f2=open('reads_matrix_filtered.txt','w')
    num=0
    num1=0

    with open('gene_count_matrix.txt','r') as f1:
        for line in f1.readlines():
            line=line.strip()
            if(num==0):

                each_part=line.split('\t')
                f2.write( 'GeneID'+ '\t')
                #for each in each_part[1:]:
                f2.write(str('\t'.join(each_part[1:])))
                print(str('\t'.join(each_part[1:])))
                f2.write('\n')

            else:


                each_part = line.split('\t')
                for each in each_part[1:]:
                    if(int(float(each))>50):
                        line1=line.replace('\t','\t')
                        f2.write(line1)
                        f2.write('\n')
                        num1+=1
                        break


            num+=1
    f2.close()
    print('total_filtered(not include head line):' + str(num1))
    print('total(include head line):'+str(num))


def main():
    args = parse_input()
    process(args.id)


if __name__ == '__main__':
    main()