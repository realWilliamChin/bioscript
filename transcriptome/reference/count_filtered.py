#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import pandas as pd
import argparse
import os
from shutil import copyfile
import sys
import numpy

import sys

# os.system("sed 's/,/\t/g' gene_count_matrix.csv > gene_count_matrix.txt")
gene_count_matrix_df = pd.read_csv('gene_count_matrix.csv')
gene_count_matrix_df.rename(columns={"gene_id": "GeneID"}, inplace=True)
if gene_count_matrix_df['GeneID'].str.contains('|').mean() > 0.6:
    if gene_count_matrix_df['GeneID'].str.contains('gene-LOC').mean() > 0.6:
        gene_count_matrix_df['GeneID'] = gene_count_matrix_df['GeneID'].str.split('|').str[0].str.split('gene-LOC').str[1]
    else:
        gene_count_matrix_df['GeneID'] = gene_count_matrix_df['GeneID'].str.split('|').str[0]

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
                if(int(each)>50):
                    line1=line.replace('\t','\t')
                    f2.write(line1)
                    f2.write('\n')
                    num1+=1
                    break


        num+=1
f2.close()
print('total_filtered(not include head line):' + str(num1))
print('total(include head line):'+str(num))