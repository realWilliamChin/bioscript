# -*- coding: UTF-8 -*-

import argparse
import os
from shutil import copyfile
import sys
import numpy

import sys
os.system("sed 's/,/\t/g' gene_fpkm_matrix.csv > gene_fpkm_matrix.txt")
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