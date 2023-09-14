# -*- coding: UTF-8 -*-
import argparse
import os
from shutil import copyfile
import sys
import numpy

import sys




os.system("sed 's/,/\t/g' gene_count_matrix.csv > gene_count_matrix.txt")

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