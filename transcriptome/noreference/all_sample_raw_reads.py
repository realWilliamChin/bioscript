# -*- coding: UTF-8 -*-
import argparse
import os
from shutil import copyfile
import sys
import numpy
#import pandas

import sys


cur_path=os.curdir
doc_list=os.listdir(cur_path)
dic_raw1={}
dic_fpkm={}


sample_list=[]
column_num=0

for each_doc in doc_list:
    if '.genes.results' in each_doc:
        column_num += 1
        sample=each_doc.split('.genes.results')[0]
        sample_list.append(sample)
#获得sample样本数和sample list



doc_num=0
for each_doc in doc_list:
    if '.genes.results' in each_doc:
        doc_num+=1
        #正在运行的文件的位置，可以定位到列数
        with open (each_doc,'r') as f1:
            row_num=0
            #当前文件读的行数
            for each_line1 in f1:
                if row_num==0:
                    row_num+=1
                    continue
                else:
                    gene_id=each_line1.split('\t')[0]
                    content1=each_line1.split('\t')[4]
                    if gene_id in dic_raw1:
                        dic_raw1.setdefault(gene_id,[]).append(content1)
                    else:
                        i=0
                        while(i<(doc_num-1)):
                            dic_raw1.setdefault(gene_id, []).append('0')
                            i+=1
                        dic_raw1.setdefault(gene_id, []).append(content1)



f2=open('gene_count_matrix.txt','w')
f2.write('gene_id')
for each in sample_list:
    f2.write('\t' + each)
f2.write('\n')
for each_key in dic_raw1:
    raw_reads='\t'.join(dic_raw1[each_key])
    f2.write(each_key+'\t'+raw_reads+'\n')







#
