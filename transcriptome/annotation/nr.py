#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/14 17:05
# Author        : William GoGo
import argparse
import os
from shutil import copyfile
import sys
import numpy
import re
import pandas as pd

######################################################################################################################
###输入
input=os.getcwd()
parser=argparse.ArgumentParser(description="nr")
parser.add_argument('--input',dest='input',default=input,help='path of nr.blast')

arg=parser.parse_args()


dic_nr={}
cur_path=arg.input

output_path=cur_path
#output_path=cur_path+os.sep+'swiss_delete_duplicate_all_info'
isExists = os.path.exists(output_path)

# 判断结果
if not isExists:

    os.makedirs(output_path)

temp=output_path+os.sep+'temp'
isExists = os.path.exists(temp)

# 判断结果
if not isExists:

    os.makedirs(temp)
dir_list=os.listdir(cur_path)
for each_doc in dir_list:
    if ('nr.blast' in each_doc):
    #if ('swiss.blast' in each_doc):
        new_name1=each_doc.replace('.blast','_uniq.blast')
        new_name2 = each_doc.replace('.blast', '_after_sort_all_info.txt')
        # new_name3 = each_doc.replace('.blast', '_nr_gene_def.txt')
        # new_name4 = each_doc.replace('.blast','_nr_TF_def.txt')
        new_name3 = each_doc.replace('nr.blast', 'nr_gene_def.txt')
        new_name4 = each_doc.replace('nr.blast','nr_TF_def.txt')
        data_frame = pd.read_csv(cur_path+os.sep+each_doc, header=None ,sep='\t')  # 获取日期数据
        print(data_frame)
        data_frame2=data_frame.sort_values(by=[0, 2],ascending=[True,False])

        data_frame2.to_csv(path_or_buf=temp+os.sep+new_name2,sep='\t',header=None, index = False)


dic_nr1={}
dic_nr2={}
f1=open(temp+os.sep+new_name2,'r')
f2=open(output_path+os.sep+new_name1,'w')
f2.write('qseqid'+'\t'+ 'sseqid'+'\t'+ 'pident'+'\t'+ 'length'+'\t'+ 'mismatch'+'\t'+ 'gapopen'+'\t' +'qstart'+'\t'+ 'qend''\t' +'sstart'+'\t' +'send' +'\t'+'evalue' +'\t'+'bitscore' +'\t'+'stitle'+'\n')
f3=open(output_path+os.sep+new_name3,'w')
f3.write('GeneID'+'\t'+'NCBI_ID'+'\t'+'NR_Def'+'\n')

f4=open(output_path+os.sep+new_name4,'w')
f4.write('GeneID'+'\t'+'NCBI_ID'+'\t'+'NR_Def'+'\n')

for each_line in f1:
    id=each_line.split('\t')[0]
    content1=each_line.split('\t',1)[-1]
    content2 = each_line.split('\t')[-1]
    if id in dic_nr1:
        continue
    else:
        dic_nr1[id]=content1
        dic_nr2[id]=content2.split(' ')[0]+'\t'+content2.split(' ',1)[-1]
for each_key in dic_nr1:
    f2.write(each_key+'\t'+dic_nr1[each_key])
for each_key in dic_nr2:
    f3.write(each_key+'\t'+dic_nr2[each_key])
    if 'transcription' in dic_nr2[each_key]:
        f4.write(each_key+'\t'+dic_nr2[each_key])

f1.close()
f2.close()
f3.close()
f4.close()












# dic_nr={}
# for each_doc in dir_list:
#     if ('nr_diamond.blast' in each_doc):
#
#
#
#
#
#
#
#
#
#
#
# input_doc='Trifolium_repens_Baisanye_unigene_CDS_nr_diamond.blast'
# f1=open(input_doc,'r')
# f2=open(input_doc.split('.blast')[0]+'_extract.txt','w')
# for each_line in f1:
#     gene_id=each_line.split('\t')[0]
#     annotation=each_line.strip().split('\t')[-1]
#     if gene_id in dic_nr:
#         dic_nr[gene_id]+=';'+annotation
#     else:
#         dic_nr[gene_id]=annotation
# for each_key in dic_nr:
#     f2.write(each_key+'\t'+dic_nr[each_key]+'\n')
#
#
#
# f1.close()
#
# f2.close()

