# -*- coding: UTF-8 -*-
import argparse
import os
from shutil import copyfile
import sys
import numpy
cur_path=os.getcwd()
doc=os.listdir(cur_path)

'''
在当前文件夹中提出特定_ID.txt文件里的序列，
需要--input 文件  文件中包含.fasta文件 _ID.txt文件


'''

parser=argparse.ArgumentParser(description="把当前文件夹中的Gene_def文件夹中的所有注释文件，匹配_ID.txt文件里的geneid，输出到_ID_all_def.txt中 ")

parser.add_argument('--input',dest='input',default=os.getcwd(),
        help='path of input file（包含.fasta文件和_ID.txt）')




arg=parser.parse_args()



cur_path=os.getcwd()



dic1={}

with open('Put_cluster.txt','r') as f1:
    num1=0
    for each_line1 in f1:
        if num1==0:
            num1+=1
        else:
            group=each_line1.strip().split('\t')[-1]
            id=each_line1.split('\t')[0]
            if group in dic1:
                dic1[group]+='\t'+id
            else:
                dic1[group]=id
dic_aba={}
with open('Put_time_normalization.txt','r') as f2:
    num2=0
    for each_line2 in f2:
        if num2==0:
            num2+=1
            title=each_line2.replace('\"','')
            continue
        else:
            id2=each_line2.strip().replace('\"','').split('\t')[0]
            content=each_line2.strip().replace('\"','').split('\t',1)[-1]
            dic_aba[id2]=content



for each in dic1:
    new_doc='group'+each+'.txt'
    f1=open(new_doc,'w')
    f1.write('\t'+title)
    id_list=dic1[each].split('\t')

    for each_id in id_list:
        if each_id not in dic_aba:
            print(each_id)
            continue
        out=each_id +'\t'+dic_aba[each_id]
        f1.write(out+'\n')
    f1.close()















