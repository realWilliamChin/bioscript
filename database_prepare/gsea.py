# -*- coding: utf-8 -*-
# @FileName : gsea.py
# @DATE ：16:14 2022/4/18
# @Author : CS:GO

import os
import re


cur_path=os.getcwd()
output_path='./'
dir_list=os.listdir(cur_path)
dic1={}

for each in dir_list:
    if ('.txt' in each)&(('GO' in each)|('KEGG' in each)):
        if os.path.exists(output_path)==False:
            os.mkdir(output_path)
        # if os.path.exists('./GMT')==True:
        #     os.chdir('./GMT')
        #     command='rm *.GMT'
        #     os.system(command)
        #     os.chdir('../')
        key1=each.replace('.txt','').replace('_clean','')
        if 'BP' in key1:
            key='GO_BP'
        if 'CC' in key1:
            key = 'GO_CC'
        if 'MF' in key1:
            key = 'GO_MF'
        if 'KEGG' in key1:
            key = 'KEGG'
        f2=open(output_path+os.sep+key1+'.GMT','w')
        dic1=key
        globals()[dic1] = {}
        dic={}


        with open(each,'r') as f1:
            num=0
            for each_line in f1:
                if num==0:
                    num+=1
                    continue
                else:
                    if 'KEGG' in  key1:
                        id=each_line.split('\t')[0].strip()
                        go_id=each_line.split('\t')[1].split(':')[0].strip()
                        content=each_line.split(':')[-1].strip()
                        if go_id in globals()[dic1]:
                            globals()[dic1][go_id] += '\t'+id
                        else:
                            globals()[dic1][go_id] = content+';;'+key+'\t'+id


                    else:
                        id=each_line.split('\t')[0].strip()
                        go_id=each_line.split('\t')[1].split('_')[0].strip()
                        content=each_line.split('_')[-1].strip()
                        if go_id in globals()[dic1]:
                            globals()[dic1][go_id] += '\t'+id
                        else:
                            globals()[dic1][go_id] = content+';;'+key+'\t'+id
        for each2 in globals()[dic1]:
            output=globals()[dic1][each2].split(';;')[0]+' '+'['+each2+']'+'\t'+globals()[dic1][each2].split(';;')[-1]+'\n'
            f2.write(output)






















