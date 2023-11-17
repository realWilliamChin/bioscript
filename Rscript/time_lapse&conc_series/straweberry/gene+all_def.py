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
def_list=cur_path+os.sep+'Gene_def'
doc_list=os.listdir(def_list)

dic1={}
deflist=[]
for each_doc in doc_list:

    if '_def.txt' in each_doc:

        dic1=each_doc.replace('_def.txt','')
        deflist.append(dic1)
        globals()[dic1]={}
        with open(def_list+os.sep+each_doc, 'r') as f:
            for each_line in f:
                geneid=each_line.split('\t')[0]
                content=each_line.strip().split('\t')[-1]
                globals()[dic1][geneid]=content
def_num=len(deflist)
doc1_list=os.listdir(cur_path)
for each_doc in doc1_list:
    if '_cluster.txt' in each_doc:
        new_name=each_doc.replace('.txt','_all_def.txt')
        f1=open(new_name,'w')

        row=0
        with open(each_doc,'r') as f2:
            for each_line in f2:
                if row==0:
                    f1.write('\t'+each_line.strip()+'\t'+'\t'.join(deflist)+'\n')
                    row+=1
                    continue
                else:
                    gene=each_line.replace('gene-','').strip().split('\t')[0]
                    annotation=''
                    for each_def in deflist:

                        if gene in globals()[each_def]:
                            annotation+='\t'+globals()[each_def][gene]
                        else:
                            annotation+='\t'+'N/A'
                    f1.write(each_line.replace('gene-','').strip()+annotation+'\n')
        f1.close()





















# #################################################################################33
# def rpsblast_annotation(input_file):
#     dic_cdd = {}
#     # cdd  numeber对应的content
#     with open('/opt/biosoft/RpsbProc-x64-linux/data/cddid.tbl', 'r') as f1:
#         for each_line in f1:
#             cdd_number = each_line.split('\t')[0]
#             content = each_line.strip().strip('\t', 1)[-1]
#             dic_cdd[cdd_number] = content
#
#     parser = argparse.ArgumentParser(description="把当前文件夹中的Gene_def文件夹中的所有注释文件，匹配_ID.txt文件里的geneid，输出到_ID_all_def.txt中 ")
#
#     parser.add_argument('--f', dest='input',
#                         help='cog/smart/tifr/...的rpsblast结果')
#
#     with open (input_file) as f2:
#         new_doc=input_file.replace('.out','_annotation.blast')
#         f3=open(new_doc,'w')
#         for each_line2 in f2:
#             cdd=each_line2.strip().split('\t')[1].replace('CDD:','')
#             if cdd in dic_cdd:
#                 f3.write(each_line2.strip()+'\t'+dic_cdd[cdd]+'\n')
#         f3.close()




