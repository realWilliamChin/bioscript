# -*- coding: UTF-8 -*-
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools'))
from load_input import load_table

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/transcriptome'))
from reorder_genetable_with_samplesdes import reindex


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', default='输入目录')
    parser.add_argument('-o', default='gene_fpkm_matrix.txt')
    parser.add_argument('-s', default='samples_described.txt')
    
    return parser.parse_args()


def main():
    args = parse_input()
    
    doc_list=os.listdir(args.i)
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
            with open (os.path.join(args.i, each_doc),'r') as f1:
                row_num=0
                #当前文件读的行数
                for each_line1 in f1:
                    if row_num==0:
                        row_num+=1
                        continue
                    else:
                        gene_id=each_line1.split('\t')[0]
                        content1=each_line1.split('\t')[6].strip()
                        dic_raw1.setdefault(gene_id,[]).append(content1)

    f2=open(args.o,'w')
    f2.write('GeneID')
    for each in sample_list:
        f2.write('\t' + each)
    f2.write('\n')
    for each_key in dic_raw1:
        raw_reads='\t'.join(dic_raw1[each_key])
        f2.write(each_key+'\t'+raw_reads+'\n')
    
    sample_lst = load_table(args.s)['sample'].tolist()
    reindex(sample_lst, args.o, args.o)
    

if __name__ == '__main__':
    main()







#
