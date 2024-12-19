# -*- coding: UTF-8 -*-
import os, sys
import argparse
from loguru import logger

def parse_input():
    parser = argparse.ArgumentParser(description=' ')

    parser.add_argument('-r', default='reads_matrix_filtered.txt')
    parser.add_argument('-i', default='gene_fpkm_matrix.txt')
    parser.add_argument('-o', default='fpkm_matrix_filtered.txt')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_input()
    reads_file = args.r
    input_file = args.i
    output_file = args.o
    
    
    num=0
    num1=0
    num3=0
    id_list=[]
    f3=open(reads_file,'r')
    for line in f3.readlines():
            id=line.split('\t')[0]
            id_list.extend([id])
            num1+=1
    logger.info('num of count_matri_filter(include head line)：'+str(num1))
    f2=open(output_file,'w')
    num=0
    with open(input_file,'r') as f1:
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
    logger.info('mapped fpkm list(not include head line):'+str(num3))
    logger.info('num of fpkm_matrix(include head line)'+str(num))