# -*- coding: UTF-8 -*-
import os, sys
import argparse
from loguru import logger



#os.system("sed 's/,/\t/g' gene_count_matrix.csv > gene_count_matrix.txt")
def parse_input():
    parser = argparse.ArgumentParser(description=' ')

    parser.add_argument('-n', type=int, default=50)
    parser.add_argument('-i', default='gene_count_matrix.txt')
    parser.add_argument('-o', default='reads_matrix_filtered.txt')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_input()
    
    greater_than = int(args.n)
    input_file = args.i
    output_file = args.o

    f2=open(output_file,'w')
    num=0
    num1=0

    with open(input_file,'r') as f1:
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
                    if(float(each)>greater_than):
                        line1=line.replace('\t','\t')
                        f2.write(line1)
                        f2.write('\n')
                        num1+=1
                        break
            num+=1
    f2.close()
    logger.info('total_filtered(not include head line):' + str(num1))
    logger.info('total(include head line):'+str(num))
