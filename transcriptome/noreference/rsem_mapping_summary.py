# -*- coding: UTF-8 -*-

import argparse
import os
from shutil import copyfile
import sys
import numpy

import sys


f1=open(r'mapping_summary.txt','w')
f1.write('Sample'+'\t'+'Total reads'+'\t'+'Mapped reads'+'\t''Unmapped reads'+'\t'+'Unique mapped reads'+'\t'+'Multiple mapped reads'+'\t'+'alignment rate'+'\n')
dirlist=os.listdir()

for dir_name in dirlist:
    if ('_mapping.txt' in dir_name):
        sample=dir_name.split('_mapping.txt')[0]
        f1.write(sample+'\t')
        with open(dir_name) as f2:
            data_list=[]
            row_num = 0
            for line in f2:
                line=line.strip()
                '''
                if row_num==0:
                    total_reads=line.split(' ')[0]
                    f1.write(total_reads+'\t')
                if row_num==3:
                    mapped_reads=line.split(' ')[0]
                    f1.write(total_reads+'\t')
                row_num+=1
                '''
                data_list.append(line.strip().split(' ')[0])


            total_reads=int(data_list[0])
            mapped_reads=int(int(data_list[3])+int(data_list[4]))
            unmapped_reads=int(data_list[2])
            unique_mapped=int(data_list[3])
            multi_mapped_reads=int(data_list[4])
            alignment_rate=data_list[-1]
            f1.write(str(total_reads)+'\t'+str(mapped_reads)+'\t'+str(unmapped_reads)+'\t'+str(unique_mapped)+'\t'+str(multi_mapped_reads)+'\t'+str(alignment_rate)+'\n')



f1.close()



