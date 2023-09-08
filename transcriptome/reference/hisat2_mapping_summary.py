#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/20 21:45
# Author        : chenshen
import os

out_doc = 'mapping_summary.txt'
f1=open(out_doc,'w')
f1.write('Sample'+'\t'+'Total reads'+'\t'+'Mapped reads'+'\t''Unmapped reads'+'\t'+'Unique mapped reads'+'\t'+'Multiple mapped reads'+'\n')

for dir_name in os.listdir():
    if '_mapping.txt' in dir_name:
        sample=dir_name.split('_mapping.txt')[0]
        f1.write(sample+'\t')
        with open(dir_name) as f2:
            data_list=[]
            row_num = 0
            for line in f2:
                if 'nohup' in line:
                    continue
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
                data_list.append(line.split(' ')[0])

            total_reads=int(data_list[0])
            mapped_reads=int(int(data_list[3])+int(data_list[4])+int(data_list[7])+(int(data_list[12])+int(data_list[13]))/2)
            unmapped_reads=total_reads-mapped_reads
            unique_mapped=int(data_list[3])+int(data_list[7])
            multi_mapped_reads=mapped_reads-unique_mapped
            f1.write(str(total_reads)+'\t'+str(mapped_reads)+'\t'+str(unmapped_reads)+'\t'+str(unique_mapped)+'\t'+str(multi_mapped_reads)+'\n')
f1.close()



