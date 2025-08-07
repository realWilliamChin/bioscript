#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/21 12:00
# Author        : William GoGo

import argparse
import pandas as pd
import os
from shutil import copyfile


def parse_input():
    args = argparse.ArgumentParser(description='')
    args.add_argument('-i', '--input', help='mapping_summary 文件目录', default='.')
    args.add_argument('-s', '--samples', help='sampels_described.txt')
    args.add_argument('-o', '--output', help='output mapping summary file (默认: alignment_report.txt)', default='alignment_report.txt')

    parsed_args = args.parse_args()
    return parsed_args


def mapping_summary(data_dir, samples_file, output_file):
    samples_data = pd.read_csv(samples_file, sep='\t', usecols=['sample'])
    open(output_file, 'w').write('Sample'+'\t'+'Total reads'+'\t'+'Mapped reads'+'\t''Unmapped reads'+'\t'+'Unique mapped reads'+'\t'+'Multiple mapped reads'+'\t'+'alignment rate'+'\n')
    for sample in samples_data['sample'].to_list():
        mapping_file = os.path.join(data_dir, f'{sample}_mapping_stat.txt')
        
        with open(output_file, 'a') as f1, open(mapping_file, 'r') as f2:
            f1.write(sample+'\t')
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


def main():
    args = parse_input()
    mapping_summary(args.input, args.samples, args.output)
    

if __name__ == '__main__':
    main()


