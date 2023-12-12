#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/20 21:45
# Author        : William GoGo, chenshen
import os
import argparse


def parse_input():
    args = argparse.ArgumentParser(description='hisat2_mapping_summary.py')
    args.add_argument('-i', '--input', help='input hisat2 mapping dir', default='.')
    args.add_argument('-o', '--output', help='output mapping summary file', default='mapping_summary.txt')
    args.add_argument('-s', '--samples', help='samples_described.txt', default='samples_described.txt')

    return args.parse_args()


def process_mapping(mapping_file):
    with open(mapping_file) as f2:
        data_list=[]
        for line in f2:
            if 'nohup' in line:
                continue
            line=line.strip()
            data_list.append(line.split(' ')[0])

        total_reads=int(data_list[0])
        mapped_reads=int(int(data_list[3])+int(data_list[4])+int(data_list[7])+(int(data_list[12])+int(data_list[13]))/2)
        unmapped_reads=total_reads-mapped_reads
        unique_mapped=int(data_list[3])+int(data_list[7])
        multi_mapped_reads=mapped_reads-unique_mapped
        overall_rate = data_list[14].split(' ')[0]
    return total_reads, mapped_reads, unmapped_reads, unique_mapped, multi_mapped_reads, overall_rate


def main():
    args = parse_input()
    samples_data = open(args.samples, 'r').readlines()
    sample_data_list = [os.sep.join([args.input, x.split('\t')[1] + '_mapping.txt']) for x in samples_data[1:] if x.strip() != '']
    
    column_name = ['Sample', 'Total reads', 'Mapped reads', 'Unmapped reads', 'Unique mapped reads', 'Multiple mapped reads', 'overall alignment rate']
    open(args.output, 'w').write('\t'.join(column_name)+'\n')
    
    for sample_data in sample_data_list:
        sample_name = os.path.basename(sample_data).replace('_mapping.txt', '')
        total_reads, mapped_reads, unmapped_reads, unique_mapped, multi_mapped_reads, overall_rate = process_mapping(sample_data)
        open(args.output, 'a').write(f'{sample_name}\t{total_reads}\t{mapped_reads}\t{unmapped_reads}\t{unique_mapped}\t{multi_mapped_reads}\t{overall_rate}\n')
    


if __name__ == '__main__':
    main()



