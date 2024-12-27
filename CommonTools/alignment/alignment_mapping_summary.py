#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/12/17 14:27
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    args = argparse.ArgumentParser(description='')
    args.add_argument('-i', '--input', help='mapping_summary 文件目录', default='.')
    args.add_argument('-s', '--samples', help='samples_described.txt', default='samples_described.txt')
    args.add_argument('--msf', help='output mapping summary file (默认: mapping_summary.txt)', default='mapping_summary.txt')

    parsed_args = args.parse_args()
    return parsed_args


def summary_hisat2_mapping(mapping_file_dir, samples_file, output_file):
    """ 对所有比对文件进行汇总，输出一个总的比对结果

    Args:
        mapping_file_dir (str): mapping summary 输入目录
        samples_file (str): samples_described.txt
        output_file (str): 输出文件
    """
    column_name = ['Sample', 'Total reads', 'Mapped reads', 'Unmapped reads', 'Unique mapped reads', 'Multiple mapped reads', 'overall alignment rate']
    open(output_file, 'w').write('\t'.join(column_name)+'\n')
    samples_data = pd.read_csv(samples_file, sep='\t', usecols=['sample'])
    for sample in samples_data['sample'].to_list():
        mapping_file = os.path.join(mapping_file_dir, f'{sample}_mapping.txt')
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
            
        open(output_file, 'a').write(f'{sample}\t{total_reads}\t{mapped_reads}\t{unmapped_reads}\t{unique_mapped}\t{multi_mapped_reads}\t{overall_rate}\n')


if __name__ == '__main__':
    args = parse_input()
    summary_hisat2_mapping(args.input, args.samples, args.msf)
