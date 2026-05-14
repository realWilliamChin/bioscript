#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/12/17 14:27
# Author        : William GoGo

import os, sys
import argparse
import pandas as pd
from loguru import logger


def parse_input():
    args = argparse.ArgumentParser()
    args.add_argument('-i', '--input', nargs='+', default=['.'],
                      help='mapping_summary 文件目录或多个mapping_stat.txt文件路径')
    args.add_argument('-s', '--samples', default='samples_described.txt',
                      help='samples_described.txt (当-i指定为目录时需要)')
    args.add_argument('-o', '--output', default='alignment_report.txt', 
                      help='output mapping summary file (默认: alignment_report.txt)')
    
    args = args.parse_args()
    return args


def summary_hisat2_mapping(input_paths, samples_file, output_file):
    """ 对所有比对文件进行汇总，输出一个总的比对结果
    解析格式如下：
    24827879 reads; of these:
      24827879 (100.00%) were paired; of these:
        12424614 (50.04%) aligned concordantly 0 times
        10625156 (42.80%) aligned concordantly exactly 1 time
        1778109 (7.16%) aligned concordantly >1 times
        ----
        12424614 pairs aligned concordantly 0 times; of these:
          1324528 (10.66%) aligned discordantly 1 time
        ----
        11100086 pairs aligned 0 times concordantly or discordantly; of these:
          22200172 mates make up the pairs; of these:
            15019088 (67.65%) aligned 0 times
            6129907 (27.61%) aligned exactly 1 time
            1051177 (4.73%) aligned >1 times
    69.75% overall alignment rate

    Args:
        input_paths (list): mapping summary 输入目录或多个mapping_stat.txt文件路径
        samples_file (str): samples_described.txt (当-i指定为目录时需要)
        output_file (str): 输出文件
    """
    column_name = ['Sample', 'Total reads', 'Mapped reads', 'Unmapped reads', 'Unique mapped reads', 'Multiple mapped reads', 'overall alignment rate']
    open(output_file, 'w').write('\t'.join(column_name)+'\n')
    
    # 判断是文件模式还是目录模式
    mapping_files = []
    if len(input_paths) == 1 and os.path.isdir(input_paths[0]):
        # 目录模式，按照samples文件处理
        mapping_file_dir = input_paths[0]
        samples_data = pd.read_csv(samples_file, sep='\t', usecols=['sample'])
        for sample in samples_data['sample'].to_list():
            # 支持两种后缀格式
            mapping_file_stat = os.path.join(mapping_file_dir, f'{sample}_mapping_stat.txt')
            mapping_file_normal = os.path.join(mapping_file_dir, f'{sample}_mapping.txt')
            if os.path.exists(mapping_file_stat):
                mapping_file = mapping_file_stat
            elif os.path.exists(mapping_file_normal):
                mapping_file = mapping_file_normal
            else:
                logger.warning(f"找不到映射文件: {mapping_file_stat} 或 {mapping_file_normal}，跳过该样本")
                continue
            mapping_files.append((sample, mapping_file))
    else:
        # 文件模式，直接处理所有输入文件
        for file_path in input_paths:
            if not os.path.exists(file_path):
                logger.warning(f"找不到文件: {file_path}，跳过")
                continue
            # 支持两种后缀格式
            if file_path.endswith('_mapping_stat.txt'):
                sample_name = os.path.basename(file_path).replace('_mapping_stat.txt', '')
            elif file_path.endswith('_mapping.txt'):
                sample_name = os.path.basename(file_path).replace('_mapping.txt', '')
            else:
                logger.warning(f"文件 {file_path} 不是_mapping.txt或_mapping_stat.txt格式，跳过")
                continue
            mapping_files.append((sample_name, file_path))
    
    # 处理所有mapping文件
    for sample, mapping_file in mapping_files:
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
    summary_hisat2_mapping(args.input, args.samples, args.output)