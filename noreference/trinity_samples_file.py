#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/05/20 21:41
# Author        : William GoGo
import os
import argparse
import pandas as pd
import re


def parse_input():
    argparser = argparse.ArgumentParser(description='This is a script to generate samples trinity file.')
    argparser.add_argument('-s', '--samples_described', default='samples_described.txt',
                           help='samples described file')
    argparser.add_argument('-i', '--input_dir', help='input dir', default='02_Pinjiedata')
    argparser.add_argument('-o', '--output', help='output file', default='samples_trinity.txt')
    argparser.add_argument('-t', '--type', choices=['all', 'max', 'custom'], default='planA',
                           help='all: all samples; planA: 每组中选文件最长的(指定-s); custom: custom samples')
    args = argparser.parse_args()
    if args.type == 'planA':
        if not args.samples_described:
            print('Please input samples described file.')
    return args


def main():
    args = parse_input()
    samples_file = args.samples_described
    # samples_file 读取为字典
    samples_df = pd.read_csv(samples_file, sep='\t', header=0)
    samples_dict = samples_df.groupby('group')['sample'].apply(list).to_dict()
    print(samples_dict)
    
    seq_file_lst = [x for x in os.listdir(args.input_dir) if '.fq' in x or '.fastq' in x or '.fasta' in x or '.fa' in x]
    seq_file_pair_lst = []
    for i in range(0, len(seq_file_lst), 2):
        for group in samples_dict.keys():
            group_lst = []
            for sample in samples_dict[group]:
                if sample in seq_file_lst[i] and re.search(r'{}[^a-zA-Z0-9]'.format(re.escape(sample)), seq_file_lst[i]):
                    file_1_abspath = os.path.abspath(os.path.join(args.input_dir, seq_file_lst[i]))
                    file_2_abspath = os.path.abspath(os.path.join(args.input_dir, seq_file_lst[i+1]))
                    group_lst.append([group, sample, file_1_abspath, file_2_abspath, os.path.getsize(file_1_abspath)])
            seq_file_pair_lst.append(group_lst)
            
    print(seq_file_pair_lst)
    # 每组中只选文件最大的
    if args.type == 'max':
        for group in seq_file_pair_lst:
            group.sort(key=lambda x: x[4], reverse=True)
            with open(args.output, 'a') as f:
                f.write(dir + '\t' + dir + '_rep1' + '\t' + group[0][2] + '\t' + group[1][3] + '\n')
    # group_samples_size = []
    # samples_paired = []
    # group_samples_paired = []
    # last_group = samples_dict.keys()[0]
    # for group, samples in samples_dict.items():
    #     if args.type == 'planA':
    #         for item in os.scandir(args.input):
    #             if samples in item.name and re.search(r'samples[^a-zA-Z0-9]', item.name):
    #                 samples_paired.append(item.name)
    #         if group == last_group:
    #             group == last_group
    #             group_samples_paired.append(samples_paired)
    #             samples_paired = []
    #         else:
            
    #     for item in os.scandir(args.input):
    #         if samples in item.name and re.search(r'samples[^a-zA-Z0-9]', item.name):
    #             samples_paired.append(item.name)

    #     samples_paired.sort()
        
                

    # paired = []
    # # 生成 samples trinity file
    # for item in os.scandir(args.input):
    #     if '_1.clean.fq' in dir:
    #         # 获取两个文件的绝对路径
    #         file1_path = item.path
    #         paired.append(file1_path)
    #     elif '_2.clean.fq' in dir:
    #         file2_path = item.path
    #         paired.append(file2_path)
    #     if len(paired) == 2:
    #         # 生成 samples trinity file
    #         with open(args.output, 'a') as f:
    #             f.write(dir + '\t' + dir + '_rep1' + '\t' + paired[0] + '\t' + paired[1] + '\n')
    #         paired = []


if __name__ == '__main__':
    main()
