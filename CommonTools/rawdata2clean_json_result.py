#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/12/12 16:30
# Author        : William GoGo
import json
import os
import argparse


def parse_input():
    args = argparse.ArgumentParser(description='rawdata2clean_json_result.py')
    args.add_argument('-s', '--samples', help='samples_described.txt')
    args.add_argument('-i', '--input', help='input rawdata dir', default='.')
    args.add_argument('-o', '--output', help='output json file', default='测序质量表.txt')
    
    return args.parse_args()


def parse_json(json_file):
    with open(json_file, 'r') as f:
        json_data = json.load(f)
    alfter_filtering = json_data['summary']['after_filtering']
    clean_reads = int(alfter_filtering['total_reads'] / 2)
    clean_base = round(float(alfter_filtering['total_bases'] / (10 ** 9)), 2)
    q20 = str(round(alfter_filtering['q20_rate'] * 100, 2)) + '%'
    q30 = str(round(alfter_filtering['q30_rate'] * 100, 2)) + '%'
    gc = str(round(alfter_filtering['gc_content'] * 100, 2)) + '%'
    return clean_reads, clean_base, q20, q30, gc


def main():
    args = parse_input()
    if args.samples:
        samples_data = open(args.samples, 'r').readlines()
        sample_name_list = [x.split('\t')[1].strip() for x in samples_data[1:] if x.strip() != '' and x.split('\t')[1].strip() != '']
        sample_data_list = [os.sep.join([args.input, x.split('\t')[1].strip() + '_fastp.json']) for x in samples_data[1:] if x.strip() != '']
    else:
        sample_name_list = [x.replace('_fastp.json', '') for x in os.listdir(args.input) if x.endswith('_fastp.json')]
        sample_data_list = [os.sep.join([args.input, x]) for x in os.listdir(args.input) if x.endswith('_fastp.json')]
    open(args.output, 'w').write('Sample\tClean_reads\tClean_base\tQ20\tQ30\tGC\n')
    sample_list = zip(sample_name_list, sample_data_list)
    for name, json_data in sample_list:
        if not os.path.exists(json_data):
            print(f'{json_data} not exists!')
            continue
        clean_reads, clean_base, q20, q30, gc = parse_json(json_data)
        # sample_name = os.path.basename(json_data).replace('_fastp.json', '')
        open(args.output, 'a').write(f'{name}\t{clean_reads}\t{clean_base}\t{q20}\t{q30}\t{gc}\n')
    
    print('\nDone!')


if __name__ == '__main__':
    main()
