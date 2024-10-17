#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2024/09/27 11:01
# Author        : William GoGo
import os, sys
import pandas as pd
import subprocess
import argparse
import fastq_read_write as fastq
import time
from loguru import logger


def parse_input():
    p = argparse.ArgumentParser(description="usage: python fastq_stat.py <fastq_file1> <fastq_file2>")
    p.add_argument(dest='fastq_files', nargs='+', help='输入 fastq 文件')
    p.add_argument('-o', '--output', default='fastq_stat_summary.csv',
                   help='输出 fastq_stat_summary.csv 文件')
    
    return p.parse_args()


def convert_bases_to_gb(total_bases):
    # 1 GB = 1,073,741,824 字节
    bytes_in_gb = 1073741824
    return total_bases / bytes_in_gb
    

def qual_stat(qstr):
    q20 = 0
    q30 = 0
    for q in qstr:
        qual = ord(q) - 33
        if qual >= 30:
            q30 += 1
            q20 += 1
        elif qual >= 20:
            q20 += 1
    return q20, q30


def q20_q30_stat(filename):
    reader = fastq.Reader(filename)
    total_count = 0
    q20_count = 0
    q30_count = 0
    while True:
        read = reader.nextRead()
        if read == None:
            break
        total_count += len(read[3])
        q20, q30 = qual_stat(read[3])
        q20_count += q20
        q30_count += q30

    total_count_gb = round(convert_bases_to_gb(total_count), 2)
    q20_percents = round(100 * float(q20_count)/float(total_count), 2)
    q30_percents = round(100 * float(q30_count)/float(total_count), 2)
    
    return {
        'total_base': total_count_gb,
        # 'q20_base': q20_count,
        # 'q30_base': q30_count,
        '%q20': q20_percents,
        '%q30': q30_percents
    }


def get_gc_stat(file_path):
    try:
        result = subprocess.run(f"seqtk fqchk {file_path} | head -n 3 | tail -n 2",
                                shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise Exception(f"Error running seqkit: {result.stderr}")

        output = result.stdout.strip().split('\n')
        
        header = output[0].split()
        values = output[1].split()
        
        g_value = float(values[4])
        c_value = float(values[3])
        gc_value = round(g_value + c_value, 2)
        
        return {r'%GC': gc_value}
    
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        return None


def get_seqkit_stat(file_path):
    try:
        result = subprocess.run(['seqkit', 'stat', file_path], capture_output=True, text=True)
        
        if result.returncode != 0:
            raise Exception(f"Error running seqkit: {result.stderr}")

        output = result.stdout.strip().split('\n')

        # 解析第二行数据
        headers = output[0].split()
        values = output[1].split()

        return {
            headers[3]: int(values[3].replace(',', '')),  # num_seqs
            headers[4]: int(values[4].replace(',', '')),  # sum_len
            headers[5]: int(values[5]),  # min_len
            headers[6]: float(values[6]),  # avg_len
            headers[7]: int(values[7])   # max_len
        }

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        return None


def main():
    args = parse_input()
    fastq_files = args.fastq_files
    all_stats = []
    for fastq_file in fastq_files:
        logger.info(f'calculating {fastq_file}')
        fastq_all_stat = {
            'File': fastq_file,
            **q20_q30_stat(fastq_file),
            **get_seqkit_stat(fastq_file),
            **get_gc_stat(fastq_file)
        }
        all_stats.append(fastq_all_stat)
    df = pd.DataFrame(all_stats)
    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    time1 = time.time()
    main()
    time2 = time.time()
    logger.info('Time used: ' + str(time2-time1))