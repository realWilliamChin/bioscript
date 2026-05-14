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
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from multithreads_task_runner import run_multithreads_tasks


def parse_input():
    p = argparse.ArgumentParser(description="usage: python fastq_stat.py <fastq_file1> <fastq_file2>")
    p.add_argument(dest='fastq_files', nargs='+', help='输入 fastq 文件')
    p.add_argument('-o', '--output', default='fastq_stat_summary.csv',
                   help='输出 fastq_stat_summary.csv 文件')
    p.add_argument('-p', '--threads', type=int, default=3,
                   help='并行进程数，默认 3')
    
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


def process_single_fastq(fastq_file):
    """处理单个 fastq 文件的统计"""
    fastq_all_stat = {
        'File': fastq_file,
        **q20_q30_stat(fastq_file),
        **get_seqkit_stat(fastq_file),
        **get_gc_stat(fastq_file)
    }
    return fastq_all_stat


def merge_pe_stats(stats_list):
    """合并双端测序 R1/R2 的统计结果"""
    if len(stats_list) < 2:
        return stats_list
    
    # 按文件名分组，识别 R1/R2 配对
    pe_groups = {}
    import re
    for stat in stats_list:
        filename = stat['File']
        basename = os.path.basename(filename)
        
        # 移除文件扩展名
        name_without_ext = re.sub(r'\.(fq|fastq)(\.gz)?$', '', basename, flags=re.IGNORECASE)
        
        # 识别 R1/R2 模式，支持:
        # _R1, _R2, .R1, .R2, _1, _2, .1, .2
        # 匹配末尾的 R1/R2 或 1/2 标识
        r1_match = re.search(r'[._](R)?1$', name_without_ext, flags=re.IGNORECASE)
        r2_match = re.search(r'[._](R)?2$', name_without_ext, flags=re.IGNORECASE)
        
        # 提取样本名
        if r1_match:
            sample_name = name_without_ext[:r1_match.start()]
            is_r1 = True
            is_r2 = False
        elif r2_match:
            sample_name = name_without_ext[:r2_match.start()]
            is_r1 = False
            is_r2 = True
        else:
            # 没有识别到 R1/R2 标记，作为单端处理
            continue
        
        if sample_name not in pe_groups:
            pe_groups[sample_name] = {'R1': None, 'R2': None}
        
        if is_r1:
            pe_groups[sample_name]['R1'] = stat
        elif is_r2:
            pe_groups[sample_name]['R2'] = stat
    
    # 合并 R1/R2 统计
    merged_stats = []
    for sample_name, pair in pe_groups.items():
        r1_stat = pair['R1']
        r2_stat = pair['R2']
        
        if r1_stat and r2_stat:
            # 合并统计
            merged = {
                'File': sample_name + '_PE',
                'total_base': round(r1_stat['total_base'] + r2_stat['total_base'], 2),
                '%q20': round((r1_stat['%q20'] + r2_stat['%q20']) / 2, 2),
                '%q30': round((r1_stat['%q30'] + r2_stat['%q30']) / 2, 2),
                'num_seqs': r1_stat['num_seqs'] + r2_stat['num_seqs'],
                'sum_len': r1_stat['sum_len'] + r2_stat['sum_len'],
                'min_len': min(r1_stat['min_len'], r2_stat['min_len']),
                'avg_len': round((r1_stat['avg_len'] + r2_stat['avg_len']) / 2, 2),
                'max_len': max(r1_stat['max_len'], r2_stat['max_len']),
                '%GC': round((r1_stat['%GC'] + r2_stat['%GC']) / 2, 2),
            }
            merged_stats.append(merged)
    
    return stats_list + merged_stats


def main():
    args = parse_input()
    fastq_files = args.fastq_files
    
    # 构建任务参数列表
    tasks_args = [(f,) for f in fastq_files]
    
    # 并行执行任务
    all_stats = run_multithreads_tasks(
        task_func=process_single_fastq,
        tasks_args=tasks_args,
        max_workers=args.threads,
        show_log=True
    )
    
    # 过滤掉 None（出错的任务）
    all_stats = [s for s in all_stats if s is not None]
    
    # 合并双端测序统计
    all_stats = merge_pe_stats(all_stats)
    
    df = pd.DataFrame(all_stats)
    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    time1 = time.time()
    main()
    time2 = time.time()
    logger.info('Time used: ' + str(time2-time1))