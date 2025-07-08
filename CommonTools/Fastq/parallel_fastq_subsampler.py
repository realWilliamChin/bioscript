#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/08 16:48
# Author        : William GoGo
"""
parallel_fastq_subsampler.py - 并行随机抽样FASTQ文件工具

功能：
1. 自动识别当前目录下的FASTQ.GZ文件对（R1和R2）
2. 随机抽取每对文件中指定比例的reads（默认50%）
3. 保持R1/R2配对一致性
4. 并行处理多个文件对以加速处理

使用示例：
$ python parallel_fastq_subsampler.py -p 0.3 -j 4 -o sampled_data

参数说明：
-p, --proportion: 抽样比例 (0.0-1.0, 默认0.5)
-j, --jobs: 并行任务数 (默认CPU核心数)
-o, --output: 输出目录 (默认当前目录)
--suffix: 输出文件后缀 (默认"_sampled")
"""

import gzip
import os
import random
import argparse
import sys
import time
import math
from multiprocessing import Pool, cpu_count
from functools import partial
from collections import defaultdict

def count_reads(filename):
    """统计FASTQ.GZ文件中的reads数量"""
    count = 0
    try:
        with gzip.open(filename, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 3:  # 每4行完成一个read
                    count += 1
    except Exception as e:
        print(f"ERROR: Failed to read {filename}: {str(e)}")
        return 0
    return count

def process_file_pair(file_pair, output_dir, sample_proportion, suffix):
    """处理单个FASTQ文件对"""
    r1_file, r2_file = file_pair
    
    # 统计总reads
    n1 = count_reads(r1_file)
    n2 = count_reads(r2_file)
    
    if n1 == 0 or n2 == 0:
        print(f"ERROR: Empty file detected in pair: {r1_file} or {r2_file}")
        return None
    
    if n1 != n2:
        print(f"WARNING: Read count mismatch in {r1_file} ({n1}) and {r2_file} ({n2})")
        min_reads = min(n1, n2)
        if min_reads == 0:
            print("ERROR: Cannot process pair with zero reads")
            return None
        print(f"Using minimum read count: {min_reads}")
        total_reads = min_reads
    else:
        total_reads = n1
    
    # 计算抽样数量
    sample_size = max(1, math.floor(total_reads * sample_proportion))
    
    # 生成随机索引
    indices = sorted(random.sample(range(total_reads), sample_size))
    
    # 准备输出路径
    os.makedirs(output_dir, exist_ok=True)
    r1_basename = os.path.basename(r1_file)
    r2_basename = os.path.basename(r2_file)
    
    # 处理输出文件名（保留可能的多个扩展名）
    if r1_basename.endswith('.fastq.gz'):
        r1_out = os.path.join(output_dir, r1_basename.replace('.fastq.gz', f'{suffix}.fastq.gz'))
    elif r1_basename.endswith('.fq.gz'):
        r1_out = os.path.join(output_dir, r1_basename.replace('.fq.gz', f'{suffix}.fq.gz'))
    else:
        r1_out = os.path.join(output_dir, f"{os.path.splitext(r1_basename)[0]}{suffix}.fastq.gz")
    
    if r2_basename.endswith('.fastq.gz'):
        r2_out = os.path.join(output_dir, r2_basename.replace('.fastq.gz', f'{suffix}.fastq.gz'))
    elif r2_basename.endswith('.fq.gz'):
        r2_out = os.path.join(output_dir, r2_basename.replace('.fq.gz', f'{suffix}.fq.gz'))
    else:
        r2_out = os.path.join(output_dir, f"{os.path.splitext(r2_basename)[0]}{suffix}.fastq.gz")
    
    # 执行抽样
    j = 0  # 已处理的有效索引计数
    current_read = 0  # 当前读取的read索引
    processed = 0
    
    try:
        with gzip.open(r1_file, 'rt') as f1, \
             gzip.open(r2_file, 'rt') as f2, \
             gzip.open(r1_out, 'wt') as out1, \
             gzip.open(r2_out, 'wt') as out2:
            
            while j < len(indices):
                # 读取R1的4行
                r1_lines = [f1.readline() for _ in range(4)]
                # 读取R2的4行
                r2_lines = [f2.readline() for _ in range(4)]
                
                # 检查文件结束
                if not all(r1_lines) or not all(r2_lines):
                    break
                    
                # 如果当前read在抽取列表中，则写入
                if current_read == indices[j]:
                    out1.writelines(r1_lines)
                    out2.writelines(r2_lines)
                    j += 1
                    processed += 1
                    
                current_read += 1
                
                # 进度显示（每处理1%显示一次）
                if current_read % max(1, total_reads // 100) == 0:
                    percent = current_read / total_reads * 100
                    print(f"  {r1_basename}: {percent:.1f}% processed", end='\r')
    
    except Exception as e:
        print(f"\nERROR processing {r1_file}: {str(e)}")
        # 删除可能不完整的输出文件
        if os.path.exists(r1_out):
            os.remove(r1_out)
        if os.path.exists(r2_out):
            os.remove(r2_out)
        return None
    
    print(f"\nProcessed: {r1_basename} and {r2_basename}")
    print(f"  Total reads: {total_reads:,}")
    print(f"  Sampled reads: {sample_size:,} ({sample_proportion*100:.1f}%)")
    print(f"  Output: {os.path.basename(r1_out)} and {os.path.basename(r2_out)}\n")
    
    return (r1_file, r2_file, total_reads, sample_size, r1_out, r2_out)

def main():
    parser = argparse.ArgumentParser(description='Parallel FASTQ Subsampler')
    parser.add_argument('-p', '--proportion', type=float, default=0.5,
                        help='抽样比例 (0.0-1.0, 默认0.5)')
    parser.add_argument('-j', '--jobs', type=int, default=cpu_count(),
                        help=f'并行任务数 (默认: {cpu_count()})')
    parser.add_argument('-o', '--output', default='.',
                        help='输出目录 (默认当前目录)')
    parser.add_argument('--suffix', default='_sampled',
                        help='输出文件后缀 (默认"_sampled")')
    args = parser.parse_args()
    
    # 验证参数
    if args.proportion <= 0 or args.proportion > 1:
        print("ERROR: Proportion must be between 0 and 1")
        sys.exit(1)
    
    if args.jobs < 1:
        args.jobs = 1
    
    start_time = time.time()
    
    print("\n" + "="*60)
    print(f"Parallel FASTQ Subsampler - Sampling {args.proportion*100:.1f}% of reads")
    print(f"Parallel jobs: {args.jobs}")
    print(f"Output directory: {args.output}")
    print(f"Output suffix: {args.suffix}")
    print("="*60 + "\n")
    
    # 获取所有FASTQ.GZ文件
    files = [f for f in os.listdir() if f.endswith(('.fastq.gz', '.fq.gz'))]
    
    if not files:
        print("ERROR: No FASTQ.GZ files found in current directory")
        sys.exit(1)
    
    # 按样本分组文件
    sample_files = defaultdict(list)
    for f in files:
        # 尝试识别常见命名模式
        if '_R1' in f:
            sample_name = f.split('_R1')[0]
            sample_files[sample_name].append(('R1', f))
        elif '_R2' in f:
            sample_name = f.split('_R2')[0]
            sample_files[sample_name].append(('R2', f))
        elif '_1' in f:
            sample_name = f.split('_1')[0]
            sample_files[sample_name].append(('R1', f))
        elif '_2' in f:
            sample_name = f.split('_2')[0]
            sample_files[sample_name].append(('R2', f))
        else:
            # 无法识别的文件
            print(f"WARNING: Could not classify file: {f}")
    
    # 创建文件对列表
    file_pairs = []
    for sample, files in sample_files.items():
        r1_files = [f[1] for f in files if f[0] == 'R1']
        r2_files = [f[1] for f in files if f[0] == 'R2']
        
        if not r1_files or not r2_files:
            print(f"WARNING: Incomplete pair for sample {sample}")
            continue
        
        # 对于有多个lane的情况，只取第一个匹配对
        if len(r1_files) == 1 and len(r2_files) == 1:
            file_pairs.append((r1_files[0], r2_files[0]))
        else:
            # 尝试匹配多个文件
            for r1 in r1_files:
                # 寻找匹配的R2
                base_name = os.path.basename(r1)
                possible_r2 = base_name.replace('_R1', '_R2').replace('_1', '_2')
                
                if possible_r2 in r2_files:
                    file_pairs.append((r1, possible_r2))
                else:
                    print(f"WARNING: Could not find matching R2 for {r1}")
    
    if not file_pairs:
        print("ERROR: No valid file pairs found")
        sys.exit(1)
    
    print(f"Found {len(file_pairs)} FASTQ file pairs to process:")
    for r1, r2 in file_pairs:
        print(f"  {r1} -> {r2}")
    print()
    
    # 设置并行处理
    process_func = partial(process_file_pair, 
                          output_dir=args.output,
                          sample_proportion=args.proportion,
                          suffix=args.suffix)
    
    with Pool(processes=args.jobs) as pool:
        results = pool.map(process_func, file_pairs)
    
    # 汇总结果
    success_count = 0
    total_input_reads = 0
    total_output_reads = 0
    
    print("\n" + "="*60)
    print("Processing Summary:")
    for res in results:
        if res:
            success_count += 1
            _, _, input_reads, output_reads, r1_out, r2_out = res
            total_input_reads += input_reads
            total_output_reads += output_reads
            print(f"  {os.path.basename(r1_out)}: {input_reads:,} -> {output_reads:,} reads")
    
    print("\nSummary Statistics:")
    print(f"  Processed pairs: {success_count}/{len(file_pairs)}")
    print(f"  Total input reads: {total_input_reads:,}")
    print(f"  Total sampled reads: {total_output_reads:,}")
    print(f"  Sampling rate: {total_output_reads/total_input_reads*100:.2f}%")
    
    elapsed = time.time() - start_time
    print(f"\nTotal time: {elapsed:.2f} seconds")
    print("="*60)

if __name__ == "__main__":
    main()