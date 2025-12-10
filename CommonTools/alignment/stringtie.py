#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/08/14 11:53
# Author        : William GoGo
import os, sys
import pandas as pd
import argparse
import subprocess
from loguru import logger
import concurrent.futures

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='使用StringTie进行转录本组装和定量分析')
    parser.add_argument('-i', '--input', required=True, help='bam 文件目录')
    parser.add_argument('-s', '--samples', required=True, help='samples_described.txt')
    parser.add_argument('-g', '--gff', required=True, help='参考基因组GFF/GTF文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    parser.add_argument('-t', '--threads', type=int, default=4, help='并行线程数，默认4')
    
    args = parser.parse_args()
    
    return args


def stringtie(sample_name, bam_file_name, gff_file, result_dir):
    """运行StringTie命令"""
    # 检查输入文件是否存在
    if not os.path.exists(bam_file_name):
        logger.error(f'BAM文件不存在: {bam_file_name}')
        return False
    
    if not os.path.exists(gff_file):
        logger.error(f'GFF文件不存在: {gff_file}')
        return False
    
    # 创建输出目录
    fpkm_dir = os.path.join(result_dir, 'fpkm')
    ballgown_dir = os.path.join(result_dir, 'ballgown', sample_name)
    os.makedirs(fpkm_dir, exist_ok=True)
    os.makedirs(ballgown_dir, exist_ok=True)

    stringtie_cmd = (
        f'stringtie -e -B '
        f'-G {gff_file} '
        f'-A {result_dir}/fpkm/{sample_name}_fpkm.txt '
        f'-o {result_dir}/ballgown/{sample_name}/{sample_name}.gtf '
        f'{bam_file_name}'
    )
    
    logger.info(f'运行StringTie命令: {stringtie_cmd}')
    
    try:
        # 使用subprocess运行命令，不使用nohup和后台运行
        result = subprocess.run(stringtie_cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info(f'样本 {sample_name} 处理完成')
            return True
        else:
            logger.error(f'样本 {sample_name} 处理失败: {result.stderr}')
            return False
            
    except Exception as e:
        logger.error(f'样本 {sample_name} 处理异常: {str(e)}')
        return False


def process_samples_parallel(samples_df, input_dir, gff_file, result_dir, max_workers=1):
    """并行处理所有样本"""
    logger.info(f'开始并行处理 {len(samples_df)} 个样本，使用 {max_workers} 个线程')
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 创建任务映射
        future_to_sample = {}
        for _, row in samples_df.iterrows():
            sample_name = row['sample']
            bam_file = os.path.join(input_dir, f"{sample_name}.bam")
            future = executor.submit(
                stringtie, 
                sample_name, 
                bam_file, 
                gff_file, 
                result_dir
            )
            future_to_sample[future] = sample_name
        
        # 收集结果
        results = []
        for future in concurrent.futures.as_completed(future_to_sample):
            sample_name = future_to_sample[future]
            try:
                success = future.result()
                results.append((sample_name, success))
            except Exception as e:
                logger.error(f'样本 {sample_name} 处理异常: {str(e)}')
                results.append((sample_name, False))
    
    # 统计结果
    successful = sum(1 for _, success in results if success)
    failed = len(results) - successful
    
    logger.info(f'处理完成: 成功 {successful} 个样本，失败 {failed} 个样本')
    
    return results


def main():
    args = parse_input()
    samples_df = load_table(args.samples)
    
    # 并行处理样本
    results = process_samples_parallel(samples_df, args.input, args.gff, args.output, args.threads)


if __name__ == '__main__':
    main()