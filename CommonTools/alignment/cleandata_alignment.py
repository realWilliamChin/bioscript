#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/20 21:45
# Author        : William GoGo
import os, sys
import pandas as pd
import argparse
import subprocess
from loguru import logger
import concurrent.futures


def parse_input():
    args = argparse.ArgumentParser()
    args.add_argument('-t', '--type', dest='alignment_type', choices=['hisat2', 'bowtie2', 'bwa', 'salmon'], help=("选择使用哪种方式比对，并确保提供必要的参数"))
    args.add_argument('-s', '--samples', required=True, help='samples_described.txt (默认: samples_described.txt)', default='samples_described.txt')
    args.add_argument('-d', '--cd', help='cleandata file dir')
    args.add_argument('-o', '--result-dir', dest='result_dir', help='输出结果文件夹')
    args.add_argument('-r', '--ref', help='reference index，运行 salmon 时只需要输入目录')
    args.add_argument('-n', '--cpu', help='单个样本比对线程数量 (默认: 4)', default=4, type=int)
    args.add_argument('-p', '--parallel', dest='parallel_num', default=3, type=int, help='同时运行的样本数量 (默认: 3)')
    args.add_argument('-m', '--msf', help='output mapping summary file (默认: alignment_report.txt)', default='alignment_report.txt')

    parsed_args = args.parse_args()
    return parsed_args


def alignment(alignment_type, sample_name, R1, R2, result_dir, ref_index, num_threads, gff_file=None):
    
        mapping_file_name = os.path.join(result_dir, f'{sample_name}_mapping.txt')
        bam_file_name = os.path.join(result_dir, f'{sample_name}.bam')
        
        if alignment_type.lower() == 'hisat2':
            alignment_cmd = f'hisat2 -x {ref_index} \
                -p {num_threads} \
                -I 200 -X 400 --fr \
                --min-intronlen 20 --max-intronlen 4000 \
                -1 {R1} \
                -2 {R2} \
                2> {mapping_file_name} | \
                samtools sort --threads {num_threads} -O BAM -o - > {bam_file_name}'
                
        elif alignment_type.lower() == 'bowtie2':
            alignment_cmd = f'bowtie2 -x {ref_index} \
                -p {num_threads} \
                -1 {R1} \
                -2 {R2} \
                2> {mapping_file_name} | \
                samtools sort --threads {num_threads} -O BAM -o - > {bam_file_name}'
                
        elif alignment_type.lower() == 'bwa':
            alignment_cmd = f'bwa mem -M \
                -t {num_threads} \
                -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tLB:WES\\tPL:Illumina" \
                {ref_index} \
                {R1} {R2} |\
                samtools sort -O bam -@ {num_threads} -o {bam_file_name}'
                
        elif alignment_type.lower() == 'salmon':
            gene_map_file = os.path.join(ref_index, 'gene_map.txt')
            salmon_output_dir = os.path.join(result_dir, f'{sample_name}')
            alignment_cmd = f'/home/data/opt/biosoft/salmon-latest_linux_x86_64/bin/salmon quant \
                -p {num_threads} \
                --validateMappings \
                -i {ref_index} \
                -l A \
                -g {gene_map_file} \
                -1 {R1} -2 {R2} \
                -o {salmon_output_dir}'
        
        rep = subprocess.run(alignment_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
        if rep.returncode != 0:
            logger.error(f'样本 {sample_name} 比对失败')
            logger.error(rep.stderr.decode())
            logger.error(rep.stdout.decode())
            
            return False
        else:
            logger.success('比对完成')
        
        if os.access(mapping_file_name, os.R_OK):
            with open(mapping_file_name) as f:
                last_line = f.readlines()[-1]
                logger.success(f'{sample_name} {last_line}')
        
        return bam_file_name


def main():
    logger.warning('比对前建议检查文件和样本是否匹配')
    args = parse_input()
    
    # 变量
    alignment_type = args.alignment_type
    ref = args.ref
    samples_file, gff_file, mapping_summary_file = args.samples, args.gff, args.msf
    cleandata_dir, result_dir = args.cd, args.result_dir
    num_threads, parallel_num = args.cpu, args.parallel_num
    
    # 创建结果文件夹
    os.makedirs(result_dir, exist_ok=True)
    
    samples_df = pd.read_csv(samples_file, sep='\t', usecols=['sample', 'R1', 'R2'])

    # 控制同时运行的样本数量
    with concurrent.futures.ProcessPoolExecutor(max_workers=parallel_num) as executor:
        # 活跃的任务列表
        active_futures = []
        # 样本索引
        sample_index = 0
        total_samples = len(samples_df)
        
        # 首先提交parallel_num个任务
        while len(active_futures) < parallel_num and sample_index < total_samples:
            row = samples_df.iloc[sample_index]
            sample_name = row['sample']
            R1 = os.path.join(cleandata_dir, row['R1'])
            R2 = os.path.join(cleandata_dir, row['R2'])
            logger.info(f'开始比对样本 {sample_name}, {R1} 和 {R2}')
            
            future = executor.submit(
                alignment, alignment_type, sample_name, R1, R2, result_dir, ref, num_threads, gff_file
            )
            active_futures.append((future, sample_name))
            sample_index += 1
        
        # 处理完成的任务并提交新任务
        while active_futures:
            # 等待任意一个任务完成
            done, not_done = concurrent.futures.wait(
                [f for f, _ in active_futures],
                return_when=concurrent.futures.FIRST_COMPLETED
            )
            
            # 处理完成的任务
            for future in done:
                # 找到对应的样本名
                for i, (f, name) in enumerate(active_futures):
                    if f == future:
                        sample_name = name
                        active_futures.pop(i)
                        break
                
                try:
                    result = future.result()
                    logger.info(f'样本 {sample_name} 比对完成，结果: {result}')
                except Exception as exc:
                    logger.error(f'样本 {sample_name} 运行时出错: {exc}')
                
                # 如果还有未处理的样本，提交新任务
                if sample_index < total_samples:
                    row = samples_df.iloc[sample_index]
                    new_sample_name = row['sample']
                    R1 = os.path.join(cleandata_dir, row['R1'])
                    R2 = os.path.join(cleandata_dir, row['R2'])
                    logger.info(f'开始比对样本 {new_sample_name}, {R1} 和 {R2}')
                    
                    new_future = executor.submit(
                        alignment, alignment_type, new_sample_name, R1, R2, result_dir, ref, num_threads, gff_file
                    )
                    active_futures.append((new_future, new_sample_name))
                    sample_index += 1


if __name__ == '__main__':
    main()