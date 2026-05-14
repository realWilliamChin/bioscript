#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/03/16 11:23
# Author        : William GoGo

import argparse
import subprocess
import os
import sys
from loguru import logger

# 导入通用多线程任务运行器
sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from multithreads_task_runner import run_multithreads_tasks

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='Calculate statistics for multiple BAM files')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input BAM files (multiple allowed)')
    parser.add_argument('-o', '--output-dir', default='.', help='Output directory, default: current directory')
    parser.add_argument('-p', '--parallel', default=10, type=int, help='Number of parallel samtools processes, default: 10')
    
    args = parser.parse_args()
    
    # 参数校验
    for bam in args.input:
        if not os.path.exists(bam):
            logger.error(f"BAM file {bam} does not exist")
            sys.exit(1)
    
    # 创建输出目录
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    return args


def check_bam_index(bam_file):
    """检查BAM索引是否存在，不存在则创建"""
    bai_file = f"{bam_file}.bai"
    if not os.path.exists(bai_file):
        logger.info(f"Index file not found for {bam_file}, creating index...")
        try:
            subprocess.run(['samtools', 'index', bam_file], check=True)
            logger.info(f"Successfully created index for {bam_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create index for {bam_file}: {e}")
            return False
    else:
        logger.debug(f"Index file {bai_file} already exists")
    return True

def parse_flagstat(flagstat_file):
    """解析samtools flagstat输出文件"""
    with open(flagstat_file, 'r') as f:
        lines = f.readlines()
    
    all_reads = int(lines[0].split()[0])
    paired_mapped = int(lines[9].split()[0])
    single_mapped = int(lines[10].split()[0])
    unmapped = all_reads - single_mapped - paired_mapped
    
    if all_reads > 0:
        align_ratio = round((single_mapped + paired_mapped) * 100.0 / all_reads, 2)
    else:
        align_ratio = 0.00
    
    return {
        'all_reads': all_reads,
        'single_mapped': single_mapped,
        'paired_mapped': paired_mapped,
        'unmapped': unmapped,
        'align_ratio': align_ratio
    }


def process_single_bam(bam_file, output_dir):
    """处理单个BAM文件：运行flagstat并返回统计数据"""
    sample_name = os.path.basename(bam_file).rsplit('.', 1)[0]
    stat_file = os.path.join(output_dir, f"{sample_name}_mapping_stat.txt")
    
    logger.info(f"Running samtools flagstat on {bam_file}")
    try:
        with open(stat_file, 'w') as f:
            subprocess.run(['samtools', 'flagstat', bam_file], stdout=f, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run flagstat on {bam_file}: {e}")
        return None
    
    stats = parse_flagstat(stat_file)
    return {
        'sample': sample_name,
        'bam_file': bam_file,
        **stats
    }


def process_single_chrom_stat(bam_file):
    """处理单个BAM文件的染色体统计"""
    sample_name = os.path.basename(bam_file).rsplit('.', 1)[0]
    logger.info(f"Calculating chromosome statistics for {bam_file}")
    try:
        # 使用samtools idxstats获取染色体统计
        result = subprocess.run(['samtools', 'idxstats', bam_file],
                              capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')
        
        sample_data = {}
        chrom_set = set()
        for line in lines:
            if line.strip():
                fields = line.split()
                if len(fields) == 4:
                    chrom, length, mapped, unmapped = fields
                    if chrom != '*':  # 跳过未比对的
                        sample_data[chrom] = {
                            'length': length,
                            'mapped': mapped,
                            'unmapped': unmapped
                        }
                        chrom_set.add(chrom)
        
        return {
            'sample': sample_name,
            'bam_file': bam_file,
            'data': sample_data,
            'chrom_set': chrom_set
        }
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to calculate chromosome statistics for {bam_file}: {e}")
        return None


def run_flagstat(bam_files, output_dir, parallel_num):
    """运行flagstat并生成汇总报告，使用多进程并行处理"""
    summary_file = os.path.join(output_dir, 'Mapping_summary.txt')
    
    # 准备所有任务的参数
    tasks_args = [(bam, output_dir) for bam in bam_files]
    
    # 使用多线程运行器并行处理
    logger.info(f"Running flagstat with {parallel_num} parallel processes...")
    results = run_multithreads_tasks(process_single_bam, tasks_args, max_workers=parallel_num)
    
    # 收集成功的结果
    summary_data = [result for result in results if result is not None]
    
    # 检查是否有失败的任务
    failed_count = len(results) - len(summary_data)
    if failed_count > 0:
        logger.error(f"{failed_count} BAM files failed to process during flagstat")
        sys.exit(1)
    
    # 写入汇总文件
    with open(summary_file, 'w') as f:
        f.write("Sample\tAllReads\tSingleMappedReads\tPairedMappedReads\tUnmappedReads\tAlignRatio(%)\n")
        for data in summary_data:
            f.write(f"{data['sample']}\t{data['all_reads']}\t{data['single_mapped']}\t"
                   f"{data['paired_mapped']}\t{data['unmapped']}\t{data['align_ratio']}%\n")
    
    logger.info(f"Mapping summary saved to {summary_file}")
    return summary_data

def chromosome_stat(bam_files, output_dir, parallel_num, flagstat_summary):
    """统计所有BAM文件每个染色体的reads数，汇总到一个表格，使用多进程并行处理"""
    chr_summary_file = os.path.join(output_dir, 'Chromosome_stat_summary.txt')
    all_chrom_data = {}
    chrom_set = set()
    
    # 准备所有任务的参数
    tasks_args = [(bam,) for bam in bam_files]
    
    # 使用多线程运行器并行处理
    logger.info(f"Running chromosome statistics with {parallel_num} parallel processes...")
    results = run_multithreads_tasks(process_single_chrom_stat, tasks_args, max_workers=parallel_num)
    
    # 收集成功的结果
    valid_results = [result for result in results if result is not None]
    
    # 检查是否有失败的任务
    failed_count = len(results) - len(valid_results)
    if failed_count > 0:
        logger.error(f"{failed_count} BAM files failed to process during chromosome statistics")
        sys.exit(1)
    
    # 收集所有染色体数据，并从flagstat结果获取总reads数
    for result in valid_results:
        sample_name = result['sample']
        all_chrom_data[sample_name] = result['data']
        chrom_set.update(result['chrom_set'])
    
    # 从flagstat结果获取每个样本的总reads数
    sample_total_reads = {}
    for summary in flagstat_summary:
        sample_total_reads[summary['sample']] = summary['all_reads']
    
    # 写入汇总文件：Sample\tChromosome\tLength\tPairedMappedReads\tSingleMappedReads\tUnmappedReads\tAlignRatio(%)
    with open(chr_summary_file, 'w') as f:
        f.write("Sample\tChromosome\tLength\tPairedMappedReads\tSingleMappedReads\tUnmappedReads\tAlignRatio(%)\n")
        for sample_name, sample_data in all_chrom_data.items():
            for chrom, data in sample_data.items():
                # 根据需求重命名和计算：
                # mapped -> PairedMappedReads (因为一对reads占据两个比对结果)
                # unmapped -> SingleMappedReads
                # UnmappedReads = 染色体长度 - PairedMappedReads - SingleMappedReads
                paired_mapped = int(data['mapped']) // 2
                single_mapped = int(data['unmapped'])
                unmapped = 0
                if unmapped < 0:
                    unmapped = 0
                    
                # AlignRatio = (PairedMappedReads + SingleMappedReads) / 该样本所有染色体总reads数
                total_all = int(sample_total_reads[sample_name])
                if total_all > 0:
                    align_ratio = round((paired_mapped + single_mapped) * 100.0 / total_all, 2)
                else:
                    align_ratio = 0.00
                
                f.write(f"{sample_name}\t{chrom}\t{data['length']}\t{paired_mapped}\t{single_mapped}\t{unmapped}\t{align_ratio}\n")
    
        logger.info(f"Chromosome summary statistics saved to {chr_summary_file}")
        logger.info(f"Total chromosomes: {len(chrom_set)}")

def process_single_coverage(bam_file, output_dir):
    """处理单个BAM文件的覆盖度统计"""
    sample_name = os.path.basename(bam_file).rsplit('.', 1)[0]
    cov_stat_file = os.path.join(output_dir, f"{sample_name}_cov_stat.txt")
    
    logger.info(f"Calculating coverage statistics for {bam_file}")
    try:
        # 运行bedtools genomecov计算每个碱基的覆盖深度，流式处理输出
        bedtools_cmd = ['bedtools', 'genomecov', '-d', '-ibam', bam_file]
        with subprocess.Popen(bedtools_cmd, stdout=subprocess.PIPE, text=True, bufsize=1) as p:
            total = 0
            cov_ge_1 = 0
            cov_ge_4 = 0
            cov_ge_10 = 0
            cov_ge_20 = 0
            cov_ge_30 = 0
            
            # 逐行流式处理，内存占用极低
            for line in p.stdout:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 3:
                    continue
                try:
                    depth = int(parts[2])
                except ValueError:
                    continue
                total += 1
                if depth >= 1:
                    cov_ge_1 += 1
                if depth >= 4:
                    cov_ge_4 += 1
                if depth >= 10:
                    cov_ge_10 += 1
                if depth >= 20:
                    cov_ge_20 += 1
                if depth >= 30:
                    cov_ge_30 += 1
            
            # 检查子进程是否正常结束
            if p.wait() != 0:
                logger.error(f"bedtools genomecov failed for {bam_file} with exit code {p.returncode}")
                return None
        
        # 计算百分比
        if total > 0:
            p1 = round(100 * cov_ge_1 / total, 2)
            p4 = round(100 * cov_ge_4 / total, 2)
            p10 = round(100 * cov_ge_10 / total, 2)
            p20 = round(100 * cov_ge_20 / total, 2)
            p30 = round(100 * cov_ge_30 / total, 2)
        else:
            p1 = p4 = p10 = p20 = p30 = 0.00
        
        cov_values = [str(p1), str(p4), str(p10), str(p20), str(p30)]
        
        # 保存到样本单独的统计文件
        with open(cov_stat_file, 'w') as f:
            f.write("cov_ge_1(%)\tcov_ge_4(%)\tcov_ge_10(%)\tcov_ge_20(%)\tcov_ge_30(%)\n")
            f.write('\t'.join(cov_values) + '\n')
        
        # 返回统计结果
        return {
            'sample': sample_name,
            'cov_ge_1': p1,
            'cov_ge_4': p4,
            'cov_ge_10': p10,
            'cov_ge_20': p20,
            'cov_ge_30': p30
        }
        
    except Exception as e:
        logger.error(f"Failed to calculate coverage statistics for {bam_file}: {str(e)}", exc_info=True)
        return None

def run_coverage_stat(bam_files, output_dir, parallel_num):
    """运行覆盖度统计并生成汇总报告，使用多进程并行处理"""
    summary_file = os.path.join(output_dir, 'Coverage_summary.txt')
    
    # 准备所有任务的参数
    tasks_args = [(bam, output_dir) for bam in bam_files]
    
    # 使用多线程运行器并行处理
    logger.info(f"Running coverage statistics with {parallel_num} parallel processes...")
    results = run_multithreads_tasks(process_single_coverage, tasks_args, max_workers=parallel_num)
    
    # 收集成功的结果
    summary_data = [result for result in results if result is not None]
    
    # 检查是否有失败的任务
    failed_count = len(results) - len(summary_data)
    if failed_count > 0:
        logger.error(f"{failed_count} BAM files failed to process during coverage statistics")
        sys.exit(1)
    
    # 写入汇总文件
    with open(summary_file, 'w') as f:
        f.write("Sample\tcov_ge_1(%)\tcov_ge_4(%)\tcov_ge_10(%)\tcov_ge_20(%)\tcov_ge_30(%)\n")
        for data in summary_data:
            f.write(f"{data['sample']}\t{data['cov_ge_1']}\t{data['cov_ge_4']}\t"
                   f"{data['cov_ge_10']}\t{data['cov_ge_20']}\t{data['cov_ge_30']}\n")
    
    logger.info(f"Coverage summary saved to {summary_file}")
    return summary_data

def main():
    args = parse_args()
    
    logger.info(f"Starting BAM statistics calculation")
    logger.info(f"Input BAM files: {len(args.input)} BAM files")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Parallel processes: {args.parallel}")
    
    # 检查并创建索引（这里索引创建比较快，单线程足够，samtools index本身也是单线程）
    logger.info("Checking BAM indexes...")
    for bam in args.input:
        if not check_bam_index(bam):
            sys.exit(1)
    
    # 运行flagstat统计，并行处理
    logger.info("Running mapping statistics...")
    flagstat_summary = run_flagstat(args.input, args.output_dir, args.parallel)
    
    # 染色体统计，并行处理
    logger.info("Running chromosome-level statistics...")
    chromosome_stat(args.input, args.output_dir, args.parallel, flagstat_summary)
    
    # 覆盖度统计，并行处理
    logger.info("Running coverage statistics...")
    run_coverage_stat(args.input, args.output_dir, args.parallel)
    
    logger.info("=" * 60)
    logger.info("BAM statistics completed successfully!")
    logger.info("=" * 60)

if __name__ == "__main__":
    main()
