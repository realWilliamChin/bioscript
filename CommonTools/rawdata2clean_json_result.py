#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/12/12 16:30
# Author        : William GoGo
import json
import os
import glob
import argparse
import subprocess


def parse_input():
    args = argparse.ArgumentParser(description='rawdata 质控合并脚本: 运行 fastp 生成 cleandata 并汇总测序质量表。'
                                               '可基于 samples_described.txt, 或在未提供样本表时自动发现 rawdata 目录下的双端 fastq 文件。')
    args.add_argument('-s', '--samples', help='samples_described.txt (含 group/sample 列; 运行 fastp 时还需 R1/R2 列)。'
                                              '不提供时改为自动扫描 --input 目录下的双端 fastq 文件', required=False)
    args.add_argument('-i', '--input', help='rawdata 目录, 自动发现模式下从此目录递归查找 *1.fastq.gz / *1.fq.gz', default='.')
    args.add_argument('-o', '--output', help='output json file', default='测序质量表.txt')
    args.add_argument('--run', help='运行模式: beirui 或 nuohuo (需要运行 fastp 时必需; 自动发现模式默认 nuohuo)', choices=['beirui', 'nuohuo'], required=False)

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


def run_fastp_beirui(f1, f2, cleandata_dir, key_name, fastp_report_dir):
    """运行 beirui 模式的 fastp 命令"""
    os.makedirs(cleandata_dir, exist_ok=True)
    os.makedirs(fastp_report_dir, exist_ok=True)
    
    output_r1 = os.path.join(cleandata_dir, f'{key_name}_R1.clean.fq.gz')
    output_r2 = os.path.join(cleandata_dir, f'{key_name}_R2.clean.fq.gz')
    json_output = os.path.join(fastp_report_dir, f'{key_name}_fastp.json')
    
    cmd = [
        'fastp',
        '-q', '5',
        '-u', '20',
        '-l', '150',
        '--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
        '--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
        '-i', f1,
        '-I', f2,
        '-o', output_r1,
        '-O', output_r2,
        '-j', json_output
    ]
    
    print(f'运行 fastp (beirui模式): {key_name}')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'错误: fastp 运行失败 - {key_name}')
        print(result.stderr)
        return False, None, None
    return True, output_r1, output_r2


def run_fastp_nuohuo(f1, f2, cleandata_dir, key_name, fastp_report_dir):
    """运行 nuohuo 模式的 fastp 命令"""
    os.makedirs(cleandata_dir, exist_ok=True)
    os.makedirs(fastp_report_dir, exist_ok=True)
    
    output_r1 = os.path.join(cleandata_dir, f'{key_name}_R1.clean.fq.gz')
    output_r2 = os.path.join(cleandata_dir, f'{key_name}_R2.clean.fq.gz')
    json_output = os.path.join(fastp_report_dir, f'{key_name}_fastp.json')
    
    cmd = [
        'fastp',
        '-G',
        '-q', '5',
        '-u', '50',
        '-A',
        '-n', '15',
        '-l', '150',
        '--overlap_diff_limit', '1',
        '--overlap_diff_percent_limit', '10',
        '-i', f1,
        '-I', f2,
        '-o', output_r1,
        '-O', output_r2,
        '-j', json_output
    ]
    
    print(f'运行 fastp (nuohuo模式): {key_name}')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'错误: fastp 运行失败 - {key_name}')
        print(result.stderr)
        return False, None, None
    return True, output_r1, output_r2




def discover_samples(input_dir):
    """自动发现 input_dir 下的双端 fastq 文件 (复刻原 rawdata2cleandata.sh 的 find 行为),
    返回 [{'group','sample','R1','R2'}, ...]。同时兼容 .fastq.gz 与 .fq.gz。"""
    r1_files = []
    for pattern in ('*1.fastq.gz', '*1.fq.gz'):
        r1_files.extend(glob.glob(os.path.join(input_dir, '**', pattern), recursive=True))
    r1_files = sorted(set(r1_files))

    samples = []
    for r1 in r1_files:
        # 推导 R2: 将文件名中的 1.fastq.gz / 1.fq.gz 替换为 2.*
        if r1.endswith('1.fastq.gz'):
            r2 = r1[:-len('1.fastq.gz')] + '2.fastq.gz'
            suffix = '1.fastq.gz'
        else:
            r2 = r1[:-len('1.fq.gz')] + '2.fq.gz'
            suffix = '1.fq.gz'
        if not os.path.exists(r2):
            print(f'警告: 未找到 {r1} 对应的 R2 文件 {r2}，跳过')
            continue
        # 样本名: 去掉路径与 _1.fastq.gz / _1.fq.gz 后缀
        base = os.path.basename(r1)
        sample = base[:-len(suffix)].rstrip('_.')
        samples.append({'group': sample, 'sample': sample, 'R1': r1, 'R2': r2})
    return samples


def run_fastp(run_mode, f1, f2, cleandata_dir, key_name, fastp_report_dir):
    """按模式分派 fastp 运行。"""
    if run_mode == 'beirui':
        return run_fastp_beirui(f1, f2, cleandata_dir, key_name, fastp_report_dir)
    return run_fastp_nuohuo(f1, f2, cleandata_dir, key_name, fastp_report_dir)


def main():
    args = parse_input()

    # 创建输出目录
    cleandata_dir = './cleandata'
    fastp_report_dir = './fastp_report'

    # 自动发现模式: 未提供样本表 -> 扫描 input 目录，必定需要运行 fastp
    if not args.samples:
        run_mode = args.run or 'nuohuo'
        print(f'未提供 --samples，进入自动发现模式 (运行模式: {run_mode})，从 {args.input} 扫描双端 fastq 文件')
        discovered = discover_samples(args.input)
        if not discovered:
            print(f'错误: 在 {args.input} 下未发现可配对的双端 fastq 文件')
            return

        processed_samples = []
        sample_name_list = []
        sample_data_list = []
        for info in discovered:
            success, new_r1, new_r2 = run_fastp(run_mode, info['R1'], info['R2'], cleandata_dir, info['sample'], fastp_report_dir)
            if not success:
                continue
            processed_samples.append({'group': info['group'], 'sample': info['sample'], 'R1': new_r1, 'R2': new_r2})
            json_file = os.path.join(fastp_report_dir, f"{info['sample']}_fastp.json")
            if os.path.exists(json_file):
                sample_name_list.append(info['sample'])
                sample_data_list.append(json_file)

        if processed_samples:
            with open('samples_described.txt', 'w') as f:
                f.write('group\tsample\tR1\tR2\n')
                for s in processed_samples:
                    f.write(f"{s['group']}\t{s['sample']}\t{s['R1']}\t{s['R2']}\n")
            print('\n已生成新的样本文件: samples_described.txt')

        write_quality_table(args.output, sample_name_list, sample_data_list)
        return

    # 检查 fastp_report_dir 是否存在
    fastp_report_exists = os.path.exists(fastp_report_dir) and os.path.isdir(fastp_report_dir)
    
    # 如果 fastp_report_dir 不存在，需要运行 fastp，此时 --run 参数必需
    if not fastp_report_exists:
        if not args.run:
            print('错误: fastp_report 目录不存在，必须指定 --run 参数 (beirui 或 nuohuo)')
            return
    
    # 读取样本文件
    with open(args.samples, 'r') as f:
        samples_data = f.readlines()
    
    if len(samples_data) < 2:
        print('错误: 样本文件格式不正确，至少需要表头行和一行数据')
        return
    
    # 解析表头
    header = samples_data[0].strip().split('\t')
    if 'group' not in header or 'sample' not in header:
        print('错误: 样本文件必须包含 group, sample 列')
        return
    
    group_idx = header.index('group')
    sample_idx = header.index('sample')
    
    # 存储处理后的样本信息
    processed_samples = []
    sample_name_list = []
    sample_data_list = []
    
    # 如果 fastp_report_dir 已存在，直接从已有的 json 文件读取
    if fastp_report_exists:
        print(f'检测到 fastp_report 目录已存在，跳过 fastp 运行，直接读取已有结果')
        
        # 处理每个样本
        for line in samples_data[1:]:
            if not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < max(group_idx, sample_idx) + 1:
                continue
            
            group = fields[group_idx].strip()
            sample = fields[sample_idx].strip()
            
            if not sample:
                continue
            
            # 检查 json 文件是否存在
            json_file = os.path.join(fastp_report_dir, f'{sample}_fastp.json')
            if os.path.exists(json_file):
                sample_name_list.append(sample)
                sample_data_list.append(json_file)
            else:
                print(f'警告: {json_file} 不存在，跳过')
    else:
        # 需要运行 fastp，检查 R1 和 R2 列
        if 'R1' not in header or 'R2' not in header:
            print('错误: 运行 fastp 时，样本文件必须包含 R1, R2 列')
            return
        
        r1_idx = header.index('R1')
        r2_idx = header.index('R2')
        
        # 处理每个样本
        for line in samples_data[1:]:
            if not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < max(group_idx, sample_idx, r1_idx, r2_idx) + 1:
                continue
            
            group = fields[group_idx].strip()
            sample = fields[sample_idx].strip()
            r1_path = fields[r1_idx].strip()
            r2_path = fields[r2_idx].strip()
            
            if not sample or not r1_path or not r2_path:
                continue
            
            # 检查输入文件是否存在
            if not os.path.exists(r1_path):
                print(f'警告: {r1_path} 不存在，跳过')
                continue
            if not os.path.exists(r2_path):
                print(f'警告: {r2_path} 不存在，跳过')
                continue
            
            # 运行 fastp
            if args.run == 'beirui':
                success, new_r1_path, new_r2_path = run_fastp_beirui(r1_path, r2_path, cleandata_dir, sample, fastp_report_dir)
            elif args.run == 'nuohuo':
                success, new_r1_path, new_r2_path = run_fastp_nuohuo(r1_path, r2_path, cleandata_dir, sample, fastp_report_dir)
            else:
                print(f'错误: 未知的运行模式 {args.run}')
                return
            
            if not success or not new_r1_path or not new_r2_path:
                continue
            
            # 保存处理后的样本信息
            processed_samples.append({
                'group': group,
                'sample': sample,
                'R1': new_r1_path,
                'R2': new_r2_path
            })
            
            # 收集质量表信息
            json_file = os.path.join(fastp_report_dir, f'{sample}_fastp.json')
            if os.path.exists(json_file):
                sample_name_list.append(sample)
                sample_data_list.append(json_file)
        
        # 生成新的 samples_described.txt（仅当运行了 fastp 时）
        if processed_samples:
            new_samples_file = 'samples_described.txt'
            with open(new_samples_file, 'w') as f:
                f.write('group\tsample\tR1\tR2\n')
                for sample_info in processed_samples:
                    f.write(f"{sample_info['group']}\t{sample_info['sample']}\t{sample_info['R1']}\t{sample_info['R2']}\n")
            print(f'\n已生成新的样本文件: {new_samples_file}')
    
    # 生成质量表
    write_quality_table(args.output, sample_name_list, sample_data_list)


def write_quality_table(output_file, sample_name_list, sample_data_list):
    """根据 fastp json 列表写出测序质量表。"""
    with open(output_file, 'w') as f:
        f.write('Sample\tClean_reads\tClean_base\tQ20\tQ30\tGC\n')
        for name, json_data in zip(sample_name_list, sample_data_list):
            if not os.path.exists(json_data):
                print(f'警告: {json_data} 不存在')
                continue
            try:
                clean_reads, clean_base, q20, q30, gc = parse_json(json_data)
                f.write(f'{name}\t{clean_reads}\t{clean_base}\t{q20}\t{q30}\t{gc}\n')
            except Exception as e:
                print(f'警告: 解析 {json_data} 失败: {e}')
                continue

    print(f'\n已生成质量表: {output_file}')
    print('\nDone!')


if __name__ == '__main__':
    main()
