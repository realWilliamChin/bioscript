#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/20 21:45
# Author        : William GoGo
import os, sys
import pandas as pd
import argparse
import subprocess
from loguru import logger


def parse_input():
    args = argparse.ArgumentParser(description='hisat2_mapping_summary.py')
    args.add_argument(dest='alignment_type', choices=['hisat2', 'bowtie2', 'bwa'], help='使用哪种方式比对')
    args.add_argument('-s', '--samples', help='samples_described.txt', default='samples_described.txt')
    args.add_argument('--gff', help='[有参流程] stringtie 所需的 gff 文件')
    args.add_argument('--md', help='hisat2 mapping dir')
    args.add_argument('--cd', help='cleandata file dir', default='./02_Cleandata')
    
    args.add_argument('--bd', help='bam output dir')
    args.add_argument('--ref', help='reference specie')
    args.add_argument('--cpu', help='运行线程数量', default=30)
    args.add_argument('--msf', help='output mapping summary file', default='mapping_summary.txt')
    
    return args.parse_args()


def check_data(samples_file, cleandata_dir):
    """检查 cleandata 文件是否存在

    Args:
        samples_file (str): samples_described.txt
        cleandata_dir (str): cleandata 文件路径
    """
    samples_df = pd.read_csv(samples_file, sep='\t', usecols=['sample', 'filename1', 'filename2'])
    for index, row in samples_df.iterrows():
        samples_name = row['sample']
        file1 = os.path.join(cleandata_dir, row['filename1'])
        file2 = os.path.join(cleandata_dir, row['filename2'])
        if not os.access(file1, os.R_OK) or not os.access(file2, os.R_OK):
            logger.error(f'文件 {file1} 或 {file2} 不存在，将跳过此样本 {samples_name} 的比对')
            continue
    logger.info('cleandata 文件检查完毕')
    
    # if not os.access(f'{ref_data}.1.ht2', os.R_OK):
    #    logger.error(f'文件 {ref_data} 不存在，将退出程序')
    #    sys.exit(1)


def alignment(alignment_type, samples_file, cleandata_dir, mapping_file_dir, bam_file_dir, ref_specie, num_threads, gff_file=None):
    """# hisat2 and bowtie2 alignment 比对

    Args:
        samples_file (str): samples_described.txt
        mapping_file_dir (str): 输出 mapping 文件路径
        bam_file_dir (str): 输出 bam 文件路径
        ref_specie (str): 参考物种的路径路径
        num_threads (str): 运行线程数量
    """
    os.makedirs(mapping_file_dir, exist_ok=True)
    os.makedirs(bam_file_dir, exist_ok=True)
    samples_df = pd.read_csv(samples_file, sep='\t', usecols=['sample', 'filename1', 'filename2'])

    for index, row in samples_df.iterrows():
        
        sample_name = row['sample']
        file1 = os.path.join(cleandata_dir, row['filename1'])
        file2 = os.path.join(cleandata_dir, row['filename2'])
        mapping_file_name = os.path.join(mapping_file_dir, f'{sample_name}_mapping.txt')
        bam_file_name = os.path.join(bam_file_dir, f'{sample_name}.bam')
        
        logger.info(f'开始比对样本 {sample_name}, {file1} 和 {file2}，输出文件 {mapping_file_name} 和 {bam_file_name}')
        
        if alignment_type == 'hisat2':
            alignment_cmd = f'hisat2 -x {ref_specie} \
                -p {num_threads} \
                -I 200 -X 400 --fr \
                --min-intronlen 20 --max-intronlen 4000 \
                -1 {file1} \
                -2 {file2} \
                2> {mapping_file_name} | \
                samtools sort --threads {num_threads} -O BAM -o - > {bam_file_name}'
        elif alignment_type == 'bowtie2':
            alignment_cmd = f'bowtie2 -x {ref_specie} \
                -p {num_threads} \
                -1 {file1} \
                -2 {file2} \
                2> {mapping_file_name} | \
                samtools sort --threads {num_threads} -O BAM -o - > {bam_file_name}'
        elif alignment_type == 'bwa':
            alignment_cmd = f'bwa mem -M \
                -t {num_threads} \
                -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tLB:WES\\tPL:Illumina" \
                {ref_specie} \
                {file1} {file2} |\
                samtools sort -O bam -@ {num_threads} -o - > {bam_file_name}'
        
        rep = subprocess.run(alignment_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
        if rep.returncode != 0:
            logger.error(f'样本 {sample_name} 比对失败')
            logger.error(rep.stderr.decode())
            logger.error(rep.stdout.decode())
        else:
            with open(mapping_file_name) as f:
                last_line = f.readlines()[-1]
                logger.success(f'样本 {sample_name} 比对成功 {last_line}')
        
        if gff_file:
            stringtie_cmd = f'nohup stringtie -e -B \
                -G {gff_file} \
                -A {bam_file_dir}/fpkm/{sample_name}_fpkm.txt \
                -o {bam_file_dir}/ballgown/{sample_name}/{sample_name}.gtf \
                {bam_file_name} > /dev/null 2>&1 &'
            logger.debug(stringtie_cmd)
            os.system(stringtie_cmd)


def summary_hisat2_mapping(mapping_file_dir, samples_file, output_file):
    """ 对所有比对文件进行汇总，输出一个总的比对结果

    Args:
        mapping_file_dir (str): mapping summary 输入目录
        samples_file (str): samples_described.txt
        output_file (str): 输出文件
    """
    column_name = ['Sample', 'Total reads', 'Mapped reads', 'Unmapped reads', 'Unique mapped reads', 'Multiple mapped reads', 'overall alignment rate']
    open(output_file, 'w').write('\t'.join(column_name)+'\n')
    samples_data = pd.read_csv(samples_file, sep='\t', usecols=['sample'])
    for sample in samples_data['sample'].to_list():
        mapping_file = os.path.join(mapping_file_dir, f'{sample}_mapping.txt')
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
    if args.cd:
        check_data(args.samples, args.cd)
        alignment(args.alignment_type, args.samples, args.cd, args.md, args.bd, args.ref, args.cpu, args.gff)
    mapping_summary_file = os.path.join(args.md, args.msf)
    summary_hisat2_mapping(args.md, args.samples, mapping_summary_file)



