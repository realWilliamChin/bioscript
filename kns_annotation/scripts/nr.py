#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/14 17:05
# Author        : William GoGo
"""
nr 注释程序
生成 nr.blast 文件，nr_gene_def.txt 文件 和 nr_TF_def.txt 文件
"""
import argparse
import os, sys
import subprocess
import pandas as pd
from loguru import logger
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/Fasta/'))
from get_sequence_from_list import get_seq_from_idlist


def parse_input():
    parser = argparse.ArgumentParser(description='输入 nr.blast 文件的路径')
    # parser.add_argument(dest='running_type', choices=['all', 'blast', 'parse_blast'], help='运行类型')
    
    parser.add_argument('-b', '--blast', type=str, help='nr.blast file')
    parser.add_argument('--basicinfo', type=str, dest='basicinfo', 
                        help='gff 类型, embl or ncbi，默认自动检测，检测失败手动输入，如果 basicinfo 文件 Gene_Def 都是 NA，怎不需要添加了')
    parser.add_argument('-p', '--prefix', help='输出文件的前缀')
    parser.add_argument('--not-annotationed', action='store_true', dest='not_annotationed',
                        help='输出没有注释上的基因信息（必须输入 --fasta 参数）')
    
    annotation = parser.add_argument_group('需要注释添加一下参数')
    annotation.add_argument('-f', '--fasta', help='输入需要去注释的 cds fasta 或 unigene fasta 文件（去重，名字精简）')
    annotation.add_argument('-t', '--threads', type=int, help='运行 nr 注释线程数量(好像不咋管用)')
    annotation.add_argument('-n', '--outfmt', help="nr 输出格式",
                            default="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,stitle")

    args = parser.parse_args()
        
    return args


def nr_annotation(fasta_file, blast_file, outfmt, num_threads):
    os.mkdir('./temp')
    anno_cmd = f'diamond blastx --db /home/data/ref_data/db/diamond_nr/diamond_nr \
        --threads ${num_threads} \
        --query {fasta_file} \
        --out {blast_file} \
        --outfmt 6 {outfmt} \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --id 30 \
        --block-size 20.0 \
        --tmpdir ./temp\
        --index-chunks 1'

    ret = subprocess.run(anno_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{fasta_file} nr 注释程序失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        logger.success(f"{fasta_file} nr 注释成功，输出到 {blast_file}")
    os.rmdir('./temp')
    return True


def nr(nr_blast_file, gene_basicinfo_file, columns, output_name):
    data_frame = pd.read_csv(nr_blast_file, sep='\t', names=columns, dtype=str)

    data_frame = data_frame.sort_values(by=['qseqid', 'bitscore'], ascending=[True, False])
    data_frame = data_frame.drop_duplicates(subset=['qseqid'], keep='first')
    data_frame.to_csv(output_name +  '_nr_uniq.blast', sep='\t', index=False)
    nr_gene_def_df = data_frame[['qseqid', 'sseqid', 'stitle']].copy()
    nr_gene_def_df.columns = ['GeneID', 'NCBI_ID', 'NR_Def']
    nr_gene_def_df['NR_Def'] = nr_gene_def_df['NR_Def'].str.split(n=1).str[1]
    print(f'注释上的基因数量是 {nr_gene_def_df.shape[0]} 个')
    
    # 有参基因中没有注释上的 gene 的 biotype 类型，如果是无参不需要此步骤
    if gene_basicinfo_file and os.path.exists(gene_basicinfo_file):
        nr_gene_def_df = nr_def_add_not_protein_coding(nr_gene_def_df, gene_basicinfo_file)
    else:
        print('没有检测到 gene_basicinfo 文件，或文件输入有误，如没有输入 -b 参数，忽略此条消息')
        
    nr_gene_def_df.to_csv(output_name + '_nr_gene_def.txt', sep='\t', index=False)
    
    nr_TF_def_df = nr_gene_def_df[nr_gene_def_df['NR_Def'].str.contains('transcription')]
    nr_TF_def_df = nr_TF_def_df.sort_values(by='NR_Def', key=lambda x: x.str.lower())
    nr_TF_def_df.to_csv(output_name + '_nr_TF_def.txt', sep='\t', index=False)


def nr_def_add_not_protein_coding(nr_gene_def_df, gene_basicinfo_file):
    gff_basicinfo_df = pd.read_csv(gene_basicinfo_file, sep='\t', usecols=['GeneID', 'Gene_Def'])
    print(f'总基因数量是 {gff_basicinfo_df.shape[0]} 个')
    
    # 没有注释上的
    gene_protein_coding_df = gff_basicinfo_df[gff_basicinfo_df['Gene_Def'].str.contains('protein_coding')].copy()
    non_annotationed_df = gene_protein_coding_df[~gene_protein_coding_df['GeneID'].isin(nr_gene_def_df['GeneID'])]
    non_annotationed_df = non_annotationed_df.rename(columns={'Gene_Def': 'NR_Def'})
    logger.info(f'protein coding 没有注释上的有 {non_annotationed_df.shape[0]} 个')
    
    # 不是 protein_coding 的
    gene_non_protein_coding_df = gff_basicinfo_df[~gff_basicinfo_df['Gene_Def'].str.contains('protein_coding')].copy()
    gene_non_protein_coding_df = gene_non_protein_coding_df.rename(columns={'Gene_Def': 'NR_Def'})
    logger.info(f'不是 protein_coding 没有注释上的有 {gene_non_protein_coding_df.shape[0]} 个')
    
    result = pd.concat([nr_gene_def_df, non_annotationed_df, gene_non_protein_coding_df], axis=0)
    result = result.drop_duplicates(subset='GeneID', keep='first')
    result = result.sort_values(by='GeneID')
    logger.info(f'加上没有注释上的和不是 protein_coding 的有 {result.shape[0]} 个')
    if gff_basicinfo_df.shape[0] == result.shape[0]:
        logger.success('\n注释结果正确！')
    
    return result


def get_not_annotationed_fasta(fasta, nr_uniq_blast_file, output_prefix):
    nr_uniq_bast_df = pd.read_csv(nr_uniq_blast_file, sep='\t', usecols=[0], skiprows=1, names=['GeneID'])
    get_seq_from_idlist(nr_uniq_bast_df, fasta, 'off', f"{output_prefix}_nr_not_annotationed.fasta")


def main():
    args = parse_input()
    
    if args.fasta and args.threads:
        outfmt = args.outfmt.replace(',', ' ')
        ret = nr_annotation(args.fasta, args.blast, outfmt, args.threads)
        if not ret:
            sys.exit(1)
    else:
        if not args.blast:
            args.blast = [x for x in os.listdir() if x.endswith('nr.blast')][0]
    columns = args.outfmt.split(',')
    nr(args.blast, args.basicinfo, columns, args.prefix)
    
    if args.not_annotationed:
        get_not_annotationed_fasta(args.fasta, args.blast, args.prefix)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()