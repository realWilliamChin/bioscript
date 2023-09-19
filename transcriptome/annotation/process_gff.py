#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/04/16 14:13
# Author        : William GoGo
import os
import re
import argparse
from nr_gff import detect_gff_type


def parse_input():
    argparser = argparse.ArgumentParser(description='处理一些 gff3 文件，生成 chromesome geneid startpos endpos strand')
    argparser.add_argument('-g', '--gff', help='gff3 file，默认是当前文件夹下的 gff3 或者 gff 文件')
    argparser.add_argument('-p', '--prefix', required=True, help='输出文件的前缀')
    argparser.add_argument('-t', '--gfftype', choices=['embl', 'ncbi', 'auto_detect'], default='auto_detect',
                           help='gff 类型, embl or ncbi，默认自动检测，检测失败手动输入')
    args = argparser.parse_args()
    
    return args


def get_all_id(gene_basicinfo_file, key_name):
    all_id_filename = key_name + '_all_gene_id.txt'
    # TODO:以后建议更换成不使用系统命令的方式
    command = f'cut -f 1 {gene_basicinfo_file} | tail -n +2 | sort -u > {all_id_filename}'
    os.system(command)
    

def get_basic_info(gff_file, output_file, gff_type):
    # ncbi 的 chromosome
    ncbi_chromosome = None
    
    with open(output_file, 'w') as gff3_w:
        gff3_w.write('GeneID\tChr_Number\tStart\tEnd\tStrand\tGene_Def\n')
    with open(gff_file, 'r') as file:
        for line in file:
            line = line.strip()
            
            if 'chromosome=' in line:
                ncbi_chromosome = re.search('chromosome=(.*?);', line).group(1)
            if line.startswith('#'):
                continue
            if not 'ID=gene' in line:
                continue
            
            columns = line.split('\t')
            
            # 根据 embl 和 ncbi 分形式提取
            if gff_type == 'embl':
                chromosome = columns[0]
                gene_id = re.search('ID=gene:(.*?);', line).group(1)
                biotype = re.search('biotype=(.*?);', line).group(1)
            elif gff_type == 'ncbi':
                chromosome = ncbi_chromosome
                gene_id = re.search('GeneID:(.*?);', line).group(1).split(',')[0]
                biotype = line.split('gene_biotype=')[1].strip().split(';')[0]
            
            start_pos = int(columns[3])
            end_pos = int(columns[4])
            strand = columns[6]
            
            with open(output_file, 'a') as gff3_w:
                gff3_w.write(f'{gene_id}\t{chromosome}\t{start_pos}\t{end_pos}\t{strand}\t{biotype}\n')
    
    
def main():
    args = parse_input()

    if args.gff:
        gff_file = args.gff
    else:
        gff_file = [x for x in os.listdir() if x.endswith('.gff3') or x.endswith('.gff') or x.endswith('.gtf')][0]

    if args.gfftype == 'auto_detect':
        args.gfftype = detect_gff_type(gff_file)
        
    gene_basicinfo_name = args.prefix + '_gene_basicinfo.txt'
    
    # 生成 gene_basicinfo 文件
    get_basic_info(gff_file, gene_basicinfo_name, args.gfftype)
    # 生成 gene_id 的文件
    get_all_id(gene_basicinfo_name, args.prefix)

if __name__ == '__main__':
    main()
