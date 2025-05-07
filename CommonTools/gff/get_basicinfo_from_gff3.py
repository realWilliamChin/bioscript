#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/04/16 14:13
# Author        : William GoGo
import os
import re
import argparse
from loguru import logger


def parse_input():
    argparser = argparse.ArgumentParser(description='处理一些 gff3 文件，生成 chromesome geneid startpos endpos strand')
    argparser.add_argument('-g', '--gff', help='gff3 file，默认是当前文件夹下的 gff3 或者 gff 文件')
    argparser.add_argument('-p', '--prefix', required=True, help='输出文件的前缀')
    argparser.add_argument('--re', default='gene-(.*?);', help='指定正则表达式，用来提取 geneid，默认为 gene-(.*?);')
    args = argparser.parse_args()
    
    return args


def get_all_id(gene_basicinfo_file, key_name):
    all_id_filename = key_name + '_all_gene_id.txt'
    # TODO:以后建议更换成不使用系统命令的方式
    command = f'cut -f 1 {gene_basicinfo_file} | tail -n +2 | sort -u > {all_id_filename}'
    os.system(command)
    

def get_basic_info(gff_file, output_file, re_pattern):
    # ncbi 的 chromosome
    ncbi_chromosome = None
    
    with open(output_file, 'w') as gff3_w:
        gff3_w.write('GeneID\tChromosome\tStart\tEnd\tStrand\tGene_type\n')
    with open(gff_file, 'r') as file:
        for line in file:
            line = line.strip()
            columns = line.split('\t')
            # 过滤空行
            if not line or len(line.split('\t')) < 3:
                continue
            if 'chromosome=' in line:
                chromosome = re.search('chromosome=(.*?);', line).group(1)
            else:
                chromosome = columns[0]
            if line.startswith('#'):
                continue
            if 'gene' not in line.split('\t')[2].lower():
                continue
            
            try:
                biotype = re.search('biotype=(.*?)(?:;|$)', line).group(1)
            except AttributeError:
                biotype = '---'
            try:
                gene_id = re.search(re_pattern, line).group(1).split(',')[0]
            except Exception:
                logger.warning('GeneID not found, save to not_found.txt')
                open('not_found.txt', 'a').write(line + '\n')
                
            
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

    gene_basicinfo_name = args.prefix + '_gene_basicinfo.txt'
    
    # 生成 gene_basicinfo 文件
    get_basic_info(gff_file, gene_basicinfo_name, args.re)
    # 生成 gene_id 的文件
    get_all_id(gene_basicinfo_name, args.prefix)
    
    logger.info('Done!')


if __name__ == '__main__':
    main()
