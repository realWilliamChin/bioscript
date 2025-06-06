#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/04/16 14:13
# Author        : William GoGo
import os, sys
import re
import argparse
from loguru import logger
import pandas as pd

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    argparser = argparse.ArgumentParser(description='处理一些 gff3 文件，生成 chromesome geneid startpos endpos strand')
    argparser.add_argument('-g', '--gff', help='gff3 file，默认是当前文件夹下的 gff3 或者 gff 文件')
    argparser.add_argument('--re', default='gene-(.*?);', help='指定正则表达式，用来提取 geneid，默认为 gene-(.*?);')
    argparser.add_argument('-o', '--output', default=os.getcwd(), help='输出文件夹，也可以+前缀，默认当前目录')
    args = argparser.parse_args()
    
    if not args.gff:
        args.gff = [x for x in os.listdir() if x.endswith('.gff3') or x.endswith('.gff')][0]
    if not os.path.isdir(args.output) and not args.output.endswith('_'):
        args.output = args.output + '_'
    
    return args


def get_basic_info(gff_file, output, re_pattern):
    # 创建字典存储数据
    gene_info_dict = {}
    
    duplicates_geneid_lines = []
    
    with open(gff_file, 'r') as file:
        for line in file:
            line = line.strip()
            columns = line.split('\t')
            # 过滤无效行
            if not line or len(line.split('\t')) < 3:
                continue
            # 过滤注释行
            if line.startswith('#'):
                continue
            
            # 找到chromosome=，提取chromosome
            if 'chromosome=' in line:
                chromosome = re.search('chromosome=(.*?);', line).group(1)
            else:
                chromosome = columns[0]
            
            # 过滤非基因行
            if 'gene' not in line.split('\t')[2].lower():
                continue
            if 'parent' in line.lower():
                continue
            
            # 提取biotype
            try:
                biotype = re.search('biotype=(.*?)(?:;|$)', line).group(1)
            except AttributeError:
                biotype = 'N/A'
            
            # 提取geneid
            try:
                gene_id = re.search(re_pattern, line).group(1).split(',')[0]
            except Exception:
                logger.warning('GeneID not found, save to not_found.txt')
                open('tmp_geneid_not_found.txt', 'a').write(line + '\n')
                continue
            
            start_pos = int(columns[3])
            end_pos = int(columns[4])
            strand = columns[6]
            
            # 检查是否有重复的gene_id
            if gene_id in gene_info_dict:
                if line not in duplicates_geneid_lines:
                    duplicates_geneid_lines.append(line)
                duplicates_geneid_lines.append(line)

            # 将数据存入字典
            gene_info_dict[gene_id] = {
                'chromosome': chromosome,
                'start': start_pos,
                'end': end_pos,
                'strand': strand,
                'biotype': biotype
            }

    df = pd.DataFrame.from_dict(gene_info_dict, orient='index')
    df.index.name = 'GeneID'
    df.reset_index(inplace=True)
    df.columns = ['GeneID', 'Chromosome', 'Start', 'End', 'Strand', 'Gene_type']
    write_output_df(df, output+'gene_basicinfo.txt', index=False)
    
    write_output_df(df[['GeneID']], output+'all_gene_id.txt', index=False)
    
    # 将重复的geneid行写入文件
    if duplicates_geneid_lines:
        with open('duplicates_geneid_lines.txt', 'w') as dup_file:
            for line in duplicates_geneid_lines:
                dup_file.write(line + '\n')


def main():
    args = parse_input()
    get_basic_info(args.gff, args.output, args.re)
    logger.success('Done!')


if __name__ == '__main__':
    main()
