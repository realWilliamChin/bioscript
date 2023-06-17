# -*- coding: UTF-8 -*-
# Created Time  : 2023/4/26 10:50
# Author        : WilliamGoGo
import os
import pandas as pd


def get_fasta_seq_from_gene_list(fasta_file, gene_list):
    """
    从 fasta 文件中提取指定基因的序列
    """
    gene_seq_dict = {}
    gene_name = ''
    with open(fasta_file, 'r') as f:
        for line in f.readlines():
            if '>' in line:
                gene_name = line.split(' ')[0].split('>')[1].strip()
                gene_seq_dict[gene_name] = ''
            else:
                gene_seq_dict[gene_name] = gene_seq_dict[gene_name] + line.strip()
    gene_list = list(set(gene_list))
    gene_seq_dict = {k: v for k, v in gene_seq_dict.items() if k in gene_list}
    return gene_seq_dict


def parse_input():
    """
    解析输入参数
    """
    import argparse
    parser = argparse.ArgumentParser(description='从 fasta 文件中提取指定基因的序列')
    parser.add_argument('-f', '--fasta', type=str, help='fasta 文件, 如果不指定, 则默认是当前目录的 .fasta 结尾的文件')
    parser.add_argument('-g', '--genelistfile', type=str, required=True, help='gene_list file 文件，必须指定')
    parser.add_argument('-n', '--colnumber', type=str, required=True, help='gene list 表的 gene list 列，从 0 开始数')
    args = parser.parse_args()
    if not args.fasta:
        args.fasta = [x for x in os.listdir() if x.endswith('.fasta')][0]
    return args.fasta, args.genelistfile, args.colnumber


def get_gene_list_from_df(gene_list_file, col_number):
    """
    从 gene_list_file 文件中提取基因名
    """
    df = pd.read_csv(gene_list_file, sep='\t', usecols=[col_number])
    gene_list = df.iloc[:, 0].tolist()
    gene_list = list(set(gene_list))
    return gene_list


if __name__ == '__main__':
    parse_input = parse_input()
    fasta_file_name = parse_input[0]
    gene_list_file = parse_input[1]
    col_number = int(parse_input[2])
    
    gene_list = get_gene_list_from_df(gene_list_file, col_number)
    gene_seq_dict = get_fasta_seq_from_gene_list(fasta_file_name, gene_list)
    with open(gene_list_file+'.fasta', 'w') as f:
        for gene_name, gene_seq in gene_seq_dict.items():
            f.write('>' + gene_name + '\n')
            f.write(gene_seq + '\n')
