#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/14 17:37
# Author        : William GoGo
import argparse
from Bio import SeqIO
from collections import defaultdict


def get_fasta_gene_number(input_file):
    gene_count = 0
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                gene_count += 1
    return gene_count


def rename_fasta(input_file, output_file, prefix='G_', start_num=1):
    """
    重命名 FASTA 文件中的基因名称，优化版。

    Args:
        input_file (str): 输入 FASTA 文件路径。
        output_file (str): 输出 FASTA 文件路径。
        prefix (str): 输入前缀
        start_num (int, optional): 新基因名称的起始编号。默认值为 1。
    """

    id_to_index = defaultdict(list)
    for i, record in enumerate(SeqIO.parse(input_file, "fasta")):
        id_to_index[record.id].append(i)

    # 处理基因 ID 重复
    for gene_id, indices in id_to_index.items():
        if len(indices) > 1:
            print(f"Warning: Gene ID {gene_id} is duplicated.")

    all_gene_count = get_fasta_gene_number(input_file)
    num_digits = len(str(all_gene_count))
    with open(output_file, 'w') as output_handle:
        for i, record in enumerate(SeqIO.parse(input_file, "fasta")):
            gene_num = start_num + i
            num_zeros = max(0, num_digits - len(str(gene_num)))
            new_id = f"{prefix}{'0'*num_zeros}{gene_num}"
            record.id = new_id
            record.description = ""
            SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rename genes in a FASTA file')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    parser.add_argument('-p', '--prefix', default='G_', help='Gene 名称前缀，默认是 G_')
    parser.add_argument('-s', '--start', type=int, default=1, help='Starting number for new gene names')
    args = parser.parse_args()

    rename_fasta(args.input, args.output, args.prefix, args.start)