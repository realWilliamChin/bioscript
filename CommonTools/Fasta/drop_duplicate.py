#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/14 17:37
# Author        : William GoGo
import argparse
from Bio import SeqIO
from loguru import logger

def parse_arguments():
    parser = argparse.ArgumentParser(description='FASTA file formatter')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    parser.add_argument('-s', '--split', default='@@@', help='序列名以哪个字符串分割，默认不分割')
    parser.add_argument('--keep', choices=['long', 'first'], default='long',
                        help='去除重复的选项，选择保留最长的，或者是第一个')
    return parser.parse_args()


def filter_duplicates(fasta_file, strsplit, keep='long'):
    sequences = {}

    # 读取 FASTA 文件，并将序列保存到字典中
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_id = record.id.split(strsplit, 1)[0]
        sequence = str(record.seq)
        length = len(sequence)
        
        if sequence_id not in sequences:
            sequences[sequence_id] = (sequence, length)
        elif keep == 'long' and length > len(sequences[sequence_id][0]):
            sequences[sequence_id] = (sequence, length)
        # 如果 keep == 'first'，则保持第一个序列不变

    # 返回保留的序列列表，包括原始的序列名和序列
    filtered_sequences = [(sequence_id, sequence) for sequence_id, (sequence, _) in sequences.items()]

    return filtered_sequences


def write_sequences(output_file, sequences):
    with open(output_file, "w") as file:
        for sequence_id, sequence in sequences:
            file.write(f">{sequence_id}\n{sequence}\n")

def main():
    args = parse_arguments()
    filtered_sequences = filter_duplicates(args.input, args.split, args.keep)
    write_sequences(args.output, filtered_sequences)
    logger.success(f"已成功去除重复序列，并保存到 {args.output} 文件中")


if __name__ == '__main__':
    main()

