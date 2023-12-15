'''
Author: WilliamGoGo realwilliamchin@outlook.com
Date: 2023-10-19 14:02:47
LastEditors: WilliamGoGo realwilliamchin@outlook.com
LastEditTime: 2023-10-19 14:02:47
FilePath: /script/fasta/drop_duplicate.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/14 17:37
# Author        : William GoGo
import argparse
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description='FASTA file formatter')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    parser.add_argument('-p', '--prefix', default='.', help='')
    return parser.parse_args()


def filter_duplicates(fasta_file, prefix):
    sequences = {}

    # 读取 FASTA 文件，并将序列保存到字典中
    for record in SeqIO.parse(fasta_file, "fasta"):
        # sequence_id = record.id
        sequence_id = record.id.rsplit(prefix, 1)[0]
        sequence = str(record.seq)
        length = len(sequence)
        
        
        # 如果序列已存在且长度更长，则更新字典中的值
        if sequence_id in sequences and length > len(sequences[sequence_id][0]):
            sequences[sequence_id] = (sequence, length)
        # 如果序列不存在，则将其添加到字典中
        elif sequence_id not in sequences:
            sequences[sequence_id] = (sequence, length)

    # 返回保留最长序列的列表，包括原始的序列名和序列
    longest_sequences = [(sequence_id, sequence) for sequence_id, (sequence, _) in sequences.items()]

    return longest_sequences


def write_sequences(output_file, sequences):
    with open(output_file, "w") as file:
        for sequence_id, sequence in sequences:
            file.write(f">{sequence_id}\n{sequence}\n")

def main():
    args = parse_arguments()
    longest_sequences = filter_duplicates(args.input, args.prefix)
    write_sequences(args.output, longest_sequences)

if __name__ == '__main__':
    main()

