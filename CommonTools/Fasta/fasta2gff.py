#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2025/05/13 17:37
# Author        : William GoGo
import os, sys
import argparse
from loguru import logger


def parse_fasta(fasta_file):
    """解析FASTA文件，生成序列ID和序列"""
    seq_id = None
    sequence = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    yield (seq_id, ''.join(sequence))
                seq_id = line[1:].split()[0]  # 取第一个空格前的部分作为ID
                sequence = []
            else:
                sequence.append(line)
        if seq_id is not None:
            yield (seq_id, ''.join(sequence))


def write_gff(fasta_file, gff_file):
    """将FASTA转换为GFF"""
    with open(gff_file, 'w') as gff:
        for seq_id, seq in parse_fasta(fasta_file):
            length = len(seq)
            # GFF格式：seqid, source, type, start, end, score, strand, phase, attributes
            gff_line = f"{seq_id}\t.\tgene\t1\t{length}\t.\t.\t.\tID={seq_id}\n"
            gff.write(gff_line)


def parse_input():
    """使用argparse解析命令行输入输出参数"""
    parser = argparse.ArgumentParser(description='将FASTA文件转换为GFF格式')
    parser.add_argument('-i', '--input', required=True, help='输入的FASTA文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出的GFF文件路径')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_input()
    write_gff(args.input, args.output)
    logger.success('Done!')