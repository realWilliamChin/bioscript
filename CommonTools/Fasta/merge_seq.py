#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2026/03/13 17:42
# Author        : William GoGo
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Merge short scaffold sequences in genomic fasta file')
    parser.add_argument('-i', '--input', required=True, help='Input genomic fasta file')
    parser.add_argument('-o', '--output', required=True, help='Output fasta file')
    parser.add_argument('-l', '--length-threshold', type=int, required=True, 
                        help='Length threshold, sequences shorter than this will be merged')
    parser.add_argument('-g', '--gap-length', type=int, default=100, 
                        help='Number of Ns to insert between merged sequences, default: 100')
    parser.add_argument('-n', '--merged-name', default='Scaffold', 
                        help='Name of the merged sequence, default: Scaffold')
    
    args = parser.parse_args()
    
    # 参数校验
    if not os.path.exists(args.input):
        logger.error(f"Input file {args.input} does not exist")
        sys.exit(1)
    
    if args.length_threshold <= 0:
        logger.error("Length threshold must be a positive integer")
        sys.exit(1)
    
    if args.gap_length < 0:
        logger.error("Gap length cannot be negative")
        sys.exit(1)
    
    return args

def main():
    args = parse_args()
    
    logger.info(f"Starting sequence merging process")
    logger.info(f"Input file: {args.input}")
    logger.info(f"Output file: {args.output}")
    logger.info(f"Length threshold: {args.length_threshold} bp")
    logger.info(f"Gap length: {args.gap_length} Ns")
    logger.info(f"Merged sequence name: {args.merged_name}")
    
    # 读取所有序列
    records = list(SeqIO.parse(args.input, 'fasta'))
    total_sequences = len(records)
    total_length = sum(len(r) for r in records)
    
    logger.info(f"Total sequences in input: {total_sequences}")
    logger.info(f"Total base pairs in input: {total_length:,} bp")
    
    # 分离长序列和短序列
    long_sequences = []
    short_sequences = []
    
    for record in records:
        if len(record.seq) >= args.length_threshold:
            long_sequences.append(record)
        else:
            short_sequences.append(record)
    
    num_long = len(long_sequences)
    num_short = len(short_sequences)
    
    logger.info(f"Long sequences (>= {args.length_threshold} bp): {num_long}")
    logger.info(f"Short sequences (< {args.length_threshold} bp): {num_short}")
    
    if num_short == 0:
        logger.warning("No short sequences found to merge, writing all sequences to output directly")
        SeqIO.write(records, args.output, 'fasta')
        logger.info(f"Output written to {args.output}")
        sys.exit(0)
    
    # 合并短序列
    logger.info(f"Merging {num_short} short sequences...")
    
    merged_parts = []
    for i, seq in enumerate(short_sequences):
        merged_parts.append(str(seq.seq))
        if i < len(short_sequences) - 1:
            merged_parts.append('N' * args.gap_length)
    
    merged_sequence = ''.join(merged_parts)
    merged_length = len(merged_sequence)
    
    # 创建合并后的序列记录
    merged_record = SeqRecord(
        Seq(merged_sequence),
        id=args.merged_name,
        # description=f"merged_{num_short}_sequences_gap_{args.gap_length}N_total_length_{merged_length}"
    )
    
    # 准备输出序列：先长序列，后合并序列
    output_sequences = long_sequences + [merged_record]
    
    # 写入输出文件
    SeqIO.write(output_sequences, args.output, 'fasta')
    
    # 统计信息
    output_total_sequences = len(output_sequences)
    output_total_length = sum(len(r) for r in output_sequences)
    
    logger.info("=" * 60)
    logger.info("Merge completed successfully!")
    logger.info(f"Output sequences: {output_total_sequences} (reduced by {total_sequences - output_total_sequences})")
    logger.info(f"Output total length: {output_total_length:,} bp")
    logger.info(f"Merged sequence contains {num_short} short sequences")
    logger.info(f"Merged sequence length: {merged_length:,} bp")
    logger.info(f"Result saved to: {args.output}")
    logger.info("=" * 60)

if __name__ == "__main__":
    main()