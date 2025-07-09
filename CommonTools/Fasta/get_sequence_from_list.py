# -*- coding: utf-8 -*-
import os, sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger
from collections import defaultdict


def parse_input():
    argparser = argparse.ArgumentParser(description='从 id.list 提取出 fasta 序列')
    argparser.add_argument('-i', '--idlist', required=True, help='id.list file, 列名 GeneID Start End Strand TargetGeneID')
    argparser.add_argument('-f', '--fasta', required=True, help='fasta 文件')
    argparser.add_argument('-o', '--output', required=True, help='输出文件名')
    argparser.add_argument('-t', '--type', default='on', help='输入类型(on/off),默认为on，则包含那些id提出来，反之则提出来不包含那些id的序列')
    args = argparser.parse_args()
    return args


def get_seq_from_idlist(id_df, fasta, save_type, output):
    # with open(idlist, 'r') as gene_file:
    #     gene_ids = gene_file.read().splitlines()
    has_start_end = 'Start' in id_df.columns and 'End' in id_df.columns
    has_strand = 'Strand' in id_df.columns
    has_targetid = 'TargetGeneID' in id_df.columns
    only_geneid = 'GeneID' in id_df.columns and 'TargetGeneID' not in id_df.columns and 'Start' not in id_df.columns and 'End' not in id_df.columns
    # 读取原始fasta文件并将序列存储到字典中（支持同名id全部输出）
    fasta_sequences = defaultdict(list)
    for record in SeqIO.parse(fasta, 'fasta'):
        fasta_sequences[record.id].append(record)

    sequences = []
    for idx, row in id_df.iterrows():
        if has_start_end:
            start = int(row['Start'])
            end = int(row['End'])
            if end < start:
                start, end = end, start
            start = start - 1
            if row['GeneID'] in fasta_sequences:
                for full_sequence in fasta_sequences[row['GeneID']]:
                    sub_sequence = full_sequence.seq[start:end]
                    if has_strand and row['Strand'] == '-':
                        sub_sequence = sub_sequence.reverse_complement()
                    if has_targetid:
                        new_record = SeqRecord(sub_sequence, id=f"{row['TargetGeneID']}", name='', description='')
                    else:
                        new_record = SeqRecord(sub_sequence, id=row['GeneID'], name='', description='')
                    sequences.append(new_record)

        elif only_geneid:
            if save_type == 'on':
                if row['GeneID'] in fasta_sequences:
                    for rec in fasta_sequences[row['GeneID']]:
                        sequences.append(rec)

        elif has_targetid:
            if row['GeneID'] in fasta_sequences:
                for rec in fasta_sequences[row['GeneID']]:
                    new_record = SeqRecord(rec.seq, id=row['TargetGeneID'], name='', description='')
                    sequences.append(new_record)
    if only_geneid and save_type == 'off':
        # 获取所有不在id_df中的GeneID
        all_gene_ids = set(fasta_sequences.keys())
        target_gene_ids = set(id_df['GeneID'])
        for gene_id in all_gene_ids - target_gene_ids:
            for rec in fasta_sequences[gene_id]:
                sequences.append(rec)

    SeqIO.write(sequences, output, 'fasta')


def main():
    args = parse_input()
    id_df = pd.read_csv(args.idlist, sep='\t')
    get_seq_from_idlist(id_df, args.fasta, args.type, args.output)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()
