import os, sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger


def parse_input():
    p = argparse.ArgumentParser(description='FASTA file formatter')
    p.add_argument(dest='fasta_files', nargs='+', help='输入 fasta 文件')
    p.add_argument('-o', '--output', required=True, help='Output FASTA file')
    return p.parse_args()


def write_sequences(output_file, sequences):
    with open(output_file, "w") as file:
        for sequence, sequence_id in sequences.items():
            file.write(f">{sequence_id}\n{sequence}\n")


def common_seq(fasta_file1, fasta_file2, common_file):
    fasta1_sequences = {}
    common_sequences = {}
    fasta1_uniq_seqs = {}
    fasta2_uniq_seqs = {}

    # 读取 FASTA 文件，并将序列保存到字典中
    for record in SeqIO.parse(fasta_file1, "fasta"):
        sequence_id = record.id
        sequence = str(record.seq)
        fasta1_sequences[sequence] = sequence_id
    
    for record in SeqIO.parse(fasta_file2, "fasta"):
        sequence_id = record.id
        sequence = str(record.seq).replace('\n', '').replace(' ', '')
        if sequence in list(fasta1_sequences.keys()):
            common_seq_id = f'{fasta1_sequences.get(sequence)}_{sequence_id}'
            common_sequences[sequence] = common_seq_id
        else:
            fasta2_uniq_seqs[sequence] = sequence_id
    
    print(len(common_sequences))
    
    for record in SeqIO.parse(fasta_file1, "fasta"):
        sequence_id = record.id
        sequence = str(record.seq).replace('\n', '').replace(' ', '')
        if sequence not in common_sequences.keys():
            fasta1_uniq_seqs[sequence] = sequence_id

    write_sequences(f'{fasta_file1}.uniq', fasta1_uniq_seqs)
    write_sequences(f'{fasta_file2}.uniq', fasta2_uniq_seqs)
    write_sequences(common_file, common_sequences)


def main():
    args = parse_input()
    fasta1, fasta2 = args.fasta_files
    common_seq(fasta1, fasta2, args.output)


if __name__ == '__main__':
    main()