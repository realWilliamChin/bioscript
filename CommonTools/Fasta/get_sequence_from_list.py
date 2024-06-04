# -*- coding: utf-8 -*-
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger


def parse_input():
    argparser = argparse.ArgumentParser(description='从 id.list 提取出 fasta 序列')
    argparser.add_argument('-i', '--idlist', required=True, help='id.list file, 列名 GeneID Start End Strand TargetGeneID')
    argparser.add_argument('-f', '--fasta', required=True, help='fasta 文件')
    argparser.add_argument('-o', '--output', required=True, help='输出文件名')
    argparser.add_argument('-t', '--type', default='on', help='输入类型(on/off),默认为on，则包含那些id提出来，反之则提出来不包含那些id的序列')
    # argparser.add_argument('--rename', action='store_true', help='输入文件idlist为两列，一列sourceid，另一列renameid')
    args = argparser.parse_args()
    return args


def get_seq_from_idlist(idlist, fasta, save_type, output):
    # with open(idlist, 'r') as gene_file:
    #     gene_ids = gene_file.read().splitlines()
    id_df = pd.read_csv(idlist, sep='\t')
    id_list = id_df['GeneID'].values.tolist()
    has_start_end = 'Start' in id_df.columns and 'End' in id_df.columns
    has_strand = 'Strand' in id_df.columns
    # 读取原始fasta文件并将序列存储到字典中
    fasta_sequences = {record.id: record for record in SeqIO.parse(fasta, 'fasta')}

    sequences = []
    for idx, row in id_df.iterrows():
        if has_start_end:
            start = int(row['Start'])
            end = int(row['End'])
            if end < start:
                start, end = end, start
            start = start - 1
            if row['GeneID'] in fasta_sequences:
                full_sequence = fasta_sequences[row['GeneID']]
                sub_sequence = full_sequence.seq[start:end]
                if has_strand and row['Strand'] == '-':
                    sub_sequence = sub_sequence.reverse_complement()
                new_record = SeqRecord(sub_sequence, id = f"{row['GeneID']}_{row['QueryName']}_{start+1}_{end}",
                                       description=f"{row['Strand']}")
                sequences.append(new_record)
        else:
            if save_type == 'on':
                if row['GeneID'] in fasta_sequences:
                    sequences.append(fasta_sequences[row['GeneID']])
            elif save_type == 'off':
                if row['GeneID'] not in fasta_sequences:
                    sequences.append(fasta_sequences[row['GeneID']])
    
    SeqIO.write(sequences, output, 'fasta')


def main():
    args = parse_input()
    get_seq_from_idlist(args.idlist, args.fasta, args.type, args.output)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()

