# -*- coding: utf-8 -*-
import argparse
from Bio import SeqIO


def parse_input():
    argparser = argparse.ArgumentParser(description='从 id.list 提取出 fasta 序列')
    argparser.add_argument('-i', '--idlist', required=True, help='id.list file')
    argparser.add_argument('-f', '--fasta', required=True, help='fasta 文件')
    argparser.add_argument('-o', '--output', required=True, help='输出文件名')
    argparser.add_argument('-t', '--type', default='on', help='输入类型(on/off),默认为on，则包含那些id提出来，反之则提出来不包含那些id的序列')
    args = argparser.parse_args()
    return args
    
    
def get_seq_from_idlist(idlist, fasta, save_type, output):
    # 读取gene ID列表文件
    with open(idlist, 'r') as gene_file:
        gene_ids = gene_file.read().splitlines()

    # 读取原始fasta文件并提取序列
    sequences = []
    for record in SeqIO.parse(fasta, 'fasta'):
        if save_type == 'on':
            if record.id in gene_ids:
                sequences.append(record)
        elif save_type == 'off':
            if record.id not in gene_ids:
                sequences.append(record)

    # 生成新的fasta文件
    SeqIO.write(sequences, output, 'fasta')


def main():
    args = parse_input()
    get_seq_from_idlist(args.idlist, args.fasta, args.type, args.output)
    

if __name__ == '__main__':
    main()

