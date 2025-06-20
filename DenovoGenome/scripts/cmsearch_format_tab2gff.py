#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/06/28 10:21
# Author        : William GoGo
# Description   : Convert cmsearch output rfam_out.tab to gff format
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from Fasta.get_sequence_from_list import get_seq_from_idlist
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(description="Convert cmsearch output rfam_out.tab to gff")
    parser.add_argument('-i', '--input', type=str, dest='input_file', help="Input file")
    parser.add_argument('-p', '--prefix', type=str, help="Output file prefix")
    parser.add_argument('-s', type=str, dest='split_line', default='------', help="Split line symbol, default is '------'")
    parser.add_argument('-f', '--fasta', type=str, dest='fasta_file', help="Genomic Fasta file")
    
    return parser.parse_args()


def cmsearch_result_format(input_file, output_prefix, split_line='------'):
    """_summary_

    Args:
        input_file (_type_): _description_
        output_prefix (_type_): _description_
        split_line (str, optional): _description_. Defaults to '------'.
    """
    title_num_list = []
    with open(input_file, 'r') as f, open(output_prefix + '_result.txt', 'a') as w:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if split_line in line.strip().replace(' ', ''):
                line = line.strip()
                title_num_list = [j for j, char in enumerate(line) if char == " "]
                title_line = lines[i-1]
                title_list = []
                for j in range(len(title_num_list)):
                    if j == 0:
                        title_list.append(title_line[:title_num_list[j]].strip())
                    else:
                        title_list.append(title_line[title_num_list[j-1]:title_num_list[j]].strip())
                title_list.append(title_line[title_num_list[-1]:].strip())
                open(output_prefix + '_result.txt', 'w').write('\t'.join(title_list) + '\n')
            elif line.startswith('#'):
                continue
            else:
                # line_list = []
                # for j in range(len(title_num_list)):
                #     if j == 0:
                #         line_list.append(line[:title_num_list[j]].strip())
                #     else:
                #         line_list.append(line[title_num_list[j-1]:title_num_list[j]].strip())
                # line_list.append(line[title_num_list[-1]:].strip())
                # w.write('\t'.join(line_list) + '\n')
                line_list = line.split()
                w.write('\t'.join(line_list) + '\n')


def main():
    args = parse_input()
    input_file = args.input_file
    output_prefix = args.prefix
    cmsearch_result_format(input_file, output_prefix)
    
    # 转成 gff3
    df = load_table(output_prefix + '_result.txt')
    # 对列 query name 中进行统计，相同的 query name 用数字表示
    counts_df = df['query name'].value_counts().reset_index()
    write_output_df(counts_df, f'{output_prefix}_counts.txt', index=False)
    
    # 修正 gff
    mask = df['seq from'] > df['seq to']
    df.loc[mask, ['seq from', 'seq to']] = df.loc[mask, ['seq to', 'seq from']].values
    
    write_output_df(df, f'{output_prefix}_result.txt', index=False)
    
    
    # 关键词
    # 1. 5S_rRNA
    # 2. 5_8S_rRNA
    # 3. SSU_rRNA_eukarya --> 18S_rRNA
    # 4. LSU_rRNA_eukarya --> 28S_rRNA
    
    rrna_df = df[df['query name'].str.contains('5S_rRNA|5_8S_rRNA|SSU_rRNA_eukarya|LSU_rRNA_eukarya')].copy()
    rrna_df['query name'] = rrna_df['query name'].replace('SSU_rRNA_eukarya', '18S_rRNA')
    rrna_df['query name'] = rrna_df['query name'].replace('LSU_rRNA_eukarya', '28S_rRNA')
    rrna_df = rrna_df.astype({'E-value': str, 'mdl from': str, 'mdl to': str, 'gc': str})
    rrna_df['attribute'] = (
        'Rfam_family_id=' + rrna_df['query name'] +
        ';Rfam_accession=' + rrna_df['accession.1'] +
        ';E-value=' + rrna_df['E-value'] +
        ';trunc=' + rrna_df['trunc'] +
        ';model_start=' + rrna_df['mdl from'] +
        ';model_end=' + rrna_df['mdl to'] +
        ';GC-context=' + rrna_df['gc']
    )
    rrna_df['phase'] = '.'
    rrna_df['source'] = 'Rfam'
    rrna_df.rename(columns={'query name': 'type'}, inplace=True)
    rrna_df = rrna_df[['#target name', 'source', 'type', 'seq from', 'seq to', 'score', 'strand', 'phase', 'attribute']]
    rrna_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute']
    write_output_df(rrna_df, f'{output_prefix}_rRNA.gff3', index=False, header=False)
    
    # 统计 rRNA 类型数量
    rrna_type_counts = rrna_df['type'].value_counts().reset_index()
    rrna_type_counts.columns = ['rRNA_type', 'count']
    write_output_df(rrna_type_counts, f'{output_prefix}_rRNA_counts.txt', index=False)
    
    rrna_df = rrna_df[['seqid', 'attribute', 'start', 'end']]
    rrna_df['attribute'] = rrna_df['attribute'].str.split(';').str[0].str.split('=').str[1]
    rrna_df.columns = ['GeneID', 'TargetGeneID', 'Start', 'End']
    write_output_df(rrna_df, f'{output_prefix}_rRNA_GeneID.txt', index=False)
    # get_seq_from_idlist(output_prefix + '_rrna_geneid.txt', args.fasta_file, 'on', output_prefix + '_rrna.fasta')
    get_seq_from_idlist(rrna_df, args.fasta_file, 'on', f'{output_prefix}_rRNA.fasta')
    
    ncrna_df = df[~df['query name'].str.contains('rRNA')].copy()
    ncrna_df = ncrna_df[~ncrna_df['query name'].str.contains('tRNA')].copy()
    ncrna_df['type'] = 'ncRNA'
    ncrna_df = ncrna_df.astype({'E-value': str, 'mdl from': str, 'mdl to': str, 'gc': str})
    ncrna_df['attribute'] = (
        'Rfam_family_id=' + ncrna_df['query name'] +
        ';Rfam_accession=' + ncrna_df['accession.1'] +
        ';E-value=' + ncrna_df['E-value'] +
        ';trunc=' + ncrna_df['trunc'] +
        ';model_start=' + ncrna_df['mdl from'] +
        ';model_end=' + ncrna_df['mdl to'] +
        ';GC-context=' + ncrna_df['gc']
    )
    ncrna_df['phase'] = '.'
    ncrna_df['source'] = 'Rfam'
    ncrna_df = ncrna_df[['#target name', 'source', 'type', 'seq from', 'seq to', 'score', 'strand', 'phase', 'attribute']]
    ncrna_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute']
    write_output_df(ncrna_df, f'{output_prefix}_ncRNA.gff3', index=False, header=False)
        
    ncrna_df = ncrna_df[['seqid', 'attribute', 'start', 'end']]
    ncrna_df['attribute'] = ncrna_df['attribute'].str.split(';').str[0].str.split('=').str[1]
    ncrna_df.columns = ['GeneID', 'TargetGeneID', 'Start', 'End']
    write_output_df(ncrna_df, f'{output_prefix}_ncRNA_GeneID.txt', index=False)
    get_seq_from_idlist(ncrna_df, args.fasta_file, 'on', f'{output_prefix}_ncRNA.fasta')
    
    logger.success('Done')
    

if __name__ == '__main__':
    main()
