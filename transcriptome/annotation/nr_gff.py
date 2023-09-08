#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/3/29 15:15
# Author        : WilliamGoGo
import re
import pandas as pd
import argparse
import sys


def parse_input():
    argparser = argparse.ArgumentParser(description='给 nr_gene_def 下面添加一些东西')
    argparser.add_argument('-n', '--nr', required=True, help='nr_gene_def file')
    argparser.add_argument('-g', '--gff', required=True, help='gff file')
    argparser.add_argument('-t', '--gfftype', choices=['embl', 'ncbi', 'auto_detect'], default='auto_detect',
                           help='gff 类型, embl or ncbi，默认自动检测，检测失败手动输入')
    # argparser.add_argument('-b', '--gffbasicinfo', required=True, help='gff basic info file')
    args = argparser.parse_args()
    return args


def detect_gff_type(gff_file):
    with open(gff_file, 'r') as f:
        gff_data = f.read()
        embl_type_count = gff_data.count("ID=gene:")
        ncbi_type_count = gff_data.count("ID=gene-")
        if embl_type_count > 0 and ncbi_type_count == 0:
            return "embl"
        elif embl_type_count == 0 and ncbi_type_count > 0:
            return "ncbi"
        else:
            print(embl_type_count, ncbi_type_count)
            print("nr_gff.py 注释检测 gff 类型失败！")
            sys.exit(1)


def main():
    args = parse_input()
    if args.gfftype == 'auto_detect':
        args.gfftype = detect_gff_type(args.gff)

    with open(args.gff, 'r') as f:
        lines = []
        for line in f:
            # if b'ensembl' in line and b'ID=gene:' in line:
            if 'ID=gene:' in line and args.gfftype == 'embl':
                geneid = re.search('ID=gene:(.*?);', line).group(1).strip()
                ensembl_type = line.split('\t')[2].strip()
                biotype = re.search('biotype=(.*?);', line).group(1).strip()
                lines.append([geneid, ensembl_type, biotype])
            elif 'ID=gene-' in line and args.gfftype == 'ncbi':
                geneid = re.search('GeneID:(.*?);', line).group(1).split(',')[0].strip()
                ensembl_type = line.split('\t')[2].strip()
                biotype = line.split('gene_biotype=')[1].strip().split(';')[0].strip()
                lines.append([geneid, ensembl_type, biotype])

    # gene_id.txt   Rat_Rattus_novergicus.mRatBN7.2_embl_gene_id_biotype.txt
    gene_id_df = pd.DataFrame(lines, columns=['GeneID', 'protein_coding', 'biotype'])
    # gene_id_df.to_csv(gff_file.replace('.gff3', 'gff_gene.txt'), sep='\t', index=False)
    gene_id_df = gene_id_df.sort_values(by='GeneID')
    # gene_id_df.drop(columns=['protein_coding', 'biotype']).to_csv(gff_file.replace('.gff3', '_gene_id.txt'), sep='\t', index=False)

    # gene_id.txt   Rat_Rattus_novergicus.mRatBN7.2_embl_gene_id_biotype.txt
    gene_protein_coding_df = gene_id_df[gene_id_df['biotype'].str.contains('protein_coding')].reset_index(drop=True)
    # gene_protein_coding_df['GeneID'].to_csv(gff_file.replace('.gff3', '_gene_protein_coding_id.txt'), sep='\t', index=False)

    # 加到 nr 后面的 'GeneID', 'NCBI_ID', 'NR_Def'
    nr_gene_def_df = pd.read_csv(args.nr, sep='\t', dtype={'GeneID': 'object', 'NCBI_ID': 'object'})

    # 没有注释上的
    non_annotationed_df = gene_protein_coding_df[~gene_protein_coding_df['GeneID'].isin(nr_gene_def_df['GeneID'])]
    non_annotationed_df = non_annotationed_df.drop(columns='protein_coding').rename(columns={'biotype': 'NR_Def'})

    # 不是 protein_coding 的
    gene_non_protein_coding_df = gene_id_df[~gene_id_df['biotype'].str.contains('protein_coding')]
    gene_non_protein_coding_df = gene_non_protein_coding_df.drop(columns='protein_coding').rename(columns={'biotype': 'NR_Def'})

    result_df = pd.concat([nr_gene_def_df, non_annotationed_df, gene_non_protein_coding_df], axis=0)
    result_df.drop_duplicates(subset='GeneID', keep='first', inplace=True)
    result_df.to_csv(args.nr, sep='\t', index=False)
    
    
if __name__ == '__main__':
    main()
