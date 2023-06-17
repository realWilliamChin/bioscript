#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/3/29 15:15
# Author        : WilliamGoGo
import os
import re
import pandas as pd
import argparse


def parse_input():
    argparser = argparse.ArgumentParser(description='给 nr_gene_def 下面添加一些东西')
    argparser.add_argument('-n', '--nr', required=True, help='nr_gene_def file')
    argparser.add_argument('-g', '--gff', help='gff3 file，默认是当前文件夹下的 gff3')
    args = argparser.parse_args()
    nr_file = args.nr
    gff_file = args.gff
    return nr_file, gff_file

def main():
    nr_file, gff_file = parse_input()
    if not gff_file:
        gff_file = [x for x in os.listdir() if 'gff3' in x][0]

    with open(gff_file, 'r') as f:
        lines = []
        for line in f:
            # if b'ensembl' in line and b'ID=gene:' in line:
            if 'ID=gene:' in line:
                geneid = re.search('ID=gene:(.*?);', line).group(1)
                ensembl_type = line.split('\t')[2]
                biotype = re.search('biotype=(.*?);', line).group(1)
                lines.append([geneid, ensembl_type, biotype])

    # gene_id.txt   Rat_Rattus_novergicus.mRatBN7.2_embl_gene_id_biotype.txt
    gene_id_df = pd.DataFrame(lines, columns=['GeneID', 'protein_coding', 'biotype'])
    # gene_id_df.to_csv(gff_file.replace('.gff3', 'gff_gene.txt'), sep='\t', index=False)
    gene_id_df = gene_id_df.sort_values(by='GeneID')
    gene_id_df.drop(columns=['protein_coding', 'biotype']).to_csv(gff_file.replace('.gff3', '_gene_id.txt'), sep='\t', index=False)

    # gene_id.txt   Rat_Rattus_novergicus.mRatBN7.2_embl_gene_id_biotype.txt
    gene_protein_coding_df = gene_id_df[gene_id_df['biotype'].str.contains('protein_coding')]
    gene_protein_coding_df['GeneID'].to_csv(gff_file.replace('.gff3', '_gene_protein_coding_id.txt'), sep='\t', index=False)

    # 加到 nr 后面的 'GeneID', 'NCBI_ID', 'NR_Def'
    nr_gene_def_df = pd.read_csv(nr_file, sep='\t')

    # 没有注释上的
    non_annotationed_df = gene_protein_coding_df[~gene_protein_coding_df['GeneID'].isin(nr_gene_def_df['GeneID'])]
    non_annotationed_df = non_annotationed_df.drop(columns='protein_coding').rename(columns={'biotype': 'NR_def'})

    # 不是 protein_coding 的
    gene_non_protein_coding_df = gene_id_df[~gene_id_df['biotype'].str.contains('protein_coding')]
    gene_non_protein_coding_df = gene_non_protein_coding_df.drop(columns='protein_coding').rename(columns={'biotype': 'NR_def'})

    result_df = pd.concat([nr_gene_def_df, non_annotationed_df, gene_non_protein_coding_df], axis=0)
    result_df = result_df.drop_duplicates(subset='GeneID', keep='first')
    result_df.to_csv(nr_file, sep='\t', index=False)
    
    
if __name__ == '__main__':
    main()
