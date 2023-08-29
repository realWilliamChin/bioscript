#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/3/29 15:15
# Author        : WilliamGoGo
import re
import pandas as pd
import argparse


def parse_input():
    argparser = argparse.ArgumentParser(description='给 nr_gene_def 下面添加一些东西')
    argparser.add_argument('-n', '--nr', required=True, help='nr_gene_def file')
    argparser.add_argument('-g', '--gff', required=True, help='gff file')
    argparser.add_argument('-t', '--gfftype', required=True, help='gff type, embl or ncbi')
    # argparser.add_argument('-b', '--gffbasicinfo', required=True, help='gff basic info file')
    args = argparser.parse_args()
    return args

def main():
    args = parse_input()
        
    # nr_gene_def_df = pd.read_csv(args.nr, sep='\t', names=['GeneID', 'NCBI_ID', 'NR_Def'])
    # gff_basic_info_df = pd.read_csv(args.gffbasicinfo, sep='\t', skiprows=1, usecols=[0, 5], names=['GeneID', 'NR_Def'])
    
    # # convert NCBI_ID to object, gff_basic_info_df GeneID to object，不然没法合并
    # nr_gene_def_df['NCBI_ID'] = nr_gene_def_df['NCBI_ID'].astype('object')
    # gff_basic_info_df['GeneID'] = gff_basic_info_df['GeneID'].astype('object')
    
    # # 把 gff_basic_info_df 中 GeneID 在 NR 中的去掉，NR 没有注释上的
    # print(gff_basic_info_df.shape)
    # gff_basic_info_df = gff_basic_info_df[~gff_basic_info_df['GeneID'].isin(nr_gene_def_df['GeneID'])]

    # print(gff_basic_info_df.shape)
    # # 把不是 protein_coding 的去掉
    # # gff_basic_info_df = gff_basic_info_df[~gff_basic_info_df['NR_Def'].str.contains('protein_coding')]
    
    # print(gff_basic_info_df.shape)
    # # basic_info 中的 Gene_Def 加到 nr_gene_def 的 NR_Def 列后面
    # nr_gene_def_df = pd.merge(nr_gene_def_df, gff_basic_info_df, on='GeneID', how='left')
    # print(nr_gene_def_df.shape)
    # print(nr_gene_def_df.tail())
    # nr_gene_def_df.to_csv('test.txt', sep='\t', index=False)
    

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
