#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/30 12:41
# Author        : William GoGo
import os
import argparse
import pandas as pd
from collections import Counter


def parse_input():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-n', '--nr_blast_file', help='nr_diamond.blast 文件')
    argparser.add_argument('-g', '--go_file', help='swiss idNo_def 文件')
    argparser.add_argument('-k', '--kegg_file', help='KEGG_clean 文件')
    argparser.add_argument('-c', '--cog_file', help='unigene.emapper.annotations 文件')
    argparser.add_argument('-i', '--cog_file_id', help='unigene.emapper.seed_orthologs 文件')
    argparser.add_argument('-s', '--swiss_file', help='swiss_gene_def 文件')
    argparser.add_argument('-o', '--output_path', help='输出文件夹')
    return argparser.parse_args()


def get_nrid_identity_evalue_spcount(nr_file):
    df = pd.read_csv(nr_file, sep='\t', header=None)
    df_uniq = df.drop_duplicates(subset=0)
    species = df_uniq[12]
    lst = list()
    for each in species:
        lst.append(each.split('[')[-1].split(']')[0])
    head_10 = Counter(lst).most_common()[0:10]
    df_species = pd.DataFrame(head_10)
    return df_uniq[0], df_uniq[2], df_uniq[10], df_species


def main():
    args = parse_input()
    nr_file = args.nr_blast_file if args.nr_blast_file else [x for x in os.listdir() if x.endswith('nr_diamond.blast')][0]
    go_file = args.go_file if args.go_file else [x for x in os.listdir() if x.endswith('idNo_def.txt')][0]
    kegg_file = args.kegg_file if args.kegg_file else [x for x in os.listdir() if x.endswith('KEGG_clean.txt')][0]
    cog_file = args.cog_file if args.cog_file else [x for x in os.listdir() if x.endswith('emapper.annotations')][0]
    cog_file_id = args.cog_file_id if args.cog_file_id else [x for x in os.listdir() if x.endswith('emapper.seed_orthologs')][0]
    swiss_file = args.swiss_file if args.swiss_file else [x for x in os.listdir() if x.endswith('swiss_gene_def.txt')][0]

    identity_evalue = get_nrid_identity_evalue_spcount(nr_file)
    # 输出文件
    identity_evalue[0].to_csv('NR_ID.list', header=False, index=False)
    identity_evalue[1].to_csv('identity.txt', header=False, index=False)
    identity_evalue[2].to_csv('evalue.txt', header=False, index=False)
    identity_evalue[3].to_csv('species_count.txt', sep='\t', header=False, index=False)
    pd.read_csv(go_file, sep='\t', header=None, error_bad_lines=False)[0].drop_duplicates().to_csv('GO_ID.list', header=False, index=False)
    pd.read_csv(kegg_file, sep='\t', header=None)[0].drop_duplicates().to_csv('KEGG_ID.list', header=False, index=False)
    pd.read_csv(cog_file_id, sep='\t', skiprows=5, header=None)[0].to_csv('COG_ID.list', header=False, index=False)
    pd.read_csv(swiss_file, sep='\t', header=None)[0].to_csv('Swiss_ID.list', header=False, index=False)

    # 定义 COG count
    cog_source_df = pd.DataFrame({
        0: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'],
        1: ['FCDCFC', 'FCDCCC', 'BCFCFC', 'FCFCDC', 'DCFCFC', 'DCECFC', 'CCFCFC', 'DCDCFC', 'DCCCFC', 'FCCCFC', 'FCDCEC', 'FCDCDC', 'ECFCAC', 'DCFCAC', '9CFCAC', 'CCCCFC', 'BCCCFC', 'E0E0E0', 'CCCCCC', 'FCFCAC', 'ACFCAC', 'FCFCBC', 'BCFCAC', '9CFC9C', 'FCFCCC', 'CCFCAC'],
        2: ['A:RNA processing and modification',
            'B:Chromatin structure and dynamics',
            'C:Energy production and conversion',
            'D:Cell cycle control, cell division, chromosome partitioning',
            'E:Amino acid transport and metabolism',
            'F:Nucleotide transport and metabolism',
            'G:Carbohydrate transport and metabolism',
            'H:Coenzyme transport and metabolism',
            'I:Lipid transport and metabolism',
            'J:Translation, ribosomal structure and biogenesis',
            'K:Transcription',
            'L:Replication, recombination and repair',
            'M:Cell wall/membrane/envelope biogenesis',
            'N:Cell motility',
            'O:Posttranslational modification, protein turnover, chaperones',
            'P:Inorganic ion transport and metabolism',
            'Q:Secondary metabolites biosynthesis, transport and catabolism',
            'R:General function prediction only',
            'S:Function unknown',
            'T:Signal transduction mechanisms',
            'U:Intracellular trafficking, secretion, and vesicular transport',
            'V:Defense mechanisms',
            'W:Extracellular structures',
            'X:Mobilome: prophages, transposons',
            'Y:Nuclear structure',
            'Z:Cytoskeleton'
        ]
    })

    # 输出 cog_count.txt
    cog_count_source_df = pd.read_csv(cog_file, sep='\t', skiprows=3)
    cog_list = Counter(''.join(list(cog_count_source_df['COG Functional cat.'].dropna()))).most_common()
    cog_count_df = pd.DataFrame(list(cog_list))
    cog_df = pd.merge(left=cog_source_df, right=cog_count_df, on=0, how='left').fillna(0)
    cog_df.to_csv('COG_count.txt', sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
