#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/01/15 16:50
# Author        : William GoGo
import os
import pandas as pd
import argparse


def parse_input():
    parser = argparse.ArgumentParser(description='Add definition to gene file')
    parser.add_argument('--ko-list', dest="ko_list", type=str, help='ko_list file')
    parser.add_argument('--kegg-genedef', dest="kegg_genedef", type=str, help='kegg_gene_def.txt')
    parser.add_argument('--kegg-pathway', dest="kegg_pathway", type=str, help='kegg.txt')
    parser.add_argument('--basicinfo', type=str, help='basicinfo.txt')
    parser.add_argument('-p', '--output_prefix', type=str, help='output file prefix')
    
    args = parser.parse_args()
    
    if os.path.exists(args.ko_list) is False:
        raise Exception('ko_list file not exists')
    if os.path.exists(args.kegg_genedef) is False:
        raise Exception('kegg_gene_def.txt not exists')
    if os.path.exists(args.kegg_pathway) is False:
        raise Exception('kegg.txt not exists')
    # if os.path.exists(args.basicinfo) is False:
    #     raise Exception('basicinfo.txt not exists')
    if args.output_prefix is None:
        args.output_prefix = 'output'
    
    return args


# def add_def(ko_file, kegg_gene_def_file, kegg_pathway_file, basicinfo, output_prefix):
def add_def(ko_file, kegg_gene_def_file, kegg_pathway_file, output_prefix):
    basicinfo_df = pd.read_csv(basicinfo, sep='\t', usecols=[0,1,2,3])
    basicinfo_df_columns = basicinfo_df.columns.tolist()
    ko_df = pd.read_csv(ko_file, sep='\t', names=['Ko_Number'])
    kegg_gene_def = pd.read_csv(kegg_gene_def_file, sep='\t')
    kegg_pathway_df = pd.read_csv(kegg_pathway_file, sep='\t', names=['GeneID', 'Ko'])
    kegg_pathway_df['Ko_Number'] = kegg_pathway_df['Ko'].str.split(":").str[0]
    
    ko_gene_df = pd.merge(ko_df, kegg_pathway_df, how='inner', on='Ko_Number')
    ko_gene_df = pd.merge(ko_gene_df, kegg_gene_def, how='left', on='GeneID')
    
    ko_gene_df = ko_gene_df.drop(columns=['Ko_Number'])
    ko_gene_df.fillna('NA', inplace=True)
    ko_gene_df = pd.merge(ko_gene_df, basicinfo_df, how='left', on='GeneID')
    output_columns_list = basicinfo_df_columns + ['Ko', 'KEGG_ID', 'Gene_shortname', 'EC_number', 'KEGG_def']
    ko_gene_df = ko_gene_df[output_columns_list]
    ko_gene_df.to_csv(output_prefix + '_ko_geneid_def.txt', sep='\t', index=False)
    
    
    ko_gene_df['regulation'] = 1
    ko_gene_df[['KEGG_ID', 'regulation']].to_csv(output_prefix + '_keggid_regulation.txt', sep='\t', index=False)
    
    # 筛选 passed_path.txt 中的 ko
    ko_def_df = pd.read_csv(ko_file, sep='\t', names=['Ko_Number'])
    passed_path_file = '/home/colddata/qinqiang/script/Rscript/pathview/passed_path.txt'
    passed_path_df = pd.read_csv(passed_path_file, sep='\t', names=['Ko_Number', 'Ko_Def'])
    ko_def_df = passed_path_df[passed_path_df['Ko_Number'].isin(ko_def_df['Ko_Number'])]
    ko_def_df = ko_def_df.drop_duplicates(subset=['Ko_Number'])
    ko_def_df.to_csv(output_prefix + '_ko_passed_path.txt', sep='\t', index=False, header=False)
    

def main():
    # kns_def_file = 'Blautia_coccoides_NCTC11035_kns_gene_def.txt'
    args = parse_input()
    # add_def(args.ko_list, args.kegg_genedef, args.kegg_pathway, args.basicinfo, args.output_prefix)
    add_def(args.ko_list, args.kegg_genedef, args.kegg_pathway, args.output_prefix)
    

if __name__ == '__main__':
    main()