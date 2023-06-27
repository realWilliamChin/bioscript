#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/27 11:00
# Author        : William GoGo
import os
import pandas as pd
import argparse


def parse_input():
    argparser = argparse.ArgumentParser(description='')
    args = argparser.parse_args()
    return args


def all_geneid_goid():
    df_list = []
    for each_doc in os.listdir():
        if 'GO_' in each_doc and '.txt' in each_doc:
            df = pd.read_csv(each_doc, sep='\t', header=None)
            df[1] = df[1].str.split('_', expand=True)[0]
            df_list.append(df)
    pd.concat(df_list).to_csv('all_geneid_goid.txt', sep='\t', header=False, index=False)


if __name__ == '__main__':
    all_geneid_goid()
    args = parse_input()
    cur_path = os.getcwd()
    result_path = cur_path

    dic_GO = {}
    f4 = open('all_geneid_goid.txt', 'r')
    for each_line4 in f4:
        if 'GeneID' in each_line4:
            continue
        gene_id = each_line4.split('\t')[0]
        Go_id = each_line4.strip().split('\t')[1]
        if gene_id not in dic_GO:
            dic_GO[gene_id] = Go_id + '\t'
        else:
            dic_GO[gene_id] += Go_id + '\t'

    for each_doc in os.listdir():
        if '_Down_ID.txt' in each_doc:
            sample = each_doc.split('_Down_ID.txt')[0]
            f5 = open(result_path + os.sep + sample + '_DE_gene_Down_Go.txt', 'w')
            f2 = open(sample + '_Down_ID.txt', 'r')
            for each_line2 in f2:
                gene_id2 = each_line2.split('\n')[0].replace('gene-', '')
                f5.write(gene_id2 + '\t')
                if gene_id2 not in dic_GO:
                    f5.write('\n')
                else:
                    f5.write(dic_GO[gene_id2] + '\n')
            f3 = open(sample + '_Up_ID.txt', 'r')
            f6 = open(result_path + os.sep + sample + '_DE_gene_Up_Go.txt', 'w')
            for each_line3 in f3:
                gene_id3 = each_line3.split('\n')[0].replace('gene-', '')
                f6.write(gene_id3 + '\t')
                if gene_id3 not in dic_GO:
                    f6.write('\n')
                else:
                    f6.write(dic_GO[gene_id3] + '\n')
            f6.close()
    
    