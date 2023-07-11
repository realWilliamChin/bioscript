#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/07/04 18:09
# Author        : William GoGo
import os
import argparse
import pandas as pd


def parse_input():
    argparser = argparse.ArgumentParser(description="target gene")
    argparser.add_argument('-i', '--input', required=True, help='input vcf file')
    argparser.add_argument('--input-format', default='vcf', choices=['vcf', 'gff'],
                           help='input file format, vcf or gff')
    argparser.add_argument('-g', '--gff', help='.gff 结尾文件或者 .gff3 结尾的文件')
    argparser.add_argument('-o', '--output', default='output.txt', help='output file')
    return argparser.parse_args()


def target_gene_for_vcf(file, gff_file):
    # Trait	Model	Threshold	Marker	Chrom	Position	Ref	Alt	Score	Effect
    df = pd.read_csv(file, sep='\t', skiprows=318, low_memory=False)
    gff_df = pd.read_csv(gff_file, sep='\t', header=None, usecols=[0, 1, 2, 3, 4, 6, 8], names=["seqid", "source", "type", "start", "end", "strand", "attributes"])
    gff_df = gff_df[gff_df['attributes'].str.contains('ID=')]
    # gff_df = gff_df[gff_df['seqid'].str.contains('chr')]
    df['start'] = df['POS'].apply(lambda x: max(x - 10000, 0))
    df['end'] = df['POS'].apply(lambda x: max(x + 10000, 0))

    lst = df.values.tolist()
    gff_lst = gff_df.values.tolist()
    # with open('covered_gene_for_indel.txt', 'a') as f:
    #     f.write('\t'.join(['Trait', 'Model', 'Threshold', 'Marker', 'Chrom', 'Position', 'Ref', 'Alt', 'Score', 'Effect', 'Target_GeneID', 'Target_Start', 'Target_End']) + '\n')

    for each_line in lst:
        each_line = [str(i) for i in each_line]
        qtl_start = int(each_line[-2])
        qtl_end = int(each_line[-1])
        each_line = each_line[:-2]
        for gff in gff_lst:
            gff = [str(i) for i in gff]
            gff_start = int(gff[3])
            gff_end = int(gff[4])
            append_lst = [gff[-1].split(';')[0].split('=')[1], str(gff_start), str(gff_end)]
            if gff_start <= qtl_start <= gff_end or gff_start <= qtl_end <= gff_end:
                with open('target_gene_for_potato_snp_MI_specific.vcf', 'a') as f:
                    f.write('\t'.join(each_line) + '\t' + '\t'.join(append_lst) + '\n')

    # failed using pandas
    qtl_df = pd.read_csv('covered_gene_for_indel.txt', sep='\t')
    nr_df = pd.read_csv('potato_nr_diamond_nr_gene_def.txt', sep='\t', skiprows=1, header=None, usecols=[0, 2], names=['GeneID', 'Def'])
    qtl_df = pd.merge(qtl_df, nr_df, how='left', left_on='Target_GeneID', right_on='GeneID')
    qtl_df.to_csv('covered_gene_for_indel.txt', sep='\t', index=False)


def target_gene_for_csv(csv_file, gff_file, output):
    """
        csv_file: Marker,Chrom,Position,Ref,Alt,samples
    """
    csv_df = pd.read_csv(csv_file)
    gff_df = pd.read_csv(gff_file, sep='\t', header=None, usecols=[0, 1, 2, 3, 4, 6, 8], names=["seqid", "source", "type", "start", "end", "strand", "attributes"])
    gff_df = gff_df[gff_df['attributes'].str.contains('ID=')]
    csv_df['start'] = csv_df['Position'].apply(lambda x: max(x - 10000, 0))
    csv_df['end'] = csv_df['Position'].apply(lambda x: max(x + 10000, 0))
    csv_df.fillna(value='0', inplace=True)
    
    lst = csv_df.values.tolist()
    gff_lst = gff_df.values.tolist()
    
    for each_line in lst:
        each_line = [str(i) for i in each_line]
        qtl_start = int(each_line[-2])
        qtl_end = int(each_line[-1])
        each_line = each_line[:-2]
        for gff in gff_lst:
            gff = [str(i) for i in gff]
            gff_start = int(gff[3])
            gff_end = int(gff[4])
            append_lst = [gff[-1].split(';')[0].split('=')[1], str(gff_start), str(gff_end)]
            if gff_start <= qtl_start <= gff_end or gff_start <= qtl_end <= gff_end:
                with open(output, 'a') as f:
                    f.write('\t'.join(each_line) + '\t' + '\t'.join(append_lst) + '\n')

def main():
    args = parse_input()
    if not args.gff:
        args.gff = [x for x in os.listdir() if x.endswith('.gff') or x.endswith('.gff3')][0]
    if not args.output:
        args.output = 'output.txt'
    target_gene_for_csv(args.input, args.gff, args.output)
    

if __name__ == '__main__':
    main()
