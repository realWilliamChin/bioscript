# -*- coding: UTF-8 -*-
# Created Time  : 2023/6/6 20:32
# Author        : WilliamGoGo
import os
import pandas as pd
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description='biogrid，')
    parser.add_argument('-f', '--fasta', required=True, help='输入要比对的 fasta 文件')
    # parser.add_argument('-b', '--blast', required=True, help='输入 blast 文件')
    parser.add_argument('-r', '--ref', required=True, help='输入参考文件')
    parser.add_argument('-o', '--pepoutput', required=True, help='输入 pep fasta 的文件名')
    parser.add_argument('-c', '--cpu', help='输入线程数', default=16)
    return parser.parse_args()


def exec_blast(fasta_file, num_threads, pep_file_name):
    blast_file_name = pep_file_name.replace('.fasta', '.blast')
    # cds 转成 pep，后续需要清理一下序列
    os.system(f'seqkit translate {fasta_file} > {pep_file_name}')
    # pep 文件清理 _frame=1
    os.system(f"sed 's/_frame=1//g' -i {pep_file_name}")
    
    for file in os.listdir('00_Database'):
        if 'phr' in file:
            db_name = file.split('.')[0]
    blast_command = f'/opt/biosoft/ncbi-blast-2.9.0+/bin/blastp -db 00_Database/{db_name} -query {pep_file_name} -out {blast_file_name} -evalue 1e-5 -num_threads {num_threads} -outfmt "6 qacc sacc qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle"'
    os.system(blast_command)
    return blast_file_name


def drop_dup(blast_file, out_file):
    blast_names = ['qacc','sacc','qcovhsp','ppos','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','stitle']
    df = pd.read_csv(blast_file, sep='\t', names=blast_names)
    df.sort_values(by=['qacc', 'qcovhsp'], ascending=[True, False], inplace=True)
    df.drop_duplicates(subset='qacc', keep='first', inplace=True)
    df.to_csv(out_file, sep='\t', index=False)


def biogrid(in_blast_file, ref_file, out_file):
    dic = {}
    with open(in_blast_file, 'r') as f:
        num1 = 0
        for each_line1 in f:
            if num1 == 0:
                num1 = 1
                title1 = each_line1
            else:
                value = each_line1.split('\t')[0]
                key = each_line1.strip().split('\t')[1]
                if key not in dic:
                    dic[key] = value
                else:
                    dic[key] += ';' + value

    # Write results to file
    f2 = open(out_file, 'w')
    f2.write(
        'BioGRID Interaction ID\tEmbl_A\tEmbl_B\tThroughput\tSource_Database\n')
    with open(ref_file, 'r') as f1:
        num = 0
        for each_line in f1:
            if num == 0:
                title_part = each_line.split('\t')
                num = 1
            else:
                gene1 = each_line.strip().split('\t')[1]
                gene2 = each_line.strip().split('\t')[2]
                id = each_line.strip().split('\t')[0]
                other = each_line.strip().split('\t')[3]
                source_database = each_line.strip().split('\t')[4]
                if gene1 in dic:
                    if ';' in dic[gene1]:
                        geneAs = dic[gene1].split(';')
                        for each in geneAs:
                            geneA = each
                            if gene2 in dic:
                                if ';' in dic[gene2]:
                                    geneBs = dic[gene2].split(';')
                                    for each in geneBs:
                                        geneB = each
                                        content = id + '\t' + geneA + '\t' + geneB + '\t' + other + '\t' + source_database + '\n'
                                        f2.write(content)
                                else:
                                    geneB = dic[gene2]
                                    content = id + '\t' + geneA + '\t' + geneB + '\t' + other + '\t' + source_database + '\n'
                                    f2.write(content)
                    else:
                        geneA = dic[gene1]
                        if gene2 in dic:
                            if ';' in dic[gene2]:
                                geneBs = dic[gene2].split(';')
                                for each in geneBs:
                                    geneB = each
                                    content = id + '\t' + geneA + '\t' + geneB + '\t' + other + '\t' + source_database + '\n'
                                    f2.write(content)
                            else:
                                geneB = dic[gene2]
                                content = id + '\t' + geneA + '\t' + geneB + '\t' + other + '\t' + source_database + '\n'
                                f2.write(content)


def main():
    args = parse_arguments()
    blast_file = exec_blast(args.fasta, args.cpu, args.pepoutput)
    
    uniq_blast = blast_file.replace('.blast', '_uniq.blast')
    result_file = 'Biogrid_PPI_relation_from_' + str(uniq_blast).replace('.blast', '') + '.txt'
    
    drop_dup(blast_file, uniq_blast)
    biogrid(uniq_blast, args.ref, result_file)

if __name__ == '__main__':
    main()
