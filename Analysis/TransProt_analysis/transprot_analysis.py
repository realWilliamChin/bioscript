#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2026/03/26 14:46
# Author        : William GoGo
import os, sys
import pandas as pd
import argparse
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script')
from CommonTools.load_input import load_table, write_output_df
from Analysis.enrich_analysis.genes_enrich import genes_enrich


# 定义一个样本描述文件
# Transcriptome \t Protein;
# GroupA-vs-GroupB \t GroupA-vs-GroupB
# GroupC-vs-GroupD \t GroupC-vs-GroupD

def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='输入文件，组描述文件, 格式为 Transcriptome \t Protein')
    parser.add_argument('-t', '--transcriptome-deg-dir', help='输入文件夹，转录组DEG文件夹')
    parser.add_argument('-p', '--protein-deg-dir', help='输入文件夹，蛋白组DEG文件夹')
    parser.add_argument('-v', '--value', help='log2FC 的值')
    parser.add_argument('-d', '--definition', help='定义注释文件，输出的结果添加注释')
    parser.add_argument('-o', '--output-dir', help='输出文件夹')
    
    # enrich_analysis argument
    parser.add_argument('-k', '--kegg-clean', help='kegg_clean.txt 文件')
    parser.add_argument('-g', '--gene-go', help='swiss 注释出来的的 gene_go.txt 文件')
    parser.add_argument('--genesymbol', help='GeneID 替换为 GeneSymbol 进行分析，输入包含 GeneID GeneSymbol 两列的文件')
    
    # pathview analysis argument
    
    
    return parser.parse_args()


def step1_commongene(group_desc_df, transcriptome_deg_dir, protein_deg_dir, output_dir):
    
    group_desc_df['Transcriptome_file'] = group_desc_df['Transcriptome'].apply(lambda x: os.path.join(transcriptome_deg_dir, f"{x}_DEG_data.txt"))
    group_desc_df['Protein_file'] = group_desc_df['Protein'].apply(lambda x: os.path.join(protein_deg_dir, f"{x}_DEG_data.txt"))
    
    # 循环 group_desc_df，每一行对应一个分组，生成对应的交集基因文件和 DEG/蛋白数据文件
    for idx, row in group_desc_df.iterrows():
        group_name = row['Transcriptome']
        transcriptome_file = row['Transcriptome_file']
        protein_file = row['Protein_file']

        # 读取 DEG 和蛋白的基因数据
        try:
            deg_df = load_table(transcriptome_file, usecols=['GeneID', 'log2FoldChange', 'pvalue'])
            deg_df.columns = ['GeneID', 'log2FC_RNA', 'pvalue_RNA']
        except Exception as e:
            logger.error(f'无法读取转录组DEG文件: {transcriptome_file} 错误: {str(e)}')
            sys.exit(1)
        try:
            protein_df = load_table(protein_file, usecols=['GeneID', 'logFC', 'pvalue'])
            protein_df.columns = ['GeneID', 'log2FC_protein', 'pvalue_protein']
        except Exception as e:
            logger.error(f'无法读取蛋白组DEG文件: {protein_file} 错误: {str(e)}')
            sys.exit(1)

        # 求交集
        common_genes = pd.merge(deg_df, protein_df, on='GeneID', how='inner')
        common_genes.fillna(0, inplace=True)
        write_output_df(common_genes, os.path.join(output_dir, f"{group_name}_common_result.txt"), index=False)


def step2_group(group_name, output_dir):
    rna_value = 1
    protein_value = 1
    
    df = load_table(os.path.join(output_dir, f"{group_name}_common_result.txt"), dtype={'GeneID': str})
    df['log2FC_protein'] = df['log2FC_protein'].apply(lambda x: x if x <= 9 else 9)
    df['log2FC_protein'] = df['log2FC_protein'].apply(lambda x: x if x >= -9 else -9)
    df['log2FC_RNA'] = df['log2FC_RNA'].apply(lambda x: x if x <= 9 else 9)
    df['log2FC_RNA'] = df['log2FC_RNA'].apply(lambda x: x if x >= -9 else -9)

    group1 = df[(df['log2FC_protein'] > protein_value) & (df['log2FC_RNA'] < -rna_value)].copy()
    group2 = df[(df['log2FC_protein'] > protein_value) & (df['log2FC_RNA'] > -rna_value) & (df['log2FC_RNA'] < rna_value)].copy()
    group3 = df[(df['log2FC_protein'] > protein_value) & (df['log2FC_RNA'] > rna_value)].copy()
    group4 = df[(df['log2FC_protein'] < protein_value) & (df['log2FC_protein'] > -protein_value) & (df['log2FC_RNA'] < -rna_value)].copy()
    group5 = df[(df['log2FC_protein'] < protein_value) & (df['log2FC_protein'] > -protein_value) & (df['log2FC_RNA'] > -rna_value) & (df['log2FC_RNA'] < rna_value)].copy()
    group6 = df[(df['log2FC_protein'] < protein_value) & (df['log2FC_protein'] > -protein_value) & (df['log2FC_RNA'] > rna_value)].copy()
    group7 = df[(df['log2FC_protein'] < -protein_value) & (df['log2FC_RNA'] < -rna_value)].copy()
    group8 = df[(df['log2FC_protein'] < -protein_value) & (df['log2FC_RNA'] > -rna_value) & (df['log2FC_RNA'] < rna_value)].copy()
    group9 = df[(df['log2FC_protein'] < -protein_value) & (df['log2FC_RNA'] > rna_value)].copy()
    
    group1['Group'] = 'group1'
    group2['Group'] = 'group2'
    group3['Group'] = 'group3'
    group4['Group'] = 'group4'
    group5['Group'] = 'group5'
    group6['Group'] = 'group6'
    group7['Group'] = 'group7'
    group8['Group'] = 'group8'
    group9['Group'] = 'group9'
    
    group3['regulation'] = 'Up'
    group7['regulation'] = 'Down'
    
    result_df = pd.concat([group1, group2, group3, group4, group5, group6, group7, group8, group9])
    result_df.fillna(0, inplace=True)
    write_output_df(result_df, os.path.join(output_dir, f"{group_name}_result_group.txt"), index=False)
    
    write_output_df(group1.drop(columns=['Group']), os.path.join(output_dir, f"{group_name}_group1.txt"), index=False)
    write_output_df(group2.drop(columns=['Group']), os.path.join(output_dir, f"{group_name}_group2.txt"), index=False)
    write_output_df(group3.drop(columns=['Group']), os.path.join(output_dir, f"{group_name}_group3.txt"), index=False)
    write_output_df(group4.drop(columns=['Group']), os.path.join(output_dir, f"{group_name}_group4.txt"), index=False)
    write_output_df(group5.drop(columns=['Group']), os.path.join(output_dir, f"{group_name}_group5.txt"), index=False)
    write_output_df(group6.drop(columns=['Group']), os.path.join(output_dir, f"{group_name}_group6.txt"), index=False)
    write_output_df(group7.drop(columns=['Group']), os.path.join(output_dir, f"{group_name}_group7.txt"), index=False)
    write_output_df(group8.drop(columns=['Group']), os.path.join(output_dir, f"{group_name}_group8.txt"), index=False)
    write_output_df(group9.drop(columns=['Group']), os.path.join(output_dir, f"{group_name}_group9.txt"), index=False)


def step3_enrich_analysis(gene_go, kegg_clean, output_dir, genesymbol=None):
    group_list = [os.path.join(output_dir, x) for x in os.listdir(output_dir) if x.endswith('txt')]
    genes_enrich(group_list, True, gene_go, kegg_clean, output_dir, genesymbol)


def step4_pathview_analysis():
    pass




if __name__ == '__main__':
    args = parse_input()
    os.makedirs(args.output_dir, exist_ok=True)

    # step1: 读取组描述文件，求转录组与蛋白组 DEG 交集
    group_desc_df = load_table(args.input)
    step1_commongene(group_desc_df, args.transcriptome_deg_dir, args.protein_deg_dir, args.output_dir)

    # step2: 对每个分组做九象限分类
    for group_name in group_desc_df['Transcriptome']:
        step2_group(group_name, args.output_dir)

    # step3: 对分组结果做 GO/KEGG 富集
    step3_enrich_analysis(args.gene_go, args.kegg_clean, args.output_dir, args.genesymbol)

    # step4: pathview 通路分析（待实现）
    # step4_pathview_analysis()