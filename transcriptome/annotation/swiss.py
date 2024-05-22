#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2/24/2023 5:30 PM
# Author        : WilliamGoGo
import os, sys
import argparse
import pandas as pd
import numpy as np
import subprocess
from loguru import logger


def parse_input():
    parser = argparse.ArgumentParser(description='swiss 注释')
    parser.add_argument('-b', '--blast', help='指定 swiss.blast 文件名称')
    
    parse_blast = parser.add_argument_group('需要解析 blast 添加以下参数')
    parse_blast.add_argument('-p', '--prefix', help='生成文件的前缀(一般是种名加属名), 默认去掉 blast 作为前缀')
    
    annotation = parser.add_argument_group('需要注释添加一下参数')
    annotation.add_argument('--cds', help='输入需要去注释的 cds 文件（去重，名字精简）')
    annotation.add_argument('-t', '--threads', default=30, type=int, help='运行 swiss 注释线程数量')
    
    args = parser.parse_args()
    return args


def swiss_annotation(fasta_file, blast_file, num_threads):  
    anno_cmd = f'/opt/biosoft/ncbi-blast-2.9.0+/bin/blastx \
        -db /home/data/ref_data/Linux_centos_databases/2019_Unprot_databases/swissprot \
        -query {fasta_file} \
        -out {blast_file} \
        -max_target_seqs 20 \
        -evalue 1e-5 \
        -num_threads {num_threads} \
        -outfmt "6 qacc sacc pident qcovs qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle"'

    ret = subprocess.run(anno_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{fasta_file} swiss 注释程序失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        logger.success(f"{fasta_file} swiss 注释成功，输出到 {blast_file}")
    return True


def _keep_goid(s):
    """
    处理表里的元素，只留下 GO number，格式 GO:xxxxxx
    """
    if not pd.isna(s):
        if 'O:' in s and ']' in s:
            return 'G' + s.split(']')[0]
        else:
            return ''


def process_go(swiss_df, ref_df):
    df_lst = []
    for each_go in ['GO_BP', 'GO_CC', 'GO_MF']:
        # 分成 bp，cc，mf 的
        df = pd.merge(left=swiss_df[['GeneID', 'GOID']], right=ref_df[['GOID', each_go]], on='GOID', how='left')
        df = df.dropna().drop(columns=['GOID'])

        # 去除掉每个基因少于 3 个的 go
        df = df[df[each_go].str.count('GO') >= 3]

        # 把后面的 GO 功能都加到第二列
        df_expand = df[each_go].str.split(';', expand=True)
        df = pd.concat([df['GeneID'], df_expand], axis=1).melt(id_vars=['GeneID'], value_name='GOID')

        df = df.drop(columns=['variable']).dropna()

        # 对后面的 GOID 格式改成 ID 在前，使用 _ 连接
        df_expand = df['GOID'].str.split('\[GO', expand=True)
        df['Def'] = df_expand[1].str.replace(':', 'GO:', regex=True).str.replace(']', '_', regex=True) + df_expand[0].str.strip()
        # 替换某些字符，在 Funrich 中会出错的
        df['Def'] = df['Def'].str.replace(', ', '_', regex=True)
        df['Def'] = df['Def'].str.replace('\'', '', regex=True)
        df['Def'] = df['Def'].str.replace('/', '_', regex=True)
        df['Def'] = df['Def'].str.replace(',', '_', regex=True)
        df = df.drop(columns=['GOID']).sort_values('GeneID')
        df_lst.append(df)
    return df_lst


def main():
    args = parse_input()
    
    if args.cds:
        ret = swiss_annotation(args.cds, args.blast, args.threads)
        if not ret:
            sys.exit(1)
    
    swiss_file = args.blast if args.blast else [x for x in os.listdir() if '_swiss.blast' in x][0]
    key_name = args.prefix + '_swiss' if args.prefix else swiss_file.replace('.blast', '')
    # 读取 swiss 参考文件和 blast 文件，并初始化
    swiss_df = pd.read_csv(swiss_file, sep='\t', usecols=[0, 1, 14, 15], names=['GeneID', 'GOID', 'bitscore', 'Swiss_Def'], low_memory=False)
    swiss_df = swiss_df.sort_values(by=['GeneID', 'bitscore'], ascending=['True', 'False'])
    swiss_df = swiss_df.drop(columns=['bitscore'])
    swiss_df = swiss_df.drop_duplicates(subset='GeneID', keep='first')
    siwss_df_expand = swiss_df['Swiss_Def'].str.split(';', expand=True)
    swiss_df['Swiss_Def'] = siwss_df_expand.iloc[:, 0]
    # 在 swiss_gene_def 中 Swiss_Def 删掉 RecName: Full=
    swiss_df['Swiss_Def'] = swiss_df['Swiss_Def'].str.replace('RecName: Full=', '')
    # 生成 _unigene_swiss_gene_def.txt
    swiss_df.to_csv(key_name + '_gene_def.txt', sep='\t', index=False, header=['GeneID', 'Swissprot_ID', 'Swiss_Def'])

    ref_file = '/home/data/ref_data/db/swiss_go_txt/Swiss_protein_go.txt'
    ref_df = pd.read_csv(ref_file, sep='\t', skiprows=1, names=['GOID', 'GO_BP', 'GO_CC', 'GO_MF'])

    # 生成 idNO_def 文件
    idNO_def_filename = key_name + '_idNo_def.txt'
    gene_go_filename = key_name + '_gene_go.txt'
    ref_df['merge_go'] = ref_df['GO_BP'] + '_' + ref_df['GO_CC'] + '_' + ref_df['GO_MF']
    ref_df['merge_go'].replace('', np.nan, regex=True, inplace=True)
    idNo_def = pd.merge(left=swiss_df.iloc[:, [0, 1]], right=ref_df.iloc[:, [0, 4]], on='GOID', how='left')
    idNo_def = idNo_def.dropna().drop(columns=['GOID'])
    idNo_def_expand = idNo_def['merge_go'].str.split('\[G', expand=True)
    idNo_def_expand = idNo_def_expand.applymap(_keep_goid)
    idNo_def = pd.concat([idNo_def['GeneID'], idNo_def_expand], axis=1).fillna('')
    idNo_def.to_csv(idNO_def_filename, sep='\t', index=False, header=False)

    with open(idNO_def_filename, 'r') as idNo_def_file:
        os.remove(idNO_def_filename)
        for line in idNo_def_file.readlines():
            line = line.replace('\t\t', '\t')
            line = line.strip() + '\n'
            with open(idNO_def_filename, 'a') as f, open(gene_go_filename, 'a') as gene_go:
                f.write(line)
                parts = line.strip().split('\t')
                g_id = parts[0]
                go_terms = parts[1:]
                for go_term in go_terms:
                    if go_term:
                        gene_go.write(f"{g_id}\t{go_term}\n")

    result_df = process_go(swiss_df, ref_df)
    # 保存 GO_BP GO_CC GO_MF 文件
    
    result_df[0].to_csv(key_name + '_GO_BP_ID.txt', sep='\t', index=False, header=False)
    result_df[1].to_csv(key_name + '_GO_CC_ID.txt', sep='\t', index=False, header=False)
    result_df[2].to_csv(key_name + '_GO_MF_ID.txt', sep='\t', index=False, header=False)
    
    logger.success('Done')


if __name__ == '__main__':
    main()

