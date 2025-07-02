#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2/24/2023 5:30 PM
# Author        : WilliamGoGo
import os, sys
import argparse
import pandas as pd
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df

if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    logger.critical("当前 Python 版本为:", sys.version)
    sys.exit(1)


def parse_input():
    parser = argparse.ArgumentParser(description='swiss 注释')
    parser.add_argument('-b', '--blast', help='指定 swiss.blast 文件名称')
    
    parse_blast = parser.add_argument_group('需要解析 blast 添加以下参数')
    parse_blast.add_argument('-o', '--output-prefix', dest='output_prefix',
                             help='生成文件的前缀(一般是种名加属名), 默认去掉 blast 作为前缀')
    
    annotation = parser.add_argument_group('需要注释添加一下参数')
    annotation.add_argument('--fasta', help='输入需要去注释的 cds fasta 或 unigene fasta 文件（去重，名字精简）')
    annotation.add_argument('-t', '--threads', default=30, type=int, help='运行 swiss 注释线程数量')
    
    args = parser.parse_args()
    
    if not args.output_prefix.endswith('_'):
        args.output_prefix = args.output_prefix + '_'
    
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


def swiss_bpccmf_barplot(bp_file, cc_file, mf_file, output):
    # 读取 BP 数据
    BP = load_table(bp_file, header=None, dtype={'GeneID': str})
    BP_count = BP[1].value_counts().reset_index()
    BP_count.columns = ["Term", "Count"]
    BP_20 = BP_count.sort_values(by="Count", ascending=False).head(20)

    # 读取 CC 数据
    CC = load_table(cc_file, header=None, dtype={'GeneID': str})
    CC_count = CC[1].value_counts().reset_index()
    CC_count.columns = ["Term", "Count"]
    CC_20 = CC_count.sort_values(by="Count", ascending=False).head(20)

    # 读取 MF 数据
    MF = load_table(mf_file, header=None, dtype={'GeneID': str})
    MF_count = MF[1].value_counts().reset_index()
    MF_count.columns = ["Term", "Count"]
    MF_20 = MF_count.sort_values(by="Count", ascending=False).head(20)

    # 合并数据
    GO_stat = pd.concat([BP_20, CC_20, MF_20], axis=0)
    GO_stat["GO"] = ["BP"] * 20 + ["CC"] * 20 + ["MF"] * 20

    # 清理 Term 列
    GO_stat["Term"] = GO_stat["Term"].str.replace(r"GO:\d+_", "", regex=True)

    # 绘制图表
    plt.figure(figsize=(15, 10))
    sns.barplot(data=GO_stat, x="Term", y="Count", hue="GO", dodge=False)
    plt.title("GO Term Count")
    plt.xlabel("Term")
    plt.ylabel("Count")
    plt.xticks(rotation=75, ha="right", va="top")
    plt.tight_layout()

    plt.savefig(output, dpi=300, bbox_inches="tight")


def main():
    args = parse_input()
    
    if args.fasta:
        ret = swiss_annotation(args.fasta, args.blast, args.threads)
        if not ret:
            sys.exit(1)
    
    swiss_file = args.blast if args.blast else [x for x in os.listdir() if '_swiss.blast' in x][0]
    # 读取 swiss 参考文件和 blast 文件，并初始化
    swiss_df = load_table(swiss_file, usecols=[0, 1, 14, 15], names=['GeneID', 'GOID', 'bitscore', 'Swiss_Def'], dtype={'GeneID': str})
    swiss_df = swiss_df.sort_values(by=['GeneID', 'bitscore'], ascending=[True, False])
    swiss_df = swiss_df.drop(columns=['bitscore'])
    swiss_df = swiss_df.drop_duplicates(subset='GeneID', keep='first')
    siwss_df_expand = swiss_df['Swiss_Def'].str.split(';', expand=True)
    swiss_df['Swiss_Def'] = siwss_df_expand.iloc[:, 0]
    # 在 swiss_gene_def 中 Swiss_Def 删掉 RecName: Full=
    swiss_df['Swiss_Def'] = swiss_df['Swiss_Def'].str.replace('RecName: Full=', '')
    # 生成 _unigene_swiss_gene_def.txt
    write_output_df(swiss_df, args.output_prefix + 'swiss_gene_def.txt', index=False, header=['GeneID', 'Swissprot_ID', 'Swiss_Def'])

    ref_file = '/home/data/ref_data/db/swiss_go_txt/Swiss_protein_go.txt'
    ref_df = load_table(ref_file, skiprows=1, names=['GOID', 'GO_BP', 'GO_CC', 'GO_MF'], dtype={'GOID': str})

    # 生成 idNO_def 文件
    idNO_def_filename = args.output_prefix + 'swiss_idNo_def.txt'
    gene_go_filename = args.output_prefix + 'swiss_gene_go.txt'
    ref_df['merge_go'] = ref_df['GO_BP'] + '_' + ref_df['GO_CC'] + '_' + ref_df['GO_MF']
    ref_df['merge_go'] = ref_df['merge_go'].replace('', np.nan, regex=True)
    idNo_def = pd.merge(left=swiss_df.iloc[:, [0, 1]], right=ref_df.iloc[:, [0, 4]], on='GOID', how='left')
    idNo_def = idNo_def.dropna().drop(columns=['GOID'])
    idNo_def_expand = idNo_def['merge_go'].str.split('\[G', expand=True)
    idNo_def_expand = idNo_def_expand.map(_keep_goid)
    idNo_def = pd.concat([idNo_def['GeneID'], idNo_def_expand], axis=1).fillna('')
    write_output_df(idNo_def, idNO_def_filename, index=False, header=False)

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
    
    bp_file = args.output_prefix + 'GO_BP_ID.txt'
    cc_file = args.output_prefix + 'GO_CC_ID.txt'
    mf_file = args.output_prefix + 'GO_MF_ID.txt'
    write_output_df(result_df[0], bp_file, index=False, header=False)
    write_output_df(result_df[1], cc_file, index=False, header=False)
    write_output_df(result_df[2], mf_file, index=False, header=False)
    
    swiss_bpccmf_barplot(bp_file, cc_file, mf_file, f'{args.output_prefix}GO_count.jpeg')
    
    logger.success('Done')


if __name__ == '__main__':
    main()

