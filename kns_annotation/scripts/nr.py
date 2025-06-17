#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/14 17:05
# Author        : William GoGo
"""
nr 注释程序
生成 nr.blast 文件，nr_gene_def.txt 文件 和 nr_TF_def.txt 文件
"""
import argparse
import os, sys
import random
import string
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from loguru import logger
sys.path.append('/home/colddata/qinqiang/script/CommonTools/Fasta/')
from get_sequence_from_list import get_seq_from_idlist
sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(description='输入 nr.blast 文件的路径')
    # parser.add_argument(dest='running_type', choices=['all', 'blast', 'parse_blast'], help='运行类型')
    
    parser.add_argument('-b', '--blast', type=str, help='nr.blast file')
    parser.add_argument('--basicinfo', type=str, dest='basicinfo', 
                        help='basicinfo 文件，如果 Gene_Def 都是 NA，则无需添加')
    parser.add_argument('-o', '--output-prefix', dest='output_prefix',
                        help='输出文件的前缀')
    parser.add_argument('--not-annotationed', action='store_true', dest='not_annotationed',
                        help='输出没有注释上的基因信息（必须输入 --fasta 参数）')
    
    annotation = parser.add_argument_group('需要注释添加一下参数')
    annotation.add_argument('-f', '--fasta', help='输入需要去注释的 cds fasta 或 unigene fasta 文件（去重，名字精简）')
    annotation.add_argument('-n', '--outfmt', help="nr 输出格式",
                            default="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,stitle")

    args = parser.parse_args()
    
    if not args.output_prefix.endswith('_'):
        args.output_prefix = args.output_prefix + '_'
        
    return args


def nr_annotation(fasta_file, blast_file, outfmt):
    random_char = ''.join(random.choices(string.ascii_letters + string.digits, k=8))
    temp_dir = f'./temp{random_char}'
    os.mkdir(temp_dir)
    anno_cmd = f'diamond blastx --db /home/data/ref_data/db/diamond_nr/diamond_nr \
        --query {fasta_file} \
        --out {blast_file} \
        --outfmt 6 {outfmt} \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --id 30 \
        --block-size 20.0 \
        --tmpdir {temp_dir} \
        --index-chunks 1'

    ret = subprocess.run(anno_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{fasta_file} nr 注释程序失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        logger.success(f"{fasta_file} nr 注释成功，输出到 {blast_file}")
    os.rmdir(temp_dir)
    return True


def nr(nr_blast_file, gene_basicinfo_file, nr_columns, output_prefix):
    data_frame = load_table(nr_blast_file, names=nr_columns, dtype={'GeneID': str})

    data_frame = data_frame.sort_values(by=['qseqid', 'bitscore'], ascending=[True, False])
    data_frame = data_frame.drop_duplicates(subset=['qseqid'], keep='first')
    write_output_df(data_frame, output_prefix +  'nr_uniq_blast.txt', index=False)
    nr_gene_def_df = data_frame[['qseqid', 'sseqid', 'stitle']].copy()
    nr_gene_def_df.columns = ['GeneID', 'NCBI_ID', 'NR_Def']
    nr_gene_def_df['NR_Def'] = nr_gene_def_df['NR_Def'].str.split(n=1).str[1]
    print(f'注释上的基因数量是 {nr_gene_def_df.shape[0]} 个')
    
    # 有参基因中没有注释上的 gene 的 biotype 类型，如果是无参不需要此步骤
    if gene_basicinfo_file and os.path.exists(gene_basicinfo_file):
        nr_gene_def_df = nr_def_add_not_protein_coding(nr_gene_def_df, gene_basicinfo_file)
    else:
        print('没有检测到 gene_basicinfo 文件，或文件输入有误，如没有输入 -b 参数，忽略此条消息')
    
    write_output_df(nr_gene_def_df, output_prefix + 'nr_gene_def.txt', index=False)
    
    nr_TF_def_df = nr_gene_def_df[nr_gene_def_df['NR_Def'].str.contains('transcription')]
    nr_TF_def_df = nr_TF_def_df.sort_values(by='NR_Def', key=lambda x: x.str.lower())
    write_output_df(nr_TF_def_df, output_prefix + 'nr_TF_def.txt', index=False)


def nr_def_add_not_protein_coding(nr_gene_def_df, gene_basicinfo_file):
    gff_basicinfo_df = load_table(gene_basicinfo_file, usecols=['GeneID', 'Gene_Def'], dtype={'GeneID': str})
    print(f'总基因数量是 {gff_basicinfo_df.shape[0]} 个')
    
    # 没有注释上的
    gene_protein_coding_df = gff_basicinfo_df[gff_basicinfo_df['Gene_Def'].str.contains('protein_coding')].copy()
    non_annotationed_df = gene_protein_coding_df[~gene_protein_coding_df['GeneID'].isin(nr_gene_def_df['GeneID'])]
    non_annotationed_df = non_annotationed_df.rename(columns={'Gene_Def': 'NR_Def'})
    logger.info(f'protein coding 没有注释上的有 {non_annotationed_df.shape[0]} 个')
    
    # 不是 protein_coding 的
    gene_non_protein_coding_df = gff_basicinfo_df[~gff_basicinfo_df['Gene_Def'].str.contains('protein_coding')].copy()
    gene_non_protein_coding_df = gene_non_protein_coding_df.rename(columns={'Gene_Def': 'NR_Def'})
    logger.info(f'不是 protein_coding 没有注释上的有 {gene_non_protein_coding_df.shape[0]} 个')
    
    result = pd.concat([nr_gene_def_df, non_annotationed_df, gene_non_protein_coding_df], axis=0)
    result = result.drop_duplicates(subset='GeneID', keep='first')
    result = result.sort_values(by='GeneID')
    logger.info(f'加上没有注释上的和不是 protein_coding 的有 {result.shape[0]} 个')
    if gff_basicinfo_df.shape[0] == result.shape[0]:
        logger.success('\n注释结果正确！')
    
    return result


def get_not_annotationed_fasta(fasta, nr_uniq_blast_file, output_prefix):
    nr_uniq_blast_df = load_table(nr_uniq_blast_file, header=0, usecols=[0], names=['GeneID'], dtype={'GeneID': str})
    get_seq_from_idlist(nr_uniq_blast_df, fasta, 'off', output_prefix + "nr_not_annotationed.fasta")


def specie_count(uniq_blast_file, output_file):
    # output_file species_count.txt"
    with open(uniq_blast_file, "r") as f:
        lines = f.readlines()[1:]

    # 提取倒数第二列（即 [][] 之间的内容）
    extracted = [line.split("[")[-1].split("]")[0] for line in lines if '[' in line]

    counts = Counter(extracted)
    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)

    with open(output_file, "w") as f:
        for item, count in sorted_counts:
            f.write(f"{item}\t{count}\n")


def spcie_count_plot(specie_count_file, output_pic_file):
    species = load_table(specie_count_file, header=None, dtype={0: str})
    species = species.loc[:9, :]
    species.columns = ["Species", "Count"]

    # 创建饼图
    plt.figure(figsize=(10, 7))
    # 绘制饼图，不显示标签
    patches = plt.pie(species['Count'], labels=None, startangle=90)[0]

    # 设置标题
    plt.title("Species distribution")

    # 隐藏坐标轴
    plt.axis('equal')  # 保证饼图为正圆形

    # 添加图例，按数据顺序排列
    plt.legend(
        patches,
        species['Species'],
        title='Species',
        loc='center left',
        bbox_to_anchor=(1, 0.5)  # 将图例放在饼图右侧
    )

    plt.tight_layout()
    plt.savefig(
        output_pic_file,
        dpi=320,
        bbox_inches='tight',  # 包含图例
        pad_inches=0.1        # 减少周围的空白
    )
    plt.close()


def identity_plot(uniq_blast_file, output_pic_file):
    # 读取数据
    id_nr = load_table(uniq_blast_file, skiprows=1, usecols=[2], dtype={0: str})
    id_data = id_nr.iloc[:, 0]

    # 创建分箱区间并统计
    n_bins = 10
    bins = pd.cut(id_data, bins=n_bins, include_lowest=True)
    id_count = bins.value_counts().reset_index()
    id_count.columns = ['Identity', 'Count']  # 重命名列
    id_count = id_count.sort_values('Identity')  # 按区间排序

    # 计算百分比标签
    total = id_count['Count'].sum()
    id_count['Percentage'] = id_count['Count'] / total * 100
    newlegend = [f"{row.Identity} ({row.Percentage:.2f}%)" 
                for _, row in id_count.iterrows()]

    # 创建饼图
    plt.figure(figsize=(10, 7))
    patches, texts, autotexts = plt.pie(
        id_count['Count'],
        labels=None,
        startangle=90,
        wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'},
        autopct=''  # 禁用默认百分比标签
    )

    # 添加自定义图例
    plt.legend(
        patches,
        newlegend,
        title='Identity',
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        fontsize=8
    )

    plt.axis('equal')
    plt.title("Identity distribution", pad=20)

    # 保存图片
    plt.savefig(
        output_pic_file,
        dpi=320,
        bbox_inches='tight',
        pad_inches=0.1
    )
    plt.close()


def main():
    args = parse_input()
    out_prefix = args.output_prefix
    
    if args.fasta:
        outfmt = args.outfmt.replace(',', ' ')
        ret = nr_annotation(args.fasta, args.blast, outfmt)
        if not ret:
            sys.exit(1)
    else:
        if not args.blast:
            args.blast = [x for x in os.listdir() if x.endswith('nr.blast')][0]
    nr_columns = args.outfmt.split(',')
    nr(args.blast, args.basicinfo, nr_columns, out_prefix)
    
    if args.not_annotationed:
        get_not_annotationed_fasta(args.fasta, args.blast, out_prefix)
    
    specie_count(out_prefix + 'nr_uniq_blast.txt', out_prefix + 'nr_specie_count.txt')
    spcie_count_plot(out_prefix + 'nr_specie_count.txt', out_prefix + 'sp_distribution_Top10.jpeg')
    identity_plot(out_prefix + 'nr_uniq_blast.txt', out_prefix + 'id_distribution.jpeg')
    
    logger.success('Done!')


if __name__ == '__main__':
    main()