#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/30 12:41
# Author        : William GoGo
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from loguru import logger


def parse_input():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-c', '--cog_file', help='.emapper.annotations 文件')
    argparser.add_argument('-o', '--output_prefix', help='输出前缀，输出包括 cog_count 和图片')
    return argparser.parse_args()


def cog_count_plot(cog_count_df, o_pic_name):
    
    palette = sns.color_palette("husl", len(cog_count_df["Category"].unique()))
    
    plt.figure(figsize=(18, 12))
    sns.barplot(data=cog_count_df, x="Category", y="Count", hue="Category", palette=palette, dodge=False)
    handles, labels = plt.gca().get_legend_handles_labels()
    legend_labels = cog_count_df["Function"].unique()
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=palette[i]) for i in range(len(legend_labels))]

    # 添加图例
    plt.legend(legend_handles, legend_labels, title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title("COG Function Count")
    plt.xlabel("Category")
    plt.ylabel("Count")
    # plt.legend(title="Category", labels=cog_count_df["Function"], bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # 保存图表
    plt.savefig(o_pic_name, dpi=320, bbox_inches="tight")


def main():
    args = parse_input()
    cog_file = args.cog_file if args.cog_file else [x for x in os.listdir() if x.endswith('emapper.annotations')][0]

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
    cog_count_source_df = pd.read_csv(cog_file, sep='\t', skiprows=4)
    cog_list = Counter(''.join(list(cog_count_source_df['COG_category'].dropna()))).most_common()
    cog_count_df = pd.DataFrame(list(cog_list))
    cog_df = pd.merge(left=cog_source_df, right=cog_count_df, on=0, how='left').fillna(0)
    cog_columns = ["Category", "Group", "Function", "Count"]
    cog_df.columns = cog_columns
    cog_df.to_csv(args.output_prefix + '_COG_count.txt', sep='\t', header=cog_columns, index=False)
    
    cog_count_plot(cog_df, args.output_prefix + '_COG_count_plot.jpeg')


if __name__ == '__main__':
    main()
