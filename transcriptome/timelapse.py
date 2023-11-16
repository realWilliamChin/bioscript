#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/31 15:17
# Author        : William GoGo
# ---------------------------------
# 配合 timelapse.r 使用
# ---------------------------------
import argparse
import pandas as pd


def parse_input():
    argparser = argparse.ArgumentParser(description="")
    argparser.add_argument('-f', '--fpkm', type=str, default="fpkm_matrix_filtered.txt",
                           help="输入 fpkm_matrix_filter.txt 文件，默认当前目录下的那个文件")
    argparser.add_argument('-s', '--sample', type=str, default="samples_described.txt",
                           help="输入当前目录下 samples_described.txt，默认当前目录下的那个文件")
    argparser.add_argument('-o', '--output', type=str, default='fpkm_matrix_filter_timelapse.txt',
                           help='输入输出文件名')
    return argparser.parse_args()


def timelapse(fpkm_file, samples_file, output_name):
    fpkm_df = pd.read_csv(fpkm_file, sep="\t")
    samples_df = pd.read_csv(samples_file, sep="\t", usecols=[0, 1])
    samples_lst = samples_df['sample'].values.tolist()
    group_lst = samples_df.copy().drop_duplicates(subset='group', keep='first')
    group_lst = group_lst['group'].values.tolist()
    
    # print(samples_lst)
    result_df = fpkm_df[['GeneID'] + samples_lst].copy()
    samples_df = samples_df.groupby('group')['sample'].apply(list).to_dict()
    for key, values in samples_df.items():
        result_df[key] = result_df[values].mean(axis=1)
    result_df.drop(columns=samples_lst, inplace=True)
    result_df = result_df.reindex(columns=[result_df.columns[0],] + group_lst)
    result_df.to_csv(output_name, sep='\t', index=False)


def main():
    args = parse_input()
    timelapse(args.fpkm, args.sample, args.output)


if __name__ == '__main__':
    main()