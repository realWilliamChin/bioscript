#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/12/14 11:27
# Author        : William GoGo
'''
处理 kraken2 的报告文件，将所有样本的报告文件合并成一个文件，并分出界门纲目科属种
'''
import os
import argparse
import pandas as pd


def parse_input():
    args = argparse.ArgumentParser(description='merge_kraken_new_report.py')
    args.add_argument('-i', '--input', help='input kraken2 report dir', default='.')
    args.add_argument('-s', '--samples', help='samples_described.txt', default='samples_described.txt')
    
    args = args.parse_args()
    # 判断 samples_described.txt 是否存在
    if not os.path.exists(args.samples):
        print(f"Error: {args.samples} not found!")
        exit(1)

    return args


def main():
    args = parse_input()
    
    samples_df = pd.read_csv(args.samples, sep='\t', usecols=[0, 1], names=['group', 'sample'], skiprows=1)
    
    file_list = [x for x in samples_df['sample'].values]
    file_list = [os.sep.join([args.input, x + '.new.report']) for x in file_list]
    print(f'{len(file_list)} kraken2 report files.')
    
    df = pd.read_csv(file_list[0], sep='\t', names=['specie', file_list[0].split(os.sep)[-1].replace('.new.report', '')])

    for f in file_list[1:]:
        print(f'Processing {f} ...')
        right_df = pd.read_csv(f, sep='\t', names=['specie', f.split(os.sep)[-1].replace('.new.report', '')])
        df = pd.merge(df, right_df, on='specie', how='outer')

    df.fillna(0, inplace=True)
    df.sort_values(by='specie', inplace=True)
    df = df.reindex(columns=['specie'].extend(file_list))
    df.set_index('specie').astype(int).to_csv('merge_kraken_report_result.txt', sep='\t')
    # 分出界门纲目科属种
    species_df = df[df['specie'].str.contains('s__')]
    genus_df = df[df['specie'].str.contains('g__') & ~df['specie'].str.contains('s__')]
    family_df = df[df['specie'].str.contains('f__') & ~df['specie'].str.contains('g__') & ~df['specie'].str.contains('s__')]
    order_df = df[df['specie'].str.contains('o__') & ~df['specie'].str.contains('f__') & ~df['specie'].str.contains('g__') & ~df['specie'].str.contains('s__')]
    class_df = df[df['specie'].str.contains('c__') & ~df['specie'].str.contains('o__') & ~df['specie'].str.contains('f__') & ~df['specie'].str.contains('g__') & ~df['specie'].str.contains('s__')]
    phylum_df = df[df['specie'].str.contains('p__') & ~df['specie'].str.contains('c__') & ~df['specie'].str.contains('o__') & ~df['specie'].str.contains('f__') & ~df['specie'].str.contains('g__') & ~df['specie'].str.contains('s__')]
    kingdom_df = df[df['specie'].str.contains('k__')& ~df['specie'].str.contains('p__') & ~df['specie'].str.contains('c__') & ~df['specie'].str.contains('o__') & ~df['specie'].str.contains('f__') & ~df['specie'].str.contains('g__') & ~df['specie'].str.contains('s__')]
    
    df_list = [species_df, genus_df, family_df, order_df, class_df, phylum_df, kingdom_df]
    name_list = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"]
    df_and_name = zip(df_list, name_list)
    for df, name in df_and_name:
        print(f'Processing {name} ...')
        df = df.set_index('specie')
        df = df.astype(int)
        df.to_csv(f'{name}.txt', sep='\t')
        df = df.T
        df['count'] = df.sum(axis=1)
        
        df = df.div(df['count'], axis=0).mul(100)
        df.drop(columns=['count']).T.to_csv(f'{name}_Summary_count.txt', sep='\t')


if __name__ == '__main__':
    main()