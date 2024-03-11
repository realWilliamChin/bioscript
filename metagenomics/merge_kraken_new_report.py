#!/usr/bin/env python3
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
    args.add_argument('-s', '--samples', help='samples_described.txt', default='samples_described.txt')
    
    args = args.parse_args()
    # 判断 samples_described.txt 是否存在
    if not os.path.exists(args.samples):
        print(f"Error: {args.samples} not found!")
        exit(1)

    return args


def filter_by_taxonomic_level(df, level):
    levels = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    level_index = levels.index(level)
    # 只保留所需级别的行，并确保不包含更高的级别
    return df[df['specie'].str.contains(level) & ~df['specie'].str.contains('|'.join(levels[level_index+1:]))]


def proportion(df):
    df = df.T
    df['count'] = df.sum(axis=1)
    
    df = df.div(df['count'], axis=0).mul(100)
    df = df.drop(columns=['count']).T
    
    return df


def barplot(df):
    df['count'] = df.sum(axis=1)
    df = df.sort_values(by='count', ascending=False)
    # 取前 15 行
    df = df.iloc[:15, :]
    # 最后一行为 others，等于 100- 每列前 15 行的和，others 放到 df 的最后一行
    for i in df.columns:
        if i == 'count':
            df.loc['others', i] = (len(df.columns) * 100) - df[i].iloc[:15].sum()
        else:
            df.loc['others', i] = 100 - df[i].iloc[:15].sum()
    
    return df


def merge_to_excel(raw_count, relative_abundance, output_file):
    """
    raw_count and proportion merge to excel, raw_count is sheet1 and proportion is sheet2
    """
    with pd.ExcelWriter(output_file) as writer:
        raw_count.to_excel(writer, sheet_name='RawCount', index=True)
        relative_abundance.to_excel(writer, sheet_name='Relative_Abundance', index=True)


def venn_summary(df_list, samples_df_grouped):
    #TODO: 做一个统计，下面这样的，然后输出到文件，文件名为 venn/summary.txt，以后是要插入到报告里
    # Taxonomy	B	D
    # Phylum	37	37
    # Class	73	73
    # Order	168	168
    # Family	383	384
    # Genus	1442	1447
    # Species	5516	5596
    pass


def main():
    args = parse_input()
    os.mkdir('barplot')
    os.mkdir('venn')
    os.mkdir('summary')
    filter_specie = 'k__Bacteria'
    
    samples_df = pd.read_csv(args.samples, sep='\t', usecols=[0, 1], names=['group', 'sample'], skiprows=1)
    
    file_list = [x + '.new.report' for x in samples_df['sample'].values]
    print(f'{len(file_list)} kraken2 report files.')
    
    df = pd.read_csv(file_list[0], sep='\t', names=['specie', file_list[0].split(os.sep)[-1].replace('.new.report', '')])

    for f in file_list[1:]:
        right_df = pd.read_csv(f, sep='\t', names=['specie', f.split(os.sep)[-1].replace('.new.report', '')])
        df = pd.merge(df, right_df, on='specie', how='outer')

    df.fillna(0, inplace=True)
    df.sort_values(by='specie', inplace=True)
    df = df.reindex(columns=['specie'].extend(file_list))
    df_bacteria = df[df['specie'].str.startswith(filter_specie)]
    df.set_index('specie').astype(int).to_csv('merge_kraken_report_result.txt', sep='\t')
    df_bacteria.set_index('specie').astype(int).to_csv('merge_kraken_report_filter_bacteria.txt', sep='\t')
    # 分出界门纲目科属种
    species_df = df[df['specie'].str.contains('s__')].rename(columns={'specie': 'Species'})
    genus_df = filter_by_taxonomic_level(df, 'g__').rename(columns={'specie': 'Genus'})
    family_df = filter_by_taxonomic_level(df, 'f__').rename(columns={'specie': 'Family'})
    order_df = filter_by_taxonomic_level(df, 'o__').rename(columns={'specie': 'Order'})
    class_df = filter_by_taxonomic_level(df, 'c__').rename(columns={'specie': 'Class'})
    phylum_df = filter_by_taxonomic_level(df, 'p__').rename(columns={'specie': 'Phylum'})
    kingdom_df = filter_by_taxonomic_level(df, 'k__').rename(columns={'specie': 'Kingdom'})
    
    df_list = [species_df, genus_df, family_df, order_df, class_df, phylum_df, kingdom_df]
    
    for df in df_list:
        df_name = df.columns[0]
        print(f'Processing {df_name} （只留下 k__Bacteria）...')
        
        df = df[df[df_name].str.contains("k__Bacteria")]
        # 去掉 " , # " 等字符, 在 R 语言中，这些字符会出错"
        df = df.replace({"'": '', '"': '', "#": ''}, regex=True)
        
        df = df.set_index(df_name)
        df = df.astype(int)
        df.to_csv(f'summary/{df_name}.txt', sep='\t')
        
        # proportion
        proportion_df = proportion(df)
        proportion_df.to_csv(f'summary/{df_name}_relative_abundance.txt', sep='\t', index=True)
        
        # raw_count and relative_abundance merge to excel, raw_count is sheet1 and proportion is sheet2
        merge_to_excel(df, proportion_df, f'summary/{df_name}.xlsx')
        
        # for barplot
        barplot_df = barplot(proportion_df)
        barplot_df.to_csv(f'barplot/{df_name}_Top15.txt', sep='\t', index=True)
        
        # for venn
        samples_df_grouped = samples_df.groupby('group')['sample']
        # 循环 samples_df_grouped
        for name, group in samples_df_grouped:
            if df_name == 'Kingdom':
                continue
            sample_list = group.to_list()
            cur_group_df = df[sample_list].copy()
            cur_group_df['count'] = cur_group_df.sum(axis=1)
            cur_group_df = cur_group_df[cur_group_df['count'] > 0]
            cur_group_df.index.to_series().to_csv(f'venn/{df_name}_{name}.txt', header=False, index=False) 


if __name__ == '__main__':
    main()