#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 5/30/2023 10:53 AM
# Author        : WilliamGoGo
import os
import pandas as pd
import argparse


def parse_input():
    """
    解析输入参数
    """
    parser = argparse.ArgumentParser(description='指定 keg 文件，对 keg 文件清理，处理拆分')
    parser.add_argument('-k', '--keg', type=str, help='指定 keg 文件，默认当前文件夹 keg 结尾的文件')
    parser.add_argument('-t', '--type', type=str, required=True, choices=['plant', 'animal'],
                        help='指定物种类型，植物=plant，动物=animal')
    parser.add_argument('-i', '--allid', type=str, help='指定 all_id 文件，用来生成 shortname.txt')
    args = parser.parse_args()
    if args.keg:
        keg_file = args.keg
    else:
        keg_file = [x for x in os.listdir() if x.endswith('.keg')][0]
    
    return keg_file, args.type, args.allid


def process_keg():
    parseinput = parse_input()
    keg_file, specie_type, all_id_file = parseinput[0], parseinput[1], parseinput[2]
    key_name = keg_file.replace('.keg', '')
    kegg_file = keg_file.replace('.keg','_KEGG_original.txt')
    # 初始处理 keg 文件
    command = 'perl /home/colddata/chen/03_transcript/annotation/kegg/kaas_parse_keg_file.pl -i {} -o {}'.format(keg_file, kegg_file)
    os.system(command)
    
    with open(kegg_file,'r') as f1:
        kegg_clean_file = key_name + '_KEGG_clean.txt'
        f2 = open(kegg_clean_file, 'w')
        
        # 植物物种过滤动物的一些注释
        if specie_type == 'plant':
            for each_line in f1:
                if 'A09160' in each_line:
                    continue
                if 'A09190' in each_line:
                    continue
                if 'A09150' in each_line:
                    if '09158:Development and regeneration' not in each_line and '09159:Environmental adaptation' not in each_line:
                        continue
                f2.write(each_line)
                
        # 其他物种暂时不用过滤，直接写入
        else:
            for each_line in f1:
                f2.write(each_line)
        f2.close()
    
    # 生成 KEGG_gene_def 文件
    kegg_gene_def_titles = ['GeneID', 'Pathway', 'Level2', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number']
    kegg_clean_df = pd.read_csv(kegg_clean_file, sep='\t', header=None, names=kegg_gene_def_titles, dtype={"GeneID": str})
    gene_def_df = kegg_clean_df.drop(columns=['Pathway', 'Level2', 'Level1'])
    gene_def_df.drop_duplicates(subset=['GeneID', 'KEGG_ID'], keep='first', inplace=True)
    print('生成 KEGG_gene_def 文件，去重前的数量:{}，去重后的数量:{}'.format(kegg_clean_df.shape[0]-1, gene_def_df.shape[0]-1))
    gene_def_df['EC_number'] = gene_def_df['Description EC_number'].str.split('[', expand=True)[1].str.replace(']', '').str.split(' ', expand=True)[0]
    gene_def_df['KEGG_def'] = gene_def_df['Description EC_number'].str.split('[', expand=True)[0].str.strip()
    gene_def_df.fillna(value='NA', inplace=True)
    # Gene_shortname 不能设置空为 NA，设置为空
    gene_def_df['Gene_shortname'] = gene_def_df['Gene_shortname'].str.split(',', expand=True)[0]
    gene_def_df['Gene_shortname'].fillna(value='', inplace=True)
    gene_def_df.drop(columns='Description EC_number', inplace=True)
    gene_def_df.to_csv(key_name + '_KEGG_gene_def.txt', sep='\t', index=False)
    
    # 如果指定了 -i allgeneid 文件，则生成 all_gene_id + KEGG gene_short_name 新文件
    if all_id_file:
        all_gene_df = pd.read_csv(all_id_file, sep='\t', names=['GeneID'], dtype={"GeneID": str})
        all_gene_df = pd.merge(left=all_gene_df, right=gene_def_df[['GeneID', 'Gene_shortname']], on='GeneID', how='left')
        all_gene_df.fillna(value='', inplace=True)
        all_gene_df.to_csv(key_name + '_shortname.txt', sep='\t', index=False)
    
    # 生成 tier2 和 tier3 文件
    tier2_name = key_name + '_KEGG_tier2.txt'
    tier3_name = key_name + '_KEGG.txt'
    
    tier2_df = kegg_clean_df.copy()
    tier2_df.drop(columns=['Pathway', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number'], inplace=True)
    tier2_df['Level2'] = tier2_df['Level2'].str.replace('\\','').str.replace(' / ', '_').str.replace('/', '_').str.replace(', ', '').str.replace(',', '_')
    tier2_df.drop_duplicates(keep='first', inplace=True)
    tier2_df.to_csv(tier2_name, sep='\t', index=False, header=None)
    
    tier3_df = kegg_clean_df.copy()
    tier3_df.drop(columns=['Level2', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number'], inplace=True)
    tier3_df['Pathway'] = tier3_df['Pathway'].str.replace('\\','').str.replace(' / ', '_').str.replace('/', '_').str.replace(', ', '').str.replace(',', '_')
    tier3_df.drop_duplicates(keep='first', inplace=True)
    tier3_df.to_csv(tier3_name, sep='\t', index=False, header=None)
    
    
if __name__ == '__main__':
    process_keg()