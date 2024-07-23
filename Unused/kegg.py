#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 5/30/2023 10:53 AM
# Author        : WilliamGoGo
import os, sys
import pandas as pd
import argparse

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/transcriptome/'))
from merge_fpkm_reads_matrix import merge_fpkm_reads


def parse_input():
    """
    解析输入参数
    """
    parser = argparse.ArgumentParser(description='指定 keg 文件，对 keg 文件清理，处理拆分')
    parser.add_argument('-k', '--keg', type=str, help='指定 keg 文件，默认当前文件夹 keg 结尾的文件')
    parser.add_argument('-t', '--type', type=str, help='指定物种类型，植物=plant, 动物=animal')
    parser.add_argument('-i', '--allid', type=str, help='指定 all_id 文件，用来生成 shortname.txt')
    
    # 这两个参数用于生成 ko03000_expression_data.txt 文件
    parser.add_argument('--fpkm', type=str,
                        help='指定 fpkm 文件, 用于对 ko03000 添加表达量生成新文件，如果不指定，则不生成 ko03000_expression_data 文件')
    parser.add_argument('--reads', type=str,
                        help='指定 reads 文件, 用于对 ko03000 添加表达量生成新文件, 如果指定需要 fpkm 和 reads 都存在')
    
    args = parser.parse_args()
    if args.keg and os.path.exists(args.keg) is True:
        pass
    else:
        args.keg = [x for x in os.listdir() if x.endswith('.keg')][0]
    
    # 检测 fpkm 和 reads 文件是否有效
    if args.fpkm or args.reads:
        if os.path.exists(args.fpkm) is False:
            raise Exception('fpkm 文件不存在')
        if os.path.exists(args.reads) is False:
            raise Exception('reads 文件不存在')
    
    return args


def process_keg(input_f, output_f):
    """
    TODO: 从 keg 文件开始解析, 顶替 perl 脚本, 使用 python 处理 keg 文件
    """
    pass


def ko03000(kegg_gene_df, kegg_tier3_df):
    ko03000_def_df_file = '/home/colddata/qinqiang/script/transcriptome/annotation/ko03000_def.txt'
    ko03000_def_df = pd.read_csv(ko03000_def_df_file, sep='\t')
    ko03000_geneid_list = kegg_tier3_df[kegg_tier3_df['Pathway'].str.contains("ko03000")]['GeneID'].tolist()
    ko03022_geneid_list = kegg_tier3_df[kegg_tier3_df['Pathway'].str.contains("ko03022")]['GeneID'].tolist()
    
    ko03000_df = kegg_gene_df[kegg_gene_df['GeneID'].isin(ko03000_geneid_list)][['GeneID', 'KEGG_ID']].copy()
    ko03022_df = kegg_gene_df[kegg_gene_df['GeneID'].isin(ko03022_geneid_list)][['GeneID', 'KEGG_ID']].copy()
    # ko03000_df = ko03000_df
    ko03000_df = pd.merge(left=ko03000_df, right=ko03000_def_df, on='KEGG_ID', how='left')
    
    return ko03000_df, ko03022_df

    
def main():
    args = parse_input()
    keg_file, specie_type, all_id_file = args.keg, args.type, args.allid
    key_name = keg_file.replace('.keg', '')
    kegg_file = keg_file.replace('.keg','_KEGG_original.txt')
    # 初始处理 keg 文件
    command = 'perl /home/colddata/chen/03_transcript/annotation/kegg/kaas_parse_keg_file.pl -i {} -o {}'.format(keg_file, kegg_file)
    os.system(command)
    
    # 使用 python 处理 keg 文件
    # process_keg(keg_file, kegg_file)
    
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
    # Gene_shortname 不能设置空为 NA，设置为空（张老师说的）
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
    
    # 生成 tier2 文件
    tier2_name = key_name + '_KEGG_tier2.txt'
    tier2_df = kegg_clean_df.copy()
    tier2_df.drop(columns=['Pathway', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number'], inplace=True)
    tier2_df['Level2'] = tier2_df['Level2'].str.replace('\\','').str.replace(' / ', '_').str.replace('/', '_').str.replace(', ', '').str.replace(',', '_')
    tier2_df.drop_duplicates(keep='first', inplace=True)
    tier2_df.to_csv(tier2_name, sep='\t', index=False, header=None)
    
    # tier3 文件，后来改名 KEGG.txt 了
    # geneid \t ko pathway
    # no header
    tier3_name = key_name + '_KEGG.txt'
    tier3_df = kegg_clean_df.copy()
    tier3_df.drop(columns=['Level2', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number'], inplace=True)
    tier3_df['Pathway'] = tier3_df['Pathway'].str.replace('\\','').str.replace(' / ', '_').str.replace('/', '_').str.replace(', ', '').str.replace(',', '_')
    tier3_df.drop_duplicates(keep='first', inplace=True)
    tier3_df.to_csv(tier3_name, sep='\t', index=False, header=None)
    
    # ko03022_basal_transcription_factor.txt (子集) from tier3
    # 比对 ko03000_def.txt 加上定义
    ko03000_df, ko03022_df = ko03000(gene_def_df, tier3_df)
    ko03000_df.to_csv(key_name + '_ko03000_transcription_factors.txt', sep='\t', index=False)
    ko03022_df.to_csv(key_name + '_ko03022_basal_transcription_factor.txt', sep='\t', index=False)
    
    # ko03000 需要增加表达量，生成新文件 (2024_02_22)
    if args.fpkm and args.reads:
        ko03000_fpkm_reads_df = merge_fpkm_reads(args.fpkm, args.reads)
        ko03000_fpkm_reads_df = pd.merge(ko03000_fpkm_reads_df, ko03000_df, on='GeneID', how='inner')
        # ko03000_fpkm_reads_df = ko03000_fpkm_reads_df[ko03000_fpkm_reads_df['GeneID'].isin(ko03000_df['GeneID'])]
        ko03000_fpkm_reads_df.to_csv(key_name + '_ko03000_expression_data_def.txt', sep='\t', index=False)

    
if __name__ == '__main__':
    main()