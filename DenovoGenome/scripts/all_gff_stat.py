#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/07/01 14:23
# Author        : William GoGo
# Description   : 合并 braker.gff3 rrna.gff3 trna.gff3 还有 ncRNA.gff3 并计算统计信息
import os, sys
import argparse
import pandas as pd
from loguru import logger
import re

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(description="合并 braker.gff3 rRNA.gff3 tRNA.gff3 还有 ncRNA.gff3 并计算统计信息")
    parser.add_argument("-b", "--braker", type=str, required=True, help="braker.gff3 文件")
    parser.add_argument("-r", "--rrna", type=str, required=True, help="rrna.gff3 文件")
    parser.add_argument("-n", "--ncrna", type=str, required=True, help="ncrna.gff3 文件")
    parser.add_argument("-t", "--trna", type=str, required=True, help="trna.gff3 文件")
    parser.add_argument("-p", "--prefix", type=str, required=True, help="输出文件前缀，可以加上目录")
    return parser.parse_args()


def get_braker_ends_id(braker_df):
    """
    从 braker_df 中筛选第三列为 'gene' 的行，提取第九列中 ID
    按从大到小排序后返回第一个值
    """
    gene_rows = braker_df[braker_df['type'] == 'gene']
    
    id_list = []
    for _, row in gene_rows.iterrows():
        attributes = row['attributes']
        match = re.search(r'ID=([^;]+)', attributes)
        if match:
            id_list.append(match.group(1))
        else:
            logger.warning(f'{row} 行不包括 Gene ID')
    
    id_list.sort(reverse=True)
    print(id_list[:5])
    if id_list:
        logger.info(f'braker lastest ID {id_list[0]}')
        return id_list[0]
    else:
        return None


def cmsearch_gff_results_add_ID(rrna_df, ncrna_df, latest_ID):
    prefix, num = latest_ID.rsplit('_', 1)
    
    for idx in rrna_df.index:
        num = str(int(num) + 1)
        current_attrs = rrna_df.loc[idx, 'attributes']
        rrna_df.loc[idx, 'attributes'] = f'ID={prefix}_{num};{current_attrs}'
    
    for idx in ncrna_df.index:
        num = str(int(num) + 1)
        current_attrs = ncrna_df.loc[idx, 'attributes']
        ncrna_df.loc[idx, 'attributes'] = f'ID={prefix}_{num};{current_attrs}'
    
    return f'{prefix}_{num}'


# trna 不改了
# def transcan_gff_results_add_ID(trna_df, latest_ID):
#     prefix, num = latest_ID.rsplit('_', 1)
    
#     for idx in trna_df.index:
#         current_attrs = trna_df.loc[idx, 'attributes']
#         if trna_df.loc[idx, 'type'] == 'gene':
#             num = str(int(num) + 1)
#             new_attrs = re.sub(r'ID=[^;]+', f'ID={prefix}_{num}', current_attrs)
#             trna_df.loc[idx, 'attributes'] = new_attrs
            
#         elif trna_df.loc[idx, 'type'] == 'tRNA':
#             new_attrs = re.sub(r'ID=[^\.]+', f'ID={prefix}_{num}', current_attrs)
#             new_attrs = re.sub(r'Parent=[^;]+', f'Parent={prefix}_{num}', new_attrs)
#             trna_df.loc[idx, 'attributes'] = new_attrs
            
#         elif trna_df.loc[idx, 'attributes'] == 'exon':
#             new_attrs = re.sub(r'ID=[^\.]+', f'ID={prefix}_{num}', current_attrs)
#             new_attrs = re.sub(r'Parent=[^\.]+', f'Parent={prefix}_{num}', new_attrs)
#             trna_df.loc[idx, 'attributes'] = new_attrs

def main():
    args = parse_input()
    braker_file, rrna_file, trna_file, ncrna_file = args.braker, args.rrna, args.trna, args.ncrna
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    braker_df = load_table(braker_file, header=None, names=gff3_columns)
    rrna_df = load_table(rrna_file, header=None, names=gff3_columns)
    trna_df = load_table(trna_file, header=None, names=gff3_columns)
    ncrna_df = load_table(ncrna_file, header=None, names=gff3_columns)
    
    latest_ID = get_braker_ends_id(braker_df)
    if latest_ID:
        latest_ID = cmsearch_gff_results_add_ID(rrna_df, ncrna_df, latest_ID)
        # transcan_gff_results_add_ID(trna_df, latest_ID)  # tRNA 不改了
    
    all_df = pd.concat([braker_df, rrna_df, trna_df, ncrna_df], ignore_index=True)
    
    write_output_df(all_df, f"{args.prefix}_all.gff3", index=False, header=True)
    
    # 统计每种类型的 gene id 个数，出一个统计表
    braker_id_counts = braker_df[braker_df['type'] == 'gene'].shape[0]
    rrna_id_counts = rrna_df.shape[0]
    trna_id_counts = trna_df[trna_df['type'] == 'gene'].shape[0]
    ncrna_id_counts = ncrna_df.shape[0]
    total_genes = int(braker_id_counts) + int(rrna_id_counts) + int(trna_id_counts) + int(ncrna_id_counts)
    
    
    stat_results = pd.DataFrame({
        'Gene type': ['protein_coding_gene', 'rRNA_gene', 'tRNA_gene', 'ncRNA_gene', 'Total genes'],
        'Gene counts': [braker_id_counts, rrna_id_counts, trna_id_counts, ncrna_id_counts, total_genes]
    })
    
    write_output_df(stat_results, f"{args.prefix}_gff_gene_type_counts.txt", index=False)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()