#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/29 17:19
# Author        : William GoGo
import os
import re
import sys
import pandas as pd
import argparse
import subprocess
import openpyxl
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/Rscript/')
from Rscript import draw_multigroup_heatmap
sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='GO_ID 的一个 list 文件, 列名包含 GO_ID')
    p.add_argument('-g', '--genego', help='swiss 注释解析出来的 gene_go 文件')
    p.add_argument('-f', '--fpkmmatrix', default='fpkm_matrix_filtered.txt',
                   help='fpkm_matrix')
    p.add_argument('-s', '--samples', default='samples_described.txt',
                   help='samples_described.txt 样本描述文件')
    p.add_argument('-o', '--outputdir', default='./', help='所有 heatmap 输出文件夹')
    
    p.add_argument('-e', '--expression-data', help='reads_fpkm_matrix_def.txt')
    
    args = p.parse_args()
    
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
    
    return args


def each_go_gene_expression(target_go_df, gene_go_df, expression_data, output_dir='./'):
    for _, row in target_go_df.iterrows():
        go_id = row['GO_ID']
        go_def = row['GO_def']
        each_go_id_df = gene_go_df[gene_go_df['GO_ID'] == go_id]
        if each_go_id_df.shape[0] < 1:
            logger.warning(f'没有 {go_id} 相关基因')
            continue
        each_go_id_gene_expression_df = pd.merge(each_go_id_df, expression_data, on='GeneID', how='inner')
        each_go_id_gene_expression_df.drop(columns=['GO_ID'], inplace=True)
        # 将任何不适合创建文件的字符（包括空格）变成 _
        go_replace_name = go_id.replace(':', '_') + '_' + re.sub(r'[\\/:*?"<>|\s]', '_', str(go_def))
        go_id_gene_expression_fn = os.path.join(output_dir, f'{go_replace_name}_gene_expression.xlsx')
        write_output_df(each_go_id_gene_expression_df, go_id_gene_expression_fn, index=False)


def each_go_gene_heatmap(go_id_list, gene_go_df, fpkm_matrix_df, samples_df, output_dir='./'):
    samples_df = samples_df[['sample', 'group']]
    samples_list = ['GeneID'] + samples_df['sample'].values.tolist()

    # 对 expression fpkm matrix 文件只保留 fpkm 值
    fpkm_matrix_df['GeneID'] = fpkm_matrix_df['GeneID'].astype(str)
    
    for go_id in go_id_list:
        each_go_id_df = gene_go_df[gene_go_df['GO_ID'] == go_id]
        if each_go_id_df.shape[0] < 3:
            logger.warning(f'{go_id} 相关基因数量小于 3 个，不对此 GO_ID 的相关基因出 Heatmap 图')
            
        logger.info(f'尝试对 {go_id} 的相关基因画 Heatmap 图，数量为 {each_go_id_df.shape[0]}')

        # 添加 fpkm
        each_go_id_gene_fpkm_df = pd.merge(each_go_id_df, fpkm_matrix_df, on='GeneID', how='inner')
        each_go_id_gene_fpkm_df.dropna(how='any', axis=0, inplace=True)
        if each_go_id_gene_fpkm_df.shape[0] <= 2:
            logger.error(f'{go_id} 相关基因表达量小于 2，跳过')
            continue
        
        go_replace_name = go_id.replace(':', '_')
        # 准备 multigroup heatmap 输入文件
        go_id_gene_heatmap_filename = os.path.join(output_dir, f'{go_replace_name}_gene_heatmap.xlsx')
        with pd.ExcelWriter(go_id_gene_heatmap_filename, engine='openpyxl') as writer:
            each_go_id_gene_fpkm_df[samples_list].to_excel(writer, sheet_name='Sheet1', index=False)
            samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
        
        # heatmap 画图
        go_id_heatmap_filename = os.path.join(output_dir, go_id_gene_heatmap_filename.replace('.xlsx', '.jpeg'))
        heatmap_result = draw_multigroup_heatmap(
            go_id_gene_heatmap_filename,
            go_id_heatmap_filename,
            other_args = '--cluster-rows'
        )
        if not heatmap_result:
            logger.error(f'{go_id} draw_multigroup_heatmap 程序失败')
    


def main():
    args = parse_input()
    go_df = load_table(args.input)
    if 'GO_def' not in go_df.columns:
        go_df['GO_def'] = go_df['GO_ID'].str.split('_').str[1]
    go_df['GO_ID'] = go_df['GO_ID'].str.split('_').str[0]
    
    genego_df = load_table(args.genego, header=None, names=['GeneID', 'GO_ID'])
    
    if args.fpkmmatrix:
        samples_df = load_table(args.samples)
        goid_list = go_df['GO_ID'].str.split("_").str[0].tolist()
        fpkm_df = load_table(args.fpkmmatrix)
        each_go_gene_heatmap(goid_list, genego_df, fpkm_df, samples_df, args.outputdir)
    
    if args.expression_data:
        expression_df = load_table(args.expression_data)
        each_go_gene_expression(go_df, genego_df, expression_df, args.outputdir)
    
    logger.success('Done')


if __name__ == '__main__':
    main()