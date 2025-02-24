#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/29 17:19
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse
import subprocess
import openpyxl
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
from Rscript import draw_multigroup_heatmap


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='GO_ID 的一个 list 文件, 列名包含 GO_ID')
    p.add_argument('-g', '--genego', help='swiss 注释解析出来的 gene_go 文件')
    p.add_argument('-f', '--fpkmmatrix', default='fpkm_matrix_filtered.txt',
                   help='fpkm_matrix')
    p.add_argument('-s', '--samples', default='samples_described.txt',
                   help='samples_described.txt 样本描述文件')
    p.add_argument('-o', '--outputdir', default='./', help='所有 heatmap 输出文件夹')
    
    args = p.parse_args()
    
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
    
    return args


def each_go_gene_heatmap(go_id_list, gene_go_df, fpkm_matrix_df, samples_df, output_dir='./'):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
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
        each_go_id_gene_fpkm_df = pd.merge(each_go_id_df, fpkm_matrix_df, on='GeneID', how='left')
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
    go_df = pd.read_csv(args.input, sep='\t')
    goid_list = go_df['GO_ID'].str.split('_').str[0].tolist()
    genego_df = pd.read_csv(args.genego, sep='\t', header=None, names=['GeneID', 'GO_ID'])
    fpkm_df = pd.read_csv(args.fpkmmatrix, sep='\t')
    samples_df = pd.read_csv(args.samples, sep='\t')
    
    each_go_gene_heatmap(goid_list, genego_df, fpkm_df, samples_df, args.outputdir)
    
    logger.success('Done')


if __name__ == '__main__':
    main()