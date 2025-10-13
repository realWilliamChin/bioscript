#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/03/28 15:38
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
from loguru import logger


sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
from Rscript import draw_pathview
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df

if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    logger.critical("当前 Python 版本为:", sys.version)
    sys.exit(1)


def parse_input():
    p = argparse.ArgumentParser()
    # 其他参数
    p.add_argument('-i', '--ko-list', dest="ko_list", type=str, required=True,
            help='[必须]ko_list file, 必须要有列名，KEGG_Pathway_ID')
    p.add_argument('-s', '--samples-info', dest='samples_info', type=str, required=True,
            help='[必须]输入样品信息文件，如果有则添加到输出文件中')
    p.add_argument('--fpkm-reads-merged', dest='fpkm_reads_merged', type=str,
                           help='GO 输入 fpkm_reads_merged.txt 文件')
    p.add_argument('--kegg-clean', dest="kegg_clean", required=True, type=str,
            help='[必须]输入 kegg 注释出来的 KEGG_clean 文件')
    p.add_argument('--enrich-dir', dest='enrich_dir', type=str,
            help='[必须]输入富集分析的文件夹')
    p.add_argument('--id-files', type=str, nargs='+',
        help='[必须]输入 module ids 文件（可多个，空格分隔）')
    p.add_argument('--ids-header', dest='ids_header', action='store_true',
                        help='所有 id 文件是否有 header(GeneID)')
    p.add_argument('--output-dir', dest='output_dir', type=str, default='./',
            help='输出目录, 默认当前目录')

    args = p.parse_args()
    
    if args.ids_header:
        args.ids_header = 0
    else:
        args.ids_header = None
    
    if isinstance(args.id_files, str):
        args.id_files = [args.id_files]

    return args


def genes_kegg_analysis(ids_header, target_ko_df, kegg_clean_df, enrich_dir, id_files, output_dir):
    kegg_clean_df['KEGG_Pathway_ID'] = kegg_clean_df['KEGG_Pathway_ID'].str.split(':').str[0]
    geneid_kid_df = kegg_clean_df[['GeneID', 'KEGG_ID']].copy()
    kegg_output_summary_df_list = []
    for module_id_file in id_files:
        module_name = os.path.basename(module_id_file).split('_')[0]
        module_df = load_table(
            module_id_file,
            dtype=str,
            usecols=[0],
            header=ids_header,
            names=['GeneID']
        )

        # 创建输出文件夹
        module_output_dir = os.path.join(output_dir, f'{module_name}')
        kegg_pathway_graph_dir = os.path.join(module_output_dir, 'KEGG_pathway_graph')
        
        for each_dir in [module_output_dir, kegg_pathway_graph_dir]:
            if not os.path.exists(each_dir):
                os.mkdir(each_dir)
            else:
                logger.warning(f"{each_dir} 文件夹已存在，输出将覆盖原文件")
        
        # 从 enrich 文件中提取 kegg 相关数据
        kegg_enrich_df = load_table(os.path.join(enrich_dir, f'{module_name}_EnrichmentKEGG.xlsx'))
        kegg_enrich_df = kegg_enrich_df.drop(columns=['Ontology'])
        # 过滤提取包含 kolist 的行
        kegg_enrich_df = pd.merge(target_ko_df, kegg_enrich_df, left_on='KEGG_Pathway_ID', right_on='ID', how='inner')
        kegg_enrich_df = kegg_enrich_df.drop(columns=['KEGG_Pathway_ID'])
        write_output_df(
            kegg_enrich_df,
            os.path.join(module_output_dir, f'{module_name}_Enrich.xlsx'),
            index=False
        )
        kegg_enrich_df['Module_ID'] = module_name
        kegg_output_summary_df_list.append(kegg_enrich_df)
        
        logger.info(f"====正在处理 {module_name}====")
              
        # R pathview 画图准备文件 regulation
        logger.info(f"正在准备 {module_name} 的 regulation 文件")
        
        
        regulation_df = pd.merge(module_df, geneid_kid_df, on='GeneID', how='left')
        regulation_df.drop_duplicates(subset='GeneID', keep='first', inplace=True)
        regulation_df = regulation_df.dropna().drop(columns=['GeneID'])
        regulation_df['regulation'] = 1
        
        write_output_df(
            regulation_df,
            os.path.join(module_output_dir, f"{module_name}_regulation.txt"),
            index=False
        )
        
        # R pathview 画图准备文件 passed path
        logger.info(f"正在筛选 {module_name} 的 passed_path.txt")
        
        # get passed_path
        passed_path_file = '/home/colddata/qinqiang/script/Analysis/pathview/passed_path.txt'
        passed_path_df = load_table(passed_path_file, names=['KEGG_Pathway_ID', 'Ko_Def'], dtype=str)
        ko_def_df = passed_path_df[passed_path_df['KEGG_Pathway_ID'].isin(target_ko_df['KEGG_Pathway_ID'])]
        ko_def_df = ko_def_df.drop_duplicates(subset=['KEGG_Pathway_ID'])
        write_output_df(
            ko_def_df,
            os.path.join(module_output_dir, f"{module_name}_ko_passed_path.txt"),
            index=False,
            header=False
        )
        
        # R pathview 画图
        logger.info(f"正在画 {module_name} 的 pathview 图")
        # 这两个文件用完需要删掉
        draw_pathview(
            os.path.join(module_output_dir, f"{module_name}_regulation.txt"),
            os.path.join(module_output_dir, f"{module_name}_ko_passed_path.txt")
        )
        
        # 对每个比较组的文件进行整理
        os.system(f"mv *.png {kegg_pathway_graph_dir}")


    summary_df = pd.concat(kegg_output_summary_df_list)
    summary_df_col_lst = summary_df.columns.tolist()
    fixed_cols = ['ID', 'Description', 'SubOntology', 'Ontology', 'pvalue', 'Module_ID']
    summary_df_col_lst = fixed_cols + [x for x in summary_df_col_lst if x not in fixed_cols]
    summary_df = summary_df[summary_df_col_lst]
    write_output_df(
        summary_df,
        os.path.join(output_dir, 'Target_KEGG_analysis_summary.xlsx'),
        index=False
    )


# 提取每个 ko 的表达量
def each_groups_ko_expression_def(ids_header, ko_list, kegg_clean_df, fpkm_reads_merged_df, id_files, output_dir):
    kegg_pathway_df = kegg_clean_df[['GeneID', 'KEGG_Pathway_ID']]
    for module_id_file in id_files:
        module_name = os.path.basename(module_id_file).split('_')[0]
        module_df = load_table(
            module_id_file,
            dtype=str,
            usecols=[0],
            header=ids_header,
            names=['GeneID']
        )
        module_expression_data_dir = os.path.join(output_dir, module_name, 'KO_expression_data')
        if not os.path.exists(module_expression_data_dir):
            os.makedirs(module_expression_data_dir)
        for ko_id in ko_list:
            ko_num_data_def_file = os.path.join(module_expression_data_dir, f"{module_name}_{ko_id}_data.csv")
            ko_num_data_def_df = module_df[module_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['KEGG_Pathway_ID'].str.contains(ko_id)]['GeneID'])]
            ko_num_data_def_df = pd.merge(ko_num_data_def_df, fpkm_reads_merged_df, on='GeneID', how='left')
            if ko_num_data_def_df.shape[0] == 0:
                logger.warning(f"{module_df} 的 {ko_id} 中相关的基因表达量表为空")
                continue
            write_output_df(ko_num_data_def_df, ko_num_data_def_file, index=False)


def main():
    args = parse_input()
    target_ko_df = load_table(args.ko_list)
    ko_list = target_ko_df['KEGG_Pathway_ID'].values.tolist()
    
    kegg_clean_df = load_table(
        args.kegg_clean,
        usecols=[0, 1, 4],
        names=['GeneID', 'KEGG_Pathway_ID', 'KEGG_ID'],
        dtype={'GeneID': str}
    )
    fpkm_reads_merged_df = load_table(args.fpkm_reads_merged)
    
    genes_kegg_analysis(args.ids_header, target_ko_df, kegg_clean_df,
                            args.enrich_dir, args.id_files, args.output_dir)
    
    each_groups_ko_expression_def(args.ids_header, ko_list, kegg_clean_df,
                                  fpkm_reads_merged_df, args.id_files, args.output_dir)
    
    logger.success('Done!')
    

if __name__ == '__main__':
    main()