#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/03/27 13:29
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


def kid_optmial(input_df):
    """根据输入的结果文件 DEG DataFrame，计算 up down 数量的选择
    up > down = 1
    up < down = -1
    up = down = 0

    Args:
        input_df (DataFrame): 每个 Dataframe

    Returns:
        dict: 包含每个KO最佳选择的字典，键为KO ID，值为最佳选择
        1表示上升趋势，-1表示下降趋势，0表示无明显趋势。
    """
    id_list = input_df.to_dict(orient='split')['data']
    dic = {}
    result_dic = {}
    for id in id_list:
        dic[id[0]] = {'up': 0, 'down': 0, 'nosig': 0}

    for each_line in id_list:
        id = each_line[0]
        reg = each_line[1]
        if str(reg) == '1':
            dic[id]['up'] += 1
        elif str(reg) == '0':
            dic[id]['nosig'] += 1
        elif str(reg) == '-1':
            dic[id]['down'] += 1

    for each_kid, kid_values in dic.items():
        if kid_values['up'] > kid_values['down']:
            result_dic[each_kid] = 1
        elif kid_values['up'] < kid_values['down']:
            result_dic[each_kid] = -1
        elif kid_values['up'] == kid_values['down']:
            result_dic[each_kid] = 0

    return result_dic


def load_kegg_pathway_map():
    """加载KEGG通路ID和名称的映射关系"""
    kegg_db_file = '/home/colddata/qinqiang/script/kns_annotation/scripts/kegg_db.txt'
    kegg_def_df = load_table(kegg_db_file, usecols=[0], dtype=str)
    kegg_def_df.drop_duplicates(inplace=True)
    kegg_def_df[['KEGG_pathway_ID', 'KEGG_pathway_def']] = kegg_def_df['KEGG_Pathway'].str.split(':', n=1, expand=True)
    kegg_def_df['KEGG_pathway_ID'] = kegg_def_df['KEGG_pathway_ID'].str.strip()
    kegg_def_df['KEGG_pathway_def'] = kegg_def_df['KEGG_pathway_def'].str.strip()
    # 替换通路名称中的空格和特殊字符为下划线
    kegg_def_df['KEGG_pathway_def'] = kegg_def_df['KEGG_pathway_def'].str.replace(r'[^\w]+', '_', regex=True)
    kegg_def_df.drop(columns=['KEGG_Pathway'], inplace=True)
    return kegg_def_df


def up_down_idlist_geneid2kid(down_df, up_df, output_prefix, kegg_clean_file, kegg_pathway_map_df):
    """读取上下调基因ID列表，把 regulation 变为 1,0,-1 对应 up,nosig,down
    生成的文件只包括 K_ID 和 regulation
    例如: K12345    0
    一个 K_ID 对应多个 GeneID, 有可能 up,down,nosig 都有
    对 K_ID 去重，哪个 regulation 多留下哪个，一样多的变为 0
    并生成一个中间文件，方便检查，同时返回所有相关的通路ID

    Args:
        down_df (DataFrame): down id list 列名为 GeneID
        up_df (_type_): up id list 列名为 GeneID
        output_prefix (_type_): 输出前缀
        kegg_clean (_type_): KEGG_clean 文件
        kegg_pathway_map_df (DataFrame): KEGG通路ID和名称映射表
    
    Returns:
        DataFrame: 所有相关的KEGG通路ID DataFrame
    """

    kegg_clean_df = load_table(kegg_clean_file, dtype=str, usecols=[0, 1, 4], names=['GeneID', 'KEGG_pathway_ID', 'KEGG_ID'])
    kegg_clean_df['KEGG_pathway_ID'] = kegg_clean_df['KEGG_pathway_ID'].str.split(':').str[0]
    
    geneid_kid_df = kegg_clean_df[['GeneID', 'KEGG_ID']].copy()
    k_pathway_id_kid_df = kegg_clean_df[['KEGG_pathway_ID', 'KEGG_ID']].copy()
    
    down_df = pd.merge(left=down_df, right=geneid_kid_df, on='GeneID', how='left')
    down_df = down_df.drop_duplicates(subset='GeneID', keep='first')
    down_df = down_df.dropna().drop(columns=['GeneID'])
    down_df['regulation'] = -1

    up_df = pd.merge(left=up_df, right=geneid_kid_df, on='GeneID', how='left')
    up_df = up_df.drop_duplicates(subset='GeneID', keep='first')
    up_df = up_df.dropna().drop(columns=['GeneID'])
    up_df['regulation'] = 1
    
    df = pd.concat([down_df, up_df])
    result_dic = kid_optmial(df)
    
    regulation_filename = output_prefix + '_regulation.txt'
    with open(regulation_filename, 'w') as f:
        f.write(f'KEGG_ID\tregulation\n')
        for each_kid, kid_values in result_dic.items():
            f.write(f'{each_kid}\t{kid_values}\n')
    
    regulation_df = load_table(regulation_filename, dtype=str)
    pre_passed_path_df = pd.merge(regulation_df, k_pathway_id_kid_df, on='KEGG_ID', how='left')
    pre_passed_path_df.drop(columns=['KEGG_ID'], inplace=True)
    pre_passed_path_df = pre_passed_path_df.drop_duplicates(subset='KEGG_pathway_ID', keep='first')
    pre_passed_path_df = pre_passed_path_df.dropna()
    
    # 合并通路名称
    pre_passed_path_df = pd.merge(pre_passed_path_df, kegg_pathway_map_df, on='KEGG_pathway_ID', how='left')
    # 确保列顺序正确：第一列KEGG_pathway_ID，第二列KEGG_pathway_def
    pre_passed_path_df = pre_passed_path_df[['KEGG_pathway_ID', 'KEGG_pathway_def']]
    pre_passed_path_df = pre_passed_path_df.dropna()
    
    write_output_df(pre_passed_path_df, output_prefix + '_all_passed_path.txt', index=False, header=False)
    return pre_passed_path_df


def parse_input():
    parser = argparse.ArgumentParser(description='KEGG Pathview 分析脚本，输入上下调基因ID列表生成pathview通路图')
    parser.add_argument('-3', '--group3', dest="group3_file", type=str, required=True, help='[必须]group3 上调基因ID列表文件，基因ID列为第一列')
    parser.add_argument('-7', '--group7', dest="group7_file", type=str, required=True, help='[必须]group7 下调基因ID列表文件，基因ID列为第一列')
    parser.add_argument('-k', '--kegg-clean', dest="kegg_clean", type=str, required=True, help='[必须]输入 kegg 注释出来的 KEGG_clean 文件')
    parser.add_argument('-o', '--output-dir', dest='output_dir', type=str, default='./', help='输出目录, 默认当前目录')
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    return args


def main():
    args = parse_input()
    
    logger.info("===== 开始 KEGG Pathview 分析 =====")
    logger.info(f"上调基因文件: {args.group3_file}")
    logger.info(f"下调基因文件: {args.group7_file}")
    logger.info(f"KEGG_clean文件: {args.kegg_clean}")
    logger.info(f"输出目录: {args.output_dir}")
    
    # 读取基因ID列表
    group3_df = load_table(args.group3_file, dtype=str, usecols=[0], names=['GeneID'])
    group7_df = load_table(args.group7_file, dtype=str, usecols=[0], names=['GeneID'])
    
    logger.info(f"读取到上调基因数量: {group3_df.shape[0]}")
    logger.info(f"读取到下调基因数量: {group7_df.shape[0]}")
    
    # 加载KEGG通路映射表
    logger.info("正在加载KEGG通路名称映射...")
    kegg_pathway_map_df = load_kegg_pathway_map()
    
    # 输出文件前缀
    output_prefix = os.path.join(args.output_dir, 'pathview_analysis')
    
    # 生成regulation文件和通路列表
    logger.info("正在生成regulation文件和通路列表...")
    all_pathway_df = up_down_idlist_geneid2kid(group7_df, group3_df, output_prefix, args.kegg_clean, kegg_pathway_map_df)
    
    logger.info(f"共找到相关KEGG通路数量: {all_pathway_df.shape[0]}")
    if all_pathway_df.shape[0] == 0:
        logger.error("未找到任何相关的KEGG通路，请检查输入文件是否正确！")
        sys.exit(1)
    
    # 调用pathview画图
    logger.info("正在生成pathview通路图...")
    draw_pathview(
        f"{output_prefix}_regulation.txt",
        f"{output_prefix}_all_passed_path.txt"
    )
    
    # 移动生成的png文件到输出目录
    logger.info("正在整理输出文件...")
    os.system(f"mv *.png {args.output_dir} 2>/dev/null")
    
    logger.success("===== KEGG Pathview 分析完成 =====")
    logger.info(f"通路图已保存到: {args.output_dir}")
    logger.info(f"regulation文件: {output_prefix}_regulation.txt")
    logger.info(f"通路列表文件: {output_prefix}_all_passed_path.txt")
    

if __name__ == '__main__':
    main()