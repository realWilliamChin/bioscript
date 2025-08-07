#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/29 17:19
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
from Rscript import draw_multigroup_heatmap
from Rscript import draw_twogroup_heatmap
from Rscript import draw_pathview
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df
sys.path.append('/home/colddata/qinqiang/script/Analysis/kegg_analysis')
from each_ko_gene_heatmap import each_ko_gene_heatmap
from get_passed_path import passed_path

if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    logger.critical("当前 Python 版本为:", sys.version)
    sys.exit(1)


def kid_optmial(input_df):
    """根据输入的 MultiDESeq 结果文件 DEG DataFrame，计算 up down 数量的选择
    up > down = 1
    up < down = -1
    up = down = 0

    Args:
        input_df (DataFrame): DEG 的每个 Dataframe

    Returns:
        dict: 包含每个孩子最佳选择的字典，键为孩子ID，值为最佳选择
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


def up_down_idlist_geneid2kid(down_df, up_df, output_prefix, kegg_clean_file):
    """可以读取 de_results_add_def.py 生成的 DEG_data.txt 文件
    把 regulation 变为 1,0,-1 对应 up,nosig,down
    生成的文件只包括 K_ID 和 regulation
    例如: K12345    0
    一个 K_ID 对应多个 GeneID, 有可能 up,down,nosig 都有
    对 K_ID 去重，哪个 regulation 多留下哪个，一样多的变为 0
    并生成一个中间文件，方便检查

    Args:
        down_df (DataFrame): down id list 列名为 GeneID
        up_df (_type_): up id list 列名为 GeneID
        output_prefix (_type_): 输出前缀
        kegg_clean (_type_): KEGG_clean 文件
    """

    kegg_clean_df = load_table(kegg_clean_file, dtype=str, usecols=[0, 1, 4], names=['GeneID', 'KEGG_Pathway_ID', 'KEGG_ID'])
    kegg_clean_df['KEGG_Pathway_ID'] = kegg_clean_df['KEGG_Pathway_ID'].str.split(':').str[0]
    
    geneid_kid_df = kegg_clean_df[['GeneID', 'KEGG_ID']].copy()
    k_pathway_id_kid_df = kegg_clean_df[['KEGG_Pathway_ID', 'KEGG_ID']].copy()
    
    down_df = pd.merge(left=down_df, right=geneid_kid_df, on='GeneID', how='left')
    down_df = down_df.drop_duplicates(subset='GeneID', keep='first')
    down_df = down_df.dropna().drop(columns=['GeneID'])
    down_df['regulation'] = -1

    up_df = pd.merge(left=up_df, right=kegg_clean_df, on='GeneID', how='left')
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
    regulation_df.drop(columns=['regulation'], inplace=True)
    pre_passed_path_df = pd.merge(regulation_df, k_pathway_id_kid_df, on='KEGG_ID', how='left')
    pre_passed_path_df.drop(columns=['KEGG_ID'], inplace=True)
    write_output_df(pre_passed_path_df, output_prefix + '_pre_passed_path.txt', index=False)


def kegg_summary(input_target_df, df_list):
    """kegg summary 汇总文件

    Args:
        input_target_df (str): 输入文件 KEGG_Pathway_ID，Ontology，Subontology ...
        df_list (list): 需要汇总的 dataframe list

    Returns:
        pd.DataFrame: 输出结果 dataframe
    """
    kegg_output_summary_df = pd.concat(df_list)
    kegg_output_summary_df.drop(columns=['Ontology'], inplace=True)
    kegg_output_summary_df.rename(columns={'ID': 'KEGG_Pathway_ID'}, inplace=True)
    if 'Ontology' in input_target_df.columns and 'SubOntology' in input_target_df.columns:
        kegg_output_summary_df = pd.merge(
            input_target_df[['KEGG_Pathway_ID', 'Ontology', 'SubOntology']],
            kegg_output_summary_df,
            how='right',
            on='KEGG_Pathway_ID'
        )
    return kegg_output_summary_df


def ko03000_deg_data_summary(input_dir):
    """每个组间 ko03000 相关目标基因有表达结果汇总

    Args:
        input_dir (str): 结果目录

    Returns:
        pd.DataFrame: 结果 dataframe
    """
    ko03000_df_list = []
    for each_dir in os.listdir(input_dir):
        if '-vs-' not in each_dir:
            continue
        ko03000_data_df = load_table(
            os.path.join(input_dir, each_dir, 'DEG_expression_data', f'{each_dir}_ko03000_DEG_data.xlsx'),
            header = 0,
            usecols = ['GeneID', 'sampleA', 'sampleB', 'pvalue', 'padj', 'regulation', 'FC', 'KEGG_Description']
        )
        ko03000_data_df = ko03000_data_df[ko03000_data_df['regulation'].str.lower() != 'nosignificant']
        ko03000_df_list.append(ko03000_data_df)
    ko03000_df = pd.concat(ko03000_df_list)
    return ko03000_df


def deg_kegg_analysis(ko_list_file, enrich_dir, deg_data_dir, samples_described_df, kegg_pathway_df,
                  kegg_clean_file, treat_color, control_color, output_dir):
    """
    DEG 数据的 KEGG 分析
    """
    input_target_ko_df = load_table(ko_list_file)
    ko_list = input_target_ko_df['KEGG_Pathway_ID'].values.tolist()
    deg_data_list = [x for x in os.listdir(deg_data_dir) if x.endswith('_DEG_data.txt')]
    if len(deg_data_list) == 0:
        logger.error(f'{deg_data_dir} 没有找到 DEG_data.txt 文件')
        sys.exit(1)
    
    kegg_output_summary_df_list = []
    for deg_data_file in deg_data_list:
        # 读取 DEG_data.txt 文件，获取比较组信息
        deg_data_df = load_table(os.path.join(deg_data_dir, deg_data_file), dtype=str)
        treat_group = deg_data_df.iloc[0, 1]
        control_group = deg_data_df.iloc[0, 2]
        compare_info = treat_group + '-vs-' + control_group
        compare_info_dir = os.path.join(output_dir, compare_info)
        
        # 转录组 KEGG 分析输出文件夹创建
        deg_expression_data_dir = os.path.join(compare_info_dir, 'DEG_expression_data')
        kegg_heatmap_dir = os.path.join(compare_info_dir, 'KEGG_heatmap')
        kegg_pathway_graph_dir = os.path.join(compare_info_dir, 'KEGG_pathway_graph')
        
        for each_dir in [compare_info_dir, deg_expression_data_dir, kegg_heatmap_dir, kegg_pathway_graph_dir]:
            if not os.path.exists(each_dir):
                os.mkdir(each_dir)
            else:
                logger.warning(f"{each_dir} 文件夹已存在，输出将覆盖原文件")
        
        # 从 enrich 文件中提取 kegg 相关数据
        up_enrich_df = load_table(os.path.join(enrich_dir, f'{compare_info}_Up_EnrichmentKEGG.xlsx'))
        down_enrich_df = load_table(os.path.join(enrich_dir, f'{compare_info}_Down_EnrichmentKEGG.xlsx'))
        # 过滤提取包含 kolist 的行
        kegg_up_enrich_df = up_enrich_df[up_enrich_df['ID'].isin(ko_list)]
        kegg_down_enrich_df = down_enrich_df[down_enrich_df['ID'].isin(ko_list)]
        write_output_df(
            kegg_up_enrich_df,
            os.path.join(compare_info_dir, f'{compare_info}_Up_Enrich.xlsx'),
            index=False
        )
        write_output_df(
            kegg_down_enrich_df,
            os.path.join(compare_info_dir, f'{compare_info}_Down_Enrich.xlsx'),
            index=False
        )
        
        # 上面 kegg_up/down_enrich_df 输出汇总到一个文件 2025——01——16
        kegg_up_enrich_df.insert(0, 'Group', compare_info)
        kegg_up_enrich_df.insert(1, 'Regulation', 'Up')
        kegg_down_enrich_df.insert(0, 'Group', compare_info)
        kegg_down_enrich_df.insert(1, 'Regulation', 'Down')
        kegg_output_summary_df_list.append(kegg_up_enrich_df)
        kegg_output_summary_df_list.append(kegg_down_enrich_df)
        
        logger.info(f"====正在处理 {compare_info}====")
        up_df = deg_data_df[deg_data_df['regulation'] == 'Up']['GeneID']
        down_df = deg_data_df[deg_data_df['regulation'] == 'Down']['GeneID']

        crt_group_samples = samples_described_df[samples_described_df['group'].isin([treat_group, control_group])]['sample'].to_list()
        fpkm_df = deg_data_df[["GeneID"] + [x + '_FPKM' for x in crt_group_samples]]
        # 把 fpkm_df 的所有列名带 _FPKM 的去掉
        fpkm_df.columns = ["GeneID"] + crt_group_samples
        
        treat_group_samples_name = samples_described_df[samples_described_df['group'].isin([treat_group,])]['sample'].to_list()
        control_group_samples_name = samples_described_df[samples_described_df['group'].isin([control_group,])]['sample'].to_list()
        logger.debug(f"crt_group_samples: {crt_group_samples}")
        
        # 画热图准备文件 Sheet2 samples group color
        heatmap_sheet2_samplegroupcolor_df = pd.DataFrame(columns=['samples', 'group', 'colors'])
        for each_sample in crt_group_samples:
            if each_sample in treat_group_samples_name:
                group = treat_group
                color = treat_color
            else:
                group = control_group
                color = control_color
            row_data = {'samples': each_sample, 'group': group, 'colors': color}
            row_df = pd.DataFrame([row_data])
            heatmap_sheet2_samplegroupcolor_df = pd.concat([heatmap_sheet2_samplegroupcolor_df, row_df], ignore_index=True)

        # 循环每个 ko number 筛选表达量画热图
        for ko_number in ko_list:
            ko_num_fpkm_expr_df = fpkm_df[fpkm_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['KEGG_Pathway'].str.contains(ko_number, na=False)]['GeneID'])].copy()
            if ko_num_fpkm_expr_df.shape[0] == 0:
                logger.warning(f"{compare_info} 的 {ko_number} 中相关的基因表达量表为空")
                continue
            # 去除掉除 GeneID 列之外全为空的行
            ko_num_fpkm_expr_df.dropna(subset=crt_group_samples, how='all', inplace=True)
            ko_num_fpkm_expr_df = ko_num_fpkm_expr_df.loc[~(ko_num_fpkm_expr_df[crt_group_samples] == 0).all(axis=1)]
            
            ko_num_fpkm_expr_df_file = os.path.join(kegg_heatmap_dir, f"{compare_info}_{ko_number}_fpkm_expression.xlsx")
            ko_num_fpkm_expr_pic_name = os.path.join(kegg_heatmap_dir, f"{compare_info}_{ko_number}_heatmap.jpeg")

            with pd.ExcelWriter(ko_num_fpkm_expr_df_file) as writer:
                ko_num_fpkm_expr_df.to_excel(writer, index=False, sheet_name='Sheet1')
                heatmap_sheet2_samplegroupcolor_df.to_excel(writer, index=False, sheet_name='Sheet2')
            if 100 >= ko_num_fpkm_expr_df.shape[0] > 1:
                logger.info(f"正在画 {compare_info} {ko_number} 相关基因表达量的热图")
                draw_twogroup_heatmap(ko_num_fpkm_expr_df_file, ko_num_fpkm_expr_pic_name)
            else:
                logger.warning(f"{compare_info} 的 {ko_number} 中相关的基因表达量表只有一行或超过 100 行，不执行画热图")
        
            # 输出每个 ko 的 deg_data 文件
            ko_num_deg_data_file = os.path.join(deg_expression_data_dir, f"{compare_info}_{ko_number}_DEG_data.xlsx")
            ko_num_deg_data_df = deg_data_df[deg_data_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['KEGG_Pathway'].str.contains(ko_number, na=False)]['GeneID'])]
            write_output_df(ko_num_deg_data_df, ko_num_deg_data_file, index=False)
            
        # R pathview 画图准备文件 regulation
        # print(f"正在准备 {compare_info} 的 regulation 文件")
        up_down_idlist_geneid2kid(down_df, up_df, os.path.join(compare_info_dir, compare_info), kegg_clean_file)
        
        # R pathview 画图准备文件 passed path
        # print(f"正在筛选 {compare_info} 的 passed_path.txt")
        compare_passed_path_filename = os.path.join(compare_info_dir, f"{compare_info}_ko_passed_path.txt")
        # pre_passed_path_name = f'{compare_info}_pre_passed_path.txt'
        # passed_path(pre_passed_path_name, compare_passed_path_filename)
        passed_path(ko_list_file, compare_passed_path_filename)
        
        # R pathview 画图
        logger.info(f"正在画 {compare_info} 的 pathview 图")
        # 这两个文件用完需要删掉
        draw_pathview(
            os.path.join(compare_info_dir, f"{compare_info}_regulation.txt"),
            os.path.join(compare_info_dir, f"{compare_info}_ko_passed_path.txt")
        )
        
        # 对每个比较组的文件进行整理
        os.system(f"mv *.png {kegg_pathway_graph_dir}")

    
    kegg_summary_df = kegg_summary(input_target_ko_df, kegg_output_summary_df_list)
    write_output_df(kegg_summary_df, os.path.join(output_dir, 'Target_KEGG_analysis_summary.xlsx'), index=False)
    
    ko03000_summary_df = ko03000_deg_data_summary(output_dir)
    write_output_df(ko03000_summary_df, os.path.join(output_dir, 'ko03000_DEG_data_summary.xlsx'), index=False)


def parse_input():
    argparser = argparse.ArgumentParser()
    # 其他参数
    argparser.add_argument('-i', dest="target_ko_file", type=str, required=True,
            help='[必须]输入目标 ko 文件，必须要有列名，KEGG_Pathway_ID')
    argparser.add_argument('-s', '--samplesinfo', type=str,
            help='[必须]输入样品信息文件，如果有则添加到输出文件中')
    argparser.add_argument('--kegg-clean', dest="kegg_clean", required=True, type=str,
            help='[必须]输入 kegg 注释出来的 KEGG_clean 文件')
    argparser.add_argument('--enrich-dir', dest='enrich_dir', type=str,
            help='[必须]输入富集分析的文件夹')
    argparser.add_argument('--output-dir', dest='output_dir', type=str, default='./',
            help='输出目录, 默认当前目录')
    argparser.add_argument('--deg-data-dir', dest='deg_data_dir', type=str, 
                        help='[必须]输入 DEG_data.txt 文件')

    argparser.add_argument('--treat-color', dest="treat_color", type=str, default="blue",
                        help='[可选]画 heatmap 时指定 Treat Group 的颜色，默认 blue')
    argparser.add_argument('--control-color', dest="control_color", type=str, default="red",
                        help='[可选]画 heatmap 时指定 Control Group 的颜色, 默认 red')
    
    args = argparser.parse_args()

    return args


def main():
    args = parse_input()
    
    samples_described_df = load_table(args.samplesinfo, usecols=[0, 1], dtype=str)
    kegg_pathway_df = load_table(args.kegg_clean, usecols=[0, 1], names=['GeneID', 'KEGG_Pathway'], dtype=str)
    
    deg_kegg_analysis(args.target_ko_file, args.enrich_dir, args.deg_data_dir, samples_described_df, kegg_pathway_df,
                args.kegg_clean, args.treat_color, args.control_color, args.output_dir)
    
    logger.success('Done!')
    

if __name__ == '__main__':
    main()