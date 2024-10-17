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

from each_ko_gene_heatmap import each_ko_gene_heatmap
from get_passed_path import passed_path

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
from Rscript import draw_multigroup_heatmap
from Rscript import draw_twogroup_heatmap
from Rscript import draw_pathview
# from Rscript import anova_analysis

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


def up_down_idlist_geneid2kid(down_df, up_df, output_prefix, kegg_clean):
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

    kegg_clean_df = pd.read_csv(kegg_clean, sep='\t', usecols=[0, 4], names=['GeneID', 'KEGG_ID'], dtype=str)
    
    down_df = pd.merge(left=down_df, right=kegg_clean_df, on='GeneID', how='left')
    down_df = down_df.drop_duplicates(subset='GeneID', keep='first')
    down_df = down_df.dropna().drop(columns=['GeneID'])
    down_df['regulation'] = -1

    up_df = pd.merge(left=up_df, right=kegg_clean_df, on='GeneID', how='left')
    up_df = up_df.drop_duplicates(subset='GeneID', keep='first')
    up_df = up_df.dropna().drop(columns=['GeneID'])
    up_df['regulation'] = 1
    
    df = pd.concat([down_df, up_df])
    result_dic = kid_optmial(df)
    
    with open(output_prefix + '_regulation.txt', 'w') as f:
        f.write(f'KEGG_ID\tregulation\n')
        for each_kid, kid_values in result_dic.items():
            f.write(f'{each_kid}\t{kid_values}\n')


# def add_def(ko_file, kegg_gene_def_file, kegg_pathway_file, output_prefix, basicinfo=None):
def other(args):
    (
        ko_file, kegg_gene_def_file, kegg_pathway_file, 
        output_prefix, basicinfo
    ) = (
        args.ko_list, args.kegg_genedef, args.kegg_pathway,
        args.output_prefix, args.basicinfo
    )
    
    ko_df = pd.read_csv(ko_file, sep='\t', dtype=str)
    kegg_gene_def = pd.read_csv(kegg_gene_def_file, sep='\t', dtype=str)
    kegg_pathway_df = pd.read_csv(kegg_pathway_file, sep='\t', names=['GeneID', 'Ko'], dtype=str)
    kegg_pathway_df['KEGG_ID'] = kegg_pathway_df['Ko'].str.split(":").str[0]
    
    ko_gene_df = pd.merge(ko_df, kegg_pathway_df, how='inner', on='KEGG_ID')
    ko_gene_df = pd.merge(ko_gene_df, kegg_gene_def, how='left', on='GeneID')
    
    ko_gene_df = ko_gene_df.drop(columns=['KEGG_ID'])
    ko_gene_df.fillna('NA', inplace=True)
    
    if basicinfo is not None:
        basicinfo_df = pd.read_csv(basicinfo, sep='\t', usecols=[0,1,2,3], dtype=str)
        basicinfo_df_columns = basicinfo_df.columns.tolist()
        ko_gene_df = pd.merge(ko_gene_df, basicinfo_df, how='left', on='GeneID')
        output_columns_list = basicinfo_df_columns + ['Ko', 'KEGG_ID', 'Gene_shortname', 'EC_number', 'KEGG_def']
    else:
        output_columns_list = ['GeneID', 'Ko', 'KEGG_ID', 'Gene_shortname', 'EC_number', 'KEGG_def']
        
    ko_gene_df = ko_gene_df[output_columns_list]
    ko_gene_df.to_csv(output_prefix + '_ko_geneid_def.txt', sep='\t', index=False)
    
    ko_gene_df['regulation'] = 1
    ko_gene_df[['KEGG_ID', 'regulation']].to_csv(output_prefix + '_keggid_regulation.txt', sep='\t', index=False)


def transcriptome(args):
    """针对转录组的 KEGG 分析

    Args:
        args (_type_): 所有的输入参数
    """
    # 循环每个 ko number 筛选表达量画热图
    ko_df = pd.read_csv(args.ko_list, sep='\t')
    ko_list = ko_df['KEGG_Pathway_ID'].values.tolist()
    samples_described_df = pd.read_csv(args.samplesinfo, sep='\t', usecols=[0, 1], dtype=str)
    kegg_pathway_df = pd.read_csv(args.kegg_clean, sep='\t', names=['GeneID', 'KEGG_Pathway'], usecols=[0, 1], dtype=str)
    fpkm_matrix_df = pd.read_csv(args.expression, sep='\t')
    if args.kogene_heatmap:
        # transcriptome_each_ko_gene_heatmap(ko_list, args.kegg_clean, args.expression, args.samplesinfo)
        each_ko_gene_heatmap(ko_list, kegg_pathway_df, fpkm_matrix_df, samples_described_df, output_dir='Each_KO_gene_heatmap')
    
    deg_data_list = os.listdir(args.deg_data_dir)
    
    
    for deg_data_file in deg_data_list:
        
        if not deg_data_file.endswith('_DEG_data.txt'):
            continue
        deg_data_file = os.path.join(args.deg_data_dir, deg_data_file)
        deg_data_df = pd.read_csv(deg_data_file, sep='\t', dtype=str)
        
        treat_group = deg_data_df.iloc[0, 1]
        control_group = deg_data_df.iloc[0, 2]
        compare_info = treat_group + '_vs_' + control_group
        compare_info_dir = compare_info
        
        deg_expression_data_dir = os.path.join(compare_info_dir, 'DEG_expression_data')
        kegg_heatmap_dir = os.path.join(compare_info_dir, 'KEGG_heatmap')
        kegg_pathway_graph_dir = os.path.join(compare_info_dir, 'KEGG_pathway_graph')
        
        for each_dir in [compare_info_dir, deg_expression_data_dir, kegg_heatmap_dir, kegg_pathway_graph_dir]:
            if not os.path.exists(each_dir):
                os.mkdir(each_dir)
            else:
                logger.warning(f"{each_dir} 文件夹已存在，输出将覆盖原文件")
        
        # deg_data_df['GeneID'].isin(gene_go_df[gene_go_df['GO_ID'].str.contains(go_id)]['GeneID'])
        # 从 enrich 文件中提取 kegg 相关数据
        enrich_data_up_file = os.path.join(args.enrich_dir, f'{compare_info}_Up_EnrichmentKEGG.xlsx')
        enrich_data_Down_file = os.path.join(args.enrich_dir, f'{compare_info}_Down_EnrichmentKEGG.xlsx')
        up_enrich_df = pd.read_excel(enrich_data_up_file, dtype=str, engine='openpyxl')
        down_enrich_df = pd.read_excel(enrich_data_Down_file, dtype=str, engine='openpyxl')
        # 提取包含 kolist 的行
        kegg_up_enrich_df = up_enrich_df[up_enrich_df['ID'].isin(ko_list)]
        kegg_down_enrich_df = down_enrich_df[down_enrich_df['ID'].isin(ko_list)]
        # 输出到 compare_info_dir
        kegg_upenrich_outfile = os.path.join(compare_info_dir, f'{compare_info}_Up_Enrich.xlsx')
        kegg_downenrich_outfile = os.path.join(compare_info_dir, f'{compare_info}_Down_Enrich.xlsx')
        kegg_up_enrich_df.to_excel(kegg_upenrich_outfile, index=False)
        kegg_down_enrich_df.to_excel(kegg_downenrich_outfile, index=False)
        
        
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
        
        heatmap_sheet2_samplegroupcolor_df = pd.DataFrame(columns=['samples', 'group', 'colors'])
        for each_sample in crt_group_samples:
            if each_sample in treat_group_samples_name:
                group = treat_group
                color = args.treat_color
            else:
                group = control_group
                color = args.control_color
            row_data = {'samples': each_sample, 'group': group, 'colors': color}
            row_df = pd.DataFrame([row_data])
            heatmap_sheet2_samplegroupcolor_df = pd.concat([heatmap_sheet2_samplegroupcolor_df, row_df], ignore_index=True)
            # heatmap_sheet2_df = heatmap_sheet2_df.append({'samples': each_sample, 'group': group, 'colors': color}, ignore_index=True)

        for ko_number in ko_list:
            ko_num_fpkm_expr_df = fpkm_df[fpkm_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['KEGG_Pathway'].str.contains(ko_number)]['GeneID'])].copy()
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
            if ko_num_fpkm_expr_df.shape[0] > 1:
                logger.info(f"正在画 {compare_info} {ko_number} 相关基因表达量的热图")
                draw_twogroup_heatmap(ko_num_fpkm_expr_df_file, ko_num_fpkm_expr_pic_name)
            else:
                logger.warning(f"{compare_info} 的 {ko_number} 中相关的基因表达量表只有一行，无法画热图")
            
            # 输出每个 ko 的 deg_data 文件
            ko_num_deg_data_file = os.path.join(deg_expression_data_dir, f"{compare_info}_{ko_number}_DEG_data.txt")
            ko_num_deg_data_df = deg_data_df[deg_data_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['KEGG_Pathway'].str.contains(ko_number)]['GeneID'])]
            ko_num_deg_data_df.to_csv(ko_num_deg_data_file, sep='\t', index=False)
            
        # R pathview 画图准备文件 regulation
        # print(f"正在准备 {compare_info} 的 regulation 文件")
        up_down_idlist_geneid2kid(down_df, up_df, compare_info, args.kegg_clean)
        
        # R pathview 画图准备文件 passed path
        # print(f"正在筛选 {compare_info} 的 passed_path.txt")
        compare_passed_path_filename = f'{compare_info}_ko_passed_path.txt'
        passed_path(args.ko_list, compare_passed_path_filename)
        
        # R pathview 画图
        logger.info(f"正在画 {compare_info} 的 pathview 图")
        # 这两个文件用完需要删掉
        draw_pathview(f"{compare_info}_regulation.txt", f"{compare_info}_ko_passed_path.txt")
        
        # 对每个比较组的文件进行整理
        # os.system(f"mv {compare_info}_*.xlsx {compare_info}")
        os.system(f"mv *.png {kegg_pathway_graph_dir}")


def get_ko_expression(ko, kegg_clean_file, expression_file):
    """对每个 ko 相关的基因生成一个 ko 相关基因的表达量表

    Args:
        ko (_type_): ko_list 文件
        kegg_clean_file (_type_): KEGG_clean.txt
        expression_file (_type_): fpkm_matrix
    """
    if not os.path.exists(ko):
        logger.critical(f"{ko} not found")
        sys.exit(1)
    if not os.path.exists(kegg_clean_file):
        logger.critical(f'{kegg_clean_file} not found')
        sys.exit(1)
    if not os.path.exists(expression_file):
        logger.critical(f"{expression_file} not found")
        sys.exit(1)
        
    ko_list = pd.read_csv(ko, sep='\t', usecols=['GeneID']).values.tolist()
    logger.info(f"正在输出每个 ko_pathway 相关的基因表达量表：{len(ko_list)} 个 ko_pathway")
    for ko_num in ko_list:
        kegg_pathway_df = pd.read_csv(kegg_clean_file, sep='\t', names=['GeneID', 'Ko'], usecols=[0, 1], dtype=str)
        kegg_pathway_df['Ko_Number'] = kegg_pathway_df['Ko'].str.split(":").str[0]
        kegg_pathway_df = kegg_pathway_df[kegg_pathway_df['Ko_Number'] == ko_num]
        expression_df = pd.read_csv(expression_file, sep='\t', dtype=str)
        expression_df = expression_df[expression_df['GeneID'].isin(kegg_pathway_df['GeneID'])]
        if expression_df.shape[0] > 0:
            expression_df.to_csv(ko_num + '_expression.txt', sep='\t', index=False)
            logger.info(f"正在输出 {ko_num} 相关的基因表达量表，有 {expression_df.shape[0]} 个基因")
        else:
            logger.warning(f"{ko_num} 相关的基因表达量表为空")


def parse_input():
    argparser = argparse.ArgumentParser()
    # 其他参数
    argparser.add_argument('--ko-list', dest="ko_list", type=str, required=True,
                           help='[必须]ko_list file, 必须要有列名，KEGG_ID 如果运行脚本类型是其他，只需要 KEGG_ID 就行')
    argparser.add_argument('--kegg-clean', dest="kegg_clean", required=True, type=str, help='[必须]输入 kegg 注释出来的 KEGG_clean 文件')
    argparser.add_argument('--kogene-heatmap', dest='kogene_heatmap', action='store_true',
                           help='[可选]draw kogene heatmap picture')
    argparser.add_argument('--expression', type=str,
                        help='[其他可选，转录组必须]输入表达量（fpkm_matrix_filtered.txt）文件')
    
    # 针对其他
    group1 = argparser.add_argument_group('针对其他的使用参数')
    group1.add_argument('--kegg-genedef', dest="kegg_genedef", type=str, help='[必须]kegg_gene_def.txt')
    group1.add_argument('--basicinfo', type=str, help='[可选]basicinfo.txt')
    group1.add_argument('--output-prefix', dest="output_prefix", type=str, help='output file prefix, 转录组不需要加这个参数')
    
    # 针对转录组，de results 结果 up down 的 id list 画 pahtview 图
    group2 = argparser.add_argument_group('针对转录组的使用参数，如果跑转录组的项目，下面必须的参数必须输入')
    
    group2.add_argument('--deg-data-dir', dest='deg_data_dir', type=str, help='[必须]输入 DEG_data.txt 文件')
    group2.add_argument('--enrich-dir', dest='enrich_dir', type=str, help='[必须]输入富集分析的文件夹')
    group2.add_argument('-s', '--samplesinfo', type=str, help='[必须]输入样品信息文件，如果有则添加到输出文件中')
    group2.add_argument('--output-dir', dest='output_dir', type=str, help='输出目录, 默认当前目录', default='.')
    
    # 画热图添加的参数
    group2.add_argument('--treat-color', dest="treat_color", type=str, default="blue",
                        help='[可选]画 heatmap 时指定 Treat Group 的颜色，默认 blue')
    group2.add_argument('--control-color', dest="control_color", type=str, default="red",
                        help='[可选]画 heatmap 时指定 Control Group 的颜色, 默认 red')
    
    args = argparser.parse_args()

    return args


def main():
    args = parse_input()
    # 转录组
    if args.deg_data_dir:
        transcriptome(args)
    else:
        other(args)
        
        # 筛选 passed_path.txt 中的 ko
        output_passed_path_filename = f"{args.output_prefix}_passed_path.txt"
        passed_path(args.ko_list, output_passed_path_filename)
        # 画 pathview 图
        logger.info("正在画 pathview 图")
        draw_pathview(args.output_prefix + '_regulation.txt', args.output_prefix + '_ko_passed_path.txt')
            
        # 出表达量表
        if args.expression and args.kegg_clean:
            get_ko_expression(args.ko_list, args.kegg_clean, args.expression)
    
    logger.success('Done!')
    

if __name__ == '__main__':
    main()