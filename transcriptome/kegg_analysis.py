#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/29 17:19
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse
import subprocess
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
from Rscript import draw_multigroup_heatmap
from Rscript import draw_twogroup_heatmap
from Rscript import draw_pathview
from Rscript import anova_analysis


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
    
    ko_df = pd.read_csv(ko_file, sep='\t', names=['KEGG_ID'], dtype=str)
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


def transcriptome_all_ko_gene_heatmap(ko_file, kegg_clean_file, expression_file, samples_file):
    """针对一些 KEGG_ID 的相关基因画出一个 heatmap 图

    Args:
        ko_file (str): 张老师给的 KEGG_ID 文件
        kegg_clean_file (str): kegg 注释出来的 KEGG_clean.txt
        fpkm_matrix_file (str): fpkm 矩阵文件
        samples_file (str): 样本描述文件，通常为 samples_described.txt 
    """
    # 对 expression fpkm matrix 文件只保留 fpkm 值
    expression_df = pd.read_csv(expression_file, sep='\t', dtype={'GeneID': str})
    expression_df = expression_df.drop_duplicates(subset='GeneID')
    expression_df_columns = expression_df.columns.tolist()
    expression_df_columns = [expression_df_columns[0]] + [x for x in expression_df_columns if x.endswith("_fpkm")]
    expression_df = expression_df[expression_df_columns]
    
    onlyfpkm_expression_columns = [x.replace("_fpkm", "") for x in expression_df_columns]
    expression_df.columns = onlyfpkm_expression_columns
    
    samples_df = pd.read_csv(samples_file, sep='\t', usecols=[0, 1])
    samples_df = samples_df[['sample', 'group']]
    
    ko_df = pd.read_csv(ko_file, sep='\t')
    ko_df = ko_df[['KEGG_ID', 'Ontology']]  # 只保留有用的两列
    
    kegg_pathway_df = pd.read_csv(kegg_clean_file, sep='\t', names=['GeneID', 'Ko'], usecols=[0, 1], dtype=str)
    kegg_pathway_df['KEGG_ID'] = kegg_pathway_df['Ko'].str.split(':').str[0]
    kegg_pathway_df = kegg_pathway_df.drop(columns=['Ko'])
    
    gene_df = pd.merge(left=kegg_pathway_df, right=ko_df, how='left', on='KEGG_ID')
    gene_df = gene_df.dropna(subset=['Ontology', 'GeneID'])
    gene_df = gene_df[gene_df['Ontology'] != 'Others']  # kegg 画图不需要 Others
    # gene_df = gene_df.drop(columns=['KEGG_ID'])
    
    # 添加 fpkm
    gene_fpkm_df = pd.merge(left=gene_df, right=expression_df, how='left', on='GeneID')
    # gene_fpkm_df = gene_fpkm_df.drop(columns=['Ontology'])
    
    # anova 计算
    anova_file_name = 'all_ko_gene_anova_p.txt'
    gene_fpkm_df[onlyfpkm_expression_columns].to_csv(anova_file_name, sep='\t', index=False)
    anova_analysis(anova_file_name, samples_file, anova_file_name)
    anova_gene_fpkm_df = pd.read_csv(anova_file_name, sep='\t')
    anova_gene_fpkm_df = anova_gene_fpkm_df[anova_gene_fpkm_df['p_value'] <= 0.05]
    anova_gene_fpkm_df.drop(columns=['p_value', 'BH_p_value'], inplace=True)
    anova_gene_fpkm_df = pd.merge(anova_gene_fpkm_df, gene_df, on='GeneID', how='left')
    anova_gene_fpkm_df = anova_gene_fpkm_df.sort_values(by=['Ontology'])
    
    
    # multigroup_heatmap 输入文件
    all_gene_ko_heatmap_filename = 'all_ko_gene_heatmap.xlsx'
    with pd.ExcelWriter(all_gene_ko_heatmap_filename, engine='openpyxl') as writer:
        anova_gene_fpkm_df[onlyfpkm_expression_columns].to_excel(writer, sheet_name="Sheet1", index=False)
        samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
        anova_gene_fpkm_df[['GeneID', 'Ontology']].to_excel(writer, sheet_name='Sheet3', index=False)
    
    draw_multigroup_heatmap(
        all_gene_ko_heatmap_filename, 
        all_gene_ko_heatmap_filename.replace('xlsx', '.jpeg'),
        other_args="--cluster-rows"
        )


def transcriptome_each_ko_gene_heatmap(kegg_id_list, kegg_clean_file, expression_file, samples_file, output_dir='.'):
    """针对一些 KEGG_ID 的相关基因，每一个 KEGG_ID 画出一个 heatmap 图

    Args:
        kegg_id_list (list): 张老师给的 KEGG_ID list
        kegg_clean_file (str): kegg 注释出来的 KEGG_clean.txt
        fpkm_matrix_file (str): fpkm 矩阵文件
        samples_file (str): 样本描述文件，通常为 samples_described.txt
    """
    samples_df = pd.read_csv(samples_file, sep='\t', usecols=[0, 1])
    samples_df = samples_df[['sample', 'group']]
    samples_list = ['GeneID'] + samples_df['sample'].values.tolist()

    # 对 expression fpkm matrix 文件只保留 fpkm 值
    expression_df = pd.read_csv(expression_file, sep='\t', dtype={'GeneID': str})
    expression_df = expression_df.drop_duplicates(subset='GeneID')
    expression_df_columns = expression_df.columns.tolist()
    expression_df_columns = [expression_df_columns[0]] + [x for x in expression_df_columns if x.endswith("_fpkm")]
    expression_df = expression_df[expression_df_columns]
    
    onlyfpkm_expression_columns = [x.replace("_fpkm", "") for x in expression_df_columns]
    expression_df.columns = onlyfpkm_expression_columns
    for sample in samples_list:
        if sample not in onlyfpkm_expression_columns:
            logger.critical(f"{sample} 未在 {expression_file} 表里找到")
            sys.exit(1)
    expression_df = expression_df[samples_list]
    
    # ko_df = pd.read_csv(ko_file, sep='\t')
    # ko_df = ko_df[['KEGG_ID',]]  # 只保留有用的两列
    
    kegg_pathway_df = pd.read_csv(kegg_clean_file, sep='\t', names=['GeneID', 'Ko'], usecols=[0, 1], dtype=str)
    kegg_pathway_df['KEGG_ID'] = kegg_pathway_df['Ko'].str.split(':').str[0]
    kegg_pathway_df = kegg_pathway_df.drop(columns=['Ko'])
    kegg_pathway_df = kegg_pathway_df.drop_duplicates(subset=['GeneID'])
    
    # gene_df = pd.merge(left=kegg_pathway_df, right=ko_df, how='left', on='KEGG_ID')
    # gene_df = gene_df.dropna(subset=['GeneID'])                                                                                                                                                                                                                                                                                                                            
        
    crt_kegg_id_dir = os.path.join(output_dir, 'All_groups_KEGG_analysis')
    if not os.path.exists(crt_kegg_id_dir):
        os.mkdir(crt_kegg_id_dir)
        
    # kegg_id_list = gene_df['KEGG_ID'].values.tolist()
    for each_kegg_id in kegg_id_list:
        kegg_pic_dir = os.path.join(crt_kegg_id_dir, 'KEGG_pathway_heatmap')
        if not os.path.exists(kegg_pic_dir):
            os.mkdir(kegg_pic_dir)
        logger.info(f'正在对 {each_kegg_id} 的相关基因画 heatmap 图，结果输出到 {kegg_pic_dir}')
        each_kegg_id_df = kegg_pathway_df[kegg_pathway_df['KEGG_ID'] == each_kegg_id]
        
        # kegg 相关的 id 小于 10 个就跳过
        if each_kegg_id_df.shape[0] < 10:
            logger.warning(f'{each_kegg_id} 相关基因数量小于 10 个，不对此画 heatmap 图')
            continue

        # 添加 fpkm
        each_kegg_id_gene_fpkm_df = pd.merge(each_kegg_id_df, expression_df, on='GeneID', how='left')
        
        # anova 计算输入文件
        anova_file_name = f'{crt_kegg_id_dir}/{each_kegg_id}_gene_anova_p.txt'
        each_kegg_id_gene_fpkm_df[samples_list].to_csv(anova_file_name, sep='\t', index=False)
        anova_result = anova_analysis(anova_file_name, samples_file, anova_file_name)
        if not anova_result:
            logger.error(f'{each_kegg_id} anova 计算结果失败，跳过执行 multigroup heatmap')
            continue
        each_kegg_id_gene_fpkm_df = pd.read_csv(anova_file_name, sep='\t')
        each_kegg_id_gene_fpkm_df = each_kegg_id_gene_fpkm_df[each_kegg_id_gene_fpkm_df['p_value'] <= 0.05]
        each_kegg_id_gene_fpkm_df.drop(columns=['p_value', 'BH_p_value'], inplace=True)
        
        # 2024_06_14 张老师：注释掉这个，不需要过滤 p 值
        # kegg 相关的 id 小于 10 个就跳过
        # if each_kegg_id_gene_fpkm_df.shape[0] < 10:
        #    logger.warning(f'{each_kegg_id} 相关基因根据 p 值 < 0.05 筛选后数量小于 10 个，不对此画 heatmap 图')
        #    continue
        
        # multigroup heatmap 输入文件
        kegg_id_gene_ko_heatmap_filename = f'{crt_kegg_id_dir}/{each_kegg_id}_ko_gene_heatmap.xlsx'
        with pd.ExcelWriter(kegg_id_gene_ko_heatmap_filename, engine='openpyxl') as writer:
            each_kegg_id_gene_fpkm_df[samples_list].to_excel(writer, sheet_name='Sheet1', index=False)
            samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
    
        # heatmap 画图
        kegg_id_heatmap_filename = os.path.join(kegg_pic_dir, f"{each_kegg_id}_ko_gene_heatmap.jpeg")
        heatmap_result = draw_multigroup_heatmap(
            kegg_id_gene_ko_heatmap_filename,
            kegg_id_heatmap_filename,
            other_args='--cluster-rows'
            )
        if not heatmap_result:
            logger.error(f'{each_kegg_id} draw_multigroup_heatmap 结果失败，跳过执行 multigroup heatmap')


def transcriptome(args):
    """针对转录组的 KEGG 分析

    Args:
        args (_type_): 所有的输入参数
    """
    output_dir = args.output_dir
    # 循环每个 ko number 筛选表达量画热图
    ko_df = pd.read_csv(args.ko_list, sep='\t')
    ko_list = ko_df['KEGG_ID'].values.tolist()
    if 'KEGG_ID'.lower() in ko_list:
        ko_list = ko_list[1:]

    transcriptome_each_ko_gene_heatmap(ko_list, args.kegg_clean, args.expression, args.samplesinfo)
    
    deg_data_list = os.listdir(args.deg_data_dir)
    kegg_pathway_df = pd.read_csv(args.kegg_clean, sep='\t', names=['GeneID', 'Ko'], usecols=[0, 1], dtype=str)
    samples_described_df = pd.read_csv(args.samplesinfo, sep='\t', usecols=[0, 1], names=["group", "sample"], dtype=str)
    for deg_data_file in deg_data_list:
        
        if not deg_data_file.endswith('_DEG_data.txt'):
            continue
        deg_data_file = os.path.join(args.deg_data_dir, deg_data_file)
        deg_data_df = pd.read_csv(deg_data_file, sep='\t', dtype=str)
        
        treat_group = deg_data_df.iloc[0, 1]
        control_group = deg_data_df.iloc[0, 2]
        compare_info = treat_group + '_vs_' + control_group
        
        compare_info_dir = compare_info
        deg_expression_data_dir = os.path.join(compare_info, 'DEG_expression_data')
        kegg_heatmap_dir = os.path.join(compare_info, 'KEGG_heatmap')
        kegg_pathway_graph_dir = os.path.join(compare_info, 'KEGG_pathway_graph')
        
        for each_dir in [compare_info_dir, deg_expression_data_dir, kegg_heatmap_dir, kegg_pathway_graph_dir]:
            if not os.path.exists(each_dir):
                os.mkdir(each_dir)
            else:
                logger.warning(f"{each_dir} 文件夹已存在，输出将覆盖原文件")
                
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
        
        heatmap_sheet2_df = pd.DataFrame(columns=['samples', 'group', 'colors'])
        for each_sample in crt_group_samples:
            if each_sample in treat_group_samples_name:
                group = treat_group
                color = args.treat_color
            else:
                group = control_group
                color = args.control_color
            row_data = {'samples': each_sample, 'group': group, 'colors': color}
            row_df = pd.DataFrame([row_data])
            heatmap_sheet2_df = pd.concat([heatmap_sheet2_df, row_df], ignore_index=True)
            # heatmap_sheet2_df = heatmap_sheet2_df.append({'samples': each_sample, 'group': group, 'colors': color}, ignore_index=True)

        for ko_number in ko_list:
            ko_num_fpkm_expr_df = fpkm_df[fpkm_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['Ko'].str.contains(ko_number)]['GeneID'])]
            
            if ko_num_fpkm_expr_df.shape[0] == 0:
                logger.warning(f"{ko_number} 在 {compare_info} 中相关的基因表达量表为空")
                continue
            
            ko_num_fpkm_expr_df_file = f"{compare_info}_{ko_number}_fpkm_expression.xlsx"
            ko_num_fpkm_expr_pic_name = os.path.join(kegg_heatmap_dir, f"{compare_info}_{ko_number}_heatmap.jpeg")
            with pd.ExcelWriter(ko_num_fpkm_expr_df_file) as writer:
                ko_num_fpkm_expr_df.to_excel(writer, index=False, sheet_name='Sheet1')
                heatmap_sheet2_df.to_excel(writer, index=False, sheet_name='Sheet2')
            if ko_num_fpkm_expr_df.shape[0] > 1:
                logger.info(f"正在画 {compare_info} {ko_number} 相关基因表达量的热图")
                draw_twogroup_heatmap(ko_num_fpkm_expr_df_file, ko_num_fpkm_expr_pic_name)
            else:
                logger.warning(f"{compare_info} 的 {ko_number} 中相关的基因表达量表只有一行，无法画热图")
            
            # 输出每个 ko 的 deg_data 文件
            ko_num_deg_data_file = os.path.join(deg_expression_data_dir, f"{compare_info}_{ko_number}_DEG_data.txt")
            ko_num_deg_data_df = deg_data_df[deg_data_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['Ko'].str.contains(ko_number)]['GeneID'])]
            ko_num_deg_data_df.to_csv(ko_num_deg_data_file, sep='\t', index=False)
            
        # R pathview 画图准备文件 regulation
        # print(f"正在准备 {compare_info} 的 regulation 文件")
        up_down_idlist_geneid2kid(down_df, up_df, compare_info, args.kegg_clean)
        
        # R pathview 画图准备文件 passed path
        # print(f"正在筛选 {compare_info} 的 passed_path.txt")
        passed_path(args.ko_list, compare_info)
        
        # R pathview 画图
        if args.draw_pathview:
            logger.info(f"正在画 {compare_info} 的 pathview 图")
            # 这两个文件用完需要删掉
            draw_pathview(f"{compare_info}_regulation.txt", f"{compare_info}_ko_passed_path.txt")
        
        # 对每个比较组的文件进行整理
        os.system(f"mv {compare_info}_*.xlsx {compare_info}")
        os.system(f"mv *.png {kegg_pathway_graph_dir}")


def passed_path(ko_file, output_prefix):
    ko_def_df = pd.read_csv(ko_file, sep='\t', names=['KEGG_ID'], dtype=str, usecols=[0], header=0)
    passed_path_file = '/home/colddata/qinqiang/script/Rscript/pathview/passed_path.txt'
    passed_path_df = pd.read_csv(passed_path_file, sep='\t', names=['KEGG_ID', 'Ko_Def'], dtype=str)
    ko_def_df = passed_path_df[passed_path_df['KEGG_ID'].isin(ko_def_df['KEGG_ID'])]
    ko_def_df = ko_def_df.drop_duplicates(subset=['KEGG_ID'])
    ko_def_df.to_csv(output_prefix + '_ko_passed_path.txt', sep='\t', index=False, header=False)


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
    argparser = argparse.ArgumentParser(description='Add definition to gene file')
    # 其他参数
    argparser.add_argument('--ko-list', dest="ko_list", type=str, required=True,
                           help='[必须]ko_list file, 必须要有列名，KEGG_ID 如果运行脚本类型是其他，只需要 KEGG_ID 就行')
    argparser.add_argument('--kegg-clean', dest="kegg_clean", required=True, type=str, help='[必须]输入 kegg 注释出来的 KEGG_clean 文件')
    argparser.add_argument('--draw-pathview', dest='draw_pathview', action='store_true',
                           help='[可选]draw pathview picture')
    argparser.add_argument('--expression', type=str,
                        help='[其他可选，转录组必须]输入表达量（fpkm_and_reads_matrix.txt）文件, 输出每个 ko_pathway 相关的基因表达量表')
    
    # 针对其他
    group1 = argparser.add_argument_group('针对其他的使用参数')
    group1.add_argument('--kegg-genedef', dest="kegg_genedef", type=str, help='[必须]kegg_gene_def.txt')
    group1.add_argument('--basicinfo', type=str, help='[可选]basicinfo.txt')
    group1.add_argument('--output-prefix', dest="output_prefix", type=str, help='output file prefix, 转录组不需要加这个参数')
    
    # 针对转录组，de results 结果 up down 的 id list 画 pahtview 图
    group2 = argparser.add_argument_group('针对转录组的使用参数，如果跑转录组的项目，下面必须的参数必须输入')
    
    group2.add_argument('--deg-data-dir', dest='deg_data_dir', type=str, help='[必须]输入 DEG_data.txt 文件')
    group2.add_argument('-s', '--samplesinfo', type=str, help='[必须]输入样品信息文件，如果有则添加到输出文件中')
    group2.add_argument('--output-dir', dest='output_dir', type=str, help='输出目录, 默认当前目录', default='.')
    
    # 画热图添加的参数    
    group2.add_argument('--heatmap', action='store_true',
                        help='[可选]整理表达量文件，输出差异基因的热图')
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
        if not os.path.exists(args.kegg_genedef):
            logger.critical(f"Error: {args.kegg_genedef} not found")
            sys.exit(1)
        if not os.path.exists(args.kegg_clean):
            logger.critical(f"Error: {args.kegg_clean} not found")
            sys.exit(1)
        
        other(args)
        
        # 筛选 passed_path.txt 中的 ko
        passed_path(args.ko_list, args.output_prefix)
        # 画 pathview 图
        if args.draw_pathview:
            logger.info("正在画 pathview 图")
            draw_pathview(args.output_prefix + '_regulation.txt', args.output_prefix + '_ko_passed_path.txt')
            
        # 出表达量表
        if args.expression and args.kegg_clean:
            get_ko_expression(args.ko_list, args.kegg_clean, args.expression)
    
    logger.success('Done!')
    

if __name__ == '__main__':
    main()