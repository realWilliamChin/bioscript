#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/03/27 14:01
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from pathlib import Path
from loguru import logger
import datetime

import gene_interaction_network_plot
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
from Rscript import enrichment_barplot
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df
from data_check import df_drop_element_side_space
from data_check import df_replace_illegal_folder_chars
from logger_config import get_logger

if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    sys.exit(1)


def parse_input():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-i', '--input', help="输入文件，文件中至少包含三列，GO_ID、Ontology、SubOntology")
    argparser.add_argument('-c', '--compare', help='输入 compare_info.txt 文件')
    argparser.add_argument('-g', '--genego', help='GOID 文件,swiss 注释出来的的 gene_go.txt')
    argparser.add_argument('-s', '--samples-described', help='输入 samples_described.txt 文件')
    argparser.add_argument('-e', '--enrich-data-dir', help='转录组 输入 enrich.r 运行出来的结果文件夹')
    argparser.add_argument('-d', '--deg-data-dir', help='转录组 输入 DEG_data.txt 文件目录')
    argparser.add_argument('-o', '--output', default=os.getcwd(), help='输出目录，默认为当前目录，输出目录如果不存在则会尝试创建')
    argparser.add_argument('-l', '--log-name', dest='log_name', help='日志文件名（不需要包含.log后缀），默认为 go_analysis + 当前时间')
    
    args = argparser.parse_args()
    
    # 配置日志
    current_time = datetime.datetime.now().strftime("%Y%m%d%H%M")
    log_name = args.log_name if args.log_name else f"go_analysis_{current_time}.log"
    global logger
    logger = get_logger(log_name)
    
    # Data check
    for f in [args.input, args.genego]:
        if not os.access(f, os.R_OK):
            logger.critical(f"输入文件 {f} 未找到或不可读")
            sys.exit(1)
    
    return args


def goid_enrich_summary(target_go_df: pd.DataFrame, comparisons_df: pd.DataFrame, enrich_data_dir: str, output: str) -> None:
    # target_go_df = target_go_df.rename(columns={'GO_ID': 'GO_Def', 'ID': 'GO_ID'})
    target_go_list = target_go_df['GO_ID'].values.tolist()
    target_go_enrich_df_list = []
    comparisons_df['comparisons'] = comparisons_df['Treat'] + '-vs-' + comparisons_df['Control']
    comparison_list = comparisons_df['comparisons'].values.tolist()
    for comparison in comparison_list:
        up_enrich_df = load_table(os.path.join(enrich_data_dir, f'{comparison}_Up_EnrichmentGO.xlsx'))
        down_enrich_df = load_table(os.path.join(enrich_data_dir, f'{comparison}_Down_EnrichmentGO.xlsx'))
        target_go_enrich_up_df = up_enrich_df[up_enrich_df['ID'].isin(target_go_list)]
        target_go_enrich_down_df = down_enrich_df[down_enrich_df['ID'].isin(target_go_list)]
        
        target_go_enrich_up_df.insert(0, 'Group', comparison)
        target_go_enrich_up_df.insert(1, 'Regulation', 'Up')
        target_go_enrich_down_df.insert(0, 'Group', comparison)
        target_go_enrich_down_df.insert(1, 'Regulation', 'Down')
        
        target_go_enrich_df_list.append(target_go_enrich_up_df)
        target_go_enrich_df_list.append(target_go_enrich_down_df)
    
    target_go_enrich_summary_df = pd.concat(target_go_enrich_df_list)
    target_go_enrich_summary_df.drop(columns=['Ontology'], inplace=True)
    target_go_enrich_summary_df.rename(columns={'ID': 'GO_ID'}, inplace=True)
    if 'Ontology' in target_go_df.columns and 'SubOntology' in target_go_df.columns:
        target_go_enrich_summary_df = pd.merge(
            target_go_df[['GO_ID', 'Ontology', 'SubOntology']],
            target_go_enrich_summary_df,
            how='right',
            on='GO_ID'
        )
    target_go_enrich_summary_df.sort_values(
        by=['Ontology', 'SubOntology', 'Group', 'Regulation'],
        ascending=[True, True, True, False],
        inplace=True
    )
    write_output_df(target_go_enrich_summary_df, output, index=False)
    

def deg_output_summary(df_list, samples_info_df):
    processed_df_list = []
    max_samples_number = samples_info_df.groupby('group').size().max()
    for df in df_list:
        treat = df['sampleA'].values.tolist()[0]
        control = df['sampleB'].values.tolist()[0]
        treat_samples = samples_info_df[samples_info_df['group'] == treat]['sample'].values.tolist()
        control_samples = samples_info_df[samples_info_df['group'] == control]['sample'].values.tolist()
        
        # 获取原始的 FPKM 列和 reads 列
        df_treat_fpkm_column = [f'{x}_FPKM' for x in treat_samples]
        df_control_fpkm_column = [f'{x}_FPKM' for x in control_samples]
        df_treat_reads_column = [f'{x}_raw_reads' for x in treat_samples]
        df_control_reads_column = [f'{x}_raw_reads' for x in control_samples]
        
        # Treat FPKM
        for i in range(1, max_samples_number + 1):
            new_treat_fpkm = f"Treat_{i}_FPKM"
            if i <= len(df_treat_fpkm_column):
                df = df.rename(columns={df_treat_fpkm_column[i-1]: new_treat_fpkm})
            else:
                df[new_treat_fpkm] = 'N/A'

        # Control FPKM
        for i in range(1, max_samples_number + 1):
            new_control_fpkm = f"Control_{i}_FPKM"
            if i <= len(df_control_fpkm_column):
                df = df.rename(columns={df_control_fpkm_column[i-1]: new_control_fpkm})
            else:
                df[new_control_fpkm] = 'N/A'

        # Treat reads
        for i in range(1, max_samples_number + 1):
            new_treat_reads = f"Treat_{i}_reads"
            if i <= len(df_treat_reads_column):
                df = df.rename(columns={df_treat_reads_column[i-1]: new_treat_reads})
            else:
                df[new_treat_reads] = 'N/A'

        # Control reads
        for i in range(1, max_samples_number + 1):
            new_control_reads = f"Control_{i}_reads"
            if i <= len(df_control_reads_column):
                df = df.rename(columns={df_control_reads_column[i-1]: new_control_reads})
            else:
                df[new_control_reads] = 'N/A'
        
        processed_df_list.append(df)
    
    output_summary_df = pd.concat(processed_df_list)
    
    return output_summary_df


def deg_go_analysis(args, target_go_df, gene_go_df, go_id_list, ontology_list):
    output_summary_df_list = []
    # 循环每个差异数据文件
    
    for enrich_data in os.listdir(args.enrich_data_dir):
        if not enrich_data.endswith("_EnrichmentGO.xlsx") or enrich_data.startswith("~"):
            continue
        
        logger.info(f'\n====正在处理 {enrich_data}====')
        compare_info = enrich_data.replace('_EnrichmentGO.xlsx', '')  # 比对名称
        
        # 文件输出目录
        compare_output_dir = os.path.join(args.output, compare_info)
        go_analysis_network_graph_dir = os.path.join(compare_output_dir, 'GO_analysis_network_graph')
        go_analysis_bar_plot_dir = os.path.join(compare_output_dir, 'GO_analysis_bar_plot')
        go_deg_expression_data_dir = os.path.join(compare_output_dir, 'DEG_expression_data')
        for each_dir in compare_output_dir, go_analysis_bar_plot_dir, go_analysis_network_graph_dir, go_deg_expression_data_dir:
            if not os.path.exists(each_dir):
                os.mkdir(each_dir)
            else:
                logger.warning(f'{each_dir} 文件夹已存在，输出文件将覆盖源文件')
                
        enrich_data_abspath = os.path.join(args.enrich_data_dir, enrich_data)
        enrich_go_df = load_table(enrich_data_abspath)
        enrich_go_df = enrich_go_df.drop(columns=['Ontology'])
        
        
        for ontology_name in ontology_list:
            logger.info(f"正在处理 {compare_info}_{ontology_name}")
            # gene_id_goid = {}
            ontology_df = target_go_df[target_go_df['Ontology'] == ontology_name].copy()
            ontology_df = pd.merge(left=ontology_df, right=enrich_go_df, how='inner', left_on='GO_ID', right_on='ID')
            ontology_df.drop(columns=['ID'], inplace=True)
            ontology_df.rename(columns={"GO_ID": "ID"}, inplace=True)
            
            if ontology_df.shape[0] == 0:
                logger.warning(f"{compare_info}_{ontology_name}, 没有筛选出任何数据，不进行画图，跳过")
                continue
            
            # gene_id_goid 中的值从 list 变为字符串
            # for each_gene in gene_id_goid.keys():
            #     gene_id_goid[each_gene] = ','.join(gene_id_goid[each_gene])
            # gene_id_goid_df = pd.DataFrame(gene_id_goid.items(), columns=['GeneID', 'GOID'])
            # gene_id_goid_df.to_csv(os.path.join(args.output, f'{compare_name}_{ontology_name}_gene_id_goid.txt'), sep='\t', index=False)
            
            # ontology_df 中 SubOntology 的类型出现次数小于 4 个的不要
            # ontology_df = ontology_df[ontology_df['SubOntology'].map(ontology_df['SubOntology'].value_counts()) >= 4]
            # ❌【修改为下面的规则】ontology_df 中 SubOntology 的类型出现次数大于 4 个的，保留 RichFactor 值最大的 4 个
            # [2024_04_18 张老师] ontology_df 中的 SubOntology 的类型按照出现次数最少的作为标准，对其他的类型进行筛选
            # [2024_04_22 潘老师] 不进行筛选了，直接画图，有多少是多少
            # min_type_count = ontology_df['SubOntology'].value_counts().to_list()[-1]
            # ontology_df = ontology_df.groupby('SubOntology').apply(lambda x: x.nlargest(min_type_count, 'RichFactor')) #.reset_index(drop=True)

            ontology_df = ontology_df.sort_values(by=['ID'])
            output_excel_name = os.path.join(compare_output_dir, f'{compare_info}_{ontology_name}.xlsx')
            write_output_df(ontology_df, output_excel_name, index=False)
            
            add_summary_df = ontology_df.copy()
            add_summary_df.insert(0, 'Group', compare_info.rsplit('_', 1)[0])
            add_summary_df.insert(1, 'Regulation', compare_info.split('_')[-1])
            # output_summary_df_list.append(add_summary_df)
            
            jpeg_file_name = os.path.join(go_analysis_bar_plot_dir, f'{compare_info}_{ontology_name}_barplot.jpeg')
            # 运行 R barplot 脚本
            if ontology_df.shape[0] > 1:
                enrichment_barplot(output_excel_name, jpeg_file_name)
            else:
                logger.warning(f"{compare_info}_{ontology_name}, 只筛选出 1 条数据，不画 barplot 图")
            
            # 画 enrichnet 图
            enrichnet_file_name = os.path.join(go_analysis_network_graph_dir, f'{compare_info}_{ontology_name}_enrichnet.png')
            enrichnet_data = {'source': [], 'target': []}
            for index, row in ontology_df.iterrows():
                subontologys = [row['SubOntology']] * len(row['geneID'].split('/'))
                genes = row['geneID'].split('/')
                enrichnet_data['source'].extend(subontologys)
                enrichnet_data['target'].extend(genes)
            enrichnet_df = pd.DataFrame(enrichnet_data)
            enrichnet_df = enrichnet_df.drop_duplicates()
            gene_interaction_network_plot.draw_enrichnetplot(enrichnet_df, enrichnet_file_name)
        
        # 对每个 GO_ID 的相关基因，对上 DEG_data.txt 出一个表
        for go_id in go_id_list:
            compare_name = compare_info.replace('_Down', '').replace('_Up', '')
            compare_regulation = compare_info.split('_')[-1]
            # 因为循环的是 enrich data 不是 DEG_data，使用这种方式找到对应的 DEG_data.txt 文件
            deg_data_name = os.path.join(args.deg_data_dir, f"{compare_name}_DEG_data.txt")
            go_id_deg_data_filename = os.path.join(go_deg_expression_data_dir, f'{compare_name}_{go_id}_DEG_data.txt').replace(':', '_')

            deg_data_df = load_table(deg_data_name, dtype={"GeneID": str})
            go_id_deg_data_df = deg_data_df[
                deg_data_df['regulation'].str.contains(compare_regulation) &
                deg_data_df['GeneID'].isin(gene_go_df[gene_go_df['GO_ID'].str.contains(go_id)]['GeneID'])
            ]
            if go_id_deg_data_df.shape[0] == 0:
                logger.warning(f'{compare_info} 的 {go_id} 未找到相关基因')
                continue
            write_output_df(go_id_deg_data_df, go_id_deg_data_filename, index=False)
            
            go_id_deg_data_df
            output_summary_df_list.append(go_id_deg_data_df)
    
    # 输出汇总
    if len(output_summary_df_list) == 0:
        logger.warning('Target GO ID 在组间中没有任何目标基因')
    else:
        samples_info_df = load_table(args.samples_described)
        output_summary_df = deg_output_summary(output_summary_df_list, samples_info_df)
        write_output_df(output_summary_df, os.path.join(args.output, 'Target_GO_analysis_summary.xlsx'), index=False)


def main():
    args = parse_input()
    target_go_df = load_table(args.input)
    comparison_df = load_table(args.compare)
    logger.info(f'正在对输入数据进行预处理，去除两边空格，替换 Ontology 非法字符')
    target_go_df = df_drop_element_side_space(target_go_df)
    target_go_df = df_replace_illegal_folder_chars(target_go_df, ['Ontology', 'SubOntology'])
    target_go_df['GO_ID'] = target_go_df['GO_ID'].str.split("_").str[0]  # 为了和 enrich.r 出来的文件的 ID 对应上，添加一列只包含 ID 的列， GO:0010111
    ontology_list = list(set(target_go_df['Ontology'].tolist()))
    go_id_list = list(set(target_go_df['GO_ID'].tolist()))
    gene_go_df = load_table(args.genego, header=None, names=['GeneID', 'GO_ID'], dtype={"GeneID": str})
    
    logger.info(f'执行 go 分析')
    goid_enrich_summary(target_go_df, comparison_df, args.enrich_data_dir, os.path.join(args.output, 'Target_GO_Enrich_data_summary.xlsx'))
    deg_go_analysis(args, target_go_df, gene_go_df, go_id_list, ontology_list)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()