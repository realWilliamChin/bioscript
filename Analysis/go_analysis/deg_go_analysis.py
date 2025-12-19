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
from goid_enrich_summary import goid_enrich_summary
from each_go_gene_heatmap import each_go_gene_heatmap, each_go_gene_expression
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
from Rscript import enrichment_barplot, smart_heatmap
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df
from data_check import df_drop_element_side_space, df_replace_illegal_folder_chars
from logger_config import get_logger

if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    sys.exit(1)


def parse_input():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-i', '--input', help="输入文件，文件中至少包含三列，GO_pathway_ID、Ontology、SubOntology")
    argparser.add_argument('-c', '--compare', help='输入 compare_info.txt 文件')
    argparser.add_argument('-g', '--genego', help='swiss 注释出来的的 gene_go.txt')
    argparser.add_argument('-s', '--samples-described', help='输入 samples_described.txt 文件')
    argparser.add_argument('-e', '--enrich-data-dir', help='转录组 输入 enrich.r 运行出来的结果文件夹')
    argparser.add_argument('-d', '--deg-data-dir', help='转录组 输入 DEG_data.txt 文件目录')
    argparser.add_argument('-f', '--fpkm-matrix', help='输入 fpkm_matrix.txt 文件(为了输出每个 GO_pathway_ID 的相关基因的 Heatmap 图)')
    argparser.add_argument('-r', '--expression-data', help='[暂时不用] 输入 reads_fpkm_matrix_def.txt 文件(为了输出每个 GO_ID 的相关基因的表达量数据文件)')
    argparser.add_argument('--genesymbol', help='使用 GeneSymbol 作为索引画热图，人和大鼠小鼠通常使用，输入包含 GeneID GeneSymbol 两列的文件')
    
    argparser.add_argument('-o', '--output', default=os.getcwd(), help='输出目录，默认为当前目录，输出目录如果不存在则会尝试创建')
    argparser.add_argument('-l', '--log-name', dest='log_name', help='日志文件名（不需要包含.log后缀），默认为 go_analysis + 当前时间')
    args = argparser.parse_args()
    
    # 配置日志
    current_time = datetime.datetime.now().strftime("%Y%m%d%H%M")
    log_name = args.log_name if args.log_name else f"go_analysis_{current_time}.log"
    global logger
    logger = get_logger(log_name)
    
    os.makedirs(args.output, exist_ok=True)
    
    return args
    

def _each_go_pathway_deg_data_summary(df_list, samples_info_df):
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


def each_go_pathway_deg_data(target_go_df, gene_go_df, compare_info_df, samples_df, deg_data_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_summary_df_list = []
    
    for _, row in compare_info_df.iterrows():
        compare_info = row['Treat'] + '-vs-' + row['Control']
        deg_data_df = load_table(os.path.join(deg_data_dir, f"{compare_info}_DEG_data.txt"))
        for go_pathway_id in target_go_df['GO_pathway_ID'].tolist():
            go_genes_df = gene_go_df[gene_go_df['GO_pathway_ID'].str.contains(go_pathway_id)]
            go_id_deg_data_df = pd.merge(deg_data_df, go_genes_df, on='GeneID', how='inner')
            
            if go_id_deg_data_df.shape[0] == 0:
                logger.warning(f'{compare_info} 的 {go_pathway_id} 未找到相关基因')
                continue
            
            output_summary_df_list.append(go_id_deg_data_df)

    # 输出汇总
    if len(output_summary_df_list) == 0:
        logger.warning('Target GO ID 在组间中没有任何目标基因')
    else:
        output_summary_df = _each_go_pathway_deg_data_summary(output_summary_df_list, samples_df)

        # 2025_12:18 张老师：不需要这个，只需要拆分出来每个 GO_pathway_ID 的
        # write_output_df(output_summary_df, os.path.join(output_dir, 'Target_GO_analysis_DEG_summary.xlsx'), index=False)

        # 循环 target_go_df 的 GO_ID，在 output_summary_df 中找到相关数据并输出
        logger.info('正在为每个 GO_pathway_ID 输出详细的 DEG 数据文件')
        for _, row in target_go_df.iterrows():
            go_pathway_id = row['GO_pathway_ID']
            go_pathway_def = row['GO_pathway_def']

            # 在 output_summary_df 中过滤出当前 GO_ID 的数据
            go_specific_df = output_summary_df[output_summary_df['GO_pathway_ID'] == go_pathway_id].copy()

            if go_specific_df.shape[0] == 0:
                logger.warning(f'GO_pathway_ID {go_pathway_id} 在 output_summary_df 中未找到相关数据，跳过')
                continue
            
            # 生成文件名
            filename = f"{go_pathway_id.replace(':', '_')}_{go_pathway_def}.txt"
            output_path = os.path.join(output_dir, filename)

            # 移除 GO_ID 列
            go_specific_df = go_specific_df.drop(columns=['GO_pathway_ID'])

            # 输出为 txt 文件，不包含 GO_ID 列
            write_output_df(go_specific_df, output_path, index=False)
            logger.info(f'已输出 GO_pathway_ID {go_pathway_id} 的数据到文件: {filename}')
        

def each_go_pathway_deg_data_heatmap(target_go_df, gene_go_df, compare_info_df, samples_df, deg_data_dir, output_dir, genesymbol_file=None):
    os.makedirs(output_dir, exist_ok=True)
    
    # 如果提供了 genesymbol 参数，读取 GeneID 和 GeneSymbol 映射文件
    genesymbol_df = None
    if genesymbol_file:
        genesymbol_df = load_table(genesymbol_file)
        # 只保留 GeneID 和 GeneSymbol 两列
        genesymbol_df = genesymbol_df[['GeneID', 'GeneSymbol']].copy()
        # 删除空值
        genesymbol_df.dropna(subset=['GeneID', 'GeneSymbol'], inplace=True)
        logger.info(f'读取到 {genesymbol_df.shape[0]} 个 GeneID-GeneSymbol 映射')
    
    for _, row in compare_info_df.iterrows():
        compare_info = row['Treat'] + '-vs-' + row['Control']
        deg_data_df = load_table(os.path.join(deg_data_dir, f"{compare_info}_DEG_data.txt"))
        os.makedirs(os.path.join(output_dir, compare_info), exist_ok=True)
        for go_pathway_id in target_go_df['GO_pathway_ID'].tolist():
            go_genes_df = gene_go_df[gene_go_df['GO_pathway_ID'].str.contains(go_pathway_id)]
            go_id_deg_data_df = pd.merge(deg_data_df, go_genes_df, on='GeneID', how='inner')
            
            if go_id_deg_data_df.shape[0] == 0:
                logger.warning(f'{compare_info} 的 {go_pathway_id} 未找到相关基因')
                continue
            elif go_id_deg_data_df.shape[0] < 2:
                logger.warning(f'{compare_info} 的 {go_pathway_id} 相关基因数量小于 2，跳过画 Heatmap 图')
                continue

            # 提取所有 _FPKM 列
            fpkm_columns = [col for col in go_id_deg_data_df.columns if col.endswith('_FPKM')]
            if len(fpkm_columns) == 0:
                logger.warning(f'{compare_info} 的 {go_pathway_id} 在 DEG_data.txt 中未找到 _FPKM 列，跳过')
                continue
            
            # 构建 data 表：包含 GeneID 和所有 FPKM 列
            data_df = go_id_deg_data_df[['GeneID'] + fpkm_columns].copy()
            # 将列名中的 _FPKM 后缀去掉，使其与 sample 名称匹配
            rename_dict = {col: col.replace('_FPKM', '') for col in fpkm_columns}
            data_df.rename(columns=rename_dict, inplace=True)
            
            # 数据验证：检查行数、列数和有效数据
            data_columns = [col for col in data_df.columns if col != 'GeneID']
            if (data_df.shape[0] < 2 or len(data_columns) == 0 or 
                data_df[data_columns].isna().all().all()):
                logger.warning(f'{compare_info} 的 {go_pathway_id} 数据不符合热图要求（行数<2、无有效列或全为NA），跳过画 Heatmap 图')
                continue
            
            # 如果提供了 genesymbol_df，则替换 GeneID 为 GeneSymbol
            if genesymbol_df is not None:
                # 合并 GeneSymbol 映射
                data_df = pd.merge(data_df, genesymbol_df, on='GeneID', how='inner')
                if data_df.shape[0] < 2:
                    logger.warning(f'{compare_info} 的 {go_pathway_id} 相关 genesymbol 数量小于 2，跳过画 Heatmap 图')
                    continue
                # 将 GeneID 列替换为 GeneSymbol
                data_df['GeneID'] = data_df['GeneSymbol']
                data_df.drop(columns=['GeneSymbol'], inplace=True)
            
            # 构建 group 和 sample 映射表
            # 从 FPKM 列名中提取 sample 名称（去掉 _FPKM 后缀）
            sample_names = [col.replace('_FPKM', '') for col in fpkm_columns]
            group_sample_df = pd.DataFrame({
                'sample': sample_names
            })
            # 从 samples_df 中获取对应的 group
            group_sample_df = pd.merge(group_sample_df, samples_df[['sample', 'group']], on='sample', how='left')
            # 确保列顺序为 sample, group
            group_sample_df = group_sample_df[['sample', 'group']]
            
            # 生成输出文件名
            go_pathway_def_row = target_go_df[target_go_df['GO_pathway_ID'] == go_pathway_id]
            if go_pathway_def_row.shape[0] > 0:
                go_pathway_def = go_pathway_def_row['GO_pathway_def'].values[0]
            else:
                go_pathway_def = go_pathway_id
            
            safe_go_id = go_pathway_id.replace(':', '_')
            # 处理非法字符
            temp_df = pd.DataFrame({'def': [go_pathway_def]})
            temp_df = df_replace_illegal_folder_chars(temp_df, ['def'])
            safe_go_def = temp_df['def'].values[0]
            output_filename = f"{compare_info}_{safe_go_id}_{safe_go_def}_heatmap.xlsx"
            output_path = os.path.join(output_dir, compare_info, output_filename)
            
            # 输出到 Excel（使用 Sheet1 和 Sheet2 作为 sheet 名称，与 smart_heatmap 兼容）
            with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                data_df.to_excel(writer, sheet_name='Sheet1', index=False)
                group_sample_df.to_excel(writer, sheet_name='Sheet2', index=False)
            
            logger.info(f'已输出 {compare_info} 的 {go_pathway_id} 热图数据到: {output_filename}')
            
            # 使用 smart_heatmap 画热图
            heatmap_filename = os.path.join(output_dir, compare_info, output_filename.replace('.xlsx', '.jpeg'))
            heatmap_result = smart_heatmap(
                output_path,
                heatmap_filename,
                annotation_col=2,
                cluster_rows=True,
                cluster_cols=False,
                scale="row"
            )
            if not heatmap_result:
                logger.error(f'{compare_info} 的 {go_pathway_id} smart_heatmap 画图失败')
            else:
                logger.info(f'已生成 {compare_info} 的 {go_pathway_id} 热图: {os.path.basename(heatmap_filename)}')
            


def deg_go_analysis(target_go_df, enrich_data_dir, output_dir):
    ontology_list = list(set(target_go_df['Ontology'].tolist()))
    # 循环每个差异数据文件
    for enrich_data in os.listdir(enrich_data_dir):
        if not enrich_data.endswith("_EnrichmentGO.xlsx") or enrich_data.startswith("~"):
            continue
        
        logger.info(f'\n====正在处理 {enrich_data}====')
        compare_info_regulation = enrich_data.replace('_EnrichmentGO.xlsx', '')  # 比对名称 + 上下调
        compare_info = compare_info_regulation.rsplit('_', 1)[0]  # 比对名称
        
        # 文件输出目录
        compare_output_dir = os.path.join(output_dir, compare_info_regulation)
        go_analysis_network_graph_dir = os.path.join(compare_output_dir, 'GO_analysis_network_graph')
        go_analysis_bar_plot_dir = os.path.join(compare_output_dir, 'GO_analysis_bar_plot')
        # go_deg_expression_data_dir = os.path.join(compare_output_dir, 'DEG_expression_data')
        dirs_to_make = [
            compare_output_dir,
            go_analysis_bar_plot_dir,
            go_analysis_network_graph_dir
        ]
        for each_dir in dirs_to_make:
            os.makedirs(each_dir, exist_ok=True)

        enrich_data_abspath = os.path.join(enrich_data_dir, enrich_data)
        enrich_go_df = load_table(enrich_data_abspath)
        enrich_go_df = enrich_go_df.drop(columns=['Ontology'])
         
        for ontology_name in ontology_list:
            logger.info(f"正在处理 {compare_info_regulation}_{ontology_name}")
            ontology_df = target_go_df[target_go_df['Ontology'] == ontology_name].copy()
            ontology_df = pd.merge(left=ontology_df, right=enrich_go_df, how='inner', left_on='GO_pathway_ID', right_on='ID')
            
            if ontology_df.shape[0] == 0:
                logger.warning(f"{compare_info_regulation}_{ontology_name}, 没有筛选出任何数据，不进行画图，跳过")
                continue
            elif ontology_df.shape[0] > 1:
                ontology_df['GO_pathway_ID'] = ontology_df['GO_pathway_ID'] + '_' + ontology_df['GO_pathway_def']
                ontology_df.drop(columns=['ID'], inplace=True)
                ontology_df.rename(columns={"GO_pathway_ID": "ID"}, inplace=True)
                ontology_df = ontology_df.sort_values(by=['ID'])
                output_excel_name = os.path.join(compare_output_dir, f'{compare_info_regulation}_{ontology_name}.xlsx')
                write_output_df(ontology_df, output_excel_name, index=False)
                
                add_summary_df = ontology_df.copy()
                add_summary_df.insert(0, 'Group', compare_info_regulation.rsplit('_', 1)[0])
                add_summary_df.insert(1, 'Regulation', compare_info_regulation.split('_')[-1])
                # output_summary_df_list.append(add_summary_df)
                
                jpeg_file_name = os.path.join(go_analysis_bar_plot_dir, f'{compare_info_regulation}_{ontology_name}_barplot.jpeg')
                # 运行 R barplot 脚本
                enrichment_barplot(output_excel_name, jpeg_file_name)
            else:
                logger.warning(f"{compare_info_regulation}_{ontology_name}, 只筛选出 1 条数据，不画 barplot 图")
            
            # 画 enrichnet 图
            enrichnet_file_name = os.path.join(go_analysis_network_graph_dir, f'{compare_info_regulation}_{ontology_name}_enrichnet.png')
            enrichnet_data = {'source': [], 'target': []}
            for index, row in ontology_df.iterrows():
                subontologys = [row['SubOntology']] * len(row['geneID'].split('/'))
                genes = row['geneID'].split('/')
                enrichnet_data['source'].extend(subontologys)
                enrichnet_data['target'].extend(genes)
            enrichnet_df = pd.DataFrame(enrichnet_data)
            enrichnet_df = enrichnet_df.drop_duplicates()
            gene_interaction_network_plot.draw_enrichnetplot(enrichnet_df, enrichnet_file_name)


def main():
    args = parse_input()
    samples_df = load_table(args.samples_described)
    comparison_df = load_table(args.compare)
    # target_go_df 预处理
    gene_go_df = load_table(args.genego, header=None, names=['GeneID', 'GO_pathway_ID', 'GO_pathway_def'], dtype={"GeneID": str})
    go_id_def_df = gene_go_df[['GO_pathway_ID', 'GO_pathway_def']].drop_duplicates()
    target_go_df = load_table(args.input)
    target_go_df = df_drop_element_side_space(target_go_df)
    target_go_df['GO_pathway_ID'] = target_go_df['GO_pathway_ID'].str.split('_').str[0]
    
    # 合并 target_go_df 和 gene_go_df
    num_before_merge = target_go_df.shape[0]
    merged_target_go_df = pd.merge(left=go_id_def_df, right=target_go_df, on='GO_pathway_ID', how='inner')
    num_after_merge = merged_target_go_df.shape[0]
    logger.info(f"GO分析: 合并前 target_go_df 行数为: {num_before_merge}, 合并后为: {num_after_merge}")
    if num_after_merge < num_before_merge:
        lost_go_df = target_go_df[~target_go_df['GO_pathway_ID'].isin(merged_target_go_df['GO_pathway_ID'])]
        lost_go_log_dir = os.path.join(args.output, 'log')
        os.makedirs(lost_go_log_dir, exist_ok=True)
        lost_go_file = os.path.join(lost_go_log_dir, 'target_go_lost_in_merge_gene_go.txt')
        write_output_df(lost_go_df, lost_go_file, index=False)
        logger.warning(f"有 {num_before_merge-num_after_merge} 个 GO_pathway_ID 在 gene_go_df 中未找到, 详细列表见 {lost_go_file}")
    target_go_df = merged_target_go_df
    target_go_df = df_replace_illegal_folder_chars(target_go_df, ['GO_pathway_def', 'Ontology', 'SubOntology'])
    
    
    logger.info(f'执行 go 分析')
    each_go_pathway_deg_data(target_go_df, gene_go_df, comparison_df, samples_df, args.deg_data_dir, os.path.join(args.output, '00_Each_GO_pathway_DEG_data'))
    each_go_pathway_deg_data_heatmap(target_go_df, gene_go_df, comparison_df, samples_df, args.deg_data_dir, os.path.join(args.output, '00_Each_GO_pathway_heatmap'), args.genesymbol)
    goid_enrich_summary(target_go_df, comparison_df, args.enrich_data_dir, os.path.join(args.output, 'Target_GO_Enrich_data_summary.xlsx'))
    deg_go_analysis(target_go_df, args.enrich_data_dir, args.output)
    
    
    if args.fpkm_matrix:
        fpkm_matrix_df = load_table(args.fpkm_matrix)
        output_dir = os.path.join(args.output, '00_Each_GO_pathway_heatmap')
        os.makedirs(output_dir, exist_ok=True)
        each_go_gene_heatmap(target_go_df, gene_go_df, fpkm_matrix_df, samples_df, output_dir, args.genesymbol)
    else:
        logger.warning('fpkm_matrix 为空，跳过输出每个 GO_ID 的相关基因的 Heatmap 图')

    if args.expression_data:
        expression_df = load_table(args.expression_data)
        output_dir = os.path.join(args.output, '00_Each_GO_pathway_expression_data')
        os.makedirs(output_dir, exist_ok=True)
        each_go_gene_expression(target_go_df, gene_go_df, expression_df, output_dir)
    else:
        logger.warning('expression_data 为空，跳过输出每个 GO_ID 的相关基因的表达量数据文件')
    
    logger.success('Done!')


if __name__ == '__main__':
    main()