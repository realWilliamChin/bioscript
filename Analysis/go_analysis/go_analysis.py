#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/03/27 14:01
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import subprocess
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Plot/'))
from Rscript import draw_barplot
import enrichnet_plt
if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    logger.critical("当前 Python 版本为:", sys.version)
    sys.exit(1)



def parse_input():
    argparser = argparse.ArgumentParser(description='')
    argparser.add_argument('-i', '--input', required=True, type=str, 
                           help="输入文件，文件中至少包含三列，GO_ID、Ontology、SubOntology")
    # argparser.add_argument('--goid', required=True, help='GOID 文件,swiss 注释出来的的 idNo_def.txt')
    argparser.add_argument('--genego', required=True, help='GOID 文件,swiss 注释出来的的 gene_go.txt')
    argparser.add_argument('--deg-data-dir', dest='deg_data_dir', type=str, required=True,
                           help='输入 DEG_data.txt 文件目录')
    argparser.add_argument('--enrich-data-dir', dest='enrich_data_dir', type=str, required=True,
                           help='输入 enrich.r 运行出来的结果文件夹')
    argparser.add_argument('-o', '--output', default=os.getcwd(),
                           help='输出目录，默认为当前目录，输出目录如果不存在则会尝试创建')
    
    args = argparser.parse_args()
    
    for f in [args.input, args.genego]:
        if not os.access(f, os.R_OK):
            logger.critical(f"输入文件 {f} 未找到或不可读")
            sys.exit(1)
    for d in [args.deg_data_dir, args.enrich_data_dir]:
        if not os.path.exists(d):
            logger.critical(f"输入目录 {d} 未找到")
            sys.exit(1)
    
    return args


def get_go_expression(ko, kegg_clean_file, expression_file):
    """对每个 go 相关的基因生成一个 go 相关基因的表达量表

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


def main():
    args = parse_input()

    if args.input.endswith('.txt'):
        df = pd.read_csv(args.input, sep='\t')
    elif args.input.endswith('.xlsx'):
        df = pd.read_excel(args.input, engine='openpyxl')
    elif args.input.endswith('.csv'):
        df = pd.read_csv(args.input)
    
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    df['ID'] = df['GO_ID'].str.split("_").str[0]  # 为了和 enrich.r 出来的文件的 ID 对应上，添加一列只包含 ID 的列， GO:0010111
    ontology_list = list(set(df['Ontology'].tolist()))
    go_id_list = list(set(df['ID'].tolist()))
    
    gene_go_df = pd.read_csv(args.genego, header=None, names=['GeneID', 'GO_ID'], sep='\t', dtype={"GeneID": str})
    
    output_summary_df_list = []
    # 循环每个差异数据文件
    for enrich_data in os.listdir(args.enrich_data_dir):
        if not enrich_data.endswith("_EnrichmentGO.xlsx") or enrich_data.startswith("~"):
            continue
        
        logger.info(f'====正在处理 {enrich_data}====')
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
                logger.warning(f'{each_dir} 文件夹已存在，输出文件将覆盖源文件（如果有）')
                
        enrich_data_abspath = os.path.join(args.enrich_data_dir, enrich_data)
        enrich_go_df = pd.read_excel(enrich_data_abspath, engine='openpyxl')
        enrich_go_df = enrich_go_df.drop(columns=['Ontology'])
        
        for ontology_name in ontology_list:
            logger.info(f"正在处理 {compare_info}_{ontology_name}")
            # gene_id_goid = {}
            ontology_df = df[df['Ontology'] == ontology_name].copy()
            ontology_df = pd.merge(left=ontology_df, right=enrich_go_df, how='inner', on='ID')
            ontology_df.drop(columns=['ID', 'GO_type'])
            
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

            ontology_df = ontology_df.sort_values(by=['GO_ID'])
            output_excel_name = os.path.join(compare_output_dir, f'{compare_info}_{ontology_name}.xlsx')
            ontology_df.to_excel(output_excel_name, index=False)
            
            add_summary_df = ontology_df.copy()
            add_summary_df.insert(0, 'Group', compare_info.rsplit('_', 1)[0])
            add_summary_df.insert(1, 'Regulation', compare_info.split('_')[-1])
            output_summary_df_list.append(add_summary_df)
            
            jpeg_file_name = os.path.join(go_analysis_bar_plot_dir, f'{compare_info}_{ontology_name}_barplot.jpeg')
            # 运行 R barplot 脚本
            if ontology_df.shape[0] > 1:
                draw_barplot(output_excel_name, jpeg_file_name)
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
            enrichnet_plt.draw_enrichnetplot(enrichnet_df, enrichnet_file_name)
        
        
        # 对每个 GO_ID 的相关基因，对上 DEG_data.txt 出一个表
        for go_id in go_id_list:
            compare_name = compare_info.replace('_Down', '').replace('_Up', '')
            compare_regulation = compare_info.split('_')[-1]
            # 因为循环的是 enrich data 不是 DEG_data，使用这种方式找到对应的 DEG_data.txt 文件
            deg_data_name = os.path.join(args.deg_data_dir, f"{compare_name}_DEG_data.txt")
            go_id_deg_data_filename = os.path.join(go_deg_expression_data_dir, f'{compare_name}_{go_id}_DEG_data.txt').replace(':', '_')
            
            logger.info(f'正在对 {compare_info} 的 {go_id} 找出相关基因')
            
            deg_data_df = pd.read_csv(deg_data_name, sep='\t', dtype={"GeneID": str})
            go_id_deg_data_df = deg_data_df[
                deg_data_df['regulation'].str.contains(compare_regulation) &
                deg_data_df['GeneID'].isin(gene_go_df[gene_go_df['GO_ID'].str.contains(go_id)]['GeneID'])
                ]
            if go_id_deg_data_df.shape[0] == 0:
                logger.warning(f'{compare_info} 的 {go_id} 未找到相关基因')
                continue
            go_id_deg_data_df.to_csv(go_id_deg_data_filename, sep='\t', index=False)
    
    output_summary_df = pd.concat(output_summary_df_list)
    output_summary_df.to_excel(os.path.join(args.output, 'Target_GO_analysis_summary.xlsx'), engine='openpyxl', index=False)
    logger.success('Done!')

if __name__ == '__main__':
    main()