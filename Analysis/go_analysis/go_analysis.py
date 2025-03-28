#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/03/27 14:01
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import subprocess
from pathlib import Path
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Plot/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from Rscript import draw_barplot
import enrichnet_plt
from load_input import load_table, write_output_df
if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    logger.critical("当前 Python 版本为:", sys.version)
    sys.exit(1)


def parse_input():
    argparser = argparse.ArgumentParser(
        description='GO 分析工具，用于分析转录组和 WGCNA 数据的 GO 富集结果',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    argparser.add_argument('-i', '--input', required=True, type=str, 
                           help="输入文件，文件中至少包含三列，GO_ID、Ontology、SubOntology")
    argparser.add_argument('--genego', required=True, help='GOID 文件,swiss 注释出来的的 gene_go.txt')
    argparser.add_argument('--enrich-data-dir', dest='enrich_data_dir', type=str,
                           help='转录组 输入 enrich.r 运行出来的结果文件夹')
    
    subparsers = argparser.add_subparsers(dest='command', help='子命令帮助')
    
    # Transcriptome 输入数据
    transcriptome_parser = subparsers.add_parser('transcriptome', 
                                                help='转录组输入数据',
                                                description='处理转录组数据的 GO 分析')
    transcriptome_parser.add_argument('--deg-data-dir', dest='deg_data_dir', type=str,
                           help='转录组 输入 DEG_data.txt 文件目录')
    
    # WGCNA 输入数据 color module geneid
    wgcna_parser = subparsers.add_parser('wgcna', 
                                        help='WGCNA 输入数据',
                                        description='处理 WGCNA 数据的 GO 分析')
    wgcna_parser.add_argument('--color-module-dir', dest='color_module_dir', type=str,
                           help='WGCNA 输入 color module geneid 文件夹')
    wgcna_parser.add_argument('--fpkm-reads-merged', dest='fpkm_reads_merged', type=str,
                           help='WGCNA 输入 fpkm_reads_merged.txt 文件')

    # GO 输入数据
    go_parser = subparsers.add_parser('go', 
                                     help='GO 输入数据',
                                     description='处理 GO 数据的分析')
    go_parser.add_argument('--cluster-id-dir', dest='cluster_id_dir', type=str,
                           help='GO 输入 cluster id 文件夹')
    go_parser.add_argument('--fpkm-reads-merged', dest='fpkm_reads_merged', type=str,
                           help='GO 输入 fpkm_reads_merged.txt 文件')
                           
    argparser.add_argument('-o', '--output', default=os.getcwd(),
                           help='输出目录，默认为当前目录，输出目录如果不存在则会尝试创建')
    
    args = argparser.parse_args()
    
    # 如果没有子命令，打印帮助信息
    if args.command is None:
        argparser.print_help()
        sys.exit(0)
    elif args.command == 'transcriptome':
        if args.deg_data_dir is None:
            logger.critical("转录组输入数据，需要指定 deg-data-dir 参数")
            sys.exit(1)
    elif args.command == 'wgcna':
        if args.color_module_dir is None:
            logger.critical("WGCNA 输入数据，需要指定 color-module-dir 参数")
            sys.exit(1)
    elif args.command == 'go':
        if args.cluster_id_dir is None:
            logger.critical("GO 输入数据，需要指定 cluster-id-dir 参数")
            sys.exit(1)
    
    # Data check
    for f in [args.input, args.genego]:
        if not os.access(f, os.R_OK):
            logger.critical(f"输入文件 {f} 未找到或不可读")
            sys.exit(1)
    
    return args


def transcriptome_go_analysis(args, df, gene_go_df, go_id_list, ontology_list):
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
                logger.warning(f'{each_dir} 文件夹已存在，输出文件将覆盖源文件')
                
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


def wgcna_go_analysis(enrich_data_dir, color_module_dir, fpkm_reads_merged_df,
                      df, gene_go_df, go_id_list, ontology_list, output_dir):
    # 循环 WGCNA color module geneid 文件夹
    output_summary_df_list = []
    for enrich_data in os.listdir(enrich_data_dir):
        if not enrich_data.endswith("_EnrichmentGO.xlsx") or enrich_data.startswith("~"):
            continue
        
        logger.info(f'====正在处理 {enrich_data}====')
        color_module_name = enrich_data.replace('_EnrichmentGO.xlsx', '')  # 比对名称
        
        # 文件输出目录
        color_module_output_dir = os.path.join(output_dir, f'{color_module_name}_color_module')
        go_analysis_network_graph_dir = os.path.join(color_module_output_dir, 'GO_analysis_network_graph')
        go_analysis_bar_plot_dir = os.path.join(color_module_output_dir, 'GO_analysis_bar_plot')
        go_color_module_expression_data_dir = os.path.join(color_module_output_dir, 'Color_module_GO_expression_data')
        for each_dir in color_module_output_dir, go_analysis_network_graph_dir, go_analysis_bar_plot_dir, go_color_module_expression_data_dir:
            if not os.path.exists(each_dir):
                os.mkdir(each_dir)
            else:
                logger.warning(f'{each_dir} 文件夹已存在，输出文件将覆盖源文件')
        
        enrich_go_df = load_table(os.path.join(enrich_data_dir, enrich_data))
        enrich_go_df = enrich_go_df.drop(columns=['Ontology'])
        
        for ontology_name in ontology_list:
            logger.info(f"正在处理 {color_module_name}_{ontology_name}")
            ontology_df = df[df['Ontology'] == ontology_name].copy()
            ontology_df = pd.merge(left=ontology_df, right=enrich_go_df, how='inner', on='ID')
            ontology_df.drop(columns=['ID', 'GO_type'])
            
            if ontology_df.shape[0] == 0:
                logger.warning(f"{color_module_name}_{ontology_name}, 没有筛选出任何数据，不进行画图，跳过")
                continue
            
            ontology_df = ontology_df.sort_values(by=['GO_ID'])
            output_excel_name = os.path.join(color_module_output_dir, f'{color_module_name}_{ontology_name}.xlsx')
            ontology_df.to_excel(output_excel_name, index=False)

            # 输出 summary 文件
            add_summary_df = ontology_df.copy()
            add_summary_df.insert(0, 'Group', color_module_name)
            output_summary_df_list.append(add_summary_df)
            
            # 运行 R barplot 脚本
            jpeg_file_name = os.path.join(go_analysis_bar_plot_dir, f'{color_module_name}_{ontology_name}_barplot.jpeg')
            if ontology_df.shape[0] > 1:
                draw_barplot(output_excel_name, jpeg_file_name)
            else:
                logger.warning(f"{color_module_name}_{ontology_name}, 只筛选出 1 条数据，不画 barplot 图")
            
            # 画 enrichnet 图
            enrichnet_file_name = os.path.join(go_analysis_network_graph_dir, f'{color_module_name}_{ontology_name}_enrichnet.png')
            enrichnet_data = {'source': [], 'target': []}
            for index, row in ontology_df.iterrows():
                subontologys = [row['SubOntology']] * len(row['geneID'].split('/'))
                genes = row['geneID'].split('/')
                enrichnet_data['source'].extend(subontologys)
                enrichnet_data['target'].extend(genes)
            enrichnet_df = pd.DataFrame(enrichnet_data)
            enrichnet_df = enrichnet_df.drop_duplicates()
            enrichnet_plt.draw_enrichnetplot(enrichnet_df, enrichnet_file_name)

        # 对每个 GO_ID 的相关基因，对上 kns.txt 出一个表
        for go_id in go_id_list:
            for color_module_file in os.listdir(color_module_dir):
                color_module_df = load_table(
                    os.path.join(color_module_dir, color_module_file),
                    header=None,
                    usecols=[0],
                    names=['GeneID'],
                    dtype={"GeneID": str}
                )
                color_module_df = pd.merge(left=color_module_df, right=gene_go_df, how='left', on='GeneID')
                color_module_df.dropna(subset=['GO_ID'], inplace=True)
                color_module_df = color_module_df[color_module_df['GO_ID'].str.contains(go_id)]
                if color_module_df.shape[0] == 0:
                    logger.warning(f"{color_module_file}_{go_id}, 在 fpkm_reads_merged 中没有找到相关基因，跳过")
                    continue
                color_module_df = pd.merge(left=color_module_df, right=fpkm_reads_merged_df, how='left', on='GeneID')
                write_output_df(
                    color_module_df,
                    os.path.join(
                        go_color_module_expression_data_dir,
                        f'{Path(color_module_file).stem}_{go_id.replace(":", "_")}_color_module_expression.csv'
                    ),
                    index=False
                )
        
    output_summary_df = pd.concat(output_summary_df_list)
    write_output_df(
        output_summary_df,
        os.path.join(output_dir, 'Target_GO_analysis_summary.xlsx')
    )


def timelapse_go_analysis(enrich_data_dir, cluster_id_dir, fpkm_reads_merged_df,
                      df, gene_go_df, go_id_list, ontology_list, output_dir):
    # 循环 WGCNA color module geneid 文件夹
    output_summary_df_list = []
    for enrich_data in os.listdir(enrich_data_dir):
        if not enrich_data.endswith("_EnrichmentGO.xlsx") or enrich_data.startswith("~"):
            continue
        
        logger.info(f'====正在处理 {enrich_data}====')
        cluster_module_name = enrich_data.replace('_EnrichmentGO.xlsx', '')  # 比对名称
        
        # 文件输出目录
        cluster_module_output_dir = os.path.join(output_dir, f'{cluster_module_name}')
        go_analysis_network_graph_dir = os.path.join(cluster_module_output_dir, 'GO_analysis_network_graph')
        go_analysis_bar_plot_dir = os.path.join(cluster_module_output_dir, 'GO_analysis_bar_plot')
        go_cluster_id_expression_data_dir = os.path.join(cluster_module_output_dir, 'Cluster_GO_expression_data')
        for each_dir in cluster_module_output_dir, go_analysis_network_graph_dir, go_analysis_bar_plot_dir, go_cluster_id_expression_data_dir:
            if not os.path.exists(each_dir):
                os.mkdir(each_dir)
            else:
                logger.warning(f'{each_dir} 文件夹已存在，输出文件将覆盖源文件')
        
        enrich_go_df = load_table(os.path.join(enrich_data_dir, enrich_data))
        enrich_go_df = enrich_go_df.drop(columns=['Ontology'])
        
        for ontology_name in ontology_list:
            logger.info(f"正在处理 {cluster_module_name}_{ontology_name}")
            ontology_df = df[df['Ontology'] == ontology_name].copy()
            ontology_df = pd.merge(left=ontology_df, right=enrich_go_df, how='inner', on='ID')
            ontology_df.drop(columns=['ID', 'GO_type'])
            
            if ontology_df.shape[0] == 0:
                logger.warning(f"{cluster_module_name}_{ontology_name}, 没有筛选出任何数据，不进行画图，跳过")
                continue
            
            ontology_df = ontology_df.sort_values(by=['GO_ID'])
            output_excel_name = os.path.join(cluster_module_output_dir, f'{cluster_module_name}_{ontology_name}.xlsx')
            ontology_df.to_excel(output_excel_name, index=False)

            # 输出 summary 文件
            add_summary_df = ontology_df.copy()
            add_summary_df.insert(0, 'Group', cluster_module_name)
            output_summary_df_list.append(add_summary_df)
            
            # 运行 R barplot 脚本
            jpeg_file_name = os.path.join(go_analysis_bar_plot_dir, f'{cluster_module_name}_{ontology_name}_barplot.jpeg')
            if ontology_df.shape[0] > 1:
                draw_barplot(output_excel_name, jpeg_file_name)
            else:
                logger.warning(f"{cluster_module_name}_{ontology_name}, 只筛选出 1 条数据，不画 barplot 图")
            
            # 画 enrichnet 图
            enrichnet_file_name = os.path.join(go_analysis_network_graph_dir, f'{cluster_module_name}_{ontology_name}_enrichnet.png')
            enrichnet_data = {'source': [], 'target': []}
            for index, row in ontology_df.iterrows():
                subontologys = [row['SubOntology']] * len(row['geneID'].split('/'))
                genes = row['geneID'].split('/')
                enrichnet_data['source'].extend(subontologys)
                enrichnet_data['target'].extend(genes)
            enrichnet_df = pd.DataFrame(enrichnet_data)
            enrichnet_df = enrichnet_df.drop_duplicates()
            enrichnet_plt.draw_enrichnetplot(enrichnet_df, enrichnet_file_name)

        # 对每个 GO_ID 的相关基因，对上 kns.txt 出一个表
        for go_id in go_id_list:
            for cluster_id_file in os.listdir(cluster_id_dir):
                cluster_id_df = load_table(
                    os.path.join(cluster_id_dir, cluster_id_file),
                    header=None,
                    usecols=[0],
                    names=['GeneID'],
                    dtype={"GeneID": str}
                )
                cluster_id_df = pd.merge(left=cluster_id_df, right=gene_go_df, how='left', on='GeneID')
                cluster_id_df.dropna(subset=['GO_ID'], inplace=True)
                cluster_id_df = cluster_id_df[cluster_id_df['GO_ID'].str.contains(go_id)]
                if cluster_id_df.shape[0] == 0:
                    logger.warning(f"{cluster_id_file}_{go_id}, 在 fpkm_reads_merged 中没有找到相关基因，跳过")
                    continue
                cluster_id_df = pd.merge(left=cluster_id_df, right=fpkm_reads_merged_df, how='left', on='GeneID')
                write_output_df(
                    cluster_id_df,
                    os.path.join(
                        go_cluster_id_expression_data_dir,
                        f'{Path(cluster_id_file).stem}_{go_id.replace(":", "_")}_expression.csv'
                    ),
                    index=False
                )
        
    output_summary_df = pd.concat(output_summary_df_list)
    write_output_df(
        output_summary_df,
        os.path.join(output_dir, 'Target_GO_analysis_summary.xlsx')
    )


def main():
    args = parse_input()
    df = load_table(args.input)

    # 数据清洗，字符串列，去除两边空格
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    df['ID'] = df['GO_ID'].str.split("_").str[0]  # 为了和 enrich.r 出来的文件的 ID 对应上，添加一列只包含 ID 的列， GO:0010111
    ontology_list = list(set(df['Ontology'].tolist()))
    go_id_list = list(set(df['ID'].tolist()))
    
    gene_go_df = load_table(args.genego, header=None, names=['GeneID', 'GO_ID'], dtype={"GeneID": str})
    
    if args.command == 'transcriptome':
        transcriptome_go_analysis(args, df, gene_go_df, go_id_list, ontology_list)
    elif args.command == 'wgcna':
        fpkm_reads_merged_df = load_table(args.fpkm_reads_merged)
        wgcna_go_analysis(args.enrich_data_dir, args.color_module_dir, fpkm_reads_merged_df,
                         df, gene_go_df, go_id_list, ontology_list, args.output)
    elif args.command == 'go':
        fpkm_reads_merged_df = load_table(args.fpkm_reads_merged)
        timelapse_go_analysis(args.enrich_data_dir, args.cluster_id_dir, fpkm_reads_merged_df,
                              df, gene_go_df, go_id_list, ontology_list, args.output)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()