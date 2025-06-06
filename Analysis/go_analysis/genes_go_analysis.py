#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/03/31 17:22
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from pathlib import Path
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Rscript/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/Plot/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from Rscript import enrichment_barplot
import enrichnet_plt
from load_input import load_table, write_output_df
if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    logger.critical("当前 Python 版本为:", sys.version)
    sys.exit(1)


def parse_input():
    p = argparse.ArgumentParser(
        description='GO 分析工具，用于分析转录组和 WGCNA 数据的 GO 富集结果',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument('-i', '--input', required=True, type=str,
                           help="输入文件，文件中至少包含三列，GO_ID、Ontology、SubOntology")
    p.add_argument('--genego', required=True, help='GOID 文件,swiss 注释出来的的 gene_go.txt')
    p.add_argument('--enrich-data-dir', dest='enrich_data_dir', type=str, required=True,
                           help='转录组 输入 enrich.r 运行出来的结果文件夹')
    p.add_argument('-o', '--output', default=os.getcwd(),
                           help='输出目录，默认为当前目录，输出目录如果不存在则会尝试创建')
    
    # 每个组的 goid 相关基因添加定义输出
    p.add_argument('--ids-dir', dest='ids_dir', type=str,
                           help='输入 cluster id 文件夹, os.basename 为 keyname')
    p.add_argument('--ids-header', dest='ids_header', action='store_true',
                           help='所有 id 文件是否有 header(GeneID)')
    p.add_argument('--fpkm-reads-merged', dest='fpkm_reads_merged', type=str,
                           help='输入 fpkm_reads_merged.txt 文件')
    
    args = p.parse_args()
    
    if args.ids_header:
        args.ids_header = 0
    else:
        args.ids_header = None
    
    return args


def geneids_go_analysis(target_go_df, enrich_data_dir, output_dir):
    ontology_list = list(set(target_go_df['Ontology'].tolist()))
    
    output_summary_df_list = []
    for enrich_data in [x for x in os.listdir(enrich_data_dir) if x.endswith("EnrichmentGO.xlsx")]:
        
        logger.info(f'====正在处理 {enrich_data}====')
        module_name = enrich_data.split('_')[0]
        
        # 文件输出目录
        module_output_dir = os.path.join(output_dir, f'{module_name}')
        go_analysis_network_graph_dir = os.path.join(module_output_dir, 'GO_analysis_network_graph')
        go_analysis_bar_plot_dir = os.path.join(module_output_dir, 'GO_analysis_bar_plot')
        for each_dir in module_output_dir, go_analysis_network_graph_dir, go_analysis_bar_plot_dir:
            if not os.path.exists(each_dir):
                os.mkdir(each_dir)
            else:
                logger.warning(f'{each_dir} Folder exists, output will overwrite, suggest delete it before running script')
        
        enrich_go_df = load_table(os.path.join(enrich_data_dir, enrich_data))
        enrich_go_df = enrich_go_df.drop(columns=['Ontology'])
        
        for ontology_name in ontology_list:
            logger.info(f"正在处理 {module_name}_{ontology_name}")
            ontology_df = target_go_df[target_go_df['Ontology'] == ontology_name].copy()
            ontology_df = pd.merge(left=ontology_df, right=enrich_go_df, how='inner', on='ID')
            
            if ontology_df.shape[0] == 0:
                logger.warning(f"{module_name}_{ontology_name}, 没有筛选出任何数据，不进行画图，跳过")
                continue
            
            ontology_df = ontology_df.sort_values(by=['ID'])
            output_excel_name = os.path.join(module_output_dir, f'{module_name}_{ontology_name}.xlsx')
            ontology_df.to_excel(output_excel_name, index=False)

            # 输出 summary 文件
            add_summary_df = ontology_df.copy()
            add_summary_df['Module_ID'] = module_name
            output_summary_df_list.append(add_summary_df)
            
            # 运行 R barplot 脚本
            jpeg_file_name = os.path.join(go_analysis_bar_plot_dir, f'{module_name}_{ontology_name}_barplot.jpeg')
            if ontology_df.shape[0] > 1:
                enrichment_barplot(output_excel_name, jpeg_file_name)
            else:
                logger.warning(f"{module_name}_{ontology_name}, 只筛选出 1 条数据，不画 barplot 图")
            
            # 画 enrichnet 图
            enrichnet_data = {'source': [], 'target': []}
            for index, row in ontology_df.iterrows():
                subontologys = [row['SubOntology']] * len(row['geneID'].split('/'))
                genes = row['geneID'].split('/')
                enrichnet_data['source'].extend(subontologys)
                enrichnet_data['target'].extend(genes)
            enrichnet_df = pd.DataFrame(enrichnet_data)
            enrichnet_df = enrichnet_df.drop_duplicates()
            enrichnet_plt.draw_enrichnetplot(
                enrichnet_df,
                os.path.join(go_analysis_network_graph_dir, f'{module_name}_{ontology_name}_enrichnet.png')
            )
        
    summary_df = pd.concat(output_summary_df_list)
    summary_df_col_lst = summary_df.columns.tolist()
    fixed_cols = ['ID', 'Description', 'SubOntology', 'Ontology', 'pvalue', 'Module_ID']
    summary_df_col_lst = fixed_cols + [x for x in summary_df_col_lst if x not in fixed_cols]
    summary_df = summary_df[summary_df_col_lst]
    write_output_df(
        summary_df,
        os.path.join(output_dir, 'Target_GO_analysis_summary.csv'),
        index=False
    )


# 对每个 GO_ID 的相关基因，对上 fpkm_reads_merged_df 出一个表
def each_groups_go_add_expression_def(ids_header, target_go_df, gene_go_df, fpkm_reads_merged_df, ids_dir, output_dir):
    go_id_list = list(set(target_go_df['ID'].tolist()))
    for module_file in os.listdir(ids_dir):
        module_name = module_file.split('_')[0]
        module_output_dir = os.path.join(output_dir, f'{module_name}')
        module_go_id_expression_data_dir = os.path.join(module_output_dir, 'GO_expression_data')
        if not os.path.exists(module_go_id_expression_data_dir):
            os.mkdir(module_go_id_expression_data_dir)
        module_df = load_table(
            os.path.join(ids_dir, module_file),
            header=ids_header,
            usecols=[0],
            names=['GeneID'],
            dtype={"GeneID": str}
        )

        module_df = pd.merge(left=module_df, right=gene_go_df, how='left', on='GeneID')
        module_df = pd.merge(
            left=module_df,
            right=target_go_df[['ID', 'SubOntology', 'Ontology']],
            left_on='GO_ID',
            right_on='ID',
            how='left'
        )
        module_df.dropna(subset=['GO_ID'], inplace=True)
        module_df.drop(columns=['ID'], inplace=True)
        for go_id in go_id_list:
            module_go_id_df = module_df[module_df['GO_ID'].str.contains(go_id)]
            if module_go_id_df.shape[0] == 0:
                logger.warning(f"{module_file}_{go_id}, 在 fpkm_reads_merged 中没有找到相关基因，跳过")
                continue
            module_go_id_df = pd.merge(left=module_go_id_df, right=fpkm_reads_merged_df, how='left', on='GeneID')
            write_output_df(
                module_go_id_df,
                os.path.join(
                    module_go_id_expression_data_dir,
                    f'{Path(module_name).stem}_{go_id.replace(":", "_")}_data_def.csv'
                ),
                index=False
            )


def main():
    args = parse_input()
    target_go_df = load_table(args.input)
    target_go_df = target_go_df[['GO_ID', 'SubOntology', 'Ontology']]
    target_go_df['GO_ID'] = target_go_df['GO_ID'].str.split("_").str[0]  # 为了和 enrich.r 出来的文件的 ID 对应上，添加一列只包含 ID 的列， GO:0010111
    target_go_df.rename(columns={'GO_ID': 'ID'}, inplace=True)
    # 数据清洗，字符串列，去除两边空格
    target_go_df = target_go_df.map(lambda x: x.strip() if isinstance(x, str) else x)
    
    gene_go_df = load_table(args.genego, header=None, names=['GeneID', 'GO_ID'], dtype={"GeneID": str})
    
    fpkm_reads_merged_df = load_table(args.fpkm_reads_merged, dtype={'GeneID': str})
    geneids_go_analysis(target_go_df, args.enrich_data_dir, args.output)
    
    each_groups_go_add_expression_def(args.ids_header, target_go_df, gene_go_df, fpkm_reads_merged_df, args.ids_dir, args.output)
    
    logger.success('Done!')


if __name__ == '__main__':
    main()