#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/03/27 14:01
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import subprocess

import enrichnet_plt


def parse_input():
    argparser = argparse.ArgumentParser(description='')
    argparser.add_argument('-i', '--input', required=True, type=str, 
                           help="输入文件，文件中至少包含三列，Ontology、Ontology、SubOntology")
    # argparser.add_argument('--goid', required=True, help='GOID 文件,swiss 注释出来的的 idNo_def.txt')
    # argparser.add_argument('--deg-data-dir', dest='deg_data_dir', type=str, required=True,
    #                        help='输入 DEG_data.txt 文件目录')
    argparser.add_argument('--enrich-data-dir', dest='enrich_data_dir', type=str, required=True,
                           help='输入 enrich.r 运行出来的结果文件夹')
    argparser.add_argument('-o', '--output', default='.',
                           help='输出目录，默认为当前目录，输出目录如果不存在则会尝试创建')
    
    args = argparser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"输入文件 {args.input} 不存在")
        sys.exit(1)
    
    return args


def draw_barplot(input_file, output_file):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/barplot/barplot.R \
        -f {input_file} \
        -o {output_file}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        print("barplot 画图失败")
        print(ret.stdout.decode())
        print(ret.stdout.decode())
    else:
        print("barplot 画图结束")


def main():
    args = parse_input()

    if args.input.endswith('.txt'):
        df = pd.read_csv(args.input, sep='\t')
    elif args.input.endswith('.xlsx'):
        df = pd.read_excel(args.input, engine='openpyxl')
    elif args.input.endswith('.csv'):
        df = pd.read_csv(args.input)
    
    df['ID'] = df['GO_ID'].str.split("_").str[0]  # 为了和 enrich.r 出来的文件的 ID 对应上，添加一列只包含 ID 的列， GO:0010111
    ontology_list = list(set(df['Ontology'].tolist()))
    
    # 循环每个差异数据文件
    for enrich_data in os.listdir(args.enrich_data_dir):
        if not enrich_data.endswith("_EnrichmentGO.xlsx") or enrich_data.startswith("~"):
            continue
        print(f'正在处理 {enrich_data}')
        enrich_data_abspath = os.path.join(args.enrich_data_dir, enrich_data)
        compare_name = enrich_data.replace('_EnrichmentGO.xlsx', '')  # 比对名称
        enrich_go_df = pd.read_excel(enrich_data_abspath, engine='openpyxl')
        enrich_go_df = enrich_go_df.drop(columns=['Ontology'])
        
        for ontology_name in ontology_list:
            print(f"正在处理 {compare_name}_{ontology_name}")
            # gene_id_goid = {}
            ontology_df = df[df['Ontology'] == ontology_name].copy()
            ontology_df = pd.merge(left=ontology_df, right=enrich_go_df, how='inner', on='ID')
            ontology_df.drop(columns=['ID', 'GO_type'])
            
            if ontology_df.shape[0] == 0:
                print(f"{compare_name}_{ontology_name}, 没有筛选出任何数据，不进行画图，跳过")
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
            output_excel_name = os.path.join(args.output, f'{compare_name}_{ontology_name}.xlsx')
            ontology_df.to_excel(output_excel_name, index=False)
            
            jpeg_file_name = os.path.join(args.output, f'{compare_name}_{ontology_name}_barplot.jpeg')
            # 运行 R barplot 脚本
            if ontology_df.shape[0] > 1:
                draw_barplot(output_excel_name, jpeg_file_name)
            else:
                print(f"{compare_name}_{ontology_name}, 只筛选出 1 条数据，不画 barplot 图")
            
            # 画 enrichnet 图
            enrichnet_file_name = os.path.join(args.output, f'{compare_name}_{ontology_name}_enrichnet.png')
            enrichnet_data = {'source': [], 'target': []}
            for index, row in ontology_df.iterrows():
                subontologys = [row['SubOntology']] * len(row['geneID'].split('/'))
                genes = row['geneID'].split('/')
                enrichnet_data['source'].extend(subontologys)
                enrichnet_data['target'].extend(genes)
            enrichnet_df = pd.DataFrame(enrichnet_data)
            enrichnet_df = enrichnet_df.drop_duplicates()
            enrichnet_plt.draw_enrichnetplot(enrichnet_df, enrichnet_file_name)
            
        # 对每个比较组的文件进行整理
        os.mkdir(os.path.join(args.output, compare_name))
        os.system(f"mv {os.path.join(args.output, compare_name)}_* {os.path.join(args.output, compare_name)}")


if __name__ == '__main__':
    main()