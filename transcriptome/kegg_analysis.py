#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/29 17:19
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse
import subprocess

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))


def kid_optmial(input_df):
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
    """
    可以读取 de_results_add_def.py 生成的 DEG_data.txt 文件
    把 regulation 变为 1,0,-1 对应 up,nosig,down
    生成的文件只包括 K_ID 和 regulation
    例如: K12345    0
    一个 K_ID 对应多个 GeneID, 有可能 up,down,nosig 都有
    对 K_ID 去重，哪个 regulation 多留下哪个，一样多的变为 0
    并生成一个中间文件，方便检查
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


def transcriptome_all_ko_gene_heatmap(ko_file, kegg_clean_file, fpkm_matrix_file, samples_file):
    # 对 expression fpkm matrix 文件只保留 fpkm 值
    expression_df = pd.read_csv(fpkm_matrix_file, sep='\t', dtype=str)
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
    gene_df = gene_df.drop_duplicates(subset='GeneID')
    
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
    
    draw_multigroup_heatmap(all_gene_ko_heatmap_filename)
    

def transcriptome(args):
    # 循环每个 ko number 筛选表达量画热图
    ko_df = pd.read_csv(args.ko_list, sep='\t')
    # 针对两种情况，一种只有 kolist，还有一种有 Ontology，就需要画另一个图了
    if 'Ontology' in ko_df.columns:
        ko_list = ko_df['KEGG_ID'].values.tolist()
        transcriptome_all_ko_gene_heatmap(args.ko_list, args.kegg_clean, args.expression, args.samplesinfo)
    else:
        ko_list = open(args.ko_list, 'r').read().strip().split('\n')
    
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
        
        print(f"\n====正在处理 {compare_info}====")
        up_df = deg_data_df[deg_data_df['regulation'] == 'Up']['GeneID']
        down_df = deg_data_df[deg_data_df['regulation'] == 'Down']['GeneID']
        
        crt_group_samples = samples_described_df[samples_described_df['group'].isin([treat_group, control_group])]['sample'].to_list()
        fpkm_df = deg_data_df[["GeneID"] + [x + '_FPKM' for x in crt_group_samples]]
        # 把 fpkm_df 的所有列名带 _FPKM 的去掉
        fpkm_df.columns = ["GeneID"] + crt_group_samples
        
        treat_group_samples_name = samples_described_df[samples_described_df['group'].isin([treat_group,])]['sample'].to_list()
        control_group_samples_name = samples_described_df[samples_described_df['group'].isin([control_group,])]['sample'].to_list()
        print(f"crt_group_samples: {crt_group_samples}")
        
        heatmap_sheet2_df = pd.DataFrame(columns=['samples', 'group', 'colors'])
        for each_sample in crt_group_samples:
            if each_sample in treat_group_samples_name:
                group = treat_group
                color = args.treat_color
            else:
                group = control_group
                color = args.control_color
            heatmap_sheet2_df = heatmap_sheet2_df.append({'samples': each_sample, 'group': group, 'colors': color}, ignore_index=True)

        for ko_number in ko_list:
            ko_num_fpkm_expr_df = fpkm_df[fpkm_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['Ko'].str.contains(ko_number)]['GeneID'])]
            
            if ko_num_fpkm_expr_df.shape[0] == 0:
                print(f"{ko_number} 在 {compare_info} 中相关的基因表达量表为空")
                continue
            
            ko_num_fpkm_expr_df_file = f"{compare_info}_{ko_number}_fpkm_expression.xlsx"
            ko_num_fpkm_expr_pic_name = f"{compare_info}_{ko_number}_heatmap.jpeg"
            with pd.ExcelWriter(ko_num_fpkm_expr_df_file) as writer:
                ko_num_fpkm_expr_df.to_excel(writer, index=False, sheet_name='Sheet1')
                heatmap_sheet2_df.to_excel(writer, index=False, sheet_name='Sheet2')
            if ko_num_fpkm_expr_df.shape[0] > 1:
                print(f"正在画 {compare_info} {ko_number} 相关基因表达量的热图")
                draw_twogroup_heatmap(ko_num_fpkm_expr_df_file, ko_num_fpkm_expr_pic_name)
            else:
                print(f"{compare_info} 的 {ko_number} 中相关的基因表达量表只有一行，无法画热图")
            
            # 输出每个 ko 的 deg_data 文件
            ko_num_deg_data_file = f"{compare_info}_{ko_number}_DEG_data.txt"
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
            print(f"正在画 {compare_info} 的 pathview 图")
            draw_pathview(f"{compare_info}_regulation.txt", f"{compare_info}_ko_passed_path.txt")
        
        # 对每个比较组的文件进行整理
        os.mkdir(compare_info)
        os.system(f"mv {compare_info}_* {compare_info}")
        os.system(f"mv *.png {compare_info}")


def passed_path(ko_file, output_prefix):
    ko_def_df = pd.read_csv(ko_file, sep='\t', names=['KEGG_ID'], dtype=str, usecols=[0], header=0)
    passed_path_file = '/home/colddata/qinqiang/script/Rscript/pathview/passed_path.txt'
    passed_path_df = pd.read_csv(passed_path_file, sep='\t', names=['KEGG_ID', 'Ko_Def'], dtype=str)
    ko_def_df = passed_path_df[passed_path_df['KEGG_ID'].isin(ko_def_df['KEGG_ID'])]
    ko_def_df = ko_def_df.drop_duplicates(subset=['KEGG_ID'])
    ko_def_df.to_csv(output_prefix + '_ko_passed_path.txt', sep='\t', index=False, header=False)


def get_ko_expression(ko, kegg_clean_file, expression_file):
    if not os.path.exists(ko):
        sys.exit(f"Error: {ko} not found")
    if not os.path.exists(kegg_clean_file):
        sys.exit(f"Error: {kegg_clean_file} not found")
    if not os.path.exists(expression_file):
        sys.exit(f"Error: {expression_file} not found")
        
    ko_list = open(ko, 'r').read().strip().split('\n')
    print(f"正在输出每个 ko_pathway 相关的基因表达量表：{len(ko_list)} 个 ko_pathway")
    for ko_num in ko_list:
        kegg_pathway_df = pd.read_csv(kegg_clean_file, sep='\t', names=['GeneID', 'Ko'], usecols=[0, 1], dtype=str)
        kegg_pathway_df['Ko_Number'] = kegg_pathway_df['Ko'].str.split(":").str[0]
        kegg_pathway_df = kegg_pathway_df[kegg_pathway_df['Ko_Number'] == ko_num]
        expression_df = pd.read_csv(expression_file, sep='\t', dtype=str)
        expression_df = expression_df[expression_df['GeneID'].isin(kegg_pathway_df['GeneID'])]
        if expression_df.shape[0] > 0:
            expression_df.to_csv(ko_num + '_expression.txt', sep='\t', index=False)
            print(f"正在输出 {ko_num} 相关的基因表达量表，有 {expression_df.shape[0]} 个基因")
        else:
            print(f"{ko_num} 相关的基因表达量表为空")


def draw_pathview(regulation, passed_path):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/pathview/pathview.R \
        -r {regulation} \
        -p {passed_path}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        print(f"{regulation} draw_pathview 画图失败")
        print(f"标准输出：{ret.stdout.decode()}")
        print(f"标准错误: {ret.stderr.decode()}")


def draw_multigroup_heatmap(datafile):
    """输入文件应该是 3 个 Sheet
    Sheet 1: fpkm value
    ID sampleA-A sampleA-B sampleA-C sampleB-A sampleB-B smapleB-C ...
    geneA fpkm_value ...
    geneB fpkm_value ...
    ...
    
    Sheet 2: sample_annotation
    sample group
    sampleA sampleA-A
    sampleA sampleA-B
    ...
    
    Sheet 3: gene_annotation
    gene Ontology
    geneA groupA
    geneB groupB
    genec gorupB
    ...
    """
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/heatmap/heatmap_multigroup.r -f {datafile}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        print(f"{datafile} multigroup_heatmap 画图失败")
        print(f"标准输出：{ret.stdout.decode()}")
        print(f"标准错误: {ret.stderr.decode()}")


def draw_twogroup_heatmap(datafile, output_file):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/heatmap/heatmap_twogroup.r -f {datafile} -o {output_file}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        print(f"{datafile} twogroup_heatmap 画图失败")
        print(f"标准输出：{ret.stdout.decode()}")
        print(f"标准错误: {ret.stderr.decode()}")


def anova_analysis(datafile, samples_file, output_file):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/anova/anova.r -f {datafile} -s {samples_file} -o {output_file}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        print(f"{datafile} twogroup_heatmap 画图失败")
        print(f"标准输出：{ret.stdout.decode()}")
        print(f"标准错误: {ret.stderr.decode()}")


def parse_input():
    argparser = argparse.ArgumentParser(description='Add definition to gene file')
    # 其他参数
    argparser.add_argument('--ko-list', dest="ko_list", type=str, required=True,
                           help='[必须]ko_list file, 必须要有列名，KEGG_ID 和 Ontology，如果运行脚本类型是其他，只需要 KEGG_ID 就行')
    argparser.add_argument('--kegg-clean', dest="kegg_clean", required=True, type=str, help='[必须]输入 kegg 注释出来的 KEGG_clean 文件')
    argparser.add_argument('--draw-pathview', dest='draw_pathview', action='store_true',
                           help='[可选]draw pathview picture')
    
    # 针对其他
    group1 = argparser.add_argument_group('针对其他的使用参数')
    group1.add_argument('--kegg-genedef', dest="kegg_genedef", type=str, help='[必须]kegg_gene_def.txt')
    group1.add_argument('--basicinfo', type=str, help='[可选]basicinfo.txt')
    group1.add_argument('-o', '--output_prefix', type=str, help='output file prefix, 转录组不需要加这个参数')
    group1.add_argument('--expression', type=str,
                        help='[可选]输入表达量（fpkm_and_reads_matrix.txt）文件, 输出每个 ko_pathway 相关的基因表达量表')
    
    # 针对转录组，de results 结果 up down 的 id list 画 pahtview 图
    group2 = argparser.add_argument_group('针对转录组的使用参数，如果跑转录组的项目，下面必须的参数必须输入')
    
    # 2024_04_25 使用 --kegg-clean 一个文件就可以了，少写一些参数
    # group2.add_argument('--kegg-clean', dest='kegg_clean', type=str, help="[必须]输入 kegg 注释出来的 'KEGG_clean 文件'")
    group2.add_argument('--deg-data-dir', dest='deg_data_dir', type=str, help='[必须]输入 DEG_data.txt 文件')
    group2.add_argument('-s', '--samplesinfo', type=str, help='[必须]输入样品信息文件，如果有则添加到输出文件中')
    
    # 添加了 deg-data 参数，下面这些无需使用
    # group2.add_argument('-c', '--compareinfo', type=str, help='[必须]输入比较信息文件，如果有则添加到输出文件中')
    # group2.add_argument('--fpkm', type=str, help='[可选]输入 fpkm 文件，如果有则添加到输出文件中')
    
    
    # 画热图添加的参数    
    group2.add_argument('--heatmap', action='store_true',
                        help='[可选]整理表达量文件，输出差异基因的热图')
    # # 2024_04_11 不用弄颜色了，自动分配了
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
            sys.exit(f"Error: {args.kegg_genedef} not found")
        if not os.path.exists(args.kegg_clean):
            sys.exit(f"Error: {args.kegg_clean} not found")
        
        other(args)
        
        # 筛选 passed_path.txt 中的 ko
        passed_path(args.ko_list, args.output_prefix)
        # 画 pathview 图
        if args.draw_pathview:
            print("正在画 pathview 图")
            draw_pathview(args.output_prefix + '_regulation.txt', args.output_prefix + '_ko_passed_path.txt')
            
        # 出表达量表
        if args.expression and args.kegg_clean:
            get_ko_expression(args.ko_list, args.kegg_clean, args.expression)
    
    print('\nDone!\n')
    

if __name__ == '__main__':
    main()