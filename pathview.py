#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/29 17:19
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse
import subprocess


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
    kegg_clean_df = pd.read_csv(kegg_clean, sep='\t', usecols=[0, 4], names=['GeneID', 'KEGG_ID'])
    
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
    
    ko_df = pd.read_csv(ko_file, sep='\t', names=['Ko_Number'])
    kegg_gene_def = pd.read_csv(kegg_gene_def_file, sep='\t')
    kegg_pathway_df = pd.read_csv(kegg_pathway_file, sep='\t', names=['GeneID', 'Ko'])
    kegg_pathway_df['Ko_Number'] = kegg_pathway_df['Ko'].str.split(":").str[0]
    
    ko_gene_df = pd.merge(ko_df, kegg_pathway_df, how='inner', on='Ko_Number')
    ko_gene_df = pd.merge(ko_gene_df, kegg_gene_def, how='left', on='GeneID')
    
    ko_gene_df = ko_gene_df.drop(columns=['Ko_Number'])
    ko_gene_df.fillna('NA', inplace=True)
    
    if basicinfo is not None:
        basicinfo_df = pd.read_csv(basicinfo, sep='\t', usecols=[0,1,2,3])
        basicinfo_df_columns = basicinfo_df.columns.tolist()
        ko_gene_df = pd.merge(ko_gene_df, basicinfo_df, how='left', on='GeneID')
        output_columns_list = basicinfo_df_columns + ['Ko', 'KEGG_ID', 'Gene_shortname', 'EC_number', 'KEGG_def']
    else:
        output_columns_list = ['GeneID', 'Ko', 'KEGG_ID', 'Gene_shortname', 'EC_number', 'KEGG_def']
        
    ko_gene_df = ko_gene_df[output_columns_list]
    ko_gene_df.to_csv(output_prefix + '_ko_geneid_def.txt', sep='\t', index=False)
    
    
    ko_gene_df['regulation'] = 1
    ko_gene_df[['KEGG_ID', 'regulation']].to_csv(output_prefix + '_keggid_regulation.txt', sep='\t', index=False)


def passed_path(ko_file, output_prefix):
    ko_def_df = pd.read_csv(ko_file, sep='\t', names=['Ko_Number'])
    passed_path_file = '/home/colddata/qinqiang/script/Rscript/pathview/passed_path.txt'
    passed_path_df = pd.read_csv(passed_path_file, sep='\t', names=['Ko_Number', 'Ko_Def'])
    ko_def_df = passed_path_df[passed_path_df['Ko_Number'].isin(ko_def_df['Ko_Number'])]
    ko_def_df = ko_def_df.drop_duplicates(subset=['Ko_Number'])
    ko_def_df.to_csv(output_prefix + '_ko_passed_path.txt', sep='\t', index=False, header=False)


def get_ko_expression(ko, kegg_pathway_file, expression_file):
    if not os.path.exists(ko):
        sys.exit(f"Error: {ko} not found")
    if not os.path.exists(kegg_pathway_file):
        sys.exit(f"Error: {kegg_pathway_file} not found")
    if not os.path.exists(expression_file):
        sys.exit(f"Error: {expression_file} not found")
        
    ko_list = open(ko, 'r').read().strip().split('\n')
    print(f"正在输出每个 ko_pathway 相关的基因表达量表：{len(ko_list)} 个 ko_pathway")
    for ko_num in ko_list:
        kegg_pathway_df = pd.read_csv(kegg_pathway_file, sep='\t', names=['GeneID', 'Ko'])
        kegg_pathway_df['Ko_Number'] = kegg_pathway_df['Ko'].str.split(":").str[0]
        kegg_pathway_df = kegg_pathway_df[kegg_pathway_df['Ko_Number'] == ko_num]
        expression_df = pd.read_csv(expression_file, sep='\t')
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
        print("pathview 画图失败", ret)
    else:
        print("pathview 画图结束")


def draw_heatmap(datafile):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/heatmap/heatmap.r -f {datafile}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        print("heatmap 画图失败", ret)
    else:
        print("heatmap 画图结束")


def parse_input():
    argparser = argparse.ArgumentParser(description='Add definition to gene file')
    # 其他参数
    argparser.add_argument('--ko-list', dest="ko_list", type=str, required=True, help='[必须]ko_list file')
    argparser.add_argument('--kegg-pathway', dest="kegg_pathway", required=True, type=str, help='[必须]kegg.txt')
    argparser.add_argument('--draw-pathview', dest='draw_pathview', action='store_true',
                           help='[可选]draw pathview picture')
    
    # 针对其他
    group1 = argparser.add_argument_group('针对其他的使用参数')
    group1.add_argument('--kegg-genedef', dest="kegg_genedef", type=str, help='[必须]kegg_gene_def.txt')
    group1.add_argument('--basicinfo', type=str, help='[可选]basicinfo.txt')
    group1.add_argument('-o', '--output_prefix', type=str, help='output file prefix, 转录组不需要加这个参数')
    group1.add_argument('--expression', type=str, help='[可选]输入表达量文件, 输出每个 ko_pathway 相关的基因表达量表，\
        如果输入了这个参数，那么必须输入 --kegg-pathway 参数，输出文件名为 [ko]_expression.txt')
    
    # 针对转录组，de results 结果 up down 的 id list 画 pahtview 图
    group2 = argparser.add_argument_group('针对转录组的使用参数，如果跑转录组的项目，下面必须的参数必须输入')
    
    group2.add_argument('--kegg-clean', dest='kegg_clean', type=str, help="[必须]输入 kegg 注释出来的 'KEGG_clean 文件'")
    group2.add_argument('--deg-data-dir', dest='deg_data_dir', type=str, help='[必须]输入 DEG_data.txt 文件')
    group2.add_argument('-s', '--samplesinfo', type=str, help='[必须]输入样品信息文件，如果有则添加到输出文件中')
    
    # 添加了 deg-data 参数，下面这些无需使用
    # group2.add_argument('-c', '--compareinfo', type=str, help='[必须]输入比较信息文件，如果有则添加到输出文件中')
    # group2.add_argument('--fpkm', type=str, help='[可选]输入 fpkm 文件，如果有则添加到输出文件中')
    
    
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
    if args.kegg_genedef and args.kegg_pathway:
        if not os.path.exists(args.kegg_genedef):
            sys.exit(f"Error: {args.kegg_genedef} not found")
        if not os.path.exists(args.kegg_pathway):
            sys.exit(f"Error: {args.kegg_pathway} not found")
        
        other(args)
        
        # 筛选 passed_path.txt 中的 ko
        passed_path(args.ko_list, args.output_prefix)
        # 画 pathview 图
        if args.draw_pathview:
            print("正在画 pathview 图")
            draw_pathview(args.output_prefix + '_regulation.txt', args.output_prefix + '_ko_passed_path.txt')
            
        # 出表达量表
        if args.expression and args.kegg_pathway:
            get_ko_expression(args.ko_list, args.kegg_pathway, args.expression)
    
    # 转录组
    elif args.deg_data_dir:
        deg_data_list = os.listdir(args.deg_data_dir)
        kegg_pathway_df = pd.read_csv(args.kegg_pathway, sep='\t', names=['GeneID', 'Ko'])
        samples_described_df = pd.read_csv(args.samplesinfo, sep='\t', usecols=[0, 1], names=["group", "sample"])
        for deg_data_file in deg_data_list:
            deg_data_file = os.path.join(args.deg_data_dir, deg_data_file)
            deg_data_df = pd.read_csv(deg_data_file, sep='\t')
            
            treat_group = deg_data_df.iloc[0, 1]
            control_group = deg_data_df.iloc[0, 2]
            compare_info = treat_group + '_vs_' + control_group
            
            print(f"正在处理 {compare_info}")
            up_df = deg_data_df[deg_data_df['regulation'] == 'Up']['GeneID']
            down_df = deg_data_df[deg_data_df['regulation'] == 'Down']['GeneID']
            
            crt_group_samples = samples_described_df[samples_described_df['group'].isin([treat_group, control_group])]['sample'].to_list()
            fpkm_df = deg_data_df[["GeneID"] + [x + '_FPKM' for x in crt_group_samples]]
            # 把 fpkm_df 的所有列名带 _FPKM 的去掉
            fpkm_df.columns = ["GeneID"] + crt_group_samples
            
            treat_group_samples_name = samples_described_df[samples_described_df['group'].isin([treat_group,])]['sample'].to_list()
            control_group_samples_name = samples_described_df[samples_described_df['group'].isin([control_group,])]['sample'].to_list()
            print(f"crt_group_samples: {crt_group_samples}")
            print(f"treat_group_samples_name: {treat_group_samples_name}")
            print(f"control_group_samples_name: {control_group_samples_name}")
            
            heatmap_sheet2_df = pd.DataFrame(columns=['samples', 'group', 'colors'])
            for each_sample in crt_group_samples:
                if each_sample in treat_group_samples_name:
                    group = treat_group
                    color = args.treat_color
                else:
                    group = control_group
                    color = args.control_color
                heatmap_sheet2_df = heatmap_sheet2_df.append({'samples': each_sample, 'group': group, 'colors': color}, ignore_index=True)

            # 循环每个 ko number 筛选表达量画热图
            ko_list = open(args.ko_list, 'r').read().strip().split('\n')
            for ko_number in ko_list:
                ko_num_fpkm_expr_df = fpkm_df[fpkm_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['Ko'].str.contains(ko_number)]['GeneID'])]
                if ko_num_fpkm_expr_df.shape[0] == 0:
                    print(f"{ko_number} 在 {compare_info} 中相关的基因表达量表为空")
                    continue
                
                ko_num_fpkm_expr_df_file = f"{compare_info}_{ko_number}_fpkm_expression.xlsx"
                with pd.ExcelWriter(ko_num_fpkm_expr_df_file) as writer:
                    ko_num_fpkm_expr_df.to_excel(writer, index=False, sheet_name='Sheet1')
                    heatmap_sheet2_df.to_excel(writer, index=False, sheet_name='Sheet2')
                if ko_num_fpkm_expr_df.shape[0] > 1:
                    print(f"正在画 {compare_info} {ko_number} 相关基因表达量的热图")
                    draw_heatmap(ko_num_fpkm_expr_df_file)
                else:
                    print(f"{compare_info} 的 {ko_number} 中相关的基因表达量表只有一行，无法画热图")
                
                # 输出每个 ko 的 deg_data 文件
                ko_num_deg_data_file = f"{compare_info}_{ko_number}_DEG_data.txt"
                ko_num_deg_data_df = deg_data_df[deg_data_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['Ko'].str.contains(ko_number)]['GeneID'])]
                ko_num_deg_data_df.to_csv(ko_num_deg_data_file, sep='\t', index=False)
                
            # R pathview 画图准备文件 regulation
            print(f"正在准备 {compare_info} 的 regulation 文件")
            up_down_idlist_geneid2kid(down_df, up_df, compare_info, args.kegg_clean)
            
            # R pathview 画图准备文件 passed path
            print(f"正在筛选 {compare_info} 的 passed_path.txt")
            passed_path(args.ko_list, compare_info)
            
            # R pathview 画图
            if args.draw_pathview:
                print(f"正在画 {compare_info} 的 pathview 图")
                draw_pathview(f"{compare_info}_regulation.txt", f"{compare_info}_ko_passed_path.txt")
            else:
                print(f"没有画 {compare_info} 的 pathview 图")
            
            # 对每个比较组的文件进行整理
            os.mkdir(compare_info)
            os.system(f"mv {compare_info}_* {compare_info}")
            os.system(f"mv *.png {compare_info}")
  
        # compare_df = pd.read_csv(args.compareinfo, sep='\t')
        # samples_df = pd.read_csv(args.samplesinfo, sep='\t', usecols=[0, 1])
        # fpkm_df = pd.read_csv(args.fpkm, sep='\t')
        # kegg_pathway_df = pd.read_csv(args.kegg_pathway, sep='\t', names=['GeneID', 'Ko'])
        
        # # 循环每个比较组
        # for each_row in compare_df.itertuples():
        #     # 准备画 heatmap 的 excel 文件
        #     compare_info = (each_row[1] + '_vs_' + each_row[2])
        #     compares_samples_df = samples_df[samples_df['group'].isin([each_row[1], each_row[2]])]
        #     compares_samples_sample_lst = compares_samples_df['sample'].tolist()
            
        #     fpkm_expression_df = fpkm_df.reindex(columns=[fpkm_df.columns[0],] + compares_samples_sample_lst)
        #     fpkm_expression_df['sum'] = fpkm_expression_df[compares_samples_sample_lst].sum(axis=1)
        #     fpkm_expression_df = fpkm_expression_df[fpkm_expression_df['sum'] > 0]
        #     fpkm_expression_df = fpkm_expression_df.drop(columns=['sum'])
            
        #     sample_info = pd.DataFrame(columns=['samples', 'group', 'colors'])
        #     for each_sample in compares_samples_df.itertuples():
        #         sample = each_sample[2]
        #         group = each_sample[1]
        #         if group == each_row[1]:
        #             color = args.treat_color
        #         else:
        #             color = args.control_color
        #         sample_info = sample_info.append({'samples': sample, 'group': group, 'colors': color}, ignore_index=True)
            
        #     # 循环每个 ko number 筛选表达量画热图
        #     ko_list = open(args.ko_list, 'r').read().strip().split('\n')
        #     for ko_number in ko_list:
        #         ko_num_fpkm_expr_df = fpkm_expression_df[fpkm_expression_df['GeneID'].isin(kegg_pathway_df[kegg_pathway_df['Ko'].str.contains(ko_number)]['GeneID'])]
        #         if ko_num_fpkm_expr_df.shape[0] == 0:
        #             print(f"{ko_number} 在 {compare_info} 中相关的基因表达量表为空")
        #             continue
        #         elif ko_num_fpkm_expr_df.shape[0] == 1:
        #             print(f"{ko_number} 在 {compare_info} 中相关的基因表达量表只有一行，无法画热图")
        #             continue
        #         ko_num_fpkm_expr_df_file = f"{compare_info}_{ko_number}_fpkm_expression.xlsx"
        #         with pd.ExcelWriter(ko_num_fpkm_expr_df_file) as writer:
        #             ko_num_fpkm_expr_df.to_excel(writer, index=False, sheet_name='Sheet1')
        #             sample_info.to_excel(writer, index=False, sheet_name='Sheet2')
        #         if args.heatmap and args.treat_color and args.control_color and args.samplesinfo and args.fpkm:
        #             print(f"正在画 {compare_info} {ko_number} 相关基因表达量的热图")
        #             draw_heatmap(ko_num_fpkm_expr_df_file)
        #         else:
        #             print(f"参数错误，没有画 {ko_number} 的热图")
                    
        #     # 定义两个文件名
        #     down_idlist_file = compare_info + '_Down_ID.txt'
        #     up_idlist_file = compare_info + '_Up_ID.txt'

        #     # R pathview 画图准备文件 regulation
        #     print(f"正在准备 {compare_info} 的 regulation 文件")
        #     up_down_idlist_geneid2kid(down_idlist_file, up_idlist_file, compare_info, args.kegg_clean)
            
        #     # R pathview 画图准备文件 passed path
        #     print(f"正在筛选 {compare_info} 的 passed_path.txt")
        #     passed_path(args.ko_list, compare_info)
            
        #     # R pathview 画图
        #     if args.draw_pathview:
        #         print(f"正在画 {compare_info} 的 pathview 图")
        #         draw_pathview(f"{compare_info}_regulation.txt", f"{compare_info}_ko_passed_path.txt")
        #     else:
        #         print(f"没有画 {compare_info} 的 pathview 图")
            
        #     # 对每个比较组的文件进行整理
        #     os.mkdir(compare_info)
        #     os.system(f"mv {compare_info}_* {compare_info}")
        #     os.system(f"mv *.png {compare_info}")
    else:
        sys.exit("Error: Please check your input parameters")
    
    
    
    
    print('\nDone!\n')
    

if __name__ == '__main__':
    main()