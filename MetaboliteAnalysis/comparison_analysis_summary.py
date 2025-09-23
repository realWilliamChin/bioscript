import os, sys
import argparse
import pandas as pd
import numpy as np
import glob
import subprocess
from functools import reduce

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='输入目录')
    parser.add_argument('-d', help='definition_dataframe')
    
    return parser.parse_args()


def prepare_kegg_data(vip_df, compare_group_name):
    """准备KEGG富集分析数据和目录"""
    
    # 创建目录结构
    ko01000_enrich_dir = '06_KEGG_ko01100_Enrich'
    kegg_all_enrich_dir = '05_KEGG_all_Enrich'
    
    os.makedirs(ko01000_enrich_dir, exist_ok=True)
    os.makedirs(kegg_all_enrich_dir, exist_ok=True)
    os.makedirs(os.path.join(ko01000_enrich_dir, compare_group_name), exist_ok=True)
    os.makedirs(os.path.join(kegg_all_enrich_dir, compare_group_name), exist_ok=True)
    
    # 设置输出路径
    enrich_output_prefix = os.path.join(ko01000_enrich_dir, compare_group_name, compare_group_name)
    kegg_enrich_output_prefix = os.path.join(kegg_all_enrich_dir, compare_group_name, compare_group_name)
    
    # 获取上调/下调的KEGG ID
    up_df = vip_df[vip_df['regulation'] == "Up"]['KEGG'].dropna()
    down_df = vip_df[vip_df['regulation'] == 'Down']['KEGG'].dropna()
    
    # 保存KEGG ID到文件
    up_df_filename = f"{enrich_output_prefix}_Up_Compound_ID.txt"
    down_df_filename = f"{enrich_output_prefix}_Down_Compound_ID.txt"
    
    with open(up_df_filename, 'w') as f:
        for kegg_id in up_df:
            f.write(f"{kegg_id}\n")
    
    with open(down_df_filename, 'w') as f:
        for kegg_id in down_df:
            f.write(f"{kegg_id}\n")
    
    return up_df, down_df, up_df_filename, down_df_filename, enrich_output_prefix, kegg_enrich_output_prefix


def run_kegg_enrich_analysis(up_df, down_df, up_df_filename, down_df_filename, enrich_output_prefix, kegg_enrich_output_prefix, definition_df=None):
    """运行KEGG富集分析并处理结果"""
    # 运行富集分析脚本
    enrich_script_path = "/home/colddata/qinqiang/script/MetaboliteAnalysis/MetaboliteEnrich/metabolite_enrich.r"
    up_df_cmd = f"/opt/biosoft/R-4.2.2/bin/Rscript {enrich_script_path} --datatable {up_df_filename} --outputprefix {enrich_output_prefix}_Up"
    down_df_cmd = f"/opt/biosoft/R-4.2.2/bin/Rscript {enrich_script_path} --datatable {down_df_filename} --outputprefix {enrich_output_prefix}_Down"
    
    subprocess.run(up_df_cmd, shell=True)
    subprocess.run(down_df_cmd, shell=True)
    
    # 运行KEGG富集分析脚本
    kegg_enrich_script_path = "/home/colddata/qinqiang/script/MetaboliteAnalysis/MetaboliteEnrich/metabolite_kegg_enrich.r"
    up_df_cmd = f"/opt/biosoft/R-4.2.2/bin/Rscript {kegg_enrich_script_path} --datatable {up_df_filename} --outputprefix {kegg_enrich_output_prefix}_Up"
    down_df_cmd = f"/opt/biosoft/R-4.2.2/bin/Rscript {kegg_enrich_script_path} --datatable {down_df_filename} --outputprefix {kegg_enrich_output_prefix}_Down"
    
    subprocess.run(up_df_cmd, shell=True)
    subprocess.run(down_df_cmd, shell=True)
    
    # 合并定义信息（如果有提供）
    if definition_df is not None and isinstance(definition_df, pd.DataFrame) and len(definition_df) > 0:
        if "KEGG" not in definition_df.columns:
            raise ValueError("definition_df 没有 KEGG 列，无法 merge")
        
        # 合并定义信息
        up_df_merged = pd.merge(pd.DataFrame({'KEGG': up_df}), definition_df, on='KEGG', how='left')
        down_df_merged = pd.merge(pd.DataFrame({'KEGG': down_df}), definition_df, on='KEGG', how='left')
        
        # 保存合并后的数据
        up_df_merged.to_csv(up_df_filename, sep='\t', index=False)
        down_df_merged.to_csv(down_df_filename, sep='\t', index=False)


def summarize_vip_and_enrich(input_dir, definition_df=None):
    # 列出 VIP 文件
    all_files = glob.glob(os.path.join(input_dir, "**", "*"), recursive=True)
    vip_files = [f for f in all_files if "VIP" in f and f.endswith('.xlsx')]
    
    class_count_list = []
    subclass_count_list = []
    
    for vip_file in vip_files:
        print(vip_file)
        compare_group_name = os.path.basename(vip_file).replace("_VIP.xlsx", "")
        
        # 读取数据
        vip_df = pd.read_excel(vip_file)
        
        # 标注上调/下调/不显著
        conditions = [
            (vip_df['VIP'] > 1) & ((vip_df['FoldChange'] > 1.2) | (vip_df['FoldChange'] < 0.8)) & (vip_df['FoldChange'] >= 1.2),
            (vip_df['VIP'] > 1) & ((vip_df['FoldChange'] > 1.2) | (vip_df['FoldChange'] < 0.8)) & (vip_df['FoldChange'] < 0.8)
        ]
        choices = ['Up', 'Down']
        vip_df['regulation'] = np.select(conditions, choices, default='NoSignificant')
        
        # 保存标注后的数据
        vip_df.to_excel(vip_file, index=False)
        
        # ======= enrich start ======== 有 KEGG 列才能做 enrich
        if "KEGG" in vip_df.columns:
            result = prepare_kegg_data(vip_df, compare_group_name)
            if result[0] is not None:  # 检查是否有KEGG数据
                up_df, down_df, up_df_filename, down_df_filename, enrich_output_prefix, kegg_enrich_output_prefix = result
                run_kegg_enrich_analysis(up_df, down_df, up_df_filename, down_df_filename, 
                                    enrich_output_prefix, kegg_enrich_output_prefix, definition_df)
        # ======= enrich end =======
        
        # 计算总数 (包含上调、下调和无显著差异)
        if len(vip_df) >= 1:
            # 处理Class统计
            if "Class" in vip_df.columns:
                # 计算总数
                class_total = vip_df.groupby('Class')['Metabolite'].count().reset_index()
                class_total = class_total[class_total['Class'] != ""]
                class_total = class_total.rename(columns={'Metabolite': 'Total'})
                
                # 计算上调数 (VIP > 1 且 FoldChange > 1.2)
                class_up = vip_df[(vip_df['VIP'] > 1) & (vip_df['FoldChange'] > 1.2)].groupby('Class')['Metabolite'].count().reset_index()
                class_up = class_up[class_up['Class'] != ""]
                class_up = class_up.rename(columns={'Metabolite': 'Up'})
                
                # 计算下调数 (VIP > 1 且 FoldChange < 0.8)
                class_down = vip_df[(vip_df['VIP'] > 1) & (vip_df['FoldChange'] < 0.8)].groupby('Class')['Metabolite'].count().reset_index()
                class_down = class_down[class_down['Class'] != ""]
                class_down = class_down.rename(columns={'Metabolite': 'Down'})
                
                # 计算比较组总数 (上调+下调)
                class_comparison_total = vip_df[(vip_df['VIP'] > 1) & ((vip_df['FoldChange'] > 1.2) | (vip_df['FoldChange'] < 0.8))].groupby('Class')['Metabolite'].count().reset_index()
                class_comparison_total = class_comparison_total[class_comparison_total['Class'] != ""]
                class_comparison_total = class_comparison_total.rename(columns={'Metabolite': 'Comparison_total'})
                
                # 合并结果
                class_result = pd.merge(class_total, class_comparison_total, on='Class', how='outer')
                class_result = pd.merge(class_result, class_up, on='Class', how='outer')
                class_result = pd.merge(class_result, class_down, on='Class', how='outer')
                class_result.columns = ['Class', 'Total', 'Comparison_total', 'Up', 'Down']
                
                # 替换NA为0
                class_result = class_result.fillna(0)
                
                # 计算上调比例和下调比例
                class_result['Up_ratio'] = class_result['Up'] / class_result['Total']
                class_result['Down_ratio'] = class_result['Down'] / class_result['Total']
                
                # 替换无穷大和NaN为0
                class_result['Up_ratio'] = class_result['Up_ratio'].replace([np.inf, -np.inf], 0).fillna(0)
                class_result['Down_ratio'] = class_result['Down_ratio'].replace([np.inf, -np.inf], 0).fillna(0)
                
                # 转换为百分比格式
                class_result['Up_ratio'] = (class_result['Up_ratio'] * 100).round(2).astype(str) + '%'
                class_result['Down_ratio'] = (class_result['Down_ratio'] * 100).round(2).astype(str) + '%'
                
                # 重新排列列的顺序
                class_result = class_result[['Class', 'Total', 'Comparison_total', 'Up', 'Down', 'Up_ratio', 'Down_ratio']]
                
                # 添加比较组名称作为列名前缀
                comp_name = os.path.splitext(os.path.basename(vip_file))[0]
                class_result.columns = ['Class', 'Total'] + [f"{comp_name}_{col}" for col in class_result.columns[2:]]
                
                class_count_list.append(class_result)
            
            # 处理SubClass统计
            if "SubClass" in vip_df.columns:
                # 计算总数
                subclass_total = vip_df.groupby('SubClass')['Metabolite'].count().reset_index()
                subclass_total = subclass_total[subclass_total['SubClass'] != ""]
                subclass_total = subclass_total.rename(columns={'Metabolite': 'Total'})
                
                # 计算上调数 (VIP > 1 且 FoldChange > 1.2)
                subclass_up = vip_df[(vip_df['VIP'] > 1) & (vip_df['FoldChange'] > 1.2)].groupby('SubClass')['Metabolite'].count().reset_index()
                subclass_up = subclass_up[subclass_up['SubClass'] != ""]
                subclass_up = subclass_up.rename(columns={'Metabolite': 'Up'})
                
                # 计算下调数 (VIP > 1 且 FoldChange < 0.8)
                subclass_down = vip_df[(vip_df['VIP'] > 1) & (vip_df['FoldChange'] < 0.8)].groupby('SubClass')['Metabolite'].count().reset_index()
                subclass_down = subclass_down[subclass_down['SubClass'] != ""]
                subclass_down = subclass_down.rename(columns={'Metabolite': 'Down'})
                
                # 计算比较组总数 (上调+下调)
                subclass_comparison_total = vip_df[(vip_df['VIP'] > 1) & ((vip_df['FoldChange'] > 1.2) | (vip_df['FoldChange'] < 0.8))].groupby('SubClass')['Metabolite'].count().reset_index()
                subclass_comparison_total = subclass_comparison_total[subclass_comparison_total['SubClass'] != ""]
                subclass_comparison_total = subclass_comparison_total.rename(columns={'Metabolite': 'Comparison_total'})
                
                # 合并结果
                subclass_result = pd.merge(subclass_total, subclass_comparison_total, on='SubClass', how='outer')
                subclass_result = pd.merge(subclass_result, subclass_up, on='SubClass', how='outer')
                subclass_result = pd.merge(subclass_result, subclass_down, on='SubClass', how='outer')
                subclass_result.columns = ['SubClass', 'Total', 'Comparison_total', 'Up', 'Down']
                
                # 替换NA为0
                subclass_result = subclass_result.fillna(0)
                
                # 计算上调比例和下调比例
                subclass_result['Up_ratio'] = subclass_result['Up'] / subclass_result['Total']
                subclass_result['Down_ratio'] = subclass_result['Down'] / subclass_result['Total']
                
                # 替换无穷大和NaN为0
                subclass_result['Up_ratio'] = subclass_result['Up_ratio'].replace([np.inf, -np.inf], 0).fillna(0)
                subclass_result['Down_ratio'] = subclass_result['Down_ratio'].replace([np.inf, -np.inf], 0).fillna(0)
                
                # 转换为百分比格式
                subclass_result['Up_ratio'] = (subclass_result['Up_ratio'] * 100).round(2).astype(str) + '%'
                subclass_result['Down_ratio'] = (subclass_result['Down_ratio'] * 100).round(2).astype(str) + '%'
                
                # 重新排列列的顺序
                subclass_result = subclass_result[['SubClass', 'Total', 'Comparison_total', 'Up', 'Down', 'Up_ratio', 'Down_ratio']]
                
                # 添加比较组名称作为列名前缀
                comp_name = os.path.splitext(os.path.basename(vip_file))[0]
                subclass_result.columns = ['SubClass', 'Total'] + [f"{comp_name}_{col}" for col in subclass_result.columns[2:]]
                
                subclass_count_list.append(subclass_result)
        else:
            print(f"警告: vip_df_all {os.path.basename(vip_file)} 没有数据")
    
    # 合并Class统计结果并导出
    if class_count_list:
        # 使用简单的合并函数，让pandas自动处理重复列名
        def simple_merge(left, right):
            return pd.merge(left, right, on='Class', how='outer', suffixes=('', '_dup'))
        
        # 使用简单合并函数
        class_count_result = reduce(simple_merge, class_count_list)
        class_count_result = class_count_result.fillna(0)
        
        # 处理重复的Total列，保留第一个Total列
        total_cols = [col for col in class_count_result.columns if col.startswith('Total')]
        if len(total_cols) > 1:
            cols_to_keep = ['Class']
            total_kept = False
            for col in class_count_result.columns:
                if col.startswith('Total') and not total_kept:
                    cols_to_keep.append(col)
                    total_kept = True
                elif not col.startswith('Total'):
                    cols_to_keep.append(col)
            class_count_result = class_count_result[cols_to_keep]
        
        # 重新排列列的顺序: Class, Total, {comparison_1}_total, {comparison_2}_total, .._1_Up, ..._2_Up, ..1_Down, .._2_Down, ...Up_ratio, ..._Down_ratio
        reordered_cols = ['Class', 'Total']
        
        # 添加所有comparison_total列
        for col in class_count_result.columns:
            if col.endswith('_Comparison_total'):
                reordered_cols.append(col)
        
        # 添加所有Up列
        for col in class_count_result.columns:
            if col.endswith('_Up') and not col.endswith('_Up_ratio'):
                reordered_cols.append(col)
        
        # 添加所有Down列
        for col in class_count_result.columns:
            if col.endswith('_Down') and not col.endswith('_Down_ratio'):
                reordered_cols.append(col)
        
        # 添加所有Up_ratio列
        for col in class_count_result.columns:
            if col.endswith('_Up_ratio'):
                reordered_cols.append(col)
        
        # 添加所有Down_ratio列
        for col in class_count_result.columns:
            if col.endswith('_Down_ratio'):
                reordered_cols.append(col)
        
        # 重新排列DataFrame
        class_count_result = class_count_result[reordered_cols]
        
        print("Class统计结果:")
        print(class_count_result.head())
        
        class_count_result.to_excel(
            os.path.join(input_dir, "Significant_compound_count_by_class.xlsx"),
            index=False
        )
    else:
        print("警告: 没有Class统计数据可导出")
    
    # 合并SubClass统计结果并导出
    if subclass_count_list:
        # 使用简单的合并函数，让pandas自动处理重复列名
        def simple_merge_subclass(left, right):
            return pd.merge(left, right, on='SubClass', how='outer', suffixes=('', '_dup'))
        
        # 使用简单合并函数
        subclass_count_result = reduce(simple_merge_subclass, subclass_count_list)
        subclass_count_result = subclass_count_result.fillna(0)
        
        # 处理重复的Total列，保留第一个Total列
        total_cols = [col for col in subclass_count_result.columns if col.startswith('Total')]
        if len(total_cols) > 1:
            # 保留第一个Total列，删除其他的Total列
            cols_to_keep = ['SubClass']
            total_kept = False
            for col in subclass_count_result.columns:
                if col.startswith('Total') and not total_kept:
                    cols_to_keep.append(col)
                    total_kept = True
                elif not col.startswith('Total'):
                    cols_to_keep.append(col)
            subclass_count_result = subclass_count_result[cols_to_keep]
        
        # 重新排列列的顺序: SubClass, Total, {comparison_1}_total, {comparison_2}_total, .._1_Up, ..._2_Up, ..1_Down, .._2_Down, ...Up_ratio, ..._Down_ratio
        reordered_cols = ['SubClass', 'Total']
        
        # 添加所有comparison_total列
        for col in subclass_count_result.columns:
            if col.endswith('_Comparison_total'):
                reordered_cols.append(col)
        
        # 添加所有Up列
        for col in subclass_count_result.columns:
            if col.endswith('_Up') and not col.endswith('_Up_ratio'):
                reordered_cols.append(col)
        
        # 添加所有Down列
        for col in subclass_count_result.columns:
            if col.endswith('_Down') and not col.endswith('_Down_ratio'):
                reordered_cols.append(col)
        
        # 添加所有Up_ratio列
        for col in subclass_count_result.columns:
            if col.endswith('_Up_ratio'):
                reordered_cols.append(col)
        
        # 添加所有Down_ratio列
        for col in subclass_count_result.columns:
            if col.endswith('_Down_ratio'):
                reordered_cols.append(col)
        
        # 重新排列DataFrame
        subclass_count_result = subclass_count_result[reordered_cols]
        
        print("SubClass统计结果:")
        print(subclass_count_result.head())
        
        subclass_count_result.to_excel(
            os.path.join(input_dir, "Significant_compound_count_by_subclass.xlsx"),
            index=False
        )
    else:
        print("警告: 没有SubClass统计数据可导出")
    
    return True

# 使用示例
# summarize_vip_and_enrich("/path/to/output/directory")

if __name__ == '__main__':
    args = parse_input()
    input_dir = args.i
    compound_def = load_table(args.d)
    summarize_vip_and_enrich(input_dir, compound_def)