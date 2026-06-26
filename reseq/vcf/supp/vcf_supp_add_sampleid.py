#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2026/05/15 17:23
# Author        : William GoGo
import os, sys
import argparse
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df


def parse_input():
    # 使用argparse解析参数
    parser = argparse.ArgumentParser(description="VCF表格添加Sample_ID列，识别所有有变异的样本")
    parser.add_argument("input_file", help="输入VCF表格文件路径")
    parser.add_argument("output_file", help="输出结果文件路径")
    
    args = parser.parse_args()
    
    return args
    

def main():
    args = parse_input()
    logger.info(f"开始处理文件: {args.input_file}")
    
    # 读取表格
    df = load_table(args.input_file)
    logger.info(f"成功读取表格，共 {len(df)} 行, {len(df.columns)} 列")
    
    # 定位列位置
    try:
        format_idx = df.columns.get_loc('FORMAT')
        supp_idx = df.columns.get_loc('SUPP')
    except KeyError as e:
        logger.error(f"表格中缺少必要列: {e}")
        sys.exit(1)
    
    # 获取样本列名
    sample_columns = df.columns[format_idx + 1 : supp_idx].tolist()
    logger.info(f"识别到 {len(sample_columns)} 个样本列")
    
    # 获取GT在FORMAT中的位置
    format_fields = df['FORMAT'].iloc[0].split(':')
    if 'GT' not in format_fields:
        logger.error("FORMAT列中没有GT字段，无法判断样本基因型")
        sys.exit(1)
    gt_idx = format_fields.index('GT')
    
    # 定义判断样本是否有变异的函数
    def get_variant_samples(row):
        variant_samples = []
        for sample in sample_columns:
            gt = row[sample].split(':')[gt_idx]
            # 排除参考纯合和缺失基因型
            if gt not in ('0/0', './.', '0|0', '.|.', '0', '.'):
                variant_samples.append(sample)
        return ';'.join(variant_samples)
    
    # 应用函数获取Sample_ID列
    logger.info("正在处理每一行，识别变异样本...")
    df.insert(supp_idx + 1, 'Sample_ID', df.apply(get_variant_samples, axis=1))
    
    # 统计信息
    variant_rows = len(df[df['Sample_ID'] != ''])
    logger.info(f"处理完成，共有 {variant_rows} 行存在变异样本")
    
    # 保存结果
    write_output_df(df, args.output_file, index=False)
    logger.info(f"结果已保存到: {args.output_file}")

if __name__ == "__main__":
    main()