#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/03/17 11:42
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df

def parse_args():
    parser = argparse.ArgumentParser(description='Manta VCF INFO字段格式化工具')
    parser.add_argument('-i', '--input', required=True, help='输入VCF文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出表格文件路径')
    return parser.parse_args()

def parse_info(info_str):
    """解析INFO字段，提取需要的信息"""
    info_dict = {}
    items = info_str.split(';')
    for item in items:
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            # 布尔字段如IMPRECISE
            info_dict[item] = 'True'
    
    # 提取需要的字段
    result = {
        'SVTYPE': info_dict.get('SVTYPE', ''),
        'SVLEN': info_dict.get('SVLEN', ''),
        'END': info_dict.get('END', ''),
        'STRANDS': info_dict.get('STRANDS', ''),
        'CIPOS': info_dict.get('CIPOS', ''),
        'CIEND': info_dict.get('CIEND', ''),
        'IMPRECISE': info_dict.get('IMPRECISE', 'False')
    }
    
    return result

def read_vcf(vcf_path):
    """读取VCF文件，跳过注释行，返回DataFrame"""
    logger.info(f'开始读取VCF文件: {vcf_path}')
    
    # 先找到表头行
    header = None
    data_lines = []
    with open(vcf_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                header = line.lstrip('#').split('\t')
            else:
                if line:
                    data_lines.append(line.split('\t'))
    
    if not header:
        logger.error('未找到VCF表头行')
        sys.exit(1)
    
    logger.info(f'VCF文件列名: {header}')
    logger.info(f'读取到 {len(data_lines)} 条变异记录')
    
    # 创建DataFrame
    df = pd.DataFrame(data_lines, columns=header)
    
    # 解析INFO字段
    logger.info('开始解析INFO字段')
    info_dicts = df['INFO'].apply(parse_info).tolist()
    info_df = pd.DataFrame(info_dicts)
    
    # 合并原数据和INFO字段
    result_df = pd.concat([df, info_df], axis=1)
    
    # 移除全空的列
    empty_cols = [col for col in info_df.columns if (info_df[col] == '').all()]
    if empty_cols:
        logger.info(f'移除全空的INFO字段: {empty_cols}')
        result_df = result_df.drop(columns=empty_cols)
    
    # 统计各SV类型
    sv_type_counts = result_df['SVTYPE'].value_counts().to_dict()
    logger.info(f'SV类型统计: {sv_type_counts}')
    
    return result_df

def main():
    args = parse_args()
    
    # 读取并处理VCF
    df = read_vcf(args.input)
    
    # 输出结果
    logger.info(f'开始输出结果到: {args.output}')
    write_output_df(df, args.output, index=False)
    logger.info('处理完成!')

if __name__ == '__main__':
    main()