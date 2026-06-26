#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/04/29 21:35
# Author        : William GoGo
import os
import sys
import re
import gzip
import argparse
import pandas as pd
from loguru import logger
from openpyxl import Workbook
from openpyxl.styles import numbers

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import write_output_df

def parse_args():
    parser = argparse.ArgumentParser(description='INFO字段提取工具 - 支持标准VCF和Annovar结果，提取SUPP/AC/AF/DP/MQ字段')
    parser.add_argument('-i', '--input', required=True, help='输入文件路径，支持标准VCF(.vcf/.vcf.gz)和Annovar输出表格(.txt/.tsv)')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径，支持.tsv/.txt和.xlsx/.xls格式')
    parser.add_argument('--info-col', default='INFO', help='INFO字段所在的列名，Annovar结果通常为Otherinfo11，默认值: INFO')
    parser.add_argument('--only-required-cols', action='store_true', help='仅保留基础列+解析字段，默认保留所有原始列，解析字段新增到最后')
    return parser.parse_args()

def write_excel(df, output_path):
    """输出DataFrame到Excel文件，自动处理列类型防止日期转换"""
    wb = Workbook()
    ws = wb.active
    ws.title = "VCF_INFO_Fields"
    
    # 写入表头
    headers = list(df.columns)
    ws.append(headers)
    
    # 定义文本格式，防止Excel自动转日期
    text_format = numbers.FORMAT_TEXT
    
    # 判断每列类型：数字列保持数值，可能被识别为日期的列设为文本
    date_like_pattern = re.compile(r'.*[/\-|].*')  # 包含/、-、|的内容可能被识别为日期
    col_types = []
    
    for col in headers:
        # 先尝试是否可以转成数字
        try:
            pd.to_numeric(df[col].dropna().replace('', pd.NA), errors='raise')
            col_types.append('numeric')
        except:
            # 判断列中是否包含可能被识别为日期的内容
            has_date_like = df[col].astype(str).str.match(date_like_pattern).any()
            col_types.append('text' if has_date_like else 'auto')
    
    # 写入数据
    for r_idx, row in enumerate(df.itertuples(index=False, name=None), start=2):
        for c_idx, (value, col_type) in enumerate(zip(row, col_types), start=1):
            if pd.isna(value) or value == '':
                cell_value = None
            else:
                cell_value = value
            cell = ws.cell(row=r_idx, column=c_idx, value=cell_value)
            if col_type == 'text':
                cell.number_format = text_format
    
    # 自动调整列宽
    for col in ws.columns:
        max_length = 0
        col_letter = col[0].column_letter
        for cell in col:
            try:
                cell_len = len(str(cell.value))
                if cell_len > max_length:
                    max_length = cell_len
            except:
                pass
        adjusted_width = min(max_length + 2, 50)  # 最大宽度50防止过宽
        ws.column_dimensions[col_letter].width = adjusted_width
    
    wb.save(output_path)
    logger.info(f'Excel格式输出完成，文件路径: {output_path}')

def read_vcf(vcf_path):
    """读取文件，支持普通VCF、gzip压缩VCF和Annovar表格格式"""
    logger.info(f'开始读取文件: {vcf_path}')
    
    # 判断是否是gzip压缩文件
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    mode = 'rt' if vcf_path.endswith('.gz') else 'r'
    
    header = None
    data_lines = []
    
    with open_func(vcf_path, mode, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                header = line.lstrip('#').split('\t')
                logger.info(f'识别为标准VCF格式，列名: {header[:8]}...')
            else:
                if line:
                    data_lines.append(line.split('\t'))
    
    if not header:
        # 没有找到VCF表头，尝试作为普通表格读取（Annovar结果）
        logger.info('未识别到标准VCF表头，尝试作为普通表格读取（Annovar结果格式）')
        with open_func(vcf_path, mode, encoding='utf-8') as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith('##')]
        if len(lines) < 2:
            logger.error('文件内容不足，无法提取表头和数据')
            sys.exit(1)
        header = lines[0].split('\t')
        data_lines = [line.split('\t') for line in lines[1:]]
        logger.info(f'读取到表头列: {header[:10]}...')
    
    logger.info(f'共读取到 {len(data_lines)} 条记录')
    return pd.DataFrame(data_lines, columns=header)

def extract_info_fields(df, info_col='INFO', only_required_cols=False):
    """从指定的INFO列中提取SUPP/AC/AF/DP/MQ字段，支持Annovar结果"""
    logger.info(f'从列[{info_col}]中提取SUPP/AC/AF/DP/MQ字段')
    
    if info_col not in df.columns:
        logger.error(f'指定的INFO列[{info_col}]不存在于文件中，现有列: {list(df.columns)}')
        sys.exit(1)
    
    # 目标字段列表
    target_fields = ['SUPP', 'AC', 'AF', 'DP', 'MQ']
    
    # 初始化字段
    for field in target_fields:
        df[field] = pd.NA
    
    # 解析每一行INFO
    for idx, info_str in enumerate(df[info_col]):
        if pd.isna(info_str) or info_str == '':
            continue
        info_dict = {}
        items = str(info_str).split(';')
        for item in items:
            if '=' in item:
                key, value = item.split('=', 1)
                info_dict[key] = value
        
        # 填充目标字段
        for field in target_fields:
            if field in info_dict:
                df.at[idx, field] = info_dict[field]
    
    # 统计字段存在情况
    for field in target_fields:
        non_na_count = df[field].notna().sum()
        logger.info(f'{field} 字段存在记录数: {non_na_count}/{len(df)}')
    
    # 适配Annovar列名映射
    annovar_col_map = {
        'Chr': 'CHROM',
        'Start': 'POS',
        'Ref': 'REF',
        'Alt': 'ALT',
        'Qual': 'QUAL'
    }
    # 重命名列
    for old_col, new_col in annovar_col_map.items():
        if old_col in df.columns and new_col not in df.columns:
            df = df.rename(columns={old_col: new_col})
    
    # 选择输出列
    if only_required_cols:
        # 仅保留基础列+解析字段
        output_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL'] + target_fields
        # 只保留存在的列
        output_columns = [col for col in output_columns if col in df.columns]
        df = df[output_columns]
        logger.info(f'已启用精简模式，仅保留 {len(output_columns)} 列')
    else:
        # 保留所有原始列，解析字段新增到最后面
        original_cols = list(df.columns)
        # 把解析字段放到最后（确保顺序正确）
        output_columns = [col for col in original_cols if col not in target_fields] + target_fields
        df = df[output_columns]
        logger.info(f'保留所有原始列共 {len(original_cols)} 个，新增 {len(target_fields)} 个解析字段到最后')
    
    return df

def main():
    args = parse_args()
    
    # 读取文件（支持VCF和Annovar表格）
    df = read_vcf(args.input)
    
    # 提取目标字段
    df = extract_info_fields(df, info_col=args.info_col, only_required_cols=args.only_required_cols)
    
    # 输出结果
    logger.info(f'开始输出结果到: {args.output}')
    if args.output.lower().endswith(('.xlsx', '.xls')):
        write_excel(df, args.output)
    else:
        write_output_df(df, args.output, index=False)
        logger.info(f'TSV格式输出完成，文件路径: {args.output}')
    
    logger.info('处理完成! 已成功提取SUPP/AC/AF/DP/MQ字段')

if __name__ == '__main__':
    main()