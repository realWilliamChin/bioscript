#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/04/21 18:11
# Author        : William GoGo
import os, sys
import re
import argparse
import pandas as pd
from loguru import logger
from openpyxl import Workbook
from openpyxl.styles import numbers

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df

def parse_args():
    parser = argparse.ArgumentParser(description='Manta VCF + SnpEff 注释格式化工具')
    parser.add_argument('-i', '--input', required=True, help='输入VCF文件路径 (Manta输出并经过SnpEff注释的VCF)')
    parser.add_argument('-o', '--output', required=True, help='输出表格文件路径，支持.tsv/.txt和.xlsx/.xls格式，Excel格式自动防止日期转换')
    parser.add_argument('--disable-parse-format', action='store_true', help='是否禁用FORMAT字段解析，默认开启解析并展开为样本单独列')
    parser.add_argument('--disable-sr-filter', action='store_true', help='是否禁用SR列过滤（默认开启：SR逗号前为0，逗号后>3）')
    parser.add_argument('--disable-gt-filter', action='store_true', help='是否禁用GT=1/1过滤（默认开启）')
    parser.add_argument('--disable-imprecise-filter', action='store_true', help='是否禁用IMPRECISE列过滤（默认开启：只保留IMPRECISE==False的记录）')
    parser.add_argument('--disable-svtype-filter', action='store_true', help='是否禁用SVTYPE过滤（默认开启：只保留DEL/INS/DUP）')
    parser.add_argument('--disable-significance-filter', action='store_true', help='是否禁用Significance过滤（默认开启：只保留high/moderate）')
    parser.add_argument('--disable-svlen-filter', action='store_true', help='是否禁用SVLEN非空过滤（默认开启）')
    parser.add_argument('--disable-parse-info', action='store_true', help='是否禁用INFO字段解析，默认开启解析提取SV和注释信息')
    parser.add_argument('--keep-original-info', action='store_true', help='是否在输出结果中保留原始INFO列，默认不保留')
    parser.add_argument('--keep-original-format', action='store_true', help='是否在输出结果中保留原始FORMAT和样本FORMAT列，默认不保留')
    
    return parser.parse_args()

def parse_info_manta(info_str):
    """解析INFO字段，提取Manta SV相关信息"""
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

def read_vcf_manta(vcf_path, parse_info=True):
    """读取VCF文件，跳过注释行，返回DataFrame并解析Manta INFO字段"""
    logger.info(f'开始读取VCF文件: {vcf_path}')
    
    # 先找到表头行
    header = None
    header_line = None  # 保存原始表头行（#CHROM开头）
    vcf_header = []     # 保存所有##开头的注释行
    data_lines = []
    original_vcf_lines = []  # 保存所有原始数据行（不包含注释和表头）
    with open(vcf_path, 'r', encoding='utf-8') as f:
        for line in f:
            stripped_line = line.strip()
            if stripped_line.startswith('##'):
                vcf_header.append(line.rstrip('\n'))
                continue
            elif stripped_line.startswith('#CHROM'):
                header_line = line.rstrip('\n')
                header = stripped_line.lstrip('#').split('\t')
            else:
                if stripped_line:
                    data_lines.append(stripped_line.split('\t'))
                    original_vcf_lines.append(line.rstrip('\n'))
    
    if not header:
        logger.error('未找到VCF表头行')
        sys.exit(1)
    
    logger.info(f'VCF文件列名: {header}')
    logger.info(f'读取到 {len(data_lines)} 条变异记录')
    
    # 创建DataFrame
    df = pd.DataFrame(data_lines, columns=header)
    
    # 记录原始样本列（FORMAT后面的列，避免后续新增列被误判为样本）
    sample_cols = []
    if 'FORMAT' in df.columns:
        sample_cols = list(df.columns[df.columns.get_loc('FORMAT')+1:])
    
    if parse_info:
        # 解析INFO字段
        logger.info('开始解析Manta INFO字段')
        info_dicts = df['INFO'].apply(parse_info_manta).tolist()
        info_df = pd.DataFrame(info_dicts)
        
        # 合并原数据和INFO字段
        result_df = pd.concat([df, info_df], axis=1)
        
        # 移除全空的列
        empty_cols = [col for col in info_df.columns if (info_df[col] == '').all()]
        if empty_cols:
            logger.info(f'移除全空的INFO字段: {empty_cols}')
            result_df = result_df.drop(columns=empty_cols)
        
        # 统计各SV类型
        if 'SVTYPE' in result_df.columns:
            sv_type_counts = result_df['SVTYPE'].value_counts().to_dict()
            logger.info(f'SV类型统计: {sv_type_counts}')
    else:
        result_df = df
        logger.info('已禁用INFO字段解析，跳过SV信息提取')
    
    return result_df, sample_cols, vcf_header, header_line, original_vcf_lines

def parse_snpeff(vcf_df):
    """解析SnpEff注释信息"""
    logger.info('开始解析SnpEff注释信息')
    
    if 'INFO' not in vcf_df.columns:
        logger.warning('VCF文件缺少INFO列，无法解析SnpEff注释')
        vcf_df['Mutation_type'] = None
        vcf_df['Significance'] = None
        return vcf_df

    # 提取EFF字段信息
    # 取出来EFF=后面的第一个(前面的词作为突变类型
    vcf_df['Mutation_type'] = vcf_df['INFO'].str.extract(r';EFF=([^(]+)\(')[0]
    # 取出来第一个(后面的第一个|前面的内容作为显著性
    vcf_df['Significance'] = vcf_df['INFO'].str.extract(r';EFF=[^()]*\(([^|]*)\|')[0]
    
    # 统计突变类型
    mut_type_counts = vcf_df['Mutation_type'].value_counts().to_dict()
    logger.info(f'突变类型统计: {mut_type_counts}')
    
    return vcf_df

def parse_format(vcf_df, sample_cols):
    """解析VCF的FORMAT字段，并将格式及每个样本列转为新列"""
    logger.info('开始解析FORMAT字段')
    
    if 'FORMAT' not in vcf_df.columns or len(sample_cols) == 0:
        logger.warning('VCF文件缺少FORMAT列或样本列，无法解析')
        return vcf_df
    if len(sample_cols)==0:
        logger.warning('VCF文件无样本列，无需解析FORMAT成新列')
        return vcf_df

    # 针对每一行，解析FORMAT及样本内容，并展开
    format_fields = vcf_df['FORMAT'].str.split(':')
    for sample in sample_cols:
        sample_values = vcf_df[sample].str.split(':')
        # 预存字典用于组装 (行号, 字段) -> 值
        out_dict = {}
        for idx, (fields, values) in enumerate(zip(format_fields, sample_values)):
            if isinstance(fields, list) and isinstance(values, list):
                for f, v in zip(fields, values):
                    out_dict.setdefault(f, {})[idx] = v
        # 新增列，列名: sample名 + "_" + format字段名
        all_fields = set().union(*format_fields.dropna())
        logger.info(f"FORMAT字段共有 {len(all_fields)} 个键: {', '.join(all_fields)}")
        for f in all_fields:
            v = pd.Series(out_dict.get(f, {}))
            v = v.reindex(vcf_df.index)  # 对齐索引，缺失值填充为NaN
            vcf_df[f"{sample}_{f}"] = v

    logger.info(f'FORMAT解析完成，共 {len(sample_cols)} 个样本，新增 {len(set().union(*format_fields.dropna())) * len(sample_cols)} 列')
    return vcf_df

def write_excel(df, output_path):
    """输出DataFrame到Excel文件，自动处理列类型防止日期转换"""
    wb = Workbook()
    ws = wb.active
    ws.title = "SV_Annotation"
    
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

def main():
    args = parse_args()
    
    # 第一步：执行Manta VCF格式化
    info_parsed = not args.disable_parse_info
    df, sample_cols, vcf_header, header_line, original_vcf_lines = read_vcf_manta(args.input, parse_info=info_parsed)
    
    # 第二步：执行SnpEff注释解析
    if info_parsed:
        df = parse_snpeff(df)
    else:
        logger.info('已禁用INFO解析，跳过SnpEff注释提取，将自动禁用INFO相关过滤规则')
        # 自动禁用所有依赖INFO字段的过滤
        args.disable_imprecise_filter = True
        args.disable_svtype_filter = True
        args.disable_significance_filter = True
        args.disable_svlen_filter = True
    
    # 第三步：解析FORMAT字段
    format_parsed = False
    if not args.disable_parse_format:
        df = parse_format(df, sample_cols)
        format_parsed = True
    else:
        logger.info('已禁用FORMAT字段解析，不展开为样本单独列，将自动跳过SR和GT过滤')
    
    # 重命名POS列为START
    df = df.rename(columns={'POS': 'START'})
    
    # 生成统一的列顺序（过滤前和过滤后都用这个顺序）
    # 1. 基础核心列
    all_columns = ['ID', 'CHROM', 'START', 'END', 'REF', 'ALT', 'QUAL']
    
    # 2. 可选保留原始INFO列
    if args.keep_original_info:
        all_columns.append('INFO')
    
    # 3. SV解析列
    sv_columns = ['IMPRECISE', 'SVTYPE', 'SVLEN', 'STRANDS', 'CIPOS', 'CIEND']
    all_columns.extend(sv_columns)
    
    # 4. 注释列
    anno_columns = ['Mutation_type', 'Significance']
    all_columns.extend(anno_columns)
    
    # 5. 可选保留原始FORMAT和样本列
    if args.keep_original_format:
        all_columns.append('FORMAT')
        all_columns.extend(sample_cols)
    
    # 6. 解析后的样本列
    sample_suffixes = ['_PR', '_PL', '_SR', '_GQ', '_GT', '_FT']
    for sample in sample_cols:
        for suffix in sample_suffixes:
            col = f'{sample}{suffix}'
            all_columns.append(col)
    
    # 补全所有列，不存在的列填NaN
    for col in all_columns:
        if col not in df.columns:
            df[col] = pd.NA
    
    # 调整列顺序为统一顺序
    df = df.reindex(columns=all_columns)
    
    # 输出过滤前的文件：输入VCF文件名替换.vcf后缀为.xlsx
    import os
    input_filename = os.path.basename(args.input)
    pre_filter_filename = re.sub(r'\.vcf$', '.xlsx', input_filename, flags=re.IGNORECASE)
    pre_filter_path = os.path.join(os.path.dirname(os.path.abspath(args.input)), pre_filter_filename)
    write_excel(df, pre_filter_path)
    logger.info(f'过滤前原始数据已输出到: {pre_filter_path}，列顺序已统一')
    
    # 第四步：过滤数据
    logger.info('开始应用过滤条件')
    initial_count = len(df)
    
    # 1. IMPRECISE过滤
    if not args.disable_imprecise_filter:
        df = df[df['IMPRECISE'] == 'False']
        logger.info(f'过滤IMPRECISE==False后剩余: {len(df)} 条记录')
    else:
        logger.info(f'已禁用IMPRECISE过滤，跳过该条件检查')
    
    # 2. SVTYPE过滤
    if not args.disable_svtype_filter:
        df = df[df['SVTYPE'].isin(['DEL', 'INS', 'DUP'])]
        logger.info(f'过滤SVTYPE为DEL/INS/DUP后剩余: {len(df)} 条记录')
    else:
        logger.info(f'已禁用SVTYPE过滤，跳过该条件检查')
    
    # 3. Significance过滤
    if not args.disable_significance_filter:
        df = df[df['Significance'].str.contains('high|moderate', case=False, na=False)]
        logger.info(f'过滤Significance为High/Moderate后剩余: {len(df)} 条记录')
    else:
        logger.info(f'已禁用Significance过滤，跳过该条件检查')
    
    # 4. SVLEN非空过滤
    if not args.disable_svlen_filter:
        df = df[df['SVLEN'].notna() & (df['SVLEN'] != '')]
        logger.info(f'过滤SVLEN非空后剩余: {len(df)} 条记录')
    else:
        logger.info(f'已禁用SVLEN过滤，跳过该条件检查')
    
    # 5. 对所有样本，过滤SR和GT条件
    for sample in sample_cols:
        # SR列：逗号前为0，逗号后>3（可通过--disable-sr-filter禁用
        sr_col = f'{sample}_SR'
        if sr_col in df.columns and not args.disable_sr_filter:
            # 分割SR值
            sr_split = df[sr_col].str.split(',', expand=True)
            # 确保有两个部分，第一个为0，第二个转换为数字>3
            valid_sr = (sr_split[0] == '0') & (pd.to_numeric(sr_split[1], errors='coerce') > 3)
            df = df[valid_sr]
            logger.info(f'样本{sample}过滤SR条件后剩余: {len(df)} 条记录')
        elif args.disable_sr_filter:
            logger.info(f'已禁用SR过滤，跳过样本{sample}的SR条件检查')
        
        # GT列：1/1 过滤（可通过--disable-gt-filter禁用）
        gt_col = f'{sample}_GT'
        if gt_col in df.columns and not args.disable_gt_filter:
            df = df[df[gt_col] == '1/1']
            logger.info(f'样本{sample}过滤GT==1/1后剩余: {len(df)} 条记录')
        elif args.disable_gt_filter:
            logger.info(f'已禁用GT过滤，跳过样本{sample}的GT条件检查')
    
    logger.info(f'过滤完成，共保留 {len(df)} 条记录 (总过滤掉 {initial_count - len(df)} 条)')
    
    # 第五步：选择需要保留的列（保持和过滤前完全一致的顺序）
    logger.info('开始选择输出列，保持和过滤前完全一致的列顺序')
    # 所有列已经在之前补全并调整好了顺序，直接使用即可
    all_cols = all_columns
    df = df[all_cols]
    logger.info(f'最终输出列共 {len(all_cols)} 个: {all_cols}')
    
    # 按照SVTYPE DEL→INS→DUP顺序排序（仅当INFO解析启用时）
    if info_parsed and 'SVTYPE' in df.columns:
        df['SVTYPE'] = pd.Categorical(df['SVTYPE'], categories=['DEL', 'INS', 'DUP'], ordered=True)
        df = df.sort_values('SVTYPE')
        logger.info('已按照SVTYPE DEL→INS→DUP顺序排序完成')
    else:
        logger.info('已禁用INFO解析，跳过SVTYPE排序，保持原始VCF顺序')
    
    # 输出结果
    logger.info(f'开始输出结果到: {args.output}')
    if args.output.lower().endswith(('.xlsx', '.xls')):
        write_excel(df, args.output)
    else:
        # 输出TSV格式
        write_output_df(df, args.output, index=False)
        logger.info(f'TSV格式输出完成，文件路径: {args.output}，若用Excel打开请通过「数据→从文本导入」选择文本列防止日期转换')
    
    # 输出过滤后的VCF文件
    output_filename = os.path.basename(args.output)
    vcf_filename = re.sub(r'\.(xlsx|xls)$', '.vcf', output_filename, flags=re.IGNORECASE)
    vcf_path = os.path.join(os.path.dirname(os.path.abspath(args.output)), vcf_filename)
    
    # 根据df的索引获取对应的原始VCF行
    filtered_vcf_lines = [original_vcf_lines[i] for i in df.index]
    
    # 写入VCF文件
    with open(vcf_path, 'w', encoding='utf-8') as f:
        for line in vcf_header:
            f.write(f'{line}\n')
        f.write(f'{header_line}\n')
        for line in filtered_vcf_lines:
            f.write(f'{line}\n')
    
    logger.info(f'过滤后VCF文件已输出到: {vcf_path}')
    logger.info('处理完成! 已输出过滤前原始数据、过滤后结果表和过滤后VCF三个文件')

if __name__ == '__main__':
    main()