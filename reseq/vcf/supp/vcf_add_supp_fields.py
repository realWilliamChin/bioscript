#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/04/27 15:22
# Author        : William GoGo
import os, sys
import re
import argparse
from loguru import logger

def parse_args():
    parser = argparse.ArgumentParser(description='VCF添加SUPP和SUPP_VEC字段工具（类似SURVIVOR合并输出格式）')
    parser.add_argument('-i', '--input', required=True, help='输入已经合并的VCF文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出VCF文件路径')
    parser.add_argument('--gt-field', default='GT', help='指定FORMAT中基因型字段名称，默认GT')
    parser.add_argument('--consider-hom-ref', action='store_true', help='是否将0/0视为有变异（默认否，仅非0/0且非./.视为有变异）')
    
    return parser.parse_args()

def is_variant_present(gt_str, consider_hom_ref=False):
    """判断基因型是否表示存在变异"""
    if not gt_str:
        return False
    # 处理缺失基因型
    if gt_str in ['./.', '.|.', '.']:
        return False
    # 处理纯合参考
    if gt_str in ['0/0', '0|0'] and not consider_hom_ref:
        return False
    # 其他情况都认为有变异（0/1, 1/1, 1|0, 1/2等）
    return True

def process_vcf(input_path, output_path, gt_field='GT', consider_hom_ref=False):
    """处理VCF文件，添加SUPP和SUPP_VEC字段"""
    logger.info(f'开始处理VCF文件: {input_path}')
    
    vcf_header = []     # 保存所有##开头的注释行
    header_line = None  # 保存原始表头行（#CHROM开头）
    data_lines = []     # 保存所有数据行
    
    # 第一步：读取VCF文件
    with open(input_path, 'r', encoding='utf-8') as f:
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
    
    if not header:
        logger.error('未找到VCF表头行')
        sys.exit(1)
    
    logger.info(f'读取到 {len(data_lines)} 条变异记录')
    
    # 识别样本列
    if 'FORMAT' not in header:
        logger.error('VCF文件缺少FORMAT列，无法处理样本基因型')
        sys.exit(1)
    
    format_idx = header.index('FORMAT')
    sample_names = header[format_idx + 1:]
    num_samples = len(sample_names)
    logger.info(f'检测到 {num_samples} 个样本: {sample_names}')
    
    # 第二步：添加INFO字段头信息（如果不存在）
    info_supp_exists = any('##INFO=<ID=SUPP,' in line for line in vcf_header)
    info_supp_vec_exists = any('##INFO=<ID=SUPP_VEC,' in line for line in vcf_header)
    
    # 找到INFO头信息的位置，插入到所有INFO字段的最后
    insert_pos = len(vcf_header)
    for i, line in enumerate(vcf_header):
        if line.startswith('##INFO='):
            insert_pos = i + 1
    
    if not info_supp_exists:
        supp_info = '##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of samples supporting the variant">'
        vcf_header.insert(insert_pos, supp_info)
        logger.info('已添加SUPP字段的头信息')
        insert_pos += 1
    
    if not info_supp_vec_exists:
        supp_vec_info = '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Binary vector of samples supporting the variant (1=present, 0=absent)">'
        vcf_header.insert(insert_pos, supp_vec_info)
        logger.info('已添加SUPP_VEC字段的头信息')
    
    # 第三步：处理每条记录
    processed_lines = []
    for parts in data_lines:
        info = parts[7]
        format_str = parts[format_idx]
        sample_values = parts[format_idx + 1:]
        
        # 解析FORMAT字段，找到GT的位置
        format_fields = format_str.split(':')
        if gt_field not in format_fields:
            logger.warning(f'FORMAT字段中未找到{gt_field}，跳过该记录')
            processed_lines.append('\t'.join(parts))
            continue
        
        gt_idx = format_fields.index(gt_field)
        
        # 遍历所有样本，判断是否有变异
        supp_vec = []
        supp_count = 0
        for sample_val in sample_values:
            sample_fields = sample_val.split(':')
            if gt_idx >= len(sample_fields):
                # 基因型字段缺失
                supp_vec.append('0')
                continue
            
            gt = sample_fields[gt_idx]
            if is_variant_present(gt, consider_hom_ref):
                supp_vec.append('1')
                supp_count += 1
            else:
                supp_vec.append('0')
        
        supp_vec_str = ''.join(supp_vec)
        
        # 更新INFO字段：移除已有的SUPP和SUPP_VEC，添加到最前面
        info_items = info.split(';')
        new_info_items = []
        for item in info_items:
            if item.startswith('SUPP=') or item.startswith('SUPP_VEC='):
                continue
            if item:  # 跳过空项
                new_info_items.append(item)
        
        # 添加新的字段到最前面
        new_info_items.insert(0, f'SUPP={supp_count}')
        new_info_items.insert(1, f'SUPP_VEC={supp_vec_str}')
        new_info = ';'.join(new_info_items)
        
        # 替换原INFO字段
        parts[7] = new_info
        processed_lines.append('\t'.join(parts))
    
    # 第四步：写入输出文件
    logger.info(f'开始写入输出文件: {output_path}')
    with open(output_path, 'w', encoding='utf-8') as f:
        for line in vcf_header:
            f.write(f'{line}\n')
        f.write(f'{header_line}\n')
        for line in processed_lines:
            f.write(f'{line}\n')
    
    logger.info(f'处理完成！共处理 {len(data_lines)} 条记录，输出文件: {output_path}')
    logger.info(f'SUPP字段范围: 0-{num_samples}，SUPP_VEC长度: {num_samples}位')

def main():
    args = parse_args()
    process_vcf(
        input_path=args.input,
        output_path=args.output,
        gt_field=args.gt_field,
        consider_hom_ref=args.consider_hom_ref
    )

if __name__ == '__main__':
    main()