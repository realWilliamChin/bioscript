#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/05/12 17:47
# Author        : William GoGo
import os
import sys
import gzip
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import write_output_df

def get_vcf_columns(vcf_path):
    """从原始VCF文件中获取列名列表"""
    logger.info(f'读取原始VCF文件列名: {vcf_path}')
    
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    mode = 'rt' if vcf_path.endswith('.gz') else 'r'
    
    with open_func(vcf_path, mode, encoding='utf-8') as f:
        for line in f:
            if line.startswith('#CHROM'):
                columns = line.lstrip('#').strip().split('\t')
                logger.info(f'VCF原始列名: {columns}')
                logger.info(f'共{len(columns)}列')
                return columns
    logger.error('未在VCF文件中找到#CHROM表头行')
    sys.exit(1)

def restore_annovar_colnames(annovar_path, original_vcf_columns, output_path):
    """恢复ANNOVAR结果的Otherinfo列名"""
    logger.info(f'读取ANNOVAR结果文件: {annovar_path}')
    
    # 读取ANNOVAR结果
    df = pd.read_csv(annovar_path, sep='\t', dtype=str, low_memory=False)
    logger.info(f'ANNOVAR结果共 {len(df.columns)} 列，{len(df)} 行')
    
    # 找到所有Otherinfo列
    otherinfo_cols = []
    otherinfo_indexes = []
    for idx, col in enumerate(df.columns):
        if col.startswith('Otherinfo'):
            otherinfo_cols.append(col)
            otherinfo_indexes.append(idx)
    
    if not otherinfo_cols:
        logger.warning('未找到Otherinfo列，文件可能已经恢复过列名或不是ANNOVAR结果')
        write_output_df(df, output_path, index=False)
        return
    
    logger.info(f'找到 {len(otherinfo_cols)} 个Otherinfo列')
    logger.info(f'VCF原始列除了前5列(CHROM/POS/ID/REF/ALT)还有 {len(original_vcf_columns)-5} 列')
    
    # ANNOVAR的Otherinfo对应VCF中从第6列开始的所有列
    # VCF列顺序: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...
    # Otherinfo1对应QUAL，Otherinfo2对应FILTER，以此类推
    vcf_cols_to_map = original_vcf_columns[5:]  # 跳过前5列，因为ANNOVAR已经单独提取为Chr/Start/End/Ref/Alt
    
    # 校验列数是否匹配
    if len(otherinfo_cols) != len(vcf_cols_to_map):
        logger.warning(f'Otherinfo列数({len(otherinfo_cols)})与VCF剩余列数({len(vcf_cols_to_map)})不匹配，只替换匹配部分')
        min_len = min(len(otherinfo_cols), len(vcf_cols_to_map))
        otherinfo_cols = otherinfo_cols[:min_len]
        vcf_cols_to_map = vcf_cols_to_map[:min_len]
    
    # 创建列名映射
    col_mapping = dict(zip(otherinfo_cols, vcf_cols_to_map))
    logger.info(f'列名映射: {col_mapping}')
    
    # 重命名列
    df = df.rename(columns=col_mapping)
    
    # 处理ANNOVAR的位置列和VCF列名统一
    position_col_map = {
        'Chr': 'CHROM',
        'Start': 'POS',
        'Ref': 'REF',
        'Alt': 'ALT'
    }
    df = df.rename(columns={k: v for k, v in position_col_map.items() if k in df.columns})
    
    # 删除End列（和POS重复）
    if 'End' in df.columns:
        df = df.drop(columns=['End'])
    
    # 输出结果
    write_output_df(df, output_path, index=False)
    logger.info(f'恢复列名完成，结果已保存到: {output_path}')
    logger.info(f'最终文件列数: {len(df.columns)}')

def main():
    import argparse
    parser = argparse.ArgumentParser(description='ANNOVAR结果Otherinfo列名恢复工具 - 从原始VCF提取列名批量替换OtherinfoXX')
    parser.add_argument('-a', '--annovar-result', required=True, help='输入ANNOVAR注释结果文件(.txt/.tsv)')
    parser.add_argument('-v', '--original-vcf', required=True, help='对应的原始VCF文件(.vcf/.vcf.gz)')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径，支持.tsv/.txt格式')
    args = parser.parse_args()
    
    # 检查文件存在
    if not os.path.exists(args.annovar_result):
        logger.error(f'ANNOVAR结果文件不存在: {args.annovar_result}')
        sys.exit(1)
    if not os.path.exists(args.original_vcf):
        logger.error(f'原始VCF文件不存在: {args.original_vcf}')
        sys.exit(1)
    
    # 获取VCF列名
    vcf_columns = get_vcf_columns(args.original_vcf)
    
    # 恢复列名
    restore_annovar_colnames(args.annovar_result, vcf_columns, args.output)
    
    logger.info('处理完成!')

if __name__ == '__main__':
    main()