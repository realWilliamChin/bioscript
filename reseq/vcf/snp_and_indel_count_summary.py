#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/11/28 10:06
# Author        : William GoGo
"""
统计多个样本的SNP、SNP_filtered、INDEL、INDEL_filtered文件中去除#开头的行数并生成汇总表
输入文件格式：每行5列（制表符分隔）：样本名称、SNP文件路径、SNP_filtered文件路径、INDEL文件路径、INDEL_filtered文件路径
或者每行4列：SNP文件路径、SNP_filtered文件路径、INDEL文件路径、INDEL_filtered文件路径（样本名称从SNP文件名提取）
"""
import os, sys
import argparse
import pandas as pd
from loguru import logger
from pathlib import Path

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools/'))
from load_input import load_table, write_output_df


def count_non_comment_lines(vcf_file):
    """
    统计VCF文件中非注释行（不以#开头）的数量
    
    Args:
        vcf_file (str): VCF文件路径
        
    Returns:
        int: 非注释行数量，文件不存在返回0
    """
    if not os.path.exists(vcf_file):
        logger.warning(f"文件不存在: {vcf_file}")
        return 0
    
    count = 0
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    count += 1
    except Exception as e:
        logger.error(f"读取文件失败 {vcf_file}: {str(e)}")
        return 0
    
    return count


def extract_sample_name_from_path(file_path):
    """
    从文件路径中提取样本名称
    例如：/path/to/Sample1_SNP.vcf -> Sample1
    
    Args:
        file_path (str): 文件路径
        
    Returns:
        str: 样本名称
    """
    file_name = Path(file_path).stem  # 获取不带扩展名的文件名
    # 移除常见的后缀
    for suffix in ['_SNP', '_SNP_filtered', '_INDEL', '_INDEL_filtered']:
        if file_name.endswith(suffix):
            return file_name[:-len(suffix)]
    return file_name


def parse_input():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='统计多个样本的SNP、SNP_filtered、INDEL、INDEL_filtered文件中去除#开头的行数并生成汇总表',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
输入文件格式（制表符分隔）：
方式1（推荐）：每行5列
样本名称\\tSNP文件路径\\tSNP_filtered文件路径\\tINDEL文件路径\\tINDEL_filtered文件路径

方式2：每行4列（样本名称从SNP文件名自动提取）
SNP文件路径\\tSNP_filtered文件路径\\tINDEL文件路径\\tINDEL_filtered文件路径

示例：
Sample1\\t/path/to/Sample1_SNP.vcf\\t/path/to/Sample1_SNP_filtered.vcf\\t/path/to/Sample1_INDEL.vcf\\t/path/to/Sample1_INDEL_filtered.vcf
Sample2\\t/path/to/Sample2_SNP.vcf\\t/path/to/Sample2_SNP_filtered.vcf\\t/path/to/Sample2_INDEL.vcf\\t/path/to/Sample2_INDEL_filtered.vcf
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='输入文件路径，包含所有样本的文件路径列表（格式见下方说明）'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='输出文件路径（支持txt, csv, xls, xlsx格式）'
    )
    
    return parser.parse_args()


def load_sample_files(input_file):
    """
    从输入文件加载样本和对应的文件路径
    
    Args:
        input_file (str): 输入文件路径
        
    Returns:
        list: 包含字典的列表，每个字典包含样本名称和4个文件路径
    """
    sample_data = []
    
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):  # 跳过空行和注释行
                continue
            
            parts = line.split('\t')
            
            # 判断是5列还是4列格式
            if len(parts) == 5:
                # 方式1：样本名称、SNP、SNP_filtered、INDEL、INDEL_filtered
                sample_name, snp_file, snp_filtered_file, indel_file, indel_filtered_file = parts
            elif len(parts) == 4:
                # 方式2：从SNP文件名提取样本名称
                snp_file, snp_filtered_file, indel_file, indel_filtered_file = parts
                sample_name = extract_sample_name_from_path(snp_file)
            else:
                logger.warning(f"第 {line_num} 行格式不正确（应有4或5列，实际{len(parts)}列），跳过: {line}")
                continue
            
            # 去除文件路径两端的空格
            snp_file = snp_file.strip()
            snp_filtered_file = snp_filtered_file.strip()
            indel_file = indel_file.strip()
            indel_filtered_file = indel_filtered_file.strip()
            
            sample_data.append({
                'sample_name': sample_name,
                'snp_file': snp_file,
                'snp_filtered_file': snp_filtered_file,
                'indel_file': indel_file,
                'indel_filtered_file': indel_filtered_file
            })
    
    return sample_data


def count_sample_variants(sample_info):
    """
    统计单个样本的4个VCF文件中的变异数量
    
    Args:
        sample_info (dict): 包含样本名称和4个文件路径的字典
        
    Returns:
        dict: 包含样本名称和4个统计值的字典
    """
    result = {
        'Sample': sample_info['sample_name'],
        'SNP': count_non_comment_lines(sample_info['snp_file']),
        'SNP_filtered': count_non_comment_lines(sample_info['snp_filtered_file']),
        'INDEL': count_non_comment_lines(sample_info['indel_file']),
        'INDEL_filtered': count_non_comment_lines(sample_info['indel_filtered_file'])
    }
    
    logger.info(f"样本 {result['Sample']}: SNP={result['SNP']}, SNP_filtered={result['SNP_filtered']}, "
                f"INDEL={result['INDEL']}, INDEL_filtered={result['INDEL_filtered']}")
    
    return result


def main():
    """主函数"""
    args = parse_input()
    
    # 加载样本文件列表
    sample_data = load_sample_files(args.input)
    
    if not sample_data:
        logger.error("未找到任何样本数据，请检查输入文件格式")
        sys.exit(1)
    
    logger.info(f"共加载 {len(sample_data)} 个样本的文件路径")
    
    # 统计每个样本的变异数量
    results = []
    for sample_info in sample_data:
        result = count_sample_variants(sample_info)
        results.append(result)
    
    # 生成汇总表
    df = pd.DataFrame(results)
    
    # 调整列顺序：Sample, SNP, SNP_filtered, INDEL, INDEL_filtered
    column_order = ['Sample', 'SNP', 'SNP_filtered', 'INDEL', 'INDEL_filtered']
    df = df[column_order]
    
    # 输出结果
    write_output_df(df, args.output, index=False)
    logger.info(f"汇总表已保存至: {args.output}")
    logger.info(f"共统计 {len(sample_data)} 个样本")


if __name__ == '__main__':
    main()
