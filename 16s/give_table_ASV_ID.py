#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/01/22 13:52
# Author        : William GoGo
"""
将 ASV 表格中的哈希值 ID 替换为 ASV1, ASV2, ASV3... 格式
并生成 ASV 编号与原始哈希值的映射文件
"""

import argparse
import os
import sys
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def process_asv_table(input_file, output_dir):
    """
    处理 ASV 表格，将哈希值替换为 ASV 编号
    
    Args:
        input_file: 输入文件路径
        output_dir: 输出目录路径
    """
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)
    
    # 定义输出文件路径
    output_table = os.path.join(output_dir, 'table_ASV.tsv')
    output_dict = os.path.join(output_dir, 'ASV_Hash_Dic.txt')
    
    # 使用 pandas 读取表格数据
    df = load_table(input_file)
    original_hashes = df[df.columns[0]].values.tolist()
    
    # 生成新的 ASV ID
    n_asvs = len(df)
    new_asv_ids = [f'ASV{i+1}' for i in range(n_asvs)]
    
    # 替换第一列
    df[df.columns[0]] = new_asv_ids
    
    write_output_df(df, output_table, index=False)
    
    # 生成并保存映射字典
    mapping_df = pd.DataFrame({
        'ASV_ID': new_asv_ids,
        'Hash_ID': original_hashes
    })
    write_output_df(mapping_df, output_dict, index=False)
    
    logger.info(f"处理完成！")
    logger.info(f"输出表格: {output_table}")
    logger.info(f"映射文件: {output_dict}")
    logger.info(f"共处理 {n_asvs} 个 ASV")


def parse_input():
    parser = argparse.ArgumentParser(description='将 ASV 表格中的哈希值 ID 替换为 ASV1, ASV2, ASV3... 格式')
    parser.add_argument('-i', '--input', dest='input', required=True, help='输入文件路径（包含 ASV 哈希值的表格文件）')
    parser.add_argument('-o', '--output-dir', required=True, help='输出目录路径')
    return parser.parse_args()


def main():
    """主函数"""
    args = parse_input()
    process_asv_table(args.input, args.output_dir)


if __name__ == '__main__':
    main()
