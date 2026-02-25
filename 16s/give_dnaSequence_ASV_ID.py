#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/01/22 13:52
# Author        : William GoGo
"""
将 FASTA 文件中的哈希值 ID 替换为 ASV ID
使用映射文件将序列 ID 从哈希值转换为 ASV1, ASV2, ASV3... 格式
"""

import argparse
import os
import sys
import pandas as pd
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table


def process_fasta_with_mapping(fasta_file, mapping_file, output_file):
    """
    处理 FASTA 文件，将哈希值 ID 替换为 ASV ID
    
    Args:
        fasta_file: 输入的 FASTA 文件路径（包含哈希值 ID 的序列文件）
        mapping_file: 映射文件路径（ASV_ID 和 Hash_ID 的对应关系）
        output_file: 输出文件路径
    """
    # 读取映射文件
    mapping_df = load_table(mapping_file, header=None)
    
    # 创建哈希值到 ASV ID 的字典（更高效的查找）
    if len(mapping_df.columns) >= 2:
        hash_to_asv = dict(zip(mapping_df.iloc[:, 1], mapping_df.iloc[:, 0]))
    else:
        raise ValueError(f"映射文件格式错误，应包含至少两列: {mapping_file}")
    
    logger.info(f"已加载 {len(hash_to_asv)} 个映射关系")
    
    # 处理 FASTA 文件
    replaced_count = 0
    not_found_count = 0
    
    with open(fasta_file, 'r', encoding='utf-8') as f_in, \
         open(output_file, 'w', encoding='utf-8') as f_out:
        
        for line in f_in:
            if line.startswith('>'):
                # 提取哈希值 ID（去除 '>' 和换行符）
                hash_id = line.strip()[1:]
                
                # 查找对应的 ASV ID
                if hash_id in hash_to_asv:
                    asv_id = hash_to_asv[hash_id]
                    f_out.write(f'>{asv_id}\n')
                    replaced_count += 1
                else:
                    # 如果找不到映射，保留原始 ID 并记录警告
                    f_out.write(line)
                    not_found_count += 1
                    logger.warning(f"未找到映射: {hash_id}")
            else:
                # 序列行直接写入
                f_out.write(line)
    
    logger.info(f"处理完成！")
    logger.info(f"输出文件: {output_file}")
    logger.info(f"成功替换: {replaced_count} 个序列 ID")
    if not_found_count > 0:
        logger.warning(f"未找到映射: {not_found_count} 个序列 ID")


def parse_input():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='将 FASTA 文件中的哈希值 ID 替换为 ASV ID')
    parser.add_argument('-i', '--input', dest='input', required=True, help='输入的 FASTA 文件路径（包含哈希值 ID 的序列文件）')
    parser.add_argument('-m', '--mapping', dest='mapping', required=True, help='映射文件路径（ASV_ID 和 Hash_ID 的对应关系）')
    parser.add_argument('-o', '--output', dest='output', required=True, help='输出文件路径')
    return parser.parse_args()


def main():
    """主函数"""
    args = parse_input()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"输入文件不存在: {args.input}")
    if not os.path.exists(args.mapping):
        raise FileNotFoundError(f"映射文件不存在: {args.mapping}")
    
    # 创建输出目录（如果不存在）
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    process_fasta_with_mapping(args.input, args.mapping, args.output)


if __name__ == '__main__':
    main()