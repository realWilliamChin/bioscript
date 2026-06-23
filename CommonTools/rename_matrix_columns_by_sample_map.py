#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/05/19 11:00
# Author        : William GoGo
import os
import sys
import argparse
import pandas as pd
from loguru import logger
from load_input import load_table, write_output_df


def parse_args():
    parser = argparse.ArgumentParser(description="根据样本描述表重命名矩阵列名")
    parser.add_argument("-s", "--sample-map", required=True, help="样本描述表文件路径，包含group、sample、old_sample三列")
    parser.add_argument("-d", "--data-matrix", required=True, help="原始矩阵文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出重命名后的矩阵文件路径")
    return parser.parse_args()


def main():
    args = parse_args()
    
    # 读取样本映射表
    logger.info(f"读取样本描述表: {args.sample_map}")
    sample_map_df = load_table(args.sample_map)
    
    # 检查必要列是否存在
    required_columns = ['group', 'sample', 'old_sample']
    missing_cols = [col for col in required_columns if col not in sample_map_df.columns]
    if missing_cols:
        logger.error(f"样本描述表缺少必要列: {', '.join(missing_cols)}")
        sys.exit(1)
    
    # 检查重复的old_sample
    duplicated_old = sample_map_df[sample_map_df['old_sample'].duplicated(keep=False)]
    if not duplicated_old.empty:
        logger.warning(f"样本描述表中存在重复的old_sample: {', '.join(duplicated_old['old_sample'].unique())}")
        logger.warning("将使用第一个出现的映射关系")
    
    # 检查重复的sample
    duplicated_sample = sample_map_df[sample_map_df['sample'].duplicated(keep=False)]
    if not duplicated_sample.empty:
        logger.warning(f"样本描述表中存在重复的sample名称: {', '.join(duplicated_sample['sample'].unique())}")
    
    # 建立old_sample到sample的映射
    sample_map = dict(zip(sample_map_df['old_sample'], sample_map_df['sample']))
    logger.info(f"成功加载 {len(sample_map)} 个样本映射关系")
    
    # 读取矩阵表
    logger.info(f"读取原始矩阵: {args.data_matrix}")
    matrix_df = load_table(args.data_matrix)
    logger.info(f"原始矩阵包含 {matrix_df.shape[0]} 行, {matrix_df.shape[1]} 列")
    
    # 找出矩阵中存在的old_sample列
    existing_old_samples = [col for col in sample_map.keys() if col in matrix_df.columns]
    missing_old_samples = [col for col in sample_map.keys() if col not in matrix_df.columns]
    
    if missing_old_samples:
        logger.warning(f"矩阵中不存在以下old_sample列: {', '.join(missing_old_samples)}")
    
    if not existing_old_samples:
        logger.error("矩阵中没有找到任何匹配的old_sample列，无法继续")
        sys.exit(1)
    
    logger.info(f"在矩阵中找到 {len(existing_old_samples)} 个匹配的old_sample列")
    
    # 构建新的列名映射
    new_columns = {}
    for col in matrix_df.columns:
        if col in sample_map:
            new_columns[col] = sample_map[col]
        else:
            new_columns[col] = col
    
    # 重命名列
    renamed_df = matrix_df.rename(columns=new_columns)
    
    # 输出结果
    logger.info(f"输出重命名后的矩阵到: {args.output}")
    write_output_df(renamed_df, args.output, index=False)
    logger.info(f"处理完成，输出矩阵包含 {renamed_df.shape[0]} 行, {renamed_df.shape[1]} 列")


if __name__ == '__main__':
    main()