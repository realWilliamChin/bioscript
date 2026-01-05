#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/01/04 17:52
# Author        : William GoGo, Yu WangYin
import os, sys
import argparse
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table


def parse_input():
    parser = argparse.ArgumentParser(description="dumx_to_join")
    parser.add_argument('-i', '--input-dir', required=True, help='输入数据目录，包含R1和R2文件')
    parser.add_argument('-o', '--output-dir', required=True, help='输出数据目录')
    parser.add_argument('-s', '--sample-info', required=True, help='样本描述表文件路径，包含 group, sample, R1, R2 列')
    return parser.parse_args()


def main():
    args = parse_input()
    
    input_dir = args.input_dir
    output_dir = args.output_dir
    sample_info_file = args.sample_info
    
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f'创建输出目录: {output_dir}')
    
    # 读取样本描述表
    logger.info(f'读取样本描述表: {sample_info_file}')
    samples_df = load_table(sample_info_file)
    
    # 检查必需的列是否存在
    required_columns = ['group', 'sample', 'R1', 'R2']
    missing_columns = [col for col in required_columns if col not in samples_df.columns]
    if missing_columns:
        logger.error(f'样本描述表缺少必需的列: {missing_columns}')
        logger.error(f'当前表的列: {list(samples_df.columns)}')
        sys.exit(1)
    
    # 遍历每一行，处理每个样本
    logger.info(f'开始处理 {len(samples_df)} 个样本')
    for idx, row in samples_df.iterrows():
        group = str(row['group']).strip()
        sample = str(row['sample']).strip()
        r1_file = str(row['R1']).strip()
        r2_file = str(row['R2']).strip()
        
        # 构建完整的文件路径
        r1_path = os.path.join(input_dir, r1_file)
        r2_path = os.path.join(input_dir, r2_file)
        
        # 构建输出文件名（使用group或sample）
        output_prefix = os.path.join(output_dir, sample)
        
        # 构建flash命令
        command = f'flash {r1_path} {r2_path} -p 33 -r 250 -f 500 -s 100 -o {output_prefix}'
        
        logger.info(f'处理样本: {sample} (group: {group})')
        logger.debug(f'执行命令: {command}')
        
        # 执行flash命令
        exit_code = os.system(command)
        if exit_code == 0:
            logger.success(f'完成样本: {sample} (group: {group})')
        else:
            logger.error(f'处理样本失败: {sample} (group: {group}), 退出码: {exit_code}')
    
    logger.success('所有样本处理完成!')


if __name__ == '__main__':
    main()
