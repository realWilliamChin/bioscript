#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/01/06 13:54
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import numpy as np
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df
sys.path.append('/home/colddata/qinqiang/script/Rscript/')
from Rscript import volcano_plot


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-files', nargs='+',
        help='输入文件，多个文件用空格分隔, 组间比较的 VIP 文件')
    return parser.parse_args()


def main():
    args = parse_input()
    input_files = args.input_files

    logger.info(f"共检测到 {len(input_files)} 个输入文件:")
    for idx, fname in enumerate(input_files, 1):
        logger.info(f"  {idx}. {fname}")
    for idx, input_file in enumerate(input_files, 1):
        logger.info(f"[{idx}/{len(input_files)}] 正在处理文件: {input_file}")

        try:
            vip_df = load_table(input_file)
            logger.success(f"文件读取成功: {input_file}")
        except Exception as e:
            logger.error(f"文件读取失败: {input_file}, 错误: {e}")
            continue

        vip_df = vip_df.sort_values(by='VIP', ascending=False)
        
        # 添加Annotation列
        first_col = vip_df.columns[0]
        vip_df['Annotation'] = ''  # 初始化空列
        # 前15行使用第一列的值
        n_rows = min(15, len(vip_df))
        vip_df.iloc[:n_rows, vip_df.columns.get_loc('Annotation')] = vip_df.iloc[:n_rows][first_col].values
        logger.info(f"已成功添加Annotation列, 前{n_rows}行将用于标注")
        vip_df['log10FoldChange'] = np.log10(vip_df['FoldChange'])
        # 保存结果
        output_file = input_file.replace('.xlsx', '_volcano_plot.txt')
        vip_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f'数据已保存至临时文件: {output_file}')
        
        # 画图
        output_png = output_file.replace('.txt', '.png')
        logger.info(f"开始绘制火山图: {output_png}")
        plot_success = volcano_plot(
            input_file=output_file,
            output_file=output_png,
            label_col='Annotation',
            x_col='log10FoldChange',
            y_col='VIP',
            bs_neg=np.log10(0.8),
            bs_pos=np.log10(1.2),
            width=10,
            height=10,
            dpi=300,
            center_zero=True
        )
        if plot_success:
            logger.success(f"火山图绘制成功: {output_png}")
        else:
            logger.error(f"火山图绘制失败: {output_png}")

        # 删除临时文件
        try:
            os.remove(output_file)
            logger.info(f'已删除临时文件: {output_file}')
        except Exception as e:
            logger.warning(f'删除临时文件失败: {output_file}, 错误: {e}')


if __name__ == '__main__':
    main()