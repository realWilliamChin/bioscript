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
    parser.add_argument('-i', '--input-files', nargs='+', help='输入文件，多个文件用空格分隔, 组间比较的 VIP 文件')
    
    parser.add_argument('--bs-neg', type=float, default=np.log10(0.8), help='左侧阈值，默认 log10(0.8)')
    parser.add_argument('--bs-pos', type=float, default=np.log10(1.2), help='右侧阈值，默认 log10(1.2)')
    parser.add_argument('--width', type=float, default=10, help='图片宽度(inches)，默认 10')
    parser.add_argument('--height', type=float, default=10, help='图片高度(inches)，默认 10')
    parser.add_argument('--dpi', type=int, default=300, help='图片 DPI，默认 300')
    parser.add_argument('--x-left', type=float, default=-2, help='x 轴左侧范围，默认 -2')
    parser.add_argument('--x-right', type=float, default=2, help='x 轴右侧范围，默认 2')
    parser.add_argument('--title', default='', help='图标题，默认空')
    parser.add_argument('--label-col', default='Annotation', help='标注列名，默认 Annotation')
    parser.add_argument('--center-zero', action='store_true', help='是否以 0 为中心对称显示 x 轴')
    parser.add_argument('--top-n', type=int, default=15, help='标注前多少个非 NoSignificant 的条目，默认 15；设置为 0 则不标注')

    return parser.parse_args()


def df_add_annotation(df, top_n=15):
    # 添加Annotation列
    df = df.sort_values(by='VIP', ascending=False)
    first_col = df.columns[0]
    df['Annotation'] = ''  # 初始化空列
    non_sig_mask = df['regulation'] != 'NoSignificant'
    annot_idx = df[non_sig_mask].head(top_n).index
    df.loc[annot_idx, 'Annotation'] = df.loc[annot_idx, first_col].values
    logger.info(f"已成功添加Annotation列, 共标注 {len(annot_idx)} 个非 NoSignificant 条目")
    
    return df


def main():
    args = parse_input()
    input_files = args.input_files

    logger.info(f"共检测到 {len(input_files)} 个输入文件:")
    for idx, fname in enumerate(input_files, 1):
        logger.info(f"  {idx}. {fname}")
    for idx, input_file in enumerate(input_files, 1):
        logger.info(f"[{idx}/{len(input_files)}] 正在处理文件: {input_file}")

        vip_df = load_table(input_file)
        
        if args.top_n > 0:
            vip_df = df_add_annotation(vip_df, top_n=args.top_n)
        
        vip_df['log10FoldChange'] = np.log10(vip_df['FoldChange'])
        output_file = input_file.replace('.xlsx', '_volcano_plot.txt')
        vip_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f'数据已保存至临时文件: {output_file}')
        
        # 画图
        output_png = output_file.replace('.txt', '.png')
        logger.info(f"开始绘制火山图: {output_png}")
        volcano_plot(
            input_file=output_file,
            output_file=output_png,
            label_col=args.label_col,
            x_col='log10FoldChange',
            y_col='VIP',
            bs_neg=args.bs_neg,
            bs_pos=args.bs_pos,
            width=args.width,
            height=args.height,
            dpi=args.dpi,
            x_left=args.x_left,
            x_right=args.x_right,
            title=args.title,
            center_zero=args.center_zero
        )

        # 删除临时文件
        os.remove(output_file)
        logger.info(f'已删除临时文件: {output_file}')


if __name__ == '__main__':
    main()