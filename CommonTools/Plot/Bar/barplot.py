#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/04/28 16:42
# Author        : William GoGo

import os
import sys
import argparse
import matplotlib.pyplot as plt
from loguru import logger

# 添加CommonTools到路径
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../'))
from CommonTools.load_input import load_table


def parse_args():
    parser = argparse.ArgumentParser(description='绘制柱状图脚本，支持表格输入或直接数字输入')
    
    # 输入参数组
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-f', '--file', help='输入表格文件路径，支持csv/tsv/xlsx/xls格式')
    input_group.add_argument('-d', '--data', help='直接输入逗号分隔的数字，例如: "10,20,30,15,25"')
    
    # 表格相关参数
    parser.add_argument('-x', '--x-col', help='表格中作为x轴的列名，仅当使用-file参数时需要')
    parser.add_argument('-y', '--y-col', help='表格中作为y轴的列名，仅当使用-file参数时需要')
    
    # 输出参数
    parser.add_argument('-o', '--output', help='输出图片路径，支持png/jpg/pdf/svg等格式，如果不提供则直接显示图片')
    
    # 图表样式参数
    parser.add_argument('-t', '--title', default='Bar Plot', help='图表标题，默认: "Bar Plot"')
    parser.add_argument('--xlabel', help='x轴标签，默认使用x列名或"类别"')
    parser.add_argument('--ylabel', help='y轴标签，默认使用y列名或"数值"')
    parser.add_argument('-c', '--color', default='#1f77b4', help='柱子颜色，默认: #1f77b4 (蓝色)')
    parser.add_argument('--edgecolor', default='black', help='柱子边框颜色，默认: black')
    parser.add_argument('--horizontal', action='store_true', help='绘制横向柱状图，默认纵向')
    parser.add_argument('--show-value', action='store_true', help='在柱子上显示数值标签')
    parser.add_argument('--rotate', type=int, default=0, help='x轴标签旋转角度，默认0度')
    parser.add_argument('--figsize', nargs=2, type=int, default=(10, 6), help='图片尺寸，默认: 10 6')
    parser.add_argument('--dpi', type=int, default=300, help='输出图片DPI，默认: 300')
    
    return parser.parse_args()


def main():
    args = parse_args()
    
    # 加载数据
    logger.info('开始加载数据...')
    if args.file:
        if not args.x_col or not args.y_col:
            logger.error('使用文件输入时必须指定-x/--x-col和-y/--y-col参数')
            sys.exit(1)
            
        df = load_table(args.file)
        logger.info('成功加载文件: {}, 共 {} 行数据'.format(args.file, len(df)))
        
        # 检查列是否存在
        if args.x_col not in df.columns:
            logger.error('x轴列名 {} 不存在于表格中，可用列: {}'.format(args.x_col, list(df.columns)))
            sys.exit(1)
        if args.y_col not in df.columns:
            logger.error('y轴列名 {} 不存在于表格中，可用列: {}'.format(args.y_col, list(df.columns)))
            sys.exit(1)
            
        x_data = df[args.x_col].astype(str).tolist()
        y_data = df[args.y_col].tolist()
        
        # 设置默认标签
        xlabel = args.xlabel if args.xlabel else args.x_col
        ylabel = args.ylabel if args.ylabel else args.y_col
        
    else:  # 直接数字输入
        try:
            y_data = [float(num.strip()) for num in args.data.split(',')]
        except ValueError as e:
            logger.error('解析数字失败: {}, 请输入正确的逗号分隔数字'.format(e))
            sys.exit(1)
            
        x_data = ['{}'.format(i+1) for i in range(len(y_data))]
        xlabel = args.xlabel if args.xlabel else '类别'
        ylabel = args.ylabel if args.ylabel else '数值'
        logger.info('成功解析数字序列，共 {} 个数值'.format(len(y_data)))
    
    # 绘制图表
    logger.info('开始绘制柱状图...')
    fig, ax = plt.subplots(figsize=args.figsize)
    
    if args.horizontal:
        bars = ax.barh(x_data, y_data, color=args.color, edgecolor=args.edgecolor)
    else:
        bars = ax.bar(x_data, y_data, color=args.color, edgecolor=args.edgecolor)
    
    # 设置标题和标签
    ax.set_title(args.title, fontsize=14, pad=20)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    
    # 旋转x轴标签
    plt.xticks(rotation=args.rotate, ha='right' if args.rotate > 0 else 'center')
    
    # 显示数值标签
    if args.show_value:
        if args.horizontal:
            for bar in bars:
                width = bar.get_width()
                ax.text(width + max(y_data)*0.01, bar.get_y() + bar.get_height()/2,
                       '{:.2f}'.format(width) if isinstance(width, float) else '{}'.format(width),
                       va='center', fontsize=10)
        else:
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2, height + max(y_data)*0.01,
                       '{:.2f}'.format(height) if isinstance(height, float) else '{}'.format(height),
                       ha='center', fontsize=10)
    
    # 调整布局
    plt.tight_layout()
    
    # 输出或显示
    if args.output:
        plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
        logger.info('图片已保存到: {}'.format(args.output))
    else:
        logger.info('显示图片...')
        plt.show()
    
    logger.info('柱状图绘制完成!')


if __name__ == '__main__':
    main()