#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/01/23
# Author        : William GoGo
"""
Alpha Rarefaction 曲线图绘制工具
用于绘制 QIIME 2 风格的 Alpha Rarefaction 曲线图，展示不同测序深度下的物种丰富度（Chao1）变化
"""

import json
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from loguru import logger


def load_jsonp_data(jsonp_file):
    """
    读取 jsonp 文件并转换为 DataFrame
    
    Args:
        jsonp_file: jsonp 文件路径
        
    Returns:
        DataFrame: 包含 rarefaction 数据的 DataFrame
    """
    logger.info(f"正在读取文件: {jsonp_file}")
    
    if not os.path.exists(jsonp_file):
        raise FileNotFoundError(f"文件不存在: {jsonp_file}")
    
    with open(jsonp_file, 'r', encoding='utf-8') as f:
        txt = f.read()
    
    # jsonp -> json（去掉 JS 包装）
    json_str = txt[txt.find("{"):txt.rfind("}") + 1]
    obj = json.loads(json_str)
    
    # 转成 DataFrame
    df = pd.DataFrame(obj["data"], columns=obj["columns"])
    
    # 重命名字段，方便使用
    df = df.rename(columns={
        "_alpha_rarefaction_depth_column_": "depth",
        "50%": "median",
        "25%": "q25",
        "75%": "q75"
    })
    
    logger.info(f"成功读取数据，包含 {len(df)} 行，{len(df['sample-id'].unique())} 个样本")
    
    return df


def plot_alpha_rarefaction(df, output_file, figsize=(16, 9), dpi=300, 
                           box_half_width=120, line_width=1.2, 
                           alpha_index="Chao1"):
    """
    绘制 Alpha Rarefaction 曲线图（QIIME 2 风格）
    
    Args:
        df: 包含 rarefaction 数据的 DataFrame
        output_file: 输出图片文件路径
        figsize: 图片尺寸，默认 (16, 9)
        dpi: 图片分辨率，默认 300
        box_half_width: boxplot 宽度的一半，默认 120
        line_width: 线条宽度，默认 1.2
        alpha_index: Alpha 多样性指数名称，默认 "Chao1"
    """
    logger.info(f"开始绘制 Alpha Rarefaction 图: {alpha_index}")
    
    # 获取样本列表和颜色映射
    samples = df["sample-id"].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(samples)))
    color_map = dict(zip(samples, colors))
    
    logger.info(f"共 {len(samples)} 个样本需要绘制")
    
    # 创建图形
    fig, ax = plt.subplots(figsize=figsize)
    
    # 绘制每个样本的曲线
    for sample_id, sub in df.groupby("sample-id"):
        sub = sub.sort_values("depth")
        c = color_map[sample_id]
        
        # 1️⃣ 中位数连线
        ax.plot(
            sub["depth"],
            sub["median"],
            color=c,
            linewidth=1,
            alpha=0.9,
            label=sample_id
        )
        
        # 2️⃣ 每个 depth 的空心 boxplot
        for _, r in sub.iterrows():
            x = r["depth"]
            
            # box（25%–75%）
            rect = Rectangle(
                (x - box_half_width, r["q25"]),
                2 * box_half_width,
                r["q75"] - r["q25"],
                fill=False,
                edgecolor=c,
                linewidth=line_width
            )
            ax.add_patch(rect)
            
            # median 横线
            ax.hlines(
                r["median"],
                x - box_half_width,
                x + box_half_width,
                color=c,
                linewidth=line_width
            )
            
            # whisker（如果有 2% / 98%）
            if "2%" in df.columns and "98%" in df.columns:
                ax.vlines(
                    x,
                    r["2%"],
                    r["98%"],
                    color=c,
                    linewidth=line_width * 0.8
                )
    
    # 设置坐标轴标签和标题
    ax.set_xlabel("Sequencing Depth")
    ax.set_ylabel(alpha_index)
    ax.set_title(f"Alpha Rarefaction ({alpha_index})")
    
    # 3️⃣ 右侧图例（颜色 ↔ sample）
    # 调整布局，为图例留出空间
    fig.subplots_adjust(right=0.75)
    
    # 添加图例
    legend = ax.legend(
        title="Sample",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
        fancybox=True,
        shadow=False,
        fontsize=8 if len(samples) > 10 else 9
    )
    
    # 保存图片
    output_dir = os.path.dirname(output_file) if os.path.dirname(output_file) else '.'
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"创建输出目录: {output_dir}")
    
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    
    logger.success(f"图片已保存: {output_file}")


def parse_input():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='绘制 Alpha Rarefaction 曲线图（QIIME 2 风格）',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-i', '--input', help='输入的 jsonp 文件路径（QIIME 2 alpha rarefaction 输出文件）')
    parser.add_argument('-o', '--output', help='输出图片文件路径')
    parser.add_argument('--alpha-index', default='Chao1', help='Alpha 多样性指数名称，默认为 Chao1')
    parser.add_argument('--figsize', nargs=2, type=float, default=[16, 9], metavar=('WIDTH', 'HEIGHT'), help='图片尺寸（宽 高），默认 16 9')
    parser.add_argument('--dpi', type=int, default=300, help='图片分辨率，默认 300')
    parser.add_argument('--box-width', type=float, default=120, help='boxplot 宽度的一半，默认 120')
    parser.add_argument('--line-width', type=float, default=1.2, help='线条宽度，默认 1.2')
    
    return parser.parse_args()


def main():
    """主函数"""
    args = parse_input()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        logger.error(f"输入文件不存在: {args.input}")
        raise FileNotFoundError(f"输入文件不存在: {args.input}")
    
    # 读取数据
    df = load_jsonp_data(args.input)
    
    # 绘制图形
    plot_alpha_rarefaction(
        df=df,
        output_file=args.output,
        figsize=tuple(args.figsize),
        dpi=args.dpi,
        box_half_width=args.box_width,
        line_width=args.line_width,
        alpha_index=args.alpha_index
    )
    
    logger.success('Done!')


if __name__ == '__main__':
    main()
