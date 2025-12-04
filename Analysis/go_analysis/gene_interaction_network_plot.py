#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/10/31 14:56
# Author        : William GoGo
import os, sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import axes

from tabulate import tabulate
from IPython.display import Markdown

import networkx as nx

if sys.version_info < (3, 10):
    print("Python 版本低于 3.10，请升级到 3.10 或更高版本！")
    print("当前 Python 版本为:", sys.version)
    sys.exit(1)

def parse_input():
    argparser = argparse.ArgumentParser(description="")
    argparser.add_argument("-i", "--input", type=str, required=True,
        help="输入文件，格式为第一列 source (SubOntology)，第二列 target (GeneID)",
    )
    argparser.add_argument('-o', dest='out_pic_name', help='输出图片名称')

    return argparser.parse_args()


def VennNetworkPlot(
    edge_data: pd.DataFrame,
    edge_width: int | float = 1.0,
    edge_style: int = 1,
    plot_title: str = "Differential Gene Interaction Network",
    source_node_size: int | float = 100,
    source_font_size: int = 10,
    target_node_size: int | float = 5,
    target_font_size: int = 10,
    show_node_margin: bool = True,
    show_node_color: bool = False,
    show_source_label: bool = True,
    show_target_label: bool = False,
    r: float = 0.5,
    k: float = None,
    ax: axes.Axes = None,
    max_label_length: int = 20,
    label_background: bool = True,
    label_padding: float = 0.1,
):
    """
    Plot venn network.

    Parameters
    ----------
    :param edge_data: 包含三列的数据框：source（源节点）、target（目标节点）和color（颜色）
    :param edge_width: 边的宽度
    :param edge_style: 边的样式（直线或曲线）
    :param plot_title: 图表标题
    :param source_node_size: 源节点的大小
    :param source_font_size: 源节点标签的字体大小
    :param target_node_size: 目标节点的大小
    :param target_font_size: 目标节点标签的字体大小
    :param show_node_margin: 当样式为曲线时是否显示节点边距
    :param show_node_color: 是否根据不同的源节点显示不同的节点颜色
    :param show_source_label: 是否显示源节点标签
    :param show_target_label: 是否显示目标节点标签
    :param r: 基于圆心(0,0)的固定节点半径
    :param k: spring布局中节点之间的最佳距离
    :param ax: 用于绘图的坐标轴对象（默认：None）
    :param max_label_length: 标签换行前的最大长度
    :param label_background: 是否显示标签背景
    :param label_padding: 标签文本周围的内边距
    :return: 坐标轴对象
    """
    # get axes object 坐标轴范围限制
    if ax is None:
        ax = plt.gca()

    # 创建 Title
    font = {"color": "k", "fontweight": "bold", "fontsize": 16}
    ax.set_title(plot_title, font)
    # create networkx graph
    # 该函数需要三个参数：源节点名称列、目标节点名称列和边权重列。这里构建单个无向网络即可。
    G = nx.from_pandas_edgelist(
        df=edge_data,  # 图的边用列表表示
        edge_attr=[
            "color"
        ],  # 如果True，将添加所有剩余的列。如果None，则不向图中添加边属性
        create_using=nx.Graph(),  # NetworkX 图构造函数，可选（默认=nx.Graph）。要创建的图表类型。如果是图形实例，则在填充之前清除。
    )

    # get fixed nodes from source column
    fixed_nodes = edge_data["source"].unique()
    # print(fixed_nodes)

    # set position of fixed nodes
    # 以(0,0)为原点,以r为半径画圆,sinx^2 + cosx^2 = 1,从(0，1)开始按分组数均分取坐标点，
    fixed_nodes_pos = {
        node: (
            r * np.sin(2 * np.pi * i / len(fixed_nodes)),
            r * np.cos(2 * np.pi * i / len(fixed_nodes)),
        )
        for i, node in enumerate(fixed_nodes)
    }

    # get node style
    node_size = [
        (
            target_node_size if u not in fixed_nodes else source_node_size
        )  # 根据顶点的id分配设置的size大小
        for u, d in G.nodes(
            data=True
        )  # G.nodes(data=True)返回的是NodeDataView对象(一个嵌套字典例如:{'A': {}, 'B': {}}),该对象不仅包含每个顶点的id属性，还包括顶点的其他属性。
    ]
    node_options = {
        "node_color": "lightgrey",
        "node_size": node_size,
    }

    # get edge style
    edge_color = [
        d.get("color", "gray") if pd.notna(d.get("color")) else "gray"
        for u, v, d in G.edges(data=True)
    ]  # G.edges(data=True)返回NodeDataView对象(一个嵌套列表例如:[('A', 'a', {"color":"red"}), ('B', 'b',{"color":"blue"})],获取边的颜色属性，如果颜色为nan则使用默认灰色

    edge_style_line = {
        "edge_color": edge_color,
        "style": "solid",
        "alpha": 1,
        "width": edge_width,
    }
    edge_style_curve = {
        "edge_color": edge_color,
        "style": "solid",
        "alpha": 1,
        "width": edge_width,
        "arrows": True,  # If True, draw arrowheads with FancyArrowPatches (bendable and stylish). If False, draw edges using LineCollection (linear and fast).
        "arrowsize": 10,
        "arrowstyle": None,  # 不加箭头
        "connectionstyle": "arc3,rad=0.3",  # 传递连接样式参数以创建圆角半径为 rad 的弧形,https://matplotlib.org/stable/gallery/userdemo/connectionstyle_demo.html
    }
    edge_options = {"line": edge_style_line, "curve": edge_style_curve}

    # get node layout
    pos = nx.spring_layout(
        G,
        pos=fixed_nodes_pos,
        fixed=fixed_nodes,
        seed=111,
        k=k,
        iterations=200,
    )  # spring layout布局使用Fruchterman-Reingold force-directed算法定位节点，k节点之间的最佳距离。如果没有，则距离设置为 1/sqrt(n)，其中 n 是节点数。增加此值可将节点移得更远, fixed和pos设置节点保持固定的节点，pos为坐标，fixed为id。不在G.nodes 中的节点将被忽略，seed是随机种子。iterations是最大迭代次数

    # plot network
    if show_node_color:
        node_color = {}
        # 为 source 节点分配颜色（过滤掉 nan 值）
        source_color_dict = dict(zip(edge_data["source"], edge_data["color"]))
        source_color_dict = {k: v for k, v in source_color_dict.items() if pd.notna(v)}
        node_color.update(source_color_dict)
        
        edge_data_test = edge_data.copy()
        # 为重复的 target 节点设置红色（交集节点）
        edge_data_test.loc[edge_data_test["target"].duplicated(), "color"] = "red"
        # 为 target 节点分配颜色（过滤掉 nan 值，未分配的节点使用默认颜色）
        target_color_dict = dict(zip(edge_data_test["target"], edge_data_test["color"]))
        target_color_dict = {k: v for k, v in target_color_dict.items() if pd.notna(v)}
        node_color.update(target_color_dict)
        
        # 为没有颜色的节点设置默认颜色（灰色）
        default_color = "lightgrey"
        node_options["node_color"] = [
            node_color.get(node, default_color) if pd.notna(node_color.get(node, default_color)) else default_color
            for node in G.nodes
        ]
    nx.draw_networkx_nodes(G, pos, **node_options)
    nx.draw_networkx_edges(G, pos, **edge_options[edge_style])

    if not show_node_margin:
        for p in ax.patches:
            p.shrinkA = p.shrinkB = 0.0  # 箭头收缩值设置为零,删除空隙

    def wrap_label(label, max_length):
        """将长标签自动换行"""
        if len(label) <= max_length:
            return label
        words = label.split()
        lines = []
        current_line = []
        current_length = 0
        
        for word in words:
            if current_length + len(word) + 1 <= max_length:
                current_line.append(word)
                current_length += len(word) + 1
            else:
                lines.append(' '.join(current_line))
                current_line = [word]
                current_length = len(word)
        
        if current_line:
            lines.append(' '.join(current_line))
        
        return '\n'.join(lines)

    def adjust_font_size(nodes_count):
        """根据节点数量自动调整字体大小"""
        if nodes_count < 50:
            return source_font_size, target_font_size
        elif nodes_count < 100:
            return max(8, source_font_size - 2), max(6, target_font_size - 2)
        else:
            return max(6, source_font_size - 4), max(4, target_font_size - 4)

    # 根据节点数量调整字体大小
    total_nodes = len(G.nodes)
    adjusted_source_font_size, adjusted_target_font_size = adjust_font_size(total_nodes)

    # 处理标签
    if show_source_label:
        source_labels = {node: wrap_label(node, max_label_length) for node in fixed_nodes}
        nx.draw_networkx_labels(
            G,
            pos,
            labels=source_labels,
            font_size=adjusted_source_font_size,
            font_color="black",
            font_weight="bold",
            ax=ax,
            bbox=dict(
                facecolor='white',
                edgecolor='none',
                alpha=0.7 if label_background else 0,
                pad=label_padding,
                boxstyle='round,pad=0.3'
            )
        )

    if show_target_label:
        target_labels = {
            node: wrap_label(node, max_label_length)
            for node in G.nodes if node not in fixed_nodes
        }
        nx.draw_networkx_labels(
            G,
            pos,
            labels=target_labels,
            font_size=adjusted_target_font_size,
            font_color="black",
            font_weight="normal",
            ax=ax,
            bbox=dict(
                facecolor='white',
                edgecolor='none',
                alpha=0.7 if label_background else 0,
                pad=label_padding,
                boxstyle='round,pad=0.3'
            )
        )

    ax.set_axis_off()

    return ax


def draw_enrichnetplot(df, output_pic_name):

    my_set_colors = [
        "#9c27b0", "#3f51b5", "#2196f3", "#00bcd4",
        "#009688", "#4caf50", "#8bc34a", "#cddc39",
        "#ffeb3b", "#ffc107", "#795548", "#607d8b",
    ]

    groups = list(set(df["source"].values.tolist()))
    # 为每个 group 分配颜色，如果 group 数量超过颜色数量，则循环使用
    groups_colors = [my_set_colors[i % len(my_set_colors)] for i in range(len(groups))]
    groups_color_dict = dict(zip(groups, groups_colors))
    
    # 分配颜色，确保所有行都有有效的颜色（使用 fillna 处理可能的缺失值）
    df = df.assign(color=lambda x: x["source"].map(groups_color_dict).fillna("#808080"))

    if df.shape[0] < 50:
        plotsize = (8, 8)
    elif df.shape[0] < 150:
        plotsize = (12, 12)
    elif df.shape[0] < 250:
        plotsize = (16, 16)
    elif df.shape[0] < 350:
        plotsize = (20, 20)
    elif df.shape[0] < 500:
        plotsize = (25, 25)
    else:
        plotsize = (30, 30)

    fig, ax = plt.subplots(figsize=plotsize)
    VennNetworkPlot(
        df,
        edge_style="line",
        source_node_size=200,
        target_node_size=100,
        show_node_color=True,
        target_font_size=6,
        show_target_label=True,
        k=0.25,
        ax=ax,
        plot_title=os.path.basename(output_pic_name).split('.')[0],
        max_label_length=12,
        label_background=False,
        label_padding=0.2
    )
    plt.savefig(output_pic_name)


def main():
    args = parse_input()

    if args.input.endswith(".txt"):
        df = pd.read_csv(args.input, sep="\t")
    elif args.input.endswith(".xlsx"):
        df = pd.read_excel(args.input)
    elif args.input.endswith(".csv"):
        df = pd.read_csv(args.input)
    else:
        print("输入文件格式不支持，请输入 txt 或 xlsx 或 csv 格式结尾的文件")
        sys.exit(1)

    df = df.rename(columns={"SubOntology": "source", "GeneID": "target"})
    
    draw_enrichnetplot(df, args.out_pic_name)


if __name__ == "__main__":
    main()
