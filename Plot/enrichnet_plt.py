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
from matplotlib.colors import ListedColormap

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
):
    """
    Plot venn network.

    Parameters
    ----------
    :param edge_data: a dataFrame with three columns: source, target, color
    :param edge_width: width of edges
    :param edge_style: line or curve
    :param source_node_size: size of source nodes
    :param source_font_size: font size of source labels
    :param target_node_size: size of target nodes
    :param target_font_size: font size of target labels
    :param show_node_margin: whether to show margin of nodes when style is curve
    :param show_node_color: whether to show node color according to different source
    :param show_source_label: whether to show source labels
    :param show_target_label: whether to show target labels
    :param r: radius of the fixed nodes based on circle center (0, 0)
    :param k: optimal distance between nodes for spring layout
    :param ax: axes object to plot on (default: None)
    :return: axes object
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
        d["color"] for u, v, d in G.edges(data=True)
    ]  # G.edges(data=True)返回NodeDataView对象(一个嵌套列表例如:[('A', 'a', {"color":"red"}), ('B', 'b',{"color":"blue"})],获取边的颜色属性

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
        node_color.update(dict(zip(edge_data["source"], edge_data["color"])))
        edge_data_test = edge_data.copy()
        edge_data_test.loc[edge_data_test["target"].duplicated(), "color"] = (
            "red"  # 交集的点单独设置颜色
        )
        node_color.update(dict(zip(edge_data_test["target"], edge_data_test["color"])))
        node_options["node_color"] = [node_color[node] for node in G.nodes]
    nx.draw_networkx_nodes(G, pos, **node_options)
    nx.draw_networkx_edges(G, pos, **edge_options[edge_style])

    if not show_node_margin:
        for p in ax.patches:
            p.shrinkA = p.shrinkB = 0.0  # 箭头收缩值设置为零,删除空隙

    # show fixed nodes or source nodes with label
    if show_source_label:
        nx.draw_networkx_labels(
            G,
            pos,
            labels=dict(zip(fixed_nodes, fixed_nodes)),
            font_size=source_font_size,
            font_color="black",
            font_weight="bold",
            ax=ax,
        )
    # show target nodes with label
    if show_target_label:
        nx.draw_networkx_labels(
            G,
            pos,
            labels={node: node for node in G.nodes if node not in fixed_nodes},
            font_size=target_font_size,
            font_color="black",
            font_weight="normal",
            ax=ax,
        )
    ax.set_axis_off()

    return ax


def draw_enrichnetplot(df, output_pic_name):

    my_set_colors = [
        "#9c27b0", "#3f51b5", "#2196f3", "#00bcd4",
        "#009688", "#4caf50", "#8bc34a", "#cddc39",
        "#ffeb3b", "#ffc107", "#795548", "#607d8b",
    ]

    my_set_cmap = ListedColormap(my_set_colors)

    groups = list(set(df["source"].values.tolist()))
    groups_colors = plt.get_cmap(my_set_cmap).colors[: len(groups)]

    df = df.assign(color=lambda x: x["source"].map(dict(zip(groups, groups_colors))))

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
        target_font_size=7,
        show_target_label=True,
        k=0.15,
        ax=ax,
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
