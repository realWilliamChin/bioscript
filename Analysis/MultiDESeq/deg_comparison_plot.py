#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2025/04/16 16:37
# Author        : WilliamGoGo
# Description   : 绘制差异表达基因(DEGs)比较统计图 (优化版)

import os
import sys
import argparse
import warnings

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table  # noqa: E402


# ---------- 模块级常量 ----------
DEFAULT_UP_COLOR = '#E64B35'      # Nature 风格红
DEFAULT_DOWN_COLOR = '#4DBBD5'    # Nature 风格青蓝
EXPECTED_COLS = ['Comparisons', 'Total DEGs', 'Up regulated', 'Down regulated']
NUMERIC_COLS = ['Total DEGs', 'Up regulated', 'Down regulated']


def parse_input():
    parser = argparse.ArgumentParser(
        description='绘制差异表达基因(DEGs)比较统计图',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-i', type=str, default='DEG_summary.txt',
                        help='输入文件 DEG_summary.txt')
    parser.add_argument('-o', type=str, default='DEGs_comparison_count_barplot.png',
                        help='输出文件 (扩展名决定输出格式, 推荐 .png/.pdf/.svg)')
    parser.add_argument('--dpi', type=int, default=300, help='图片 DPI')
    parser.add_argument('--up-color', dest='up_color', type=str,
                        default=DEFAULT_UP_COLOR, help='Up regulated 柱子颜色')
    parser.add_argument('--down-color', dest='down_color', type=str,
                        default=DEFAULT_DOWN_COLOR, help='Down regulated 柱子颜色')
    parser.add_argument('--orient', type=str, default='v', choices=['v', 'h'],
                        help='柱状图方向: v=竖直, h=水平')
    parser.add_argument('--sort', type=str, default='none',
                        choices=['none', 'asc', 'desc'],
                        help='按 Total DEGs 排序: none/asc/desc')
    parser.add_argument('--interactive', action='store_true',
                        help='额外生成交互式 HTML 图表 (Plotly)')
    args = parser.parse_args()
    return args


# ---------- 数据预处理 ----------
def validate_and_prepare(deg_summary_df, sort='none'):
    """校验输入、列数检查、数值转换（带 NaN 警告）、排序"""
    if deg_summary_df is None or deg_summary_df.empty:
        raise ValueError('输入数据为空, 请检查 DEG_summary 文件内容')

    if deg_summary_df.shape[1] != len(EXPECTED_COLS):
        raise ValueError(
            f'输入文件列数应为 {len(EXPECTED_COLS)} ({EXPECTED_COLS}), '
            f'实际为 {deg_summary_df.shape[1]} ({list(deg_summary_df.columns)})'
        )

    deg_summary_df = deg_summary_df.copy()
    deg_summary_df.columns = EXPECTED_COLS

    # 数值转换 + NaN 检测
    for col in NUMERIC_COLS:
        converted = pd.to_numeric(deg_summary_df[col], errors='coerce')
        nan_mask = converted.isna()
        if nan_mask.any():
            bad_rows = deg_summary_df.loc[nan_mask, 'Comparisons'].tolist()
            warnings.warn(
                f'列 "{col}" 中存在无法解析为数值的项, '
                f'已自动填充为 0; 涉及比较组: {bad_rows}'
            )
        deg_summary_df[col] = converted.fillna(0).astype(int)

    # 排序
    if sort == 'asc':
        deg_summary_df = deg_summary_df.sort_values(
            'Total DEGs', ascending=True).reset_index(drop=True)
    elif sort == 'desc':
        deg_summary_df = deg_summary_df.sort_values(
            'Total DEGs', ascending=False).reset_index(drop=True)

    return deg_summary_df


def _to_long_format(deg_summary_df):
    """转换为绘图用的长表格式"""
    return pd.melt(
        deg_summary_df,
        id_vars=['Comparisons'],
        value_vars=['Up regulated', 'Down regulated'],
        var_name='Regulation',
        value_name='Count',
    )


# ---------- 静态图（matplotlib + seaborn） ----------
def deg_summary_plot(deg_summary_df, output_file, *,
                     up_color=DEFAULT_UP_COLOR,
                     down_color=DEFAULT_DOWN_COLOR,
                     orient='v', dpi=300):
    """绘制静态 DEGs 柱状图, 支持竖直 (v) 与水平 (h) 两种方向

    Legend 统一放置在图外右侧.
    """
    sns.set_theme(style='whitegrid', font_scale=1.0)

    plot_df = _to_long_format(deg_summary_df)
    palette = {'Up regulated': up_color, 'Down regulated': down_color}
    n_groups = deg_summary_df.shape[0]
    max_count = max(plot_df['Count'].max(), 1)

    if orient == 'v':
        # 竖直柱状图: 宽度随组数变化, 高度固定; 多预留宽度给图外 legend
        fig_width = max(5.0, min(20.0, n_groups * 1.2 + 1.0)) + 1.5
        fig_height = 6.0
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        sns.barplot(
            data=plot_df, x='Comparisons', y='Count',
            hue='Regulation', palette=palette, ax=ax,
        )

        # 数值标签 (智能错开): 比较同一比较组下 Up/Down 高度差, 较小者向右偏移
        _annotate_vertical_bars(ax, plot_df)

        ax.set_xlabel('')
        ax.set_ylabel('Number of Genes')
        ax.set_ylim(0, max_count * 1.15)

        # x 轴标签旋转 (避免弃用警告: 使用 tick_params + 手动 ha)
        ax.tick_params(axis='x', rotation=45)
        for label in ax.get_xticklabels():
            label.set_ha('right')
            label.set_rotation_mode('anchor')

    else:  # orient == 'h'
        # 水平条形图: 高度随组数, 宽度固定
        fig_width = 9.0 + 1.5  # 多预留 legend 宽度
        fig_height = max(4.0, n_groups * 0.5 + 1.5)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        sns.barplot(
            data=plot_df, y='Comparisons', x='Count',
            hue='Regulation', palette=palette, ax=ax, orient='h',
        )

        _annotate_horizontal_bars(ax, plot_df)

        ax.set_ylabel('')
        ax.set_xlabel('Number of Genes')
        ax.set_xlim(0, max_count * 1.15)

    # 通用设置
    ax.set_title('Differentially Expressed Genes (DEGs)',
                 fontsize=14, fontweight='bold')
    sns.despine(ax=ax, top=True, right=True)

    # Legend 统一放置在图外右侧
    ax.legend(
        title=None,
        loc='upper left',
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
        frameon=False,
    )

    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def _annotate_vertical_bars(ax, plot_df):
    """为竖直柱状图添加数值标签, 当同组 Up/Down 高度接近时智能错开偏移"""
    # 计算高度阈值: 用于判断是否会重叠 (粗略经验值)
    max_count = max(plot_df['Count'].max(), 1)
    close_threshold = max_count * 0.05

    # 按 Comparisons 分组 Up/Down 的高度
    pivot = plot_df.pivot(index='Comparisons', columns='Regulation',
                          values='Count')

    for p in ax.patches:
        height = int(p.get_height())
        if height <= 0:
            continue

        x_center = p.get_x() + p.get_width() / 2
        # 默认偏移
        x_offset = 0
        y_offset = 3

        # 判断当前 patch 属于哪个 Regulation (通过颜色或位置)
        # seaborn barplot 中, hue 相同的 patch 颜色一致; 此处简化处理:
        # 当 Up/Down 高度接近时, 都用较小水平偏移避免标签横向重叠
        try:
            comparison = ax.get_xticklabels()[
                int(round(p.get_x() + p.get_width() / 2 - p.get_width() / 2))
            ].get_text()
            if comparison in pivot.index:
                up = pivot.loc[comparison, 'Up regulated']
                down = pivot.loc[comparison, 'Down regulated']
                if abs(up - down) < close_threshold and min(up, down) > 0:
                    # 高度接近: 根据柱子位置略微左右偏移
                    bar_center = p.get_x() + p.get_width() / 2
                    group_center = round(bar_center)
                    if bar_center < group_center:
                        x_offset = -6
                    else:
                        x_offset = 6
        except (IndexError, KeyError, ValueError):
            pass

        ax.annotate(
            f'{height}',
            (x_center, height),
            ha='center', va='bottom', fontsize=9,
            xytext=(x_offset, y_offset), textcoords='offset points',
        )


def _annotate_horizontal_bars(ax, plot_df):
    """为水平条形图添加数值标签 (放在柱子右侧)"""
    for p in ax.patches:
        width = int(p.get_width())
        if width <= 0:
            continue
        y_center = p.get_y() + p.get_height() / 2
        ax.annotate(
            f'{width}',
            (width, y_center),
            ha='left', va='center', fontsize=9,
            xytext=(3, 0), textcoords='offset points',
        )


# ---------- 交互式图表 (Plotly, 独立函数) ----------
def deg_summary_interactive_plot(deg_summary_df, output_file, *,
                                  up_color=DEFAULT_UP_COLOR,
                                  down_color=DEFAULT_DOWN_COLOR,
                                  orient='v'):
    """生成交互式 HTML 图表 (Plotly)"""
    try:
        import plotly.graph_objects as go
    except ImportError:
        warnings.warn(
            'plotly 未安装, 无法生成交互式图表; 请运行: pip install plotly'
        )
        return

    comparisons = deg_summary_df['Comparisons'].tolist()
    up_values = deg_summary_df['Up regulated'].tolist()
    down_values = deg_summary_df['Down regulated'].tolist()
    total_values = deg_summary_df['Total DEGs'].tolist()

    hover_tpl = (
        '<b>%{x}</b><br>'
        if orient == 'v'
        else '<b>%{y}</b><br>'
    )
    hover_tpl += (
        '%{fullData.name}: %{customdata[0]}<br>'
        'Total DEGs: %{customdata[1]}<extra></extra>'
    )

    fig = go.Figure()

    if orient == 'v':
        fig.add_trace(go.Bar(
            name='Up regulated', x=comparisons, y=up_values,
            marker_color=up_color,
            customdata=list(zip(up_values, total_values)),
            hovertemplate=hover_tpl,
            text=up_values, textposition='outside',
        ))
        fig.add_trace(go.Bar(
            name='Down regulated', x=comparisons, y=down_values,
            marker_color=down_color,
            customdata=list(zip(down_values, total_values)),
            hovertemplate=hover_tpl,
            text=down_values, textposition='outside',
        ))
        fig.update_layout(
            xaxis_title='', yaxis_title='Number of Genes',
            xaxis_tickangle=-45,
        )
    else:
        fig.add_trace(go.Bar(
            name='Up regulated', y=comparisons, x=up_values,
            marker_color=up_color, orientation='h',
            customdata=list(zip(up_values, total_values)),
            hovertemplate=hover_tpl,
            text=up_values, textposition='outside',
        ))
        fig.add_trace(go.Bar(
            name='Down regulated', y=comparisons, x=down_values,
            marker_color=down_color, orientation='h',
            customdata=list(zip(down_values, total_values)),
            hovertemplate=hover_tpl,
            text=down_values, textposition='outside',
        ))
        fig.update_layout(
            xaxis_title='Number of Genes', yaxis_title='',
            yaxis_autorange='reversed',
        )

    fig.update_layout(
        title=dict(
            text='<b>Differentially Expressed Genes (DEGs)</b>',
            x=0.5, xanchor='center',
        ),
        barmode='group',
        template='simple_white',
        legend=dict(
            x=1.02, y=1.0, xanchor='left', yanchor='top',
            title=None,
        ),
        margin=dict(l=80, r=160, t=80, b=80),
    )

    fig.write_html(output_file)


# ---------- 入口 ----------
def main():
    args = parse_input()

    if not os.path.isfile(args.i):
        sys.exit(f'[ERROR] 输入文件不存在: {args.i}')

    deg_summary_df = load_table(args.i, comment='#', skipinitialspace=True)
    deg_summary_df = validate_and_prepare(deg_summary_df, sort=args.sort)

    deg_summary_plot(
        deg_summary_df, args.o,
        up_color=args.up_color, down_color=args.down_color,
        orient=args.orient, dpi=args.dpi,
    )

    if args.interactive:
        html_file = os.path.splitext(args.o)[0] + '.html'
        deg_summary_interactive_plot(
            deg_summary_df, html_file,
            up_color=args.up_color, down_color=args.down_color,
            orient=args.orient,
        )


if __name__ == '__main__':
    main()