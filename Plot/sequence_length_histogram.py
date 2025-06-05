import os, sys
import argparse
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from loguru import logger


def parse_input():
    parser = argparse.ArgumentParser(description='Calculate genome statistics')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to genome file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file, 默认是 input 文件去除后缀改为 jpeg')
    parser.add_argument('-r', '--range', type=int, dest='range_num', default=100, help='区间数')
    parser.add_argument('--xlabel', default="Sequence Length")
    parser.add_argument('--ylabel', default="Count")
    parser.add_argument('--title', default='Sequence Length Distribution')
    
    args = parser.parse_args()
    
    if not args.output:
        args.output = args.input.split('.')[0] + '.jpeg'
    
    return args


def calculate_genome_statistics(genome_file, range_num):
    sequence_lengths = []
    for record in SeqIO.parse(genome_file, "fasta"):
        sequence_lengths.append(len(record.seq))

    length_counts = {}
    for length in sequence_lengths:
        length_category = length // range_num
        length_counts[length_category] = length_counts.get(length_category, 0) + 1

    return length_counts


def plot_histogram(length_counts, output_file, range_num, xlabel, ylabel, plot_title):
    x_values = [key * range_num for key in sorted(length_counts.keys())]  # Range
    y_values = [length_counts[key] for key in sorted(length_counts.keys())]  # Counts
    
    # 将 x_values 转换为区间字符串
    range_labels = [f"{x}-{x + range_num}" for x in x_values]
    
    # 创建 DataFrame
    stat_file = '.'.join(output_file.split('.')[:-1]) + '_stat.csv'
    df = pd.DataFrame({'Range': range_labels, 'Counts': y_values})
    df.to_csv(stat_file, index=False)
    

    plt.style.use("classic")  # 使用更现代的样式
    plt.figure(figsize=(12, 6))  # 调整画布大小
    # x_indices = range(len(x_values))
    
    # 动态设置条形宽度（根据区间大小自动适配）
    bar_width = range_num * 0.8  # 80% 的区间宽度，避免重叠
    
    bars = plt.bar(
        x_values,
        y_values,
        width=bar_width,
        align='center',  # 条形居中在区间
        edgecolor='black',
        linewidth=0.5,
        color='#4C72B0',  # 自定义颜色
        alpha=0.8  # 透明度
    )
    
    # 添加数据标签（仅在高度大于一定值时显示）
    for bar in bars:
        height = bar.get_height()
        if height > max(y_values) * 0.05:  # 仅显示高于 5% 最大值的标签
            plt.text(
                bar.get_x() + bar.get_width() / 2,
                height + 0.5,
                f'{height}',
                ha='center',
                va='bottom',
                fontsize=2.5,
                # rotation=45  # 旋转数据标签
            )
    
    plt.xlabel(xlabel, fontsize=12, labelpad=10)
    plt.ylabel(ylabel, fontsize=12, labelpad=10)
    plt.title(plot_title, fontsize=14, pad=20)
    
    plt.xlim(-0.5, len(x_values) - 0.5)
    
    # 动态调整 x 轴标签显示间隔
    if len(x_values) > 30:  # 如果数据点过多，间隔显示标签
        plt.xticks(x_values[::3], range_labels[::3], rotation=45, ha='right', fontsize=5)
    elif len(x_values) > 60:
        plt.xticks(x_values[::6], range_labels[::6], rotation=45, ha='right', fontsize=5)
    elif len(x_values) > 100:
        plt.xticks(x_values[::9], range_labels[::9], rotation=45, ha='right', fontsize=5)
    else:
        plt.xticks(x_values, range_labels, rotation=45, ha='right', fontsize=10)
    plt.yticks(fontsize=10)
    
    # 添加网格和边框控制
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    
    # 自动调整布局防止标签重叠
    plt.tight_layout()
    # 在每个条形上添加次数标签
    # for x, y in zip(x_values, y_values):
    #     plt.text(x + 50, y + 0.05, str(y), ha='center', va='bottom')

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close


def main():
    args = parse_input()
    length_counts = calculate_genome_statistics(args.input, args.range_num)
    plot_histogram(length_counts, args.output, args.range_num, args.xlabel, args.ylabel, args.title)
    logger.success(f'Genome statistics calculated and histogram saved to {args.output}')   
    

if __name__ == '__main__':
    main()
