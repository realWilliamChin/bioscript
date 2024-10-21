from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse


def parse_input():
    parser = argparse.ArgumentParser(description='Calculate genome statistics')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to genome file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file, 默认是 input 文件去除后缀改为 jpeg')
    parser.add_argument('-r', '--range', type=int, dest='range_num', default=100, help='区间数')
    parser.add_argument('--xlabel', default="Sequence Length")
    parser.add_argument('--ylabel', default="Number")
    parser.add_argument('--title', default='Protein Length Distribution')
    
    if not parser.parse_args().output:
        parser.parse_args().output = parser.parse_args().input.split('.')[0] + '.jpeg'
    
    return parser.parse_args()


def calculate_genome_statistics(genome_file, range_num=100):
    sequence_lengths = []
    for record in SeqIO.parse(genome_file, "fasta"):
        sequence_lengths.append(len(record.seq))

    length_counts = {}
    for length in sequence_lengths:
        length_category = length // range_num
        length_counts[length_category] = length_counts.get(length_category, 0) + 1

    return length_counts


def plot_histogram(length_counts, output_file, range_num, xlabel, ylabel, plot_title):
    x_values = [key * range_num for key in sorted(length_counts.keys())]
    y_values = [length_counts[key] for key in sorted(length_counts.keys())]

    print(x_values)
    print(y_values)

    plt.bar(x_values, y_values, width=100, align='edge', edgecolor='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(plot_title)
    # 在每个条形上添加次数标签
    # for x, y in zip(x_values, y_values):
    #     plt.text(x + 50, y + 0.05, str(y), ha='center', va='bottom')

    plt.savefig(output_file)
    

if __name__ == '__main__':
    args = parse_input()
    length_counts = calculate_genome_statistics(args.input, args.range_num)
    plot_histogram(length_counts, args.output, args.range_num, args.xlabel, args.ylabel, args.title)
    print('Genome statistics calculated and histogram saved to', args.output)
