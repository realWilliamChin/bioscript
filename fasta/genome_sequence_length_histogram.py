from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse


def parse_input():
    argparser = argparse.ArgumentParser(description='Calculate genome statistics')
    argparser.add_argument('-i', '--input', type=str, required=True, help='Path to genome file')
    argparser.add_argument('-o', '--output', type=str, help='Path to output file, 默认是 input 文件去除后缀改为 jpeg')
    
    if not argparser.parse_args().output:
        argparser.parse_args().output = argparser.parse_args().input.split('.')[0] + '.jpeg'
    
    return argparser.parse_args()


def calculate_genome_statistics(genome_file):
    sequence_lengths = []
    for record in SeqIO.parse(genome_file, "fasta"):
        sequence_lengths.append(len(record.seq))

    length_counts = {}
    for length in sequence_lengths:
        length_category = length // 100
        length_counts[length_category] = length_counts.get(length_category, 0) + 1

    return length_counts


def plot_histogram(length_counts, output_file):
    x_values = [key * 100 for key in sorted(length_counts.keys())]
    y_values = [length_counts[key] for key in sorted(length_counts.keys())]

    plt.bar(x_values, y_values, width=100, align='edge', edgecolor='black')
    plt.xlabel('Protein Length')
    plt.ylabel('Number')
    plt.title('Protein Length Distribution')
    # 在每个条形上添加次数标签
    # for x, y in zip(x_values, y_values):
    #     plt.text(x + 50, y + 0.05, str(y), ha='center', va='bottom')

    plt.savefig(output_file)
    

if __name__ == '__main__':
    args = parse_input()
    length_counts = calculate_genome_statistics(args.input)
    plot_histogram(length_counts, args.output)
    print('Genome statistics calculated and histogram saved to', args.output)
