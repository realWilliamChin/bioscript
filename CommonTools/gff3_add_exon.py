import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument('-i', type=str, dest='input', help='输入文件')
parser.add_argument('-o', type=str, dest='output', help='输出文件')

args = parser.parse_args()

with open(args.input, 'r', encoding='utf-8') as input_f, open(args.output, 'a+', encoding='utf-8') as output_f:
    for line in input_f.readlines():
        print(line)
        if 'gene' in line:
            output_f.write(line)
            output_f.write(line.replace('gene', 'exon').replace('ID=', 'ID=exon-'))
        else:
            output_f.write(line)