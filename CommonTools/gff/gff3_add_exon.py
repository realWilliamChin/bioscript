import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument('-i', type=str, dest='input', help='输入文件')
parser.add_argument('-o', type=str, dest='output', help='输出文件')

args = parser.parse_args()

with open(args.input, 'r', encoding='utf-8') as input_f, open(args.output, 'a+', encoding='utf-8') as output_f:
    for line in input_f.readlines():
        # print(line)
        if 'gene' in line:
            gene_id = line.split('ID=')[1].split(';')[0].strip()
            output_f.write(line)
            output_f.write(line.replace('gene', 'mRNA').replace('ID=', 'ID=mRNA-').strip() + f';Parent={gene_id}\n')
            output_f.write(line.replace('gene', 'CDS').replace('ID=', 'ID=CDS-').strip() + f';Parent={gene_id}\n')
            output_f.write(line.replace('gene', 'exon').replace('ID=', 'ID=exon-').strip() + f';Parent={gene_id}\n')
        else:
            output_f.write(line)