import os
import argparse
import pandas as pd


def parse_input():
    args = argparse.ArgumentParser(description='')
    args.add_argument('-i', '--input', required=True, help='输入 all_gene_id.txt 文件')
    args.add_argument('-r', '--ref', help='参考文件')
    args.add_argument('-o', '--output', help='输出文件')
    return args.parse_args()


def embl2genesymbol(input_file, ref_file, output_file):
    input_df = pd.read_csv(input_file, sep='\t', names=['GeneID'])
    embl2genesymbol = pd.read_csv(ref_file, sep='\t', skiprows=1, names=['GeneID', 'GeneSymbol'])
    result_df = pd.merge(left=embl2genesymbol, right=input_df, on='GeneID', how='inner')
    # result_df.drop(columns=['GeneID'], inplace=True)
    result_df.drop_duplicates(subset='GeneID', inplace=True)
    result_df.to_csv(output_file, sep='\t', index=False)


def main():
    args = parse_input()
    if not args.output:
        args.output = args.input.split(os.sep)[-1].split('_all')[0] + '_embl2genesymbol.txt'
        
    embl2genesymbol(args.input, args.ref, args.output)


if __name__ == '__main__':
    main()