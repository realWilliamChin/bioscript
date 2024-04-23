import os
import argparse
import pandas as pd


def parse_input():
    args = argparse.ArgumentParser(description='')
    args.add_argument('-i', '--input_dir', required=True, help='输入文件夹')
    args.add_argument('-r', '--ref', help='参考文件')
    return args.parse_args()
    

def embl2genesymbol(input_dir, ref_file):
    for each_file in os.listdir(input_dir):
        input_df = pd.read_csv(each_file, sep='\t')
        embl2genesymbol = pd.read_csv(ref_file, sep='\t', skiprows=1, names=['GeneID', 'GeneSymbol'])
        result_df = pd.merge(left=embl2genesymbol, right=input_df, on='GeneID', how='inner')
        result_df.drop(columns=['GeneID'], inplace=True)
        output_filename = each_file.replace('.txt', '_genesymbol.txt')
        result_df.drop_duplicates(subset='GeneSymbol', inplace=True)
        result_df.to_csv(input_dir + '/' + output_filename, sep='\t', index=False)


def main():
    args = parse_input()
    embl2genesymbol(args.input_dir, args.ref)


if __name__ == '__main__':
    main()