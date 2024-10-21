import argparse
import pandas as pd


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default='rfam_out.tab', help='输入 rfam_out.tab 文件')
    parser.add_argument('-o', '--output', default='rfam.txt', help='输出文件')
    
    return parser.parse_args()



def main():
    args = parse_input()
    
    df = pd.read_csv(args.input, comment='#', delim_whitespace=True, header=None, names=[
        'TargetName', 'Accession', 'QueryName', 'Accession2', 'Mdl', 'MdlFrom', 'MdlTo',
        'SeqFrom', 'SeqTo', 'Strand', 'Trunc', 'Pass', 'GC', 'Bias', 'Score', 'E-value',
        'Inc', 'DescriptionOfTarget'
    ])

    df.to_csv(args.output, sep='\t', index=False)
    
    filter_df = df[['TargetName', 'QueryName', 'SeqFrom', 'SeqTo', 'Strand']]
    filter_df.columns = ['GeneID', 'QueryName', 'Start', 'End', 'Strand']
    filter_df.to_csv('rfam_processed.txt', sep='\t', index=False)
    

if __name__ == '__main__':
    main()