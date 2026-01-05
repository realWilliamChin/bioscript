import argparse
import os, sys
import pandas as pd
from pathlib import Path
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df


def parse_arguments():
    parser = argparse.ArgumentParser(description="创建manifest和group_info文件用于16S数据")
    parser.add_argument('-s', '--samples-info', default='samples_described.txt', help='样本描述表文件路径，包含 group, sample 列')
    parser.add_argument('-i', '--input-dir', default='01_Singledata', help='输入数据目录')
    parser.add_argument('--manifest', default='manifest.csv', help='输出文件路径')
    parser.add_argument('--group-info', default='group_info.tsv', help='输出文件路径')
    args = parser.parse_args()
    
    return args


def main():
    args = parse_arguments()
    samples_df = load_table(args.samples_info, usecols=[0, 1], header=0, names=['group', 'sample-id'])
    manifest_df = samples_df[['sample-id']].copy()
    manifest_df['absolute-filepath'] = manifest_df['sample-id'].apply(lambda x: os.path.join(os.getcwd(), args.input_dir, x + '.extendedFrags.fastq'))
    manifest_df['direction'] = 'forward'
    write_output_df(manifest_df, args.manifest, index=False)
    
    group_info_df = samples_df[['sample-id', 'group']].copy()
    group_info_df.columns = ['sample-id', 'group']
    write_output_df(group_info_df, args.group_info, index=False)


if __name__ == '__main__':
    main()