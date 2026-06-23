#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/06/20 14:41
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
import subprocess
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', help='输入 .out.dm 文件')
    p.add_argument('-o', help='输出文件')
    
    # 过滤参数
    p.add_argument('-e', '--evalue-max', type=float, default=1e-17, help='最高 evalue 值')
    p.add_argument('-c', '--coverage-min', type=float, default=0.45, help='最低 coverage 值')
    
    args = p.parse_args()
    
    return args


def parse_hmmscan_file(input_file: str) -> pd.DataFrame:
    """解析hmmscan输出文件"""
    tmp_file = 'tmp_hmmscan_parser.txt'
    cmd = [
        'bash', '/home/colddata/qinqiang/script/DenovoGenome/scripts/hmmscan-parser.sh', input_file
    ]
    with open(tmp_file, 'w') as fout:
        subprocess.run(cmd, stdout=fout, check=True)
    source_column_name = [
        'Target_name', 'HMM_target_length', 'Transcript_ID', 'HMM_query_length',
        'i-Evalue', 'HMM_start', 'HMM_end', 'Query_start', 'Query_end',
        'Accession_ID', 'Description', 'Coverage'
    ]
    filtered_name = [
        'Transcript_ID', 'Accession_ID', 'Target_name', 'HMM_target_length',
        'HMM_query_length', 'HMM_start', 'HMM_end', 'Query_start', 'Query_end',
        'Description', 'i-Evalue', 'Coverage'
    ]
    hmmscan_df = load_table(tmp_file, header=None, names=source_column_name)
    hmmscan_df = hmmscan_df[filtered_name]
    os.remove(tmp_file)
    return hmmscan_df


def filter_hmmscan(df: pd.DataFrame, evalue_max: float, coverage_min: float) -> None:
    df['i-Evalue'] = pd.to_numeric(df['i-Evalue'], errors='coerce')
    df['Coverage'] = pd.to_numeric(df['Coverage'], errors='coerce')

    filtered = df[
        (df['i-Evalue'] < evalue_max) & 
        (df['Coverage'] > coverage_min) &
        df['i-Evalue'].notna() &
        df['Coverage'].notna()
    ]
    
    return filtered


def main():
    args = parse_input()
    hmmscan_df = parse_hmmscan_file(args.i)
    filterd_hmmscan_df = filter_hmmscan(hmmscan_df, args.evalue_max, args.coverage_min)
    filterd_hmmscan_df = filterd_hmmscan_df.drop(columns=['i-Evalue', 'Coverage'])
    write_output_df(filterd_hmmscan_df, args.o, index=False)
    

if __name__ == "__main__":
    main()
