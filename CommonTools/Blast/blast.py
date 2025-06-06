#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2025/05/08 15:39
# Author        : William GoGo

import os
import sys
import argparse
import pandas as pd
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(description='BLAST 比对工具')
    
    # 必需参数
    parser.add_argument('--blast-type', dest='blast_type', type=str, required=True, choices=['blastn', 'blastx', 'blastp'],
                      help='BLAST 比对类型: blastn, blastx, blastp')
    parser.add_argument('--query', type=str, required=True,
                      help='查询序列文件路径')
    parser.add_argument('--db', type=str, required=True,
                      help='数据库文件路径')
    parser.add_argument('--out', type=str, required=True,
                      help='输出文件路径')
    
    # 可选参数
    parser.add_argument('--threads', type=int, default=32,
                      help='线程数 (默认: 32)')
    parser.add_argument('--evalue', type=float, default=1e-5,
                      help='E-value 阈值 (默认: 1e-5)')
    parser.add_argument('--outfmt', type=str, 
                      default='6 qacc sacc qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids',
                      help='输出格式 (默认: 6 qacc sacc qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids)')
    parser.add_argument('--sort-by', dest='sort_by', type=str, default='staxids',
                      help='排序字段 (默认: staxids)')
    parser.add_argument('--keyword', type=str, default=None,
                      help='关键词过滤 (例如: "unknown" 或 "DNA polymerase III subunit")')
    
    return parser.parse_args()


def run_blast(args) -> None:
    """
    执行 BLAST 比对
    """
    cmd = [
        args.blast_type,
        '-db', args.db,
        '-query', args.query,
        '-out', args.out,
        '-evalue', str(args.evalue),
        '-num_threads', str(args.threads),
        '-outfmt', args.outfmt
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"BLAST 比对完成，结果保存在: {args.out}")
    except subprocess.CalledProcessError as e:
        print(f"BLAST 比对失败: {e}")
        sys.exit(1)


def process_results(args) -> None:
    """
    处理 BLAST 结果
    """
    try:
        # 读取 BLAST 结果
        df = pd.read_csv(args.out, sep='\t', header=None)
        
        # 根据 outfmt 设置列名
        columns = args.outfmt.split()[1:]
        df.columns = columns
        
        # 关键词过滤
        if args.keyword:
            df = df[df['stitle'].str.contains(args.keyword, case=False, na=False)]
        
        # 排序
        if args.sort_by in df.columns:
            df = df.sort_values(by=args.sort_by)
        
        # 保存处理后的结果
        output_file = f"{os.path.splitext(args.out)[0]}_processed.txt"
        df.to_csv(output_file, sep='\t', index=False)
        print(f"处理后的结果保存在: {output_file}")
        
    except Exception as e:
        print(f"处理结果时出错: {e}")
        sys.exit(1)


def main():
    args = parse_arguments()
    run_blast(args)
    process_results(args)


if __name__ == '__main__':
    main()
    