#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/08/05 10:55
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger

from load_input import load_table, write_output_df

def parse_input():
    p = argparse.ArgumentParser(description='GeneID <-> GeneSymbol 映射工具')
    p.add_argument('-i', '--input', required=True, help='输入文件，需包含 GeneID 或 GeneSymbol 列')
    p.add_argument('-r', '--ref', required=True, help='参考映射文件，需包含 GeneID GeneSymbol 两列')
    p.add_argument('-o', '--output', help='输出文件')
    p.add_argument('--def', action='store_true', default=False, dest='include_ref_cols',
                   help='输出文件是否包含参考文件的其他列')
    p.add_argument('--mode', choices=['id2symbol', 'symbol2id'], default='id2symbol', help='映射方向')
    p.add_argument('--action', choices=['add', 'replace'], default='add', help='添加新列或替换原列')
    p.add_argument('--dup', choices=['sum', 'max', 'min', 'first'], default='first', help='重复值处理方式')
    
    args = p.parse_args()
    
    return args


def handle_duplicates(df, group_col, agg_method):
    if agg_method == 'first':
        return df.drop_duplicates(subset=group_col, keep='first')
    elif agg_method == 'sum':
        return df.groupby(group_col, as_index=False).sum(numeric_only=True)
    elif agg_method == 'max':
        return df.groupby(group_col, as_index=False).max(numeric_only=True)
    elif agg_method == 'min':
        return df.groupby(group_col, as_index=False).min(numeric_only=True)
    else:
        return df


def map_id_symbol(input_file, ref_file, output_file, mode, action, dup, include_ref_cols):
    input_df = load_table(input_file)
    ref_df = load_table(ref_file)
    ref_df.columns = [col.strip() for col in ref_df.columns]
    if 'GeneID' not in ref_df.columns or 'GeneSymbol' not in ref_df.columns:
        raise ValueError('参考文件必须包含 GeneID 和 GeneSymbol 两列')

    if mode == 'id2symbol':
        if 'GeneID' not in input_df.columns:
            raise ValueError('输入文件需包含 GeneID 列')
        ref_cols = ref_df.columns if include_ref_cols else ['GeneID', 'GeneSymbol']
        merged = pd.merge(input_df, ref_df[ref_cols], on='GeneID', how='left')
        if action == 'replace':
            merged.drop('GeneID', axis=1, inplace=True)
        out_col = 'GeneSymbol' if action == 'replace' else None
    else:  # symbol2id
        if 'GeneSymbol' not in input_df.columns:
            raise ValueError('输入文件需包含 GeneSymbol 列')
        ref_cols = ref_df.columns if include_ref_cols else ['GeneID', 'GeneSymbol']
        merged = pd.merge(input_df, ref_df[ref_cols], on='GeneSymbol', how='left')
        if action == 'replace':
            merged.drop('GeneSymbol', axis=1, inplace=True)
        out_col = 'GeneID' if action == 'replace' else None

    # 处理重复
    if out_col:
        dup_count = merged.duplicated(subset=out_col).sum()
        if dup_count > 0:
            logger.info(f"检测到 {dup_count} 个重复的 {out_col}，按 '{dup}' 方式处理重复值。")
            merged = handle_duplicates(merged, out_col, dup)
        else:
            logger.info(f"未检测到 {out_col} 列重复，无需聚合处理。")
    write_output_df(merged, output_file, index=False)

def main():
    args = parse_input()
    if not args.output:
        base = os.path.splitext(os.path.basename(args.input))[0]
        if args.mode == 'id2symbol':
            args.output = base + '_with_symbol.txt' if args.action == 'add' else base + '_symbol.txt'
        else:
            args.output = base + '_with_id.txt' if args.action == 'add' else base + '_id.txt'
    map_id_symbol(args.input, args.ref, args.output, args.mode, args.action, args.dup, args.include_ref_cols)


if __name__ == '__main__':
    main()