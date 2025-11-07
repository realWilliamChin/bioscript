#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/07/03 16:50
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd
from loguru import logger
from typing import Optional


sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='对各种表添加基因定义文件')
    parser.add_argument('-i', '--input', help='输入文件，可根据文件后缀格式类型读取，txt是tab分隔')
    parser.add_argument('--input-header', help="默认 GeneID 添加定义，如有其他列名，写列名，如没列名，输入列位置，从 0 开始数，列名不可以是数字")
    parser.add_argument('--defi', help='输入 def 文件')
    parser.add_argument('-c', '--specify-column', help='指定 def 根据哪一列进行合并，针对需要对 GeneSymbol 列进行合并优化')
    parser.add_argument('-k', '--kegg', help='KEGG_gene_def 文件，会添加 KEGG_ID 和 KEGG_Shortname 列')
    parser.add_argument('-n', '--nr', help='NR_gene_def 文件')
    parser.add_argument('-s', '--swiss', help='Swiss_gene_def 文件')
    parser.add_argument('-o', '--output', default='output.txt', help='输出文件')
    parser.add_argument('--merge-how', default='left', help='合并方式，默认 left')
    parser.add_argument('--dup-col', default='left', choices=['left', 'right', 'drop_all', 'keep_both'],
                        help='合并之后如果有重复的列，选择保留方式，left 保留输入文件的，right 保留 defi 文件的，drop_all 不保留任何一个，keep_both 保留所有')
    parser.add_argument('--no-remove-duplicates', action='store_true', help='不删除重复数据（默认会删除重复数据）')
    parser.add_argument('--na-fill', default='N/A', help='缺失值填充值（字符串列），默认 N/A')
    
    args = parser.parse_args()

    if not args.input_header and not args.specify_column:
        args.input_header = args.specify_column = 'GeneID'
    elif not args.input_header:
        args.input_header = args.specify_column
    elif not args.specify_column:
        args.specify_column = 'GeneID'
    
    # 处理去重参数：默认删除重复，使用 --no-remove-duplicates 来关闭
    args.remove_duplicates = not args.no_remove_duplicates
    
    return args


def add_def(
    file_df: pd.DataFrame,
    index_column: str = 'GeneID',
    kegg_file: Optional[str] = None,
    nr_file: Optional[str] = None,
    swiss_file: Optional[str] = None,
    def_file: Optional[str] = None,
    merge_how: str = 'left',
    remove_duplicates: bool = True,
    dup_col: str = 'left',
    na_fill: str = 'N/A'
    ) -> pd.DataFrame:
    """对输入表添加基因定义"""
    result_df = file_df.copy()
    source_shape = file_df.shape[0]

    if nr_file:
        nr_df = load_table(nr_file, header=0, names=['GeneID', 'NR_ID', 'NR_Def1'], dtype={'GeneID': str})
        # 没有 NCBI ID 直接加会丢失。先 fillna，再去掉 NA::
        nr_df.fillna(value=na_fill, inplace=True)
        nr_df['NR_ID_Des'] = nr_df['NR_ID'] + '::' + nr_df['NR_Def1']
        nr_df['NR_ID_Des'] = nr_df['NR_ID_Des'].str.replace(f'{na_fill}::', '', regex=False)
        # ID 列
        nr_id_df = nr_df[nr_df['NR_ID_Des'].str.contains('::')].copy()
        nr_id_df.drop(columns=['NR_Def1', 'NR_ID_Des'], inplace=True)
        result_df = pd.merge(left=result_df, right=nr_id_df, on='GeneID', how='left')
        # ID_Des 列
        nr_df.drop(columns=['NR_ID', 'NR_Def1'], inplace=True)
        result_df = pd.merge(left=result_df, right=nr_df, on='GeneID', how='left')

    if swiss_file:
        swiss_df = load_table(swiss_file, header=0, names=['GeneID', 'Swiss_ID', 'Swiss_Def'], dtype={'GeneID': str})
        swiss_df.fillna(value=na_fill, inplace=True)
        swiss_df['Swiss_ID_Des'] = swiss_df['Swiss_ID'] + '::' + swiss_df['Swiss_Def']
        swiss_df['Swiss_ID_Des'] = swiss_df['Swiss_ID_Des'].str.replace(f'{na_fill}::', '', regex=False)
        # ID 列
        swiss_id_df = swiss_df[swiss_df['Swiss_ID_Des'].str.contains('::')].copy()
        swiss_id_df.drop(columns=['Swiss_Def', 'Swiss_ID_Des'], inplace=True)
        result_df = pd.merge(left=result_df, right=swiss_id_df, on='GeneID', how='left')
        # ID_Des 列
        swiss_df.drop(columns=['Swiss_ID', 'Swiss_Def'], inplace=True)
        result_df = pd.merge(left=result_df, right=swiss_df, on='GeneID', how='left')
        
    if kegg_file:
        kegg_df = load_table(kegg_file, header=0, dtype={'GeneID': str},
                             names=['GeneID', 'KEGG_ID', 'KEGG_Shortname', 'EC_Number', 'KEGG_Description'])
        kegg_df.fillna(value=na_fill, inplace=True)
        result_df = pd.merge(left=result_df, right=kegg_df, on='GeneID', how='left')

    if def_file:
        def_df = load_table(def_file, dtype={index_column: 'str'})
        duplicate_cols = [col for col in def_df.columns if col in result_df.columns and col != index_column]
        if duplicate_cols:
            if dup_col == 'left':
                # 保留左侧（result_df）的列，删除右侧（def_df）的重复列
                def_df = def_df.drop(columns=duplicate_cols)
            elif dup_col == 'right':
                # 保留右侧（def_df）的列，删除左侧（result_df）的重复列
                result_df = result_df.drop(columns=duplicate_cols)
            elif dup_col == 'drop_all':
                # 都不保留，删除所有重复列
                def_df = def_df.drop(columns=duplicate_cols)
                result_df = result_df.drop(columns=duplicate_cols)
            # dup_col == 'keep_both' 时不需要做任何处理，pandas会自动添加后缀
        result_df = pd.merge(left=result_df, right=def_df, on=index_column, how=merge_how, suffixes=('_input', '_def') if dup_col == 'keep_both' else ('', ''))

    diff_col = list(set(result_df.columns) - set(file_df.columns))
    result_df[diff_col] = result_df[diff_col].fillna(value=na_fill)
    result_shape = result_df.shape[0]
    
    if source_shape != result_shape:
        logger.info(f"原表行数{source_shape}, 结果表行数{result_shape}, 可能输入文件指定合并列有重复，或注释文件有重复")
        if remove_duplicates:
            result_df.drop_duplicates(subset=index_column, inplace=True)
            logger.info('已自动去重处理')

    return result_df


def main() -> None:
    args = parse_input()
    for_merge_column = args.specify_column
    
    try:
        args.input_header = int(args.input_header)
    except ValueError:
        pass

    if isinstance(args.input_header, str):
        df = load_table(args.input, dtype={for_merge_column: str})
    elif isinstance(args.input_header, int):
        df = load_table(args.input, header=None)
    else:
        logger.error('输入的 input_header 有误，请检查')
        exit(1)

    if args.input_header != for_merge_column:
        df.rename(columns={args.input_header: for_merge_column}, inplace=True)
    logger.debug(f'{df.columns}')
    logger.debug(f'{df.dtypes}')
    
    result_df = add_def(df, for_merge_column, args.kegg, args.nr, args.swiss, args.defi, args.merge_how, args.remove_duplicates, args.dup_col, args.na_fill)
    
    if args.output.endswith('.xlsx'):
        result_df = result_df.rename(columns={for_merge_column: args.input_header})
        write_output_df(result_df, args.output, index=False)
    if type(args.input_header) == str:
        result_df = result_df.rename(columns={for_merge_column: args.input_header})
        write_output_df(result_df, args.output, index=False)
    elif type(args.input_header) == int:
        write_output_df(result_df, args.output, index=False, header=True)
        
    logger.success('Done!')


if __name__ == '__main__':
    main()
