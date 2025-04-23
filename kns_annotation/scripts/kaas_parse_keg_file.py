#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/04/18 15:42
# Author        : William GoGo
# 解析kaas注释结果文件 *.keg文件
import os, sys
import argparse
import re
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools'))
from load_input import load_table, write_output_df

def parse_input():
    parser = argparse.ArgumentParser(description='解析KEGG注释文件并生成格式化输出')
    parser.add_argument('-i', '--input', type=str, required=True, help='输入KEGG文件路径')
    parser.add_argument('-o', '--output', type=str, required=True, help='输出文件路径')
    return parser.parse_args()


def parse_kegg_file(input_file, output_file):
    """
    解析KEGG注释文件并生成格式化输出
    
    Args:
        input_file (str): 输入KEGG文件路径
        output_file (str): 输出文件路径
    """
    current_pathway = ""
    current_category = ""
    current_metabolism = ""
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write("KEGG_Pathway_ID\tSubcategory\tCategory\tKO_ID\tGene_Symbol\tEnzyme_Description\n")
        
        for line in infile:
            line = line.strip()
            if not line:
                continue
                
            # 匹配A级分类（代谢类别）
            a_match = re.match(r'^A(\d+)\s+(.+)$', line)
            if a_match:
                current_metabolism = f"A{a_match.group(1)}:{a_match.group(2)}"
                continue
                
            # 匹配B级分类（类别）
            b_match = re.match(r'^B\s+(\d+)\s+(.+)$', line)
            if b_match:
                current_category = f"{b_match.group(1)}:{b_match.group(2)}"
                continue
                
            # 匹配C级分类（通路）
            c_match = re.match(r'^C\s+\d+\s+(.+)\s+\[PATH:([a-zA-Z]+)(\d+)\]', line)
            if c_match:
                current_pathway = f"ko{c_match.group(3)}:{c_match.group(1)}"
                continue
                
            # 匹配D级分类（基因）
            d_match = re.match(r'^D\s+(\S+)\s+(.+?)\s+(K\d+)\s+(.+?)(?:\s*\[EC:|$)', line)
            if d_match:
                gene_id = d_match.group(1)
                gene_symbol = d_match.group(4).split(';')[0].strip()
                ko_id = d_match.group(3)
                enzyme_desc = d_match.group(4).strip()
                
                outfile.write(f"{current_pathway}\t{current_category}\t{current_metabolism}\t{ko_id}\t{gene_symbol}\t{enzyme_desc}\n")
    
    # 整个文件去重，解析完会有重复的，因为不需要那个 GeneID 了
    df = load_table(output_file)
    df = df.drop_duplicates()
    write_output_df(df, output_file, index=False)


def main():
    args = parse_input()
    try:
        parse_kegg_file(args.input, args.output)
        print(f"成功解析KEGG文件并生成输出: {args.output}")
    except Exception as e:
        print(f"处理文件时出错: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main() 