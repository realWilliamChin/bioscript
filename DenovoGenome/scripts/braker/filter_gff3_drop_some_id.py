#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/01/21 16:52
# Author        : William GoGo
import sys
from time import perf_counter

def should_keep_line(line, search_set, search_strings):
    """判断是否应该保留这一行"""
    if line.startswith('#') or not line.strip():
        return True
    
    parts = line.strip().split('\t')
    if len(parts) < 9:
        return True
    
    feature_type = parts[2]
    attributes = parts[8]

    # 提取 ID
    id_idx = attributes.find("ID=")
    if id_idx == -1:
        return True
    # ID=xxx; 取出xxx
    end_idx = attributes.find(";", id_idx)
    line_id = attributes[id_idx+3 : end_idx if end_idx != -1 else None]

    if feature_type == "gene":
        # 基因 ID 去掉.后面的部分去匹配，避免g1匹配到g11
        gene_id = line_id.split('.')[0] if '.' in line_id else line_id
        # 精确匹配，避免g1匹配到g11
        for s in search_strings:
            s_gene_id = s.split('.')[0] if '.' in s else s
            if gene_id == s_gene_id:
                return True
        return False
    else:
        # 其他特征允许.t后缀匹配
        for s in search_strings:
            if s in line_id:
                return True
        return False

def filter_large_file(input_file, output_file, search_set, search_strings):
    with open(output_file, 'w', encoding='utf-8') as out_f:
        with open(input_file, 'r', encoding='utf-8', errors='replace') as in_f:
            for line in in_f:
                if should_keep_line(line, search_set, search_strings):
                    out_f.write(line)

def main():
    if len(sys.argv) != 4:
        print("Usage: python fast_filter.py <list_file> <large_file> <output_file>")
        sys.exit(1)

    list_file, large_file, output_file = sys.argv[1:4]

    print("Loading search IDs...")
    start_time = perf_counter()
    with open(list_file, 'r', encoding='utf-8', errors='replace') as f:
        search_strings = [line.strip() for line in f if line.strip()]
        search_set = set(search_strings)   # 用集合加速 gene 匹配
    load_time = perf_counter() - start_time
    print(f"Loaded {len(search_strings)} IDs in {load_time:.4f} seconds")

    print(f"Processing {large_file}...")
    start_time = perf_counter()
    filter_large_file(large_file, output_file, search_set, search_strings)
    process_time = perf_counter() - start_time

    print(f"Completed in {process_time:.4f} seconds")
    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    main()
