#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/01/21 16:52
# Author        : William GoGo
import re
import sys
from time import perf_counter

def build_pattern(search_strings):
    """构建包含所有搜索字符串的正则表达式模式"""
    escaped_strings = (re.escape(s.strip()) for s in search_strings if s.strip())
    pattern_str = r'|'.join(escaped_strings)
    return re.compile(pattern_str)

def filter_large_file(input_file, output_file, pattern):
    """从大文件中过滤匹配模式的行，保持原始顺序"""
    with open(output_file, 'w', encoding='utf-8') as out_f:
        with open(input_file, 'r', encoding='utf-8', errors='replace') as in_f:
            for line in in_f:
                if pattern.search(line):  # 使用search()而不是match()进行部分匹配
                    out_f.write(line)

def main():
    if len(sys.argv) != 4:
        print("Usage: python fast_filter.py <list_file> <large_file> <output_file>")
        sys.exit(1)

    list_file, large_file, output_file = sys.argv[1:4]

    # 步骤1：读取搜索字符串并构建高效的正则表达式
    print("Building search pattern...")
    start_time = perf_counter()
    
    with open(list_file, 'r', encoding='utf-8', errors='replace') as f:
        pattern = build_pattern(f)
    
    compile_time = perf_counter() - start_time
    print(f"Pattern compiled in {compile_time:.4f} seconds")
    print(f"Using regex pattern with {pattern.pattern.count('|') + 1} search terms")

    # 步骤2：流式处理大文件
    print(f"Processing {large_file}...")
    start_time = perf_counter()
    
    filter_large_file(large_file, output_file, pattern)
    
    process_time = perf_counter() - start_time
    print(f"Completed in {process_time:.4f} seconds")
    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    main()