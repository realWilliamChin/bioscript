#!/usr/bin/env/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/01/18 10:40
# Author        : William GoGo
import json
import argparse
import pandas as pd
from bs4 import BeautifulSoup


def parse_input():
    args = argparse.ArgumentParser(description='json2tsv.py')
    args.add_argument('-k', '--ko_number', help='输入ko_number，例如 ko03000')
    args.add_argument('-o', '--output', help='输出文件')
    args = args.parse_args()

    return args


def parse_json_to_table(data, path=[]):
    table = []
    current_path = path + [data['values']]

    if not data['children']:  # 如果没有子节点
        table.append(current_path)
    else:
        for child in data['children']:
            table.extend(parse_json_to_table(child, current_path))

    return table

def fill_na(table):
    max_length = max(len(row) for row in table)
    return [row + ['NA'] * (max_length - len(row)) for row in table]


def remove_json_field(json_dict, field):
    if isinstance(json_dict, dict):
        # if field in json_dict:
        #     print(f"Removing field '{field}' from {json_dict}")
        json_dict.pop(field, None)
        for key, value in list(json_dict.items()):  # 使用list来避免在迭代时修改字典
            remove_json_field(value, field)
    elif isinstance(json_dict, list):
        for item in json_dict:
            remove_json_field(item, field)
    
    return json_dict


def get_output_file():
        

def main():
    args = parse_input()

    json_file = f'{args.ko_number}.json'
    
    with open(json_file, 'r') as f:
        json_dict = json.load(f)
        
        root_data = json_dict['root']
        root_data = remove_json_field(root_data, 'expanded')
        root_data = remove_json_field(root_data, 'isRoot')

        table = parse_json_to_table(root_data)
        table = fill_na(table)
        df = pd.DataFrame(table)
        
        df = df.apply(lambda x: list(x.str)[0] if x.str.contains('[') else x)
        
        df.to_csv('output.txt', sep='\t', index=False, header=False)
        

if __name__ == '__main__':
    main()
                
                
        
        