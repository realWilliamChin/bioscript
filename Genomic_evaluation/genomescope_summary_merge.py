#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/18 15:34
# Author        : William GoGo
import os, sys
import argparse
import pandas as pd

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(description='Merge GenomeScope summary files into one xlsx file.')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input GenomeScope summary files (space separated)')
    parser.add_argument('--genome-scope-version', type=int, help='genome scope version')
    parser.add_argument('-o', '--output', required=True, help='Output xlsx file')
    args = parser.parse_args()
    return args


def read_summary_file(filepath, genome_scope_version):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Common properties for both versions
    common_properties = ['Genome Haploid Length', 'Genome Repeat Length', 'Genome Unique Length', 'Model Fit', 'Read Error Rate']
    
    # Extract common data
    def extract_range_data(property_name):
        line = [line.replace(property_name, '').replace('bp', '') for line in lines if line.startswith(property_name)][0]
        values = line.split()
        min_val = values[0].replace(',', '') + (' bp' if 'Length' in property_name else '')
        max_val = values[1].replace(',', '') + (' bp' if 'Length' in property_name else '')
        return min_val, max_val
    
    # Extract common ranges
    common_min_data = []
    common_max_data = []
    for prop in common_properties:
        min_val, max_val = extract_range_data(prop)
        common_min_data.append(min_val)
        common_max_data.append(max_val)
    
    if genome_scope_version == 1:
        # Version 1 specific: Heterozygosity
        hetero_line = [line.replace('Heterozygosity', '').replace('bp', '') for line in lines if line.startswith('Heterozy')][0]
        hetero_values = hetero_line.split()
        
        property = ['Heterozygosity'] + common_properties
        min_data = [hetero_values[0]] + common_min_data
        max_data = [hetero_values[1]] + common_max_data
        
    elif genome_scope_version == 2:
        # Version 2 specific: Homozygous and Heterozygous
        homo_line = [line.replace('Homozygous (aa)', '') for line in lines if line.startswith('Homozygous')][0]
        hetero_line = [line.replace('Heterozygous (ab)', '') for line in lines if line.startswith('Heterozygous')][0]
        
        homo_values = homo_line.split()
        hetero_values = hetero_line.split()
        
        property = ['Homozygous (aa)', 'Heterozygous (ab)'] + common_properties
        min_data = [homo_values[0], hetero_values[0]] + common_min_data
        max_data = [homo_values[1], hetero_values[1]] + common_max_data
    
    df = pd.DataFrame({
        'property': property,
        'min': min_data,
        'max': max_data
    })
    return df
            

def merge_summaries(file_list, genome_scope_version):
    merged_df = None
    for file in file_list:
        df = read_summary_file(file, genome_scope_version)
        df.rename(columns={
            'min': f"{file}_min",
            'max': f"{file}_max"
        }, inplace=True)
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='property', how='outer')
    return merged_df


def write_to_excel(df, output_file):
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False)
        # # 先写表格
        # df.to_excel(writer, index=False, startrow=3)
        # # 写前三行注释
        # worksheet = writer.sheets['Sheet1']
        # for i, line in enumerate(headers):
        #     worksheet.write(i, 0, line.strip())


def main():
    args = parse_input()
    input_files = args.input
    output_file = args.output
    merged_df = merge_summaries(input_files, args.genome_scope_version)
    write_to_excel(merged_df, output_file)


if __name__ == '__main__':
    main()
    
