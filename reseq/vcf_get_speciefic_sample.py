#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/30 23:46
# Author        : William GoGo
import os
import argparse
import pandas as pd
from vcf_common import skip_rows


def parse_input():
    argparser = argparse.ArgumentParser(description='Get specific sample from vcf file')
    argparser.add_argument('-v', '--vcf', help='input vcf file')
    argparser.add_argument('-l', '--sample_list', help='sample list, 格式为一个文件一行一个样本名，或者直接输入,分隔的样本名')
    argparser.add_argument('-o', '--output', default='specific_sample.vcf', help='output vcf file')
    # argparser.add_argument('--save_type')
    return argparser.parse_args()


# def get_specific_sample(vcf_file, sample_lst, output_file='specific_sample.vcf'):
    with open(vcf_file, 'r') as f:
        for line in f:
            # 写入注释行
            if line.startswith('##'):
                with open(output_file, 'a') as out_f:
                    out_f.write(line)
            # 处理标题行
            elif line.startswith('#CHROM'):
                column_name_lst = line.strip().split('\t')
                column_fixed_lst = column_name_lst[:9]
                specific_index_lst = [column_name_lst.index(x) for x in sample_lst]
                with open(output_file, 'a') as out_f:
                    out_f.write('\t'.join(column_fixed_lst + sample_lst) + '\n')
                with open('error.txt', 'a') as error_f:
                    error_f.write(line + '\n')
            # 处理数据行
            else:
                element_lst = line.strip().split('\t')
                # 分为固定列（前9列），特定样本列（sample_lst），非特定样本列（除了前9列和特定样本列的其他列）
                fixed_element_lst = element_lst[:9]

                # 获取特定样本的元素
                specific_element_lst = []
                not_specific_element_lst = []
                for elem_index, element in enumerate(element_lst):
                    if elem_index in specific_index_lst:
                        specific_element_lst.append(element)
                    elif elem_index > 8:
                        not_specific_element_lst.append(element)

                # 计算特定样本的有效位点数，是 0 则会跳过
                specific_count_valid = len(specific_element_lst) - specific_element_lst.count('./././.')
                # 计算除了特定样本外，其他样本的有效位点数，不是 0 则会跳过
                not_specific_count_valid = len(not_specific_element_lst) - not_specific_element_lst.count('./././.')

                if not_specific_count_valid == 0 and specific_count_valid >= 1:
                    if len(not_specific_element_lst) != 85:
                        with open('error.txt', 'a') as error_f:
                            error_f.write('\t'.join(element_lst) + '\n')
                    with open(output_file, 'a') as out_f:
                        out_f.write('\t'.join(fixed_element_lst + specific_element_lst) + '\n')


def get_specific_sample(vcf_file, sample_lst, output_file, vcf_type='gatk'):
    skiprows = skip_rows(vcf_file, '##')
    if vcf_type == 'vcftools':
        vcf_df = pd.read_csv(vcf_file, sep='\t', skiprows=skiprows, low_memory=False)
        vcf_df['specific_valid_count'] = (vcf_df[sample_lst] != './././.').astype(int).sum(axis=1)
        
        # for not_specific_valid_count
        vcf_df['all_valid_count'] = (vcf_df != './././.').astype(int).sum(axis=1)
        vcf_df['not_specific_valid_count'] = vcf_df['all_valid_count'] - 10 - vcf_df['specific_valid_count']
        result_df = vcf_df[(vcf_df['not_specific_valid_count'] == 0) & (vcf_df['specific_valid_count'] >= 1)]

        result_df = result_df[result_df.columns.values.tolist()[:9] + sample_lst]
        print(f"{output_file}——rows, columns:{result_df.shape}")
        result_df.to_csv(output_file, sep='\t', index=False)
    elif vcf_type == 'gatk':
        vcf_df = pd.read_csv(vcf_file, sep='\t', skiprows=skiprows, low_memory=False)
        df_sample_lst = vcf_df.columns.values.tolist()[9:]
        for sample in df_sample_lst:
            vcf_df[sample + '_count'] = vcf_df[sample].str.split(':').str[0]
        pass
        


# TODO: 合并 vcf_sample_reorder.py 和 vcf_get_speciefic_sample.py
def reorder(vcf, sample_lst, output_file):
    # read headers startwith '#' and not startwith '##' as a list
    with open(vcf, 'r') as f, open(output_file, 'a') as out_f:
        column_name_lst = []
        for line in f:
            if line[:2] == '##':
                out_f.write(line)
            elif line[:1] == '#':
                column_name_lst = line.strip().split('\t')
                break
    vcf_df = pd.read_csv(vcf, sep='\t', comment='#', header=None, names=column_name_lst, low_memory=False)
    vcf_sample_lst = vcf_df.columns.tolist()[:9] + sample_lst
    result_df = vcf_df[vcf_sample_lst]
    result_df.to_csv(output_file, sep='\t', index=False, mode='a')


def main():
    args = parse_input()
    print(f'reading file {args.vcf}')
    # 检测 sample_list 是文件还是直接输入的
    if os.path.isfile(args.sample_list):
        with open(args.sample_list, 'r') as f:
            sample_lst = f.read().strip().split('\n')
    else:
        sample_lst = args.sample_list.split(',')
    print(f'getting sample {sample_lst}')
    get_specific_sample(args.vcf, sample_lst, args.output)


if '__main__' == __name__:
    main()
