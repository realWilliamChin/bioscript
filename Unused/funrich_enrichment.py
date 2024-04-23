#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created Time  : 2023/10/17 17:19
# Author        : William GoGo
import os
import argparse
import pandas as pd
import math
from scipy.stats import hypergeom
from decimal import Decimal, getcontext


def parse_input():
    args = argparse.ArgumentParser(description='Funrich enrichment result')
    # 创建背景数据库
    args.add_argument('--all-geneid', type=str, dest='all_geneid',
                        help='输入 all-geneid 文件')
    args.add_argument('--go-bp', type=str, dest='go_bp',
                        help='输入 swiss_GO_BP_ID.txt 文件')
    args.add_argument('--go-cc', type=str, dest='go_cc',
                        help='输入 swiss_GO_CC_ID.txt 文件')
    args.add_argument('--go-mf', type=str, dest='go_mf',
                        help='输入 swiss_GO_MF_ID.txt 文件')
    args.add_argument('--kegg', type=str, dest='kegg',
                        help='输入 swiss_KEGG_ID.txt 文件')
    
    # 输入待分析的所有基因文件的目录，或者单个文件
    args.add_argument('-i', '--geneid', type=str, dest='geneid',
                        help='输入待分析的基因文件')
    args.add_argument('-d', '--dir', type=str, dest='dir',
                        help='输入待分析的所有基因文件的目录')
    
    # 判断所有的输入的文件是否存在
    args = args.parse_args()
    if args.all_geneid and not os.path.exists(args.all_geneid):
        print('输入的 all-geneid 文件不存在，请检查！')
        exit()
    elif args.go_bp and not os.path.exists(args.go_bp):
        print('输入的 GO_BP 文件不存在，请检查！')
        exit()
    elif args.go_cc and not os.path.exists(args.go_cc):
        print('输入的 GO_CC 文件不存在，请检查！')
        exit()
    elif args.go_mf and not os.path.exists(args.go_mf):
        print('输入的 GO_MF 文件不存在，请检查！')
        exit()
    elif args.kegg and not os.path.exists(args.kegg):
        print('输入的 KEGG 文件不存在，请检查！')
        exit()
    elif args.geneid and not os.path.exists(args.geneid):
        print('输入的 geneid 文件不存在，请检查！')
        exit()
    elif args.dir and not os.path.exists(args.dir):
        print('输入的 dir 目录不存在，请检查！')
        exit()
    
    return args


def fold_enrichment(N, K, n, k):
    Fold_Correction_Factor = 0.01
    if n == 0:
        result2 = 1.0
    elif k == 0:
        result2 = 0.01
    else:
        Ratio_In_Sample = (k + Fold_Correction_Factor) / (n + Fold_Correction_Factor)
        expectedratio = (K + Fold_Correction_Factor) / (N + Fold_Correction_Factor)
        result = Ratio_In_Sample / expectedratio
        if result < 0.0:
            result2 = 0.0
        else:
            result2 = result
    return result2


def log_factorial(n):
    if n <= 1:
        return 0
    else:
        return sum(math.log(i) for i in range(1, n+1))

def Calculate_PValueVAN_Fix(N, K, n, k):
    if n > N:
        n = N
    if k > K:
        k = K
    if N > 89999:
        return 1.0
    resultArray = [0.0] * (n+1)
    result = 0.0
    if k == 0:
        return 1.0
    elif k > K:
        return 0.0
    elif n + K - N <= 0:
        p = log_factorial(N - K) + log_factorial(N - n) - log_factorial(N) - log_factorial(N - n - K)
        resultArray[0] = p
        for i in range(1, min(n, K) + 1):
            resultArray[i] = resultArray[i - 1] + math.log(n - i + 1) + math.log(K - i + 1) - math.log(N - n - K + i) - math.log(i)
            if i >= k:
                result += math.exp(resultArray[i])
        return result
    else:
        k_off = n + K - N
        p = log_factorial(K) + log_factorial(n) - log_factorial(k_off) - log_factorial(N)
        for i in range(1, n + 1):
            if i < k_off:
                resultArray[i] = float('-inf')
            elif i == k_off:
                resultArray[i] = p
            elif i < K:
                resultArray[i] = resultArray[i - 1] + math.log(n - i + 1) + math.log(K - i + 1) - math.log(N - n - K + i) - math.log(i)
            else:
                resultArray[i] = float('-inf')
            if i >= k:
                result += math.exp(resultArray[i])
        return result


def bonferroni_sort_count(file_lst):
    """
    对 Funrich 的 report.txt GO 和 KEGG 的结果进行统计，统计 Bonferroni method 小于 0.05 的个数
    """
    result_lst = []
    for report_file in file_lst:
        for pathway in ['GO_BP', 'GO_CC', 'GO_MF', 'KEGG']:
            df = pd.read_excel(report_file, sheet_name=pathway, skiprows=7, header=1, usecols=[0, 6]).dropna(subset=pathway)
            file_name = report_file.split('_repo')[0] + '_' + pathway
            all_mapped = df.shape[0]
            lt005_count = df[df['Bonferroni method'] < 0.05].shape[0]
            result_lst.append([file_name, all_mapped, lt005_count])
    result_df = pd.DataFrame(result_lst, columns=['file_name', 'all_mapped', 'lt005_count'])
    result_df.to_csv('summary.txt', sep='\t', index=False, header=False)


def bonferroni_sort(funrich_file):
    file_xlsx = pd.ExcelFile(funrich_file, engine='openpyxl')
    sheet_names = file_xlsx.sheet_names
    with pd.ExcelWriter("result/" + funrich_file, engine='xlsxwriter') as writer:
        for each_sheet in sheet_names:
            df = pd.read_excel(file_xlsx, each_sheet, skiprows=8, engine='openpyxl')
            
            df.to_excel(excel_writer=writer, sheet_name=each_sheet, index=False) #
    


if __name__ == '__main__':
    b = fold_enrichment(6413, 131, 154, 2)
    c = Calculate_PValueVAN_Fix(6413, 131, 154, 2)
    print(b, c)