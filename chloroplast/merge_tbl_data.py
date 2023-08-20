#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/14 23:04
# Author        : William GoGo
import os
import argparse
import pandas as pd


def parse_input():
    argparser = argparse.ArgumentParser(description="")
    argparser.add_argument('-i', '--input_dir', default='./', help="input dir")
    argparser.add_argument('-o', '--output', default='all_tbl.txt', help="output file name")
    return argparser.parse_args()
    

def merge_all_tbl(input_dir, output_name):
    df = pd.DataFrame({
        "type": ["Retroelements", "SINEs", "Penelope", "LINEs", "CRE/SLACS", "L2/CR1/Rex",
                "R1/LOA/Jockey", "R2/R4/NeSL", "RTE/Bov-B", "L1/CIN4", "LTR elementss",
                "BEL/Pao", "Ty1/Copia", "Gypsy/DIRS1", "Retroviral", "DNA transposons",
                "hobo-Activator", "Tc1-IS630-Pogo", "En-Spm", "MuDR-IS905", "PiggyBac",
                "Tourist/Harbinger", "Other (Mirage, P-element, Transib)", "Rolling-circles",
                "Unclassified", "Total interspersed repeats", "Small RNA", "Satellites",
                "Simple repeats", "Low complexity"]
    })

    for each_file in os.listdir(input_dir):
        if '.tbl' not in each_file:
            continue
        specie_name = each_file.split('_')[0]
        with open(each_file, 'r') as f:
            # # 跳过前面的 9 行
            for _ in range(10):
                f.readline()
            # 读取中间的 38 行
            lines = []
            for _ in range(38):
                lines.append(f.readline())
            
            num_of_elements = []
            percentage_of_seq = []
            
            # 选取需要的信息
            for line in lines:
                if line.strip() == '':
                    continue
                elif line.strip().startswith('P-element'):
                    continue
                elif line.strip().startswith('=='):
                    break
                else:
                    element = line[23:26].strip()
                    percent = line.split(' ')[-2] + '%'
                    if element == 'ats':
                        element = ''
                    num_of_elements.append(element)
                    percentage_of_seq.append(percent)
                    
            df[specie_name + '_number'] = pd.Series(num_of_elements)
            df[specie_name + '_percent'] = pd.Series(percentage_of_seq)

    df.to_csv(output_name, sep='\t', index=False)
    

def main():
    args = parse_input()
    merge_all_tbl(args.input_dir, args.output)


if __name__ == '__main__':
    main()
