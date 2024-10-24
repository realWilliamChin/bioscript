# -*- coding: utf-8 -*-
# @FileName : gsea.py
# @DATE ：16:14 2022/4/18
# @Author : CS:GO, WilliamGOGO

import os
import argparse

def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('--gobp', type=str, required=True, help='GO_BP')
    p.add_argument('--gocc', type=str, required=True, help='GO_CC')
    p.add_argument('--gomf', type=str, required=True, help='GO_MF')
    p.add_argument('--kegg-tier2', dest='kegg_tier2', type=str, required=True, help='KEGG_tier2.txt')
    p.add_argument('--kegg-tier3', dest='kegg_tier3', type=str, required=True, help='KEGG.txt')
    p.add_argument('--output-dir', type=str, default='./GSEA', help='output dir')
    
    args = p.parse_args()
    
    return args


def gsea(gobp, gocc, gomf, kegg_tier2, kegg_tier3, output_path):
    dir_list=[gobp, gocc, gomf, kegg_tier2, kegg_tier3]
    dic1={}

    for each in dir_list:
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        key1=each.replace('.txt','').replace('_clean','')
        if 'BP' in key1:
            key='GO_BP'
        if 'CC' in key1:
            key = 'GO_CC'
        if 'MF' in key1:
            key = 'GO_MF'
        if 'KEGG' in key1:
            key = 'KEGG'
        
        if each == dir_list[-1]:
            tier3_name = kegg_tier2.replace('.txt', '').replace('_clean', '').replace('tier2', 'tier3')
            f2=open(os.path.join(output_path, os.path.basename(tier3_name)+'.GMT'),'w')
        else:
            f2=open(os.path.join(output_path, os.path.basename(key1)+'.GMT'),'w')
        dic1=key
        globals()[dic1] = {}
        dic={}

        with open(each,'r') as f1:
            num=0
            for each_line in f1:
                if num==0:
                    num+=1
                    continue
                else:
                    if 'KEGG' in  key1:
                        id=each_line.split('\t')[0].strip()
                        go_id=each_line.split('\t')[1].split(':')[0].strip()
                        content=each_line.split(':')[-1].strip()
                        if go_id in globals()[dic1]:
                            globals()[dic1][go_id] += '\t'+id
                        else:
                            globals()[dic1][go_id] = content+';;'+key+'\t'+id


                    else:
                        id=each_line.split('\t')[0].strip()
                        go_id=each_line.split('\t')[1].split('_')[0].strip()
                        content=each_line.split('_')[-1].strip()
                        if go_id in globals()[dic1]:
                            globals()[dic1][go_id] += '\t'+id
                        else:
                            globals()[dic1][go_id] = content+';;'+key+'\t'+id
        for each2 in globals()[dic1]:
            output=globals()[dic1][each2].split(';;')[0]+' '+'['+each2+']'+'\t'+globals()[dic1][each2].split(';;')[-1]+'\n'
            f2.write(output)


def main():
    args = parse_input()
    gsea(args.gobp, args.gocc, args.gomf, args.kegg_tier2, args.kegg_tier3, args.output_dir)


if __name__ == '__main__':
    main()





















