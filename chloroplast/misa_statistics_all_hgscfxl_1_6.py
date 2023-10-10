# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 15:43:16 2020

@author: fish

使用方法：在叶绿体fasta文件夹下直接运行此脚本  环境：train  conda base
"""
import pandas as pd, os
import argparse


def parse_input():
    argparser = argparse.argparse.ArgumentParser(description='')
    
    
path_db = r'/home/data/db/misa/'
path_dir = os.getcwd()
file_list = os.listdir()
for i in file_list:
    if '.fasta' in i:
        os.system('cp ./'+i+' '+path_db)
    
        os.system('perl /home/data/db/misa/split_genome_misa.pl')


os.system('rm /home/data/db/misa/*.fasta*')



fileline = os.listdir()
t = []#序列名称
t1 = []#单核苷酸重复数目
t2 = []#二核苷酸重复数目
t3 = []#三核苷酸重复数目
t4 = []#四核苷酸重复数目
t5 = []#五核苷酸重复数目
t6 = []#六核苷酸重复数目
for i in fileline:
    
    if '.fasta.statistics' in i:
        print(i)
        t.append(i.strip().split('.fasta.statistics')[0])
        f = open(path_dir +os.sep+ i , 'r+')
        for ii in f :
            if ii.strip().startswith('1'):
                t1.append(ii.strip().split()[1])
            if ii.strip().startswith('2'):
                t2.append(ii.strip().split()[1])
            if ii.strip().startswith('3'):
                t3.append(ii.strip().split()[1])
            if ii.strip().startswith('4'):
                t4.append(ii.strip().split()[1])
            if ii.strip().startswith('5'):
                t5.append(ii.strip().split()[1])
            if ii.strip().startswith('6'):
                t6.append(ii.strip().split()[1])
        f.close()
        if len(t) == len(t1) :
            a = 1
        else:
            t1.append('0')
        if len(t) == len(t2) :
            a = 1
        else:
            t2.append('0')
        if len(t) == len(t3) :
            a = 1
        else:
            t3.append('0')
        if len(t) == len(t4) :
            a = 1
        else:
            t4.append('0')
        if len(t) == len(t5) :
            a = 1
        else:
            t5.append('0')
        if len(t) == len(t6) :
            a = 1
        else:
            t6.append('0')


df = pd.DataFrame()
df['name'] = t
df['Mononucleotide repeats'] = t1
df['Dinucleotide repeats'] = t2
df['Trinucleotide repeats'] = t3
df['Tetranucleotide repeats'] = t4
df['Pentanucleotide repeats'] = t5
df['Hexanucleotide repeats'] = t6
df.to_csv('./misa_all_rosa_35_1_species.csv' , index = 1 , header = 1)