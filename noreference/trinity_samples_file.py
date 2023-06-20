#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/05/20 21:41
# Author        : William GoGo
import os

paired = []
# 生成 samples trinity file
for dir in os.listdir('02_Pinjiedata'):
    if os.path.isdir('02_Pinjiedata' + os.sep + dir):
        for file in os.listdir('02_Pinjiedata' + os.sep + dir):
            if '_1.clean.fq' in file:
                # 获取两个文件的绝对路径
                file1_path = os.path.abspath('02_Pinjiedata' + os.sep + dir + os.sep + file)
            elif '_2.clean.fq' in file:
                file2_path = os.path.abspath('02_Pinjiedata' + os.sep + dir + os.sep + file)
        # 生成 samples trinity file
        with open('samples_trinity.txt', 'a') as f:
            f.write(dir + '\t' + dir + '_rep1' + '\t' + file1_path + '\t' + file2_path + '\n')
    else:
        if '_1.clean.fq' in dir:
            # 获取两个文件的绝对路径
            file1_path = os.path.abspath('02_Pinjiedata' + os.sep + dir)
            paired.append(file1_path)
        elif '_2.clean.fq' in dir:
            file2_path = os.path.abspath('02_Pinjiedata' + os.sep + dir)
            paired.append(file2_path)
        if len(paired) == 2:
            # 生成 samples trinity file
            with open('samples_trinity.txt', 'a') as f:
                f.write(dir + '\t' + dir + '_rep1' + '\t' + paired[0] + '\t' + paired[1] + '\n')
            paired = []
