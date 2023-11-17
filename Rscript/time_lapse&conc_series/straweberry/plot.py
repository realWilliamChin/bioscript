# -*- coding: UTF-8 -*-
import argparse
import os
from shutil import copyfile
import sys
import numpy
cur_path=os.getcwd()
doc=os.listdir(cur_path)

'''
在当前文件夹中提出特定_ID.txt文件里的序列，
需要--input 文件  文件中包含.fasta文件 _ID.txt文件


'''
import pylab




import matplotlib.pyplot as plt
ls2=[]
with open('group1.txt','r') as f1:
    num1=0
    for each_line1 in f1:
        if num1==0:
            num1+=1
        else:
            content=each_line1.strip().split('\t')[1:]
            ls2.append(content)
for each in ls2:

    plt.plot(each)

plt.show()















