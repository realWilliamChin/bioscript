import argparse
import os
from shutil import copyfile
import sys
from matplotlib import pyplot as plt

parser=argparse.ArgumentParser(description="dumx_to_join")
parser.add_argument('--inputPath',dest='inputPath',default='./04_quality_parement/length_stat.csv',
        help='path of input file')

parser.add_argument('--outputPath',dest='outputPath',default='./04_quality_parement/',
        help='path of output file')

parser.add_argument('--xlabel',dest='xlabel',default='Length of Reads',
        help='name of xlabel')

parser.add_argument('--ylabel',dest='ylabel',default='Number of Reads',
        help='name of ylabel')

parser.add_argument('--title',dest='title',default='Sequence Length Distribution Statistics Histogram',
        help='name of title')


arg=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputfilePath=arg.inputPath
outputfilePath=arg.outputPath
x_label=arg.xlabel
y_label=arg.ylabel
title=arg.title
if not os.path.exists(outputfilePath):
    os.makedirs(outputfilePath)

f1=open(inputfilePath,'r+')
DonotFirst=False
x_list=[]
y_list=[]
for line in f1:
    if DonotFirst==False:
        DonotFirst=True
        continue
    else:
        L=line.split('\n')[0].split('\t')
        x_list.append(L[0])
        y_list.append(int(L[0]))

figsize = len(x_list)//3,9
figure, ax = plt.subplots(figsize=figsize)

def autolabel(rects):

    for rect in rects:

        height = rect.get_height()

        #plt.text(rect.get_x()+rect.get_width(), 1.03*height,'%s' % int(height))
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.xticks(rotation=90)

autolabel(plt.bar(range(len(y_list)), y_list, tick_label=x_list))

plt.savefig(outputfilePath+os.sep+'length_stat.jpg')

