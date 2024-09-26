import argparse
import os
from shutil import copyfile
import sys
import numpy as np
import matplotlib.pyplot as plt
import pygal
import time
import cairosvg

parser=argparse.ArgumentParser(description="stat table and output stacked histogram from tax in different tax")
parser.add_argument('--inputPath',dest='inputPath',default='./07_Taxonomy/taxa_bar',
        help='path of input file')

parser.add_argument('--outputPath',dest='outputPath',default='./07_Taxonomy/taxa_bar_output',
        help='path of output file')

parser.add_argument('--xlabel',dest='xlabel',default='Groups',
        help='name of xlabel')

parser.add_argument('--ylabel',dest='ylabel',default='Tax',
        help='name of ylabel')

parser.add_argument('--title',dest='title',default='Statistics Stacked Histogram',
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

class_list=['group','Domain','Phylum','Class','Order','Family','Genus']
filelist=os.listdir(inputfilePath)
for file in filelist:
    print(file)
    if '.csv' not in file:
        continue
    if file.split('.csv')[0][-1] not in '123456':
        continue
    filename=class_list[int(file.split('.csv')[0][-1])]
    f1=open(inputfilePath+os.sep+file,'r+')
    f2=open(outputfilePath+os.sep+filename+'.csv','w+')
    rawlist=[]
    for line in f1:
        rawlist.append(line.split('\n')[0].split(','))
    rawlist_T=list(map(list,zip(*rawlist)))

    garph = pygal.StackedBar()  # 创建图（叠加柱状图）
    garph.x_labels = rawlist_T[0][1:]  # 设置 X 轴刻度
    for line in rawlist_T[1:]:
        # 添加数据

        if 'class,' in ','.join(line):
            continue
        print(line)
        print(1)
        garph.add(line[0], [int(i.split('.')[0]) for i in line[1:]])
        print(2)
        print(file)
    for line in rawlist_T[:-1]:
        f2.write(','.join(line)+'\n')

    f1.close()
    f2.close()
    # garph.y_label_rotation = 45           # 设置Y轴的标签旋转多少度
    # print(3)
    # garph.x_label_rotation = 90
    # print(4)
    # garph.title = filename+' '+title  # 设置图标题
    # garph.x_title = x_label  # 设置 X 轴标题
    # garph.y_title = y_label  # 设置 Y 轴标题
    # print(5)
    # garph.legend_at_bottom = True  # 设置图例位置（下面）
    #
    # garph.margin = 35  # 设置页边距（margin_bottom、margin_top、margin_left、margin_right）
    # print(6)
    # garph.show_x_guides = False  # 显示X轴的网格线
    # garph.show_y_guides = False  # 显示Y轴的网格线
    # print(7)
    # garph.render_to_file(outputfilePath + os.sep + filename + '-level Stacked Histogram.svg')
    # print(8)
    # svg_path=outputfilePath + os.sep + filename + '-level Stacked Histogram.svg'
    # print(9)
    # png_path=outputfilePath + os.sep + filename + '-level Stacked Histogram.png'
    # print(file)
    # print('finish')
    # cairosvg.svg2png(url=svg_path, write_to=png_path)

