import argparse
import os
from shutil import copyfile
import sys
import numpy as np
import matplotlib.pyplot as plt
import pygal

import cairosvg

parser=argparse.ArgumentParser(description="dumx_to_join")
parser.add_argument('--inputPath',dest='inputPath',default='./06_convertResult/table_ASV_w_tax.txt',
        help='path of input file')

parser.add_argument('--outputPath',dest='outputPath',default='./06_convertResult/',
        help='path of output file')

parser.add_argument('--xlabel',dest='xlabel',default='Groups',
        help='name of xlabel')

parser.add_argument('--ylabel',dest='ylabel',default='ASV in Tax',
        help='name of ylabel')

parser.add_argument('--title',dest='title',default='Species Taxonomy Annotation Result Statistics Stacked Histogram',
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
List=[]
DonotSecond=False
for line in f1:
    if DonotFirst==False:
        DonotFirst=True
        continue
    elif DonotSecond==False:
        DonotSecond=True
        List.append(line.split('\t')[:-1])
        continue
    else:
        L=line.split('\n')[0].split('\t')
        tax=L[-1]
        taxline=tax.split(' ')
        tax=taxline[-1][0]
        if len(taxline)==1 and 'unc' in taxline[-1].lower():
            tax='unclassified'
        newline=[L[0]]
        for i in L[1:-1]:
            if float(i)>0:
                newline.append(tax)
            else:
                newline.append(i)
        List.append(newline)

List=list(map(list,zip(*List)))

x_list=[]
y_list_domain=[]
y_list_phylum=[]
y_list_class=[]
y_list_order=[]
y_list_family=[]
y_list_genus=[]
y_list_species=[]
y_list_unclassified=[]

for line in List[1:]:
    domain = 0
    phylum = 0
    Class = 0
    order = 0
    family = 0
    genus = 0
    species = 0
    unclassified = 0

    x_list.append(line[0])
    #print(x_list)
    for item in line[1:]:
        if 'd' in item:
            domain+=1
        elif 'p' in item:
            phylum+=1
        elif 'c' in item:
            Class+=1
        elif 'o' in item:
            order+=1
        elif 'f' in item:
            family+=1
        elif 'g' in item:
            genus+=1
        elif 's' in item:
            species+=1
        elif 'unclassified' in item:
            unclassified+=1

    y_list_domain.append(domain)
    y_list_phylum.append(phylum)
    y_list_class.append(Class)
    y_list_order.append(order)
    y_list_family.append(family)
    y_list_genus.append(genus)
    y_list_species.append(species)
    y_list_unclassified.append(unclassified)

f2=open(outputfilePath+os.sep+'物种分类学注释结果统计表.csv','w+')
output=','.join(['group','domain','phylum','Class','order','family','genus','species','unclassified'])+'\n'
f2.write(output)
times=0
while times < len(x_list):
    f2.write(','.join([x_list[times],str(y_list_domain[times]),str(y_list_phylum[times]),str(y_list_class[times]),str(y_list_order[times]),str(y_list_family[times]),str(y_list_genus[times]),str(y_list_species[times]),str(y_list_unclassified[times])]))
    f2.write('\n')
    times+=1


######################################################


all=[y_list_domain,y_list_phylum,y_list_class,y_list_order,y_list_family,y_list_genus,y_list_species,y_list_unclassified]
all=list(map(list,zip(*all)))
new=[]
for line in all:
    line=[float(i) for i in line]
    newline=[i/sum(line) for i in line]
    new.append(newline)
new=list(map(list,zip(*new)))

garph = pygal.StackedBar()    # 创建图（叠加柱状图）

# 添加数据
garph.add('Domain', new[0])
garph.add('Phylum', new[1])
garph.add('Class', new[2])
garph.add('Order', new[3])
garph.add('Family', new[4])
garph.add('Genus', new[5])
garph.add('Species', new[6])
garph.add('Unclassified', new[7])


garph.x_labels = x_list            # 设置 X 轴刻度
#garph.y_label_rotation = 45           # 设置Y轴的标签旋转多少度
garph.x_label_rotation = 90

garph.title = title  # 设置图标题
garph.x_title = x_label                 # 设置 X 轴标题
garph.y_title = 'Percentage of '+y_label            # 设置 Y 轴标题

garph.legend_at_bottom = True         # 设置图例位置（下面）

garph.margin = 35    # 设置页边距（margin_bottom、margin_top、margin_left、margin_right）

garph.show_x_guides = False    # 显示X轴的网格线
garph.show_y_guides = False    # 显示Y轴的网格线

garph.render_to_file(outputfilePath+os.sep+'物种分类学注释结果百分比堆积柱状图.svg')
svg_path=outputfilePath+os.sep+'物种分类学注释结果百分比堆积柱状图.svg'
png_path=outputfilePath+os.sep+'物种分类学注释结果百分比堆积柱状图.png'
cairosvg.svg2png(url=svg_path, write_to=png_path)

################################################################################################

garph = pygal.StackedBar()    # 创建图（叠加柱状图）

# 添加数据
garph.add('Domain', y_list_domain)
garph.add('Phylum', y_list_phylum)
garph.add('Class', y_list_class)
garph.add('Order', y_list_order)
garph.add('Family', y_list_family)
garph.add('Genus', y_list_genus)
garph.add('Species', y_list_species)
garph.add('Unclassified', y_list_unclassified)


garph.x_labels = x_list            # 设置 X 轴刻度
#garph.y_label_rotation = 45           # 设置Y轴的标签旋转多少度
garph.x_label_rotation = 90

garph.title = title  # 设置图标题
garph.x_title = x_label                 # 设置 X 轴标题
garph.y_title = 'Number of '+y_label            # 设置 Y 轴标题

garph.legend_at_bottom = True         # 设置图例位置（下面）

garph.margin = 35    # 设置页边距（margin_bottom、margin_top、margin_left、margin_right）

garph.show_x_guides = False    # 显示X轴的网格线
garph.show_y_guides = False    # 显示Y轴的网格线

garph.render_to_file(outputfilePath+os.sep+'物种分类学注释结果堆积柱状图.svg')

svg_path=outputfilePath+os.sep+'物种分类学注释结果堆积柱状图.svg'
png_path=outputfilePath+os.sep+'物种分类学注释结果堆积柱状图.png'

cairosvg.svg2png(url=svg_path, write_to=png_path)





