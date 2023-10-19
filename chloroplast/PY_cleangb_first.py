# -*- coding: UTF-8 -*-

import argparse
import os
from shutil import copyfile
import sys
parser=argparse.ArgumentParser(description="clean .gb file at the first time")

parser.add_argument('--input',dest='input',default='input_data',
        help='path of the data')

parser.add_argument('--date_out',dest='data_out',default='output_data1',
        help='path of output data')

arg=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)


inputPath=arg.input
outputPath=arg.data_out
filelist=os.listdir(inputPath)



#判断fragment的去留
'''
def CouldUseFragment():
    sdas
    if sadasdsa:
        return True
    else:
        return False
'''
#每一个part第二行，必须是 /gene
def ComplateGenePart(GenPart):
    global file
    global geneNum
    #GenPart是一个大块(list格式)，继续细分研究
    output=[]
    steps=0
    #geneNum+=1
    if geneNum>=1000:
        Num=str(geneNum)
    else:
        Num=(3-len(str(geneNum)))*'0'+str(geneNum)
    locus_tag=file.split('GenBank')[0]+Num
    protein_id=locus_tag+'p'
    #'/codon_start=1'
    #'/transl_table=11'
    add=1
    intron_num=1
    exon_num=1
    while steps+add<len(GenPart):
        line=GenPart[steps+add]
        #output.append(line)
        #继续分块，第一行必然是gene
        if '                     ' not in line or steps+add==len(GenPart)-1:
            newlist=GenPart[steps:steps+add]
            if steps+add==len(GenPart)-1:
                newlist=GenPart[steps:]
            steps += add
            add = 1
            titleline = newlist[0]
            if '     gene            ' in titleline:
                #保留坐标，名字，tag，
                output+=gene(newlist,locus_tag)
            elif '     CDS             ' in titleline:
                #需要检测坐标行数
                #保留坐标，名字，conda,table,tag,id
                output+=CDS(newlist,locus_tag,protein_id)

            elif '     intron          ' in titleline or '     exon            ' in titleline:
                #保留坐标，名字，number
                if '     intron          ' in titleline:
                    output+=ON(newlist,intron_num)
                    intron_num+=1
                else:
                    output += ON(newlist,exon_num)
                    exon_num+=1
            elif '     tRNA            ' in titleline or '     rRNA            ' in titleline:
                #保留坐标，名字，tag
                output+=RNA(newlist,locus_tag)
            elif '     repeat_region   ' in titleline:
                #保留坐标/info/annotator/rpt_type
                output+=repeat(newlist)
        else:
            add+=1
    return output


def repeat(list):
    output=[list[0]]
    for line in list[1]:
        if '                     /info' in line:
            output.append(line)
        elif '                     /annotator' in line:
            output.append(line)
        elif '                     /rpt_type' in line:
            output.append(line)
    return output


def gene(list,tag):
    output=[]
    for i in list[:2]:
        output.append(i.split('\n')[0])
    output.append('                     '+'/locus_tag="'+tag+'"')
    return output

def CDS(list,tag,id):
    output=[]
    times=1

    a = list[0].strip('\n')
    while times<len(list):
        if '/gene' not in list[times]:
            b=list[times].strip('\n')
            c = b.split(' ')[-1]
            a+=c
            times+=1
        else:
            break

    output.append(a)
    allterm=''.join(list)
    if '                     /gene="rps12"' in allterm:
        output.append('                     /trans_splicing')
    for line in list:
        if '/gene=' in line:
            output.append(line.strip('\n'))
    #output.append(list[times].split('\n')[0])
    output.append('                     /codon_start=1')
    output.append('                     /transl_table=11')
    output.append('                     /locus_tag="'+tag+'"')
    output.append('                     /protein_id="'+id+'"')

    return output

def ON(list,number):
    output=[]
    titile=list[0].split('\n')[0]
    if 'complement(join(' in titile:
        newline=titile.split('join(')
        Rpart=newline[1].split(',')[0]
        titile=newline[0]+Rpart+')'
    output.append(titile)
    for line in list[1:]:
        if '/gene' in line:
            b=line

    c='                     /number='+str(number)
    output.append(c)
    output.append(b)
    return output

def RNA(list,tag):
    output=[]
    output.append(list[0].split('\n')[0])
    output.append(list[1].split('\n')[0])
    output.append('                     /locus_tag="' + tag + '"')
    return output




for file in filelist:
    if 'gb' in file:
        f1=open(inputPath+os.sep+file,'r+')
        f1line=[]

        f2 = open(outputPath + os.sep + file.split('GenBank')[0]+'Clean.gb','w+')

        output=[]
        for i in f1:
            f1line.append(i.split('\n')[0])

        Top=[]
        Tail=[]
        times=0
        UseTop=True
        #区分三大块
        while times<len(f1line):
            line=f1line[times]
            if UseTop and '     gene            ' in line:
                Top=f1line[:times]
                UseTop=False
                Start=times
            elif 'ORIGIN' in line:
                Tail=f1line[times:]
                End=times
                break
            times+=1

        Part=f1line[Start:End]
        times=0
        add=1
        geneNum=0
        Work=[]
        #细分每一个part
        while times+add<len(Part):
            line=Part[times+add]

            if '     gene            ' in line or '     misc_feature    ' in line or times+add==len(Part)-1:
                onePart=Part[times:times+add]
                if times+add==len(Part)-1:
                    onePart = Part[times:]
                times+=add
                add=1
                if '     misc_feature    ' in onePart[0]:
                    #print('#########################')
                    #print('\n'.join(onePart))
                    continue

                else:
                    all=''.join(onePart)
                    if 'rps19-fragment' in all:
                        #print(onePart)
                        continue
                    elif 'fragment' in all and 'rRNA' in all and 'coverage' in all:
                        list=all.split(', ')
                        for i in list:
                            if 'coverage' in i:
                                new=i
                        perc=new.split('coverage')[-1][1:-1]
                        if float(perc) <50:
                            #print(onePart)
                            continue
                        else:
                            geneNum += 1
                            Work += ComplateGenePart(onePart)
                    else:
                        geneNum+=1
                        Work+=ComplateGenePart(onePart)
            else:
                add+=1

        output=Top+Work+Tail















        for i in output:
            f2.write(str(i)+'\n')

        f1.close()
        f2.close()
