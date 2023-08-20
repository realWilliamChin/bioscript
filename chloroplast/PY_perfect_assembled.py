#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/08/05 19:07
# Author        : William GoGo
import argparse
import os
from shutil import copyfile
import sys


def parse_input():
    argparser = argparse.ArgumentParser(description="perfect assembled data")
    argparser.add_argument('--input1', dest='input1', default='./05_ct_sequence_right_circle',
                        help='path of the fasta')
    argparser.add_argument('--input2', dest='input2', default='./05_first_GB_file',
                        help='path of the .gb')
    argparser.add_argument('--date_out', dest='data_out', default='./06_psba_motified_sequence',
                        help='path of output data')
    return argparser.parse_args()


def main():
    arg = parse_input()

    inputPath1 = arg.input1
    inputPath2 = arg.input2
    outputPath = arg.data_out
    filelist1 = os.listdir(inputPath1)
    filelist2 = os.listdir(inputPath2)

    if outputPath == './06_psba_motified_sequence':
        if os.path.exists(outputPath) == False:
            os.mkdir(outputPath)

    for seqfile in filelist1:
        if '.fa' in seqfile:
            seqname = seqfile.split('.fa')[0]
            for gbfile in filelist2:

                if seqname in gbfile:
                    print(seqfile)
                    print(gbfile)
                    f1=open(inputPath1+os.sep+seqfile,'r+')
                    f2=open(inputPath2+os.sep+gbfile,'r+')
                    f3=open(outputPath+os.sep+seqname+'.fasta','w+')
                    gblist=[]
                    times=0
                    for line in f2:
                        gblist.append(line)
                        if '                     /gene="psbA"' in line:
                            break
                        times+=1
                    while times>=0:
                        if '     gene            complement(' in gblist[times]:
                            posline=gblist[times]
                            break
                        else:
                            times-=1

                    pos=posline.split('     gene            complement(')[-1]
                    pos=pos.split('join(')[-1]
                    start=int(pos.split('..')[0])-5
                    seqlist=[]
                    for line in f1:
                        if '>' in line:
                            name = line
                            sequence=''
                        else:
                            sequence += line.strip()
                    output=(sequence[start:]+sequence[:start])
                    #for line in output:
                    f3.write(name+output+'\n')


    print('finish')



















