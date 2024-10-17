import argparse
import os
from shutil import copyfile
import sys
parser=argparse.ArgumentParser(description="add group info title to lefse table")
parser.add_argument('--input1',dest='input1',default='./00_database/group_info.tsv',
        help='path of input table with table name')

parser.add_argument('--input2',dest='input2',default='./16_LEfse/collapse.frequency.table.lefse.tsv',
        help='path of tax file with tax name')

parser.add_argument('--outputPath',dest='outputPath',default='./16_LEfse/collapse.frequency.table.lefse_v2.tsv',
        help='path of output file')


arg=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputfilePath1=arg.input1
inputfilePath2=arg.input2
outputfilePath=arg.outputPath


infofile=open(inputfilePath1,'r+')
tablefile=open(inputfilePath2,'r+')
outputfile=open(outputfilePath,'w+')

infor_table=[]
for line in infofile:
        line=line.split('\n')[0].split('\t')
        infor_table.append(line)

infor_table=list(map(list,zip(*infor_table)))
for line in infor_table:
    outputfile.write('\t'.join(line)+'\n')

for line in tablefile:
    if 'Group' in line or '|' not in line or 'd__Bacteria' not in line:
        continue
    outputfile.write(line)
