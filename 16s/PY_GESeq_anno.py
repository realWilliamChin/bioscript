import argparse
import os
from shutil import copyfile
import sys
parser=argparse.ArgumentParser(description="GESeq_anno")
parser.add_argument('--input1',dest='input1',default='./17_DESeq_result/table_filtered_depth.txt.WA_vs_WD.DESeq2.DE_results',
        help='path of input table with table name')

parser.add_argument('--input2',dest='input2',default='./07_Taxonomy/taxonomy.tsv',
        help='path of tax file with tax name')

parser.add_argument('--outputPath',dest='outputPath',default='../17_DESeq_result/GESeq_anno.txt',
        help='path of output file')





arg=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputfilePath1=arg.input1
inputfilePath2=arg.input2
outputfilePath=arg.outputPath


Dicfile=open(inputfilePath2,'r+')
tablefile=open(inputfilePath1,'r+')
ASV_list=[]
Anno_list=[]
for line in Dicfile:
    if '#Feature' in line:
        continue
    line=line.split('\n')[0].split('\t')
    ASV_list.append(line[0])
    Anno_list.append(line[1])

outputfile=open(outputfilePath,'w+')
for line in tablefile:
    if 'sampleA' in line:
        line=line.split('\n')[0]+'\tAnno\n'
    else:
        ID=line.split('\t')[0]
        line=line.split('\n')[0]+'\t'+Anno_list[ASV_list.index(ID)]+'\n'

    outputfile.write(line)



