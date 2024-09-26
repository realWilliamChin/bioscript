import argparse
import os
from shutil import copyfile
import sys
parser=argparse.ArgumentParser(description="creat_manifest_txt:")
parser.add_argument('--dataset_dir',dest='dataset_dir',default='input_data1',
        help='path of the data with data name')
'''
parser.add_argument('--continue_train',dest='continue_train',type=bool,default=False,
        help='if continue training, load the latest model: 1: true, 0: false')
'''
parser.add_argument('--date_out',dest='data_out',default='output_data1',
        help='path of output data with data name')

arg=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputfilePath=arg.dataset_dir
outputfilePath=arg.data_out
Dic=open(r'./06_convertResult/ASV_Hash_Dic.txt','r+')
f1=open(inputfilePath,'r+')
f2=open(outputfilePath,'w+')

ASVlist=[]
Hashlist=[]
for line in Dic:
    newline=line.split('\n')[0].split('\t')
    ASVlist.append(newline[0])
    Hashlist.append(newline[-1])

for line in f1:
    if '>' in line:
        HashID=line.split('\n')[0].split('>')[-1]
        f2.write('>'+ASVlist[Hashlist.index(HashID)]+'\n')
    else:
        f2.write(line)


f1.close()
f2.close()
Dic.close()