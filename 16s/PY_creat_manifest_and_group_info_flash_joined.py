import argparse
import os
from shutil import copyfile
import sys
parser=argparse.ArgumentParser(description="creat_manifest_txt:")
parser.add_argument('--dataset_dir',dest='dataset_dir',default='input_data1',
        help='path of the data')
'''
parser.add_argument('--continue_train',dest='continue_train',type=bool,default=False,
        help='if continue training, load the latest model: 1: true, 0: false')
'''
parser.add_argument('--date_out',dest='data_out',default='output_data1',
        help='path of output data')
arg=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)


inputPath=arg.dataset_dir
outputPath=arg.data_out
filelist=os.listdir(inputPath)
f1=open(outputPath+os.sep+'manifest.txt','w+')
f2=open(outputPath+os.sep+'group_info.tsv','w+')
f1.write('sample-id,absolute-filepath,direction'+'\n')
f2.write('sample-id'+'\t'+'class'+'\n')
filelist.sort()
path=os.getcwd()
for file in filelist:
    oldneed=file.split('_')
    need=[]
    for i in oldneed:
        need+=i.split('.')
    name1=need[0]

    times = -1
    while name1[times].isdecimal():
        times -= 1
    if name1[times]=='-':
        times-=1
    classname = name1[:times + 1]


    f1.write(name1+','+path+inputPath.split('.')[-1]+os.sep+file+','+'forward'+'\n')
    f2.write(name1 + '\t' + classname + '\n')




f1.close()
f2.close()