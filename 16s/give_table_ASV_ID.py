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
        help='path of output data')

arg=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputfilePath=arg.dataset_dir
outputPath=arg.data_out

f1=open(inputfilePath,'r+')
f2=open(outputPath+os.sep+'table_ASV.tsv','w+')
f3=open(outputPath+os.sep+'ASV_Hash_Dic.txt','w+')

ASVlist=[]
Hashlist=[]
times=1
for line in f1:
    if '# Constructed from biom file' in line:
        f2.write(line)
        continue
    elif 'ASV ID' in line:
        f2.write(line)
        continue
    else:
        newline=line.split('\t')
        ASVlist.append('ASV'+str(times))
        Hashlist.append(newline[0])
        outputline=['ASV'+str(times)]+newline[1:]
        times+=1
        f2.write('\t'.join(outputline))

times=0
while times<len(ASVlist):
    line='\t'.join([ASVlist[times],Hashlist[times]])
    f3.write(line+'\n')
    times+=1


f1.close()
f2.close()
f3.close()