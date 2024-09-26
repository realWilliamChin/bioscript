import argparse
import os
from shutil import copyfile
import sys
parser=argparse.ArgumentParser(description="get alpha rare sheet")
parser.add_argument('--inputPath',dest='inputPath',default='./12_alpha_beta_diversity/Alpha_rare',
        help='path of input file')

parser.add_argument('--outputPath',dest='outputPath',default='./12_alpha_beta_diversity',
        help='path of output file')



arg=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputfilePath=arg.inputPath
outputfilePath=arg.outputPath
if inputfilePath=='input_data1':
    os.system('echo need to input inputfilePath')
    sys.exit()
if not os.path.exists(outputfilePath):
    os.makedirs(outputfilePath)

filelist=os.listdir(inputfilePath)
output=[]
IsFirst=True

for file in filelist:
    if '.csv' in file:
        f1=open(inputfilePath+os.sep+file,'r+')
        rawlist = []
        for line in f1:
            rawlist.append(line.split('\n')[0].split(','))
        rawlist_T=list(map(list,zip(*rawlist)))

        if IsFirst:
            IsFirst=False
            output.append(rawlist_T[0])

        newlist=[file[:-4]]
        for line in rawlist[1:]:
            numline=[float(i) for i in line[-11:-1]]
            newlist.append(str(sum(numline)/len(numline)))
        output.append(newlist)
        f1.close()

output_T=list(map(list,zip(*output)))
f2=open(outputfilePath+os.sep+'Alpha_stat.csv','w+')

for line in output_T:
    line=','.join(line)+'\n'
    f2.write(line)