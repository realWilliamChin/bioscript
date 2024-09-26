import argparse
import os
from shutil import copyfile
import sys


parser=argparse.ArgumentParser(description="dumx_to_join")
parser.add_argument('--inputPath',dest='inputPath',default='input_data1',
        help='path of input file')

parser.add_argument('--outputPath',dest='outputPath',default='NA',
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
if outputfilePath=='NA':
    outputfilePath='./01_1_single_rawdata'
    if not os.path.exists(outputfilePath):
        os.makedirs(outputfilePath)
    os.system('echo your outputfilePath is '+outputfilePath)



    command=arg.type+' -db '+Path+' -query '+inputfilePath+' -out '+outputfilePath+os.sep+outputname+' -max_target_seqs '+arg.max_target_seqs+' -evalue 1e-5 -num_threads '+arg.threads+' -outfmt "6 qacc sacc pident qcovs qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxid ssciname" '
    os.system(command)
    os.system('echo running '+arg.type+' '+database)


def R(Rname):
    if 'R' in Rname:
        return Rname
    else:
        return 'R'+Rname

def splice(filename):
    outputlist=[]
    for i in filename.split('.'):
        for y in i.split('_'):
            for x in y .split('-'):
                outputlist.append(x)
    return outputlist

shortnamelist=[]
fullnamelist=[]

filelist=os.listdir(inputfilePath)
for filename in filelist:
    partlist=splice(filename)
    times=-1
    if partlist[times]=='gz':
        tail='.'.join(partlist[-2:])
        times=-3
    else:
        tail=partlist[times]
        times=-2

    Rname=R(partlist[times])
    Rstation=times
    times-=1

    while True:
        if times==0-len(partlist):
            shortname='-'.join(partlist[times:Rstation])
            shortnamelist.append(shortname)
            fullnamelist.append(filename)
            break
        if partlist[times].isdecimal() or partlist[times]=='16S' or partlist[times]=='16s' or partlist[times]=='18S' or partlist[times]=='18s' or partlist[times]=='ITS':
            times-=1
            continue
        else:
            shortname='-'.join(partlist[times:Rstation])
            shortnamelist.append(shortname)
            fullnamelist.append(filename)
            break

used=[]
times=0
while times<len(shortnamelist):
    groupA=shortnamelist[times]

    if groupA not in used:
        add = 1
        while times+add<len(shortnamelist):
            groupB=shortnamelist[times+add]
            if groupA == groupB:
                used.append(groupA)
                R1=inputfilePath+os.sep+fullnamelist[times]
                R2=inputfilePath+os.sep+fullnamelist[times+add]
                command='flash '+R1+' '+R2+' -p 33 -r 250 -f 500 -s 100 -o '+outputfilePath+os.sep+groupA
                os.system(command)
                os.system('echo running ' +groupA)
                add+=1
                break
            else:
                add+=1

    times+=1
