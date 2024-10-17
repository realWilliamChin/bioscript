import argparse
import os
from shutil import copyfile
import sys
parser=argparse.ArgumentParser(description="get alpha rare sheet")
parser.add_argument('--inputPath',dest='inputPath',default='.12_alpha_beta_diversity/metrics',
        help='path of input file')

parser.add_argument('--outputPath',dest='outputPath',default='./12_alpha_beta_diversity/beta_diversity',
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
for file in filelist:
    if '.qzv' in file:
        outputPath=outputfilePath+os.sep+file[:-4]
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
        command='qiime tools export '+'--input-path '+inputfilePath+os.sep+file+' '+'--output-path '+outputPath
        os.system(command)
print('all qzv files were transformed')