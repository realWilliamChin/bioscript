
import argparse
import os
from shutil import copyfile
import sys

parser = argparse.ArgumentParser(description="assemble data")

parser.add_argument('--input1', dest='input1', default='./02_cleandata',
                    help='path of the data unsplited')
parser.add_argument('--input2', dest='input2', default='./03_sampled_data',
                    help='path of the data splited')
parser.add_argument('--date_out', dest='data_out', default='./04_assembled_data',
                    help='path of output data')
parser.add_argument('--kmer', dest='kmer', default='65,75,85,95,105',
                    help='the kmer you want')
parser.add_argument('--Ref_Seq', dest='Ref_Seq', default='embplant_pt',
                    help='which Ref_Seq you want')

arg = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
kmer = arg.kmer
Ref_Seq = arg.Ref_Seq
inputPath1 = arg.input1
inputPath2 = arg.input2
outputPath = arg.data_out
filelist1 = os.listdir(inputPath1)
times = 0
if outputPath == './04_assembled_data':
    if os.path.exists(outputPath) == False:
        os.mkdir(outputPath)

generalPath = outputPath + os.sep + 'general'
if os.path.exists(generalPath) == False:
    os.mkdir(generalPath)

for file in filelist1:
    if '_clean' not in file:
        continue
    fileline = file.split('_clean')

    # if '.R1.' in file:
    #     line = file.split('.R1.')
    #     R2file = line[0] + '.R2.' + line[-1]
    if '_1.' in file:
        line = file.split('_1_')
        R2file = line[0] + '_2_' + line[-1]
        os.system('echo running ' + line[0])
        times += 1

        command = 'get_organelle_from_reads.py -t 8' + \
                  ' -1 ' + inputPath1 + os.sep + file + \
                  ' -2 ' + inputPath1 + os.sep + R2file + \
                  ' -o ' + generalPath + os.sep + line[0] + \
                  ' -R 15 -k ' + kmer + ' -F ' + Ref_Seq
        if times % 10 == 0:
            command = command
        else:
            command += ' &'
        os.system(command)

if os.path.exists(inputPath2) == False:
    print('No splice file path, only general data will be assembled')
    sys.exit()

filelist2 = os.listdir(inputPath2)
if len(filelist2) == 0:
    print('No splice file in this path, only general data will be assembled')
    sys.exit()

for Splice in filelist2:
    newoutputPath = outputPath + os.sep + Splice
    if os.path.exists(newoutputPath) == False:
        os.mkdir(newoutputPath)
    sfilelist = os.listdir(inputPath2 + os.sep + Splice)
    for file in sfilelist:
        fileline = file.split('_clean')
        # if '.R1.' in file:
        #     line = file.split('.R1.')
        #     R2file = line[0] + '.R2.' + line[-1]
        if '_1_' in file:
            line = file.split('_1_')
            R2file = line[0] + '_2_' + line[-1]
            os.system('echo running ' + line[0])
            times += 1

            command = 'get_organelle_from_reads.py –t 8' + \
                      ' -1 ' + inputPath2 + os.sep + Splice + os.sep + file + \
                      ' -2 ' + inputPath2 + os.sep + Splice + os.sep + R2file + \
                      ' -o ' + newoutputPath + os.sep + line[0] + \
                      ' -R 15 -k ' + kmer + ' -F ' + Ref_Seq

            if times % 10 == 0:
                command = command
            else:
                command += ' &'

            os.system(command)
