# -*- coding: UTF-8 -*-

import argparse
import os
from shutil import copyfile
import sys
import numpy

def parse_args():
    parser = argparse.ArgumentParser(description='生成mapping统计汇总文件')
    parser.add_argument('--input-dir', '-i', default='.', 
                      help='输入目录路径，包含mapping.txt文件')
    parser.add_argument('--output-file', '-o', default='mapping_summary.txt',
                      help='输出文件路径')
    return parser.parse_args()

def main():
    args = parse_args()
    
    out_doc = args.output_file
    f1 = open(out_doc, 'w')
    f1.write('Sample'+'\t'+'Total reads'+'\t'+'Mapped reads'+'\t''Unmapped reads'+'\t'+'Unique mapped reads'+'\t'+'Multiple mapped reads'+'\n')
    
    dirlist = os.listdir(args.input_dir)
    
    for dir_name in dirlist:
        if ('_mapping.txt' in dir_name):
            sample = dir_name.split('_mapping.txt')[0]
            f1.write(sample+'\t')
            with open(os.path.join(args.input_dir, dir_name)) as f2:
                data_list = []
                row_num = 0
                for line in f2:
                    if 'nohup' in line:
                        continue
                    line = line.strip()
                    data_list.append(line.split(' ')[0])

                total_reads = int(data_list[0])
                mapped_reads = int(int(data_list[3])+int(data_list[4])+int(data_list[7])+(int(data_list[12])+int(data_list[13]))/2)
                unmapped_reads = total_reads-mapped_reads
                unique_mapped = int(data_list[3])+int(data_list[7])
                multi_mapped_reads = mapped_reads-unique_mapped
                f1.write(str(total_reads)+'\t'+str(mapped_reads)+'\t'+str(unmapped_reads)+'\t'+str(unique_mapped)+'\t'+str(multi_mapped_reads)+'\n')

    f1.close()

if __name__ == '__main__':
    main()



