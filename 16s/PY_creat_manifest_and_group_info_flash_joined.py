import argparse
import os
import sys
from pathlib import Path
from loguru import logger

def parse_arguments():
    parser = argparse.ArgumentParser(description="Create manifest and group info files for 16S data")
    parser.add_argument('--dataset_dir', dest='dataset_dir', default='input_data1',
            help='Path of the input data directory')
    parser.add_argument('--data_out', dest='data_out', default='output_data1',
            help='Path of output data directory')
    args = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    return args

def extract_class_name(name):
    if not name:
        logger.error("Empty sample name encountered")
        return "unknown"
    
    try:
        times = -1
        while times >= -len(name) and name[times].isdecimal():
            times -= 1
        if times >= -len(name) and name[times] == '-':
            times -= 1
        return name[:times + 1]
    except Exception as e:
        logger.error(f"Error processing sample name '{name}': {str(e)}")
        return "unknown"

def main():
    args = parse_arguments()
    
    input_path = Path(args.dataset_dir)
    output_path = Path(args.data_out)
    
    # 确保输出目录存在
    output_path.mkdir(parents=True, exist_ok=True)
    
    if not input_path.exists():
        logger.error(f"Input directory {input_path} does not exist")
        sys.exit(1)
    
    try:
        file_list = sorted(input_path.iterdir())
        current_path = Path.cwd()
        
        manifest_path = output_path / 'manifest.txt'
        group_info_path = output_path / 'group_info.tsv'
        
        with open(manifest_path, 'w') as f1, open(group_info_path, 'w') as f2:
            f1.write('sample-id,absolute-filepath,direction\n')
            f2.write('sample-id\tclass\n')
            
            for file_path in file_list:
                if not file_path.is_file():
                    continue
                    
                file_name = file_path.name
                logger.info(f"Processing file: {file_name}")
                
                try:
                    name_parts = file_name.split('_')
                    name_parts = [part.split('.')[0] for part in name_parts]
                    sample_id = name_parts[0]
                    
                    class_name = extract_class_name(sample_id)
                    
                    # 构建绝对路径
                    abs_path = current_path / input_path / file_name
                    
                    f1.write(f'{sample_id},{abs_path},forward\n')
                    f2.write(f'{sample_id}\t{class_name}\n')
                except Exception as e:
                    logger.error(f"Error processing file {file_name}: {str(e)}")
                    continue
                
        logger.info(f"Successfully created manifest and group info files in {output_path}")
        
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()