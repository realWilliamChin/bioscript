import os
import sys
import argparse
import pandas as pd
import pygal
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(
        description="stat table and output stacked histogram from tax in different tax"
    )
    parser.add_argument('-i', '--input-dir', help='输入目录')
    parser.add_argument('-o', '--output-dir', help='输出目录')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    return args


def main():
    args = parse_input()

    input_dir = args.input_dir
    output_dir = args.output_dir

    # 分类层级名称映射：level-1 对应 Domain，level-2 对应 Phylum，以此类推
    class_list = ['group', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    # 处理 level-1.csv 到 level-7.csv
    for level in range(1, 8):
        input_filename = f'level-{level}.csv'
        input_path = os.path.join(input_dir, input_filename)
        
        # 如果文件不存在，跳过
        if not os.path.exists(input_path):
            logger.warning(f'{input_filename} not found')
            continue

        tax_level_name = class_list[level]
        output_csv_path = os.path.join(output_dir, f'{tax_level_name}.csv')

        # 使用 pandas 读取 CSV 文件并转置
        df = load_table(input_path)
        df_T = df.T  # 转置
        df_T = df_T.fillna('')  # 将 NaN 值填充为空字符串
        
        # 写出一个新的 csv（保持原逻辑，仅去掉最后一行）
        write_output_df(df_T[:-1], output_csv_path, header=False)



if __name__ == '__main__':
    main()

