#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2024/5/28 14:46
# Author        : WilliamGoGo
import os, sys
import openpyxl
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from loguru import logger


def parse_input():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--input', type=str, help='输入质谱数据的 excel 文件')
    # parser.add_argument('-o', '--output', type=str, help='输出目录')
    
    return parser.parse_args()


def main():
    args = parse_input()
    input_name = args.input
    
    xlsx = pd.ExcelFile(input_name)
    sheet_names = xlsx.sheet_names
    
    logger.info(f'数据表总共 {len(sheet_names)} 个 sheet, {sheet_names}')
    
    for sheet in sheet_names:
        os.mkdir(sheet)
        data_df = pd.read_excel(input_name, engine='openpyxl', skiprows=2, sheet_name=sheet)
        data_df = data_df.drop(columns=['原始编号'])
        
        # 读取第一行，第三列的 wavelength 值
        wavelength = pd.read_excel(args.input, engine='openpyxl', header=None)
        wavelength = wavelength.iloc[0, 2]

        data_df_columns = data_df.columns.to_list()
        data_df_columns = [data_df_columns[i:i+2] for i in range(0, len(data_df_columns), 2)]

        max_x = 0
        max_y = 0
        for each_sample in data_df_columns:
            sample_name = each_sample[0]
            sample_df = data_df[each_sample].copy()
            sample_df.columns = sample_df.iloc[0]
            sample_df.drop([0], inplace=True)
            x = sample_df['time(min)'].max()
            y = sample_df['relative intension'].max()
            
            if x > max_x:
                max_x = x
            if y > max_y:
                max_y = y
        
        logger.info(f'{max_x}, {max_y}')

        for each_sample in data_df_columns:
            sample_name = each_sample[0]
            logger.info(f"processing {sheet} {sample_name}")
            sample_df = data_df[each_sample].copy()
            sample_df.columns = sample_df.iloc[0]
            sample_df.drop([0], inplace=True)
            x_width = int(max_x * 2)
            y_height = int(max_y / 100000)
            plt.figure(figsize=(2 + x_width, 2 + y_height))
            plt.xticks(np.arange(0, max_x + 1, step=0.5))
            plt.yticks(np.arange(0, max_y + 100000, step=50000))
            plt.xlim(-0.1, max_x + 1)
            plt.ylim(-23000, max_y + 50000)
            plt.plot(sample_df['time(min)'], sample_df['relative intension'], alpha=0.5, linewidth=1, label='acc', color='black')
            plt.xlabel('Time (min)')
            plt.ylabel('Relative Intensity')
            plt.title(f'Chromagram ({wavelength} nm)')
            plt.savefig(f"{sheet}/{sample_name}_{input_name.replace('.xlsx', '')}_Chromagram.jpeg")


if __name__ == '__main__':
    main()