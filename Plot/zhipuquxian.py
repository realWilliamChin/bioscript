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
from matplotlib import ticker
from loguru import logger


def parse_input():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--input', type=str, help='输入质谱数据的 excel 文件')
    # parser.add_argument('-o', '--output', type=str, help='输出目录')
    
    return parser.parse_args()


def sci_notation(x, pos):  
    # 将浮点数转换为科学计数法的字符串  
    a, b = '{:.0e}'.format(x).split('e')  
    b = int(b)  
    # return r'${} \times 10^{{{}}}$'.format(a, b)
    if b < 10:
        b = f'0{b}'
    if str(a) == str(0):
        return str(0)
    else:
        return f'{a}E+{b}'


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
        wavelength = pd.read_excel(args.input, engine='openpyxl', header=None, sheet_name=sheet)
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
            
            # plt.figure(figsize=(3 + x_width, 5 + y_height))
            plt.figure(figsize=(15, 10))
            plt.xticks(np.arange(0, max_x + 1, step=0.5))
            plt.yticks(np.arange(0, (max_y + 100000), step=(max_y / 10)))
            plt.xlim(-(max_x / 20), max_x + (max_x / 20))
            plt.ylim(-(max_y / 12), max_y + (max_y / 12))
            plt.plot(sample_df['time(min)'], sample_df['relative intension'], alpha=0.5, linewidth=1, label='acc', color='black')
            # plt.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
            plt.xlabel('Time (min)')
            plt.ylabel('Relative Intensity')
            plt.title(f'Chromagram ({wavelength} nm)')
            
            # 创建一个FuncFormatter实例，并将其应用于y轴的刻度  
            formatter = ticker.FuncFormatter(sci_notation)  
            plt.gca().yaxis.set_major_formatter(formatter) 
            
            plt.savefig(f"{sheet}/{sample_name}_{input_name.replace('.xlsx', '')}_Chromagram.jpeg")


if __name__ == '__main__':
    main()