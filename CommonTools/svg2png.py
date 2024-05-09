#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/4/28 14:46
# Author        : WilliamGoGo
import os
import argparse
import cairosvg

def convert_svg_to_png(input_folder, output_folder):
    # 确保输出文件夹存在
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 获取输入文件夹中所有的 SVG 文件
    svg_files = [f for f in os.listdir(input_folder) if f.endswith('.svg')]

    # 遍历所有的 SVG 文件并进行转换
    for svg_file in svg_files:
        input_path = os.path.join(input_folder, svg_file)
        output_path = os.path.join(output_folder, svg_file.replace('.svg', '.png'))
        cairosvg.svg2png(url=input_path, write_to=output_path)
        print(f"Converted {svg_file} to PNG")

    print("Conversion complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert SVG files to PNG files.')
    parser.add_argument('input_folder', type=str, help='Input folder containing SVG files')
    parser.add_argument('output_folder', type=str, help='Output folder to save PNG files')

    args = parser.parse_args()
    input_folder = args.input_folder
    output_folder = args.output_folder

    convert_svg_to_png(input_folder, output_folder)
