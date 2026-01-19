#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/04/25 10:33
# Author        : William GoGo
import os
import sys
import pandas as pd
import argparse
import subprocess
from loguru import logger


def _build_rscript_cmd(script_path, **kwargs):
    """构建通用的 Rscript 命令，统一处理参数解析逻辑。"""
    cmd_parts = [f'Rscript {script_path}']
    for key, value in kwargs.items():
        if value is None:
            continue
        # 跳过空字符串（避免传递给 R 脚本时出现问题）
        if isinstance(value, str):
            if not value.strip():
                continue
            # 字符串参数，如果包含空格则加引号
            if ' ' in value:
                arg_val = f'"{value}"'
            else:
                arg_val = value
        # 逻辑值参数转换为 R 格式
        elif isinstance(value, bool):
            arg_val = str(value).upper()
        else:
            arg_val = value

        cmd_parts.append(f'--{key} {arg_val}')

    return ' '.join(cmd_parts)


def enrichment_barplot(input_file, output_file):
    cmd = f"Rscript /home/colddata/qinqiang/script/Analysis/enrich_analysis/enrichment_barplot.r \
        -f {input_file} \
        -o {output_file}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{input_file} enrichment_barplot 画图失败")
        logger.info(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        return True


def draw_pathview(regulation, passed_path):
    cmd = f"Rscript /home/colddata/qinqiang/script/Analysis/pathview/pathview.R \
        -r {regulation} \
        -p {passed_path}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{regulation} draw_pathview 画图失败")
        logger.info(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        
        return False
    else:
        return True
        

def draw_multigroup_heatmap(datafile, output, other_args=''):
    """# 输入文件应该是 3 个 Sheet
    Sheet 1: fpkm value
    ID sampleA-A sampleA-B sampleA-C sampleB-A sampleB-B smapleB-C ...
    geneA fpkm_value ...
    geneB fpkm_value ...
    ...
    
    Sheet 2: sample_annotation
    sample group
    sampleA sampleA-A
    sampleA sampleA-B
    ...
    第三个 Sheet 可以不输入
    Sheet 3: gene_annotation
    gene Ontology
    geneA groupA
    geneB groupB
    genec gorupB
    ...
    """
    if not datafile.endswith('xlsx'):
        logger.error(f'文件不是 xlsx 格式 \n 参考yi下说明 {draw_multigroup_heatmap.__doc__}')
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/heatmap/heatmap_multigroup.r -f {datafile} -o {output} {other_args}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{datafile} multigroup_heatmap 画图失败")
        logger.info(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        
        return False
    else:
        return True


def draw_twogroup_heatmap(datafile, output_file):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/heatmap/heatmap_twogroup.r -f {datafile} -o {output_file}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{datafile} twogroup_heatmap 画图失败")
        logger.info(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        
        return False
    else:
        return True


def anova_analysis(datafile, samples_file, output_file):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/anova/anova.r -f {datafile} -s {samples_file} -o {output_file}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{datafile} anova.r 执行失败")
        logger.info(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        
        return False
    else:
        return True
    
    
def enrich_analysis(input_file, genego_file, keggclean_file, output_dir):
    script_path = '/home/colddata/qinqiang/script/Analysis/enrich_analysis/enrich.r'
    cmd = f'Rscript {script_path} \
        --inputidfile {input_file} \
        --genego {genego_file} \
        --keggclean {keggclean_file} \
        --outputdir {output_dir}'
    logger.debug(f'运行命令 {cmd}')
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"执行失败")
        logger.info(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        
        return False
    else:
        return True


def smart_heatmap(input_file, output_file=None, **kwargs):
    """调用 heatmap_1.r
    
    参数:
        input_file: 输入数据文件（txt/csv/tsv/xlsx/xls），首列为行名（必需）
        output_file: 输出图文件名，自动根据后缀选择格式（pdf/png/tiff/svg/jpg）
        other_args: 额外的命令行参数（字符串），例如 '--no-cluster-rows' 或 '--cluster_rows FALSE'
                    注意：此参数会直接追加到命令末尾，用于向后兼容
        **kwargs: 其他参数，支持 heatmap_1.r 的所有命令行参数，例如：
            - input_sheet: Excel Sheet 名或序号
            - annotation_row: 行注释文件
            - annotation_col: 列注释文件
            - main: 主标题
            - scale: 缩放方式（none/row/column）
            - cluster_rows: 是否聚类行（True/False）
            - cluster_cols: 是否聚类列（True/False）
            - show_rownames: 显示行名（True/False）
            - show_colnames: 显示列名（True/False）
            - fontsize_row: 行名字体大小
            - fontsize_col: 列名字体大小
            - cellwidth: 单元格宽度
            - cellheight: 单元格高度
            - treeheight_row: 行树高度
            - treeheight_col: 列树高度
            - legend: 显示图例（True/False）
            - dpi: 输出 DPI
            - autoset: 启用自动规格设置（True/False）
    
    返回:
        bool: 成功返回 True，失败返回 False
    """
    script_path = '/home/colddata/qinqiang/script/Plot/Heatmap/heatmap_1.r'

    # 构建命令（统一参数处理）
    cmd = _build_rscript_cmd(script_path, input=input_file, output=output_file, **kwargs)
    logger.debug(f'运行命令: {cmd}')
    
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{input_file} smart_heatmap 画图失败")
        logger.debug(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        return True


def volcano_plot(input_file, output_file, **kwargs):
    script_path = '/home/colddata/qinqiang/script/Plot/Volcano/volcano_1.r'

    # 构建命令（统一参数处理）
    cmd = _build_rscript_cmd(script_path, input=input_file, output=output_file, **kwargs)

    # 运行命令
    logger.debug(f'运行命令: {cmd}')
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{input_file} volcano_plot 画图失败")
        logger.debug(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        return True