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


def draw_barplot(input_file, output_file):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/barplot/barplot.R \
        -f {input_file} \
        -o {output_file}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{input_file} draw_barplot 画图失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        return True


def draw_pathview(regulation, passed_path):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/pathview/pathview.R \
        -r {regulation} \
        -p {passed_path}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{regulation} draw_pathview 画图失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
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
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        
        return False
    else:
        return True


def draw_twogroup_heatmap(datafile, output_file):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/heatmap/heatmap_twogroup.r -f {datafile} -o {output_file}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{datafile} twogroup_heatmap 画图失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        
        return False
    else:
        return True


def anova_analysis(datafile, samples_file, output_file):
    cmd = f"Rscript /home/colddata/qinqiang/script/Rscript/anova/anova.r -f {datafile} -s {samples_file} -o {output_file}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"{datafile} anova.r 执行失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        
        return False
    else:
        return True
    
    
def enrich_analysis(input_file, genego_file, keggclean_file, output_dir):
    script_path = '/home/colddata/qinqiang/script/Rscript/enrich_analysis/enrich.r'
    cmd = f'Rscript {script_path} \
        --inputidfile {input_file} \
        --genego {genego_file} \
        --keggclean {keggclean_file} \
        --outputdir {output_dir}'
    logger.info(f'运行命令 {cmd}')
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"执行失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        
        return False
    else:
        return True