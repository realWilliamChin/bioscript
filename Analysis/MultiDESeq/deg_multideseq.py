#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/4/28 14:46
# Author        : WilliamGoGo
import os, sys
import argparse
import subprocess
import pandas as pd
import openpyxl
from loguru import logger

from deg_comparison_plot import deg_summary_plot
sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from check_SampDesAndCompInfo import check_sample_comp
from load_input import load_table, write_output_df
sys.path.append('/home/colddata/qinqiang/script/Analysis/enrich_analysis')
from deg_enrich import deg_enrich
from deg_enrich_data_merge import deg_enrich_data_merge


if sys.version_info < (3, 10):
    logger.critical("Python 版本低于 3.10，请使用 conda 激活 python310 环境运行程序！")
    logger.critical("当前 Python 版本为:", sys.version)
    sys.exit(1)


def parse_input():
    p = argparse.ArgumentParser(description='输入 kegg, nr, swiss file 的路径')
    p.add_argument('--run-type', dest='run_type', type=str, choices=['deseq', 'limma_deseq'], default='deseq',
                   help='运行类型: deseq 或 limma_deseq（不需要输入 reads_matrix）')
    p.add_argument('--degvalue', type=float, help='deg value FC 值小于多少')
    p.add_argument('--kns', type=str, help='输入 kns_def.txt')
    p.add_argument('--genego', type=str, help='gene_go swiss 注释出来的文件')
    p.add_argument('--keggclean', type=str, help='KEGG_clean.txt kegg 注释出来的文件')
    p.add_argument('-o', '--output-dir', help='分析结果输出目录，默认当前目录', default=os.getcwd())
    
    # 下面输入基本默认即可
    p.add_argument('--samples', type=str, default='samples_described.txt', help='默认 samples_described.txt')
    p.add_argument('--compare', type=str, default='compare_info.txt', help='默认 compare_info.txt')
    p.add_argument('--fpkm', type=str, default='fpkm_matrix_filtered.txt', help='默认 fpkm_matrix_filtered.txt')
    p.add_argument('--reads', type=str, default='reads_matrix_filtered.txt', help='默认 reads_matrix_filtered.txt')
    p.add_argument('--filter-col', type=str, choices=['pvalue', 'padj'], default='padj', help='默认 padj')
    p.add_argument('--filter-value', type=float, default=0.05, help='默认 0.05')
    

    args = p.parse_args()

    return args


def deseq(fpkm_file, reads_file, samples_file, compare_file, filter_col, filter_value, deg_value, output_dir):
    logger.info(f'检查 {samples_file} 和 {compare_file} 文件是否符合要求')
    rep = check_sample_comp(samples_file, compare_file)
    if rep == 1:
        logger.critical(f'样本描述文件 {samples_file} 和比较文件 {compare_file} 不符合要求')
        sys.exit(1)
        
    cmd = f"/home/data/opt/biosoft/R-422/bin/Rscript /home/colddata/qinqiang/script/Analysis/MultiDESeq/multiple_samples_DESeq2.r \
        --degvalue {deg_value} \
        --fpkm {fpkm_file} \
        --reads {reads_file} \
        --samples {samples_file} \
        --compare {compare_file} \
        --filtercol {filter_col} \
        --filtervalue {filter_value} \
        --outputdir {output_dir}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"multiple_samples_DESeq2 程序运行失败")
        logger.info(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        return True


def deseq_limma(fpkm_file, samples_file, compare_file, filter_col, filter_value, deg_value, output_dir):
    logger.info(f'检查 {samples_file} 和 {compare_file} 文件是否符合要求')
    rep = check_sample_comp(samples_file, compare_file)
    if rep == 1:
        logger.critical(f'样本描述文件 {samples_file} 和比较文件 {compare_file} 不符合要求')
        sys.exit(1)
        
    cmd = f"/home/data/opt/biosoft/R-422/bin/Rscript /home/colddata/qinqiang/script/Analysis/MultiDESeq/limma.R \
        --degvalue {deg_value} \
        --fpkm {fpkm_file} \
        --samples {samples_file} \
        --compare {compare_file} \
        --filtercol {filter_col} \
        --filtervalue {filter_value} \
        --outputdir {output_dir}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"limma 程序运行失败")
        logger.info(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        return True


def deg_enrich_distribution(work_dir):
    cur_dir = os.getcwd()
    os.chdir(work_dir)
    logger.info(f'运行 enrich_distribution_plot 脚本中')
    cmd = f"/home/data/opt/biosoft/R-422/bin/Rscript /home/colddata/qinqiang/script/Analysis/enrich_analysis/enrich_distribution_plot.r"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os.chdir(cur_dir)
    if ret.returncode != 0:
        logger.error(f"enrich_distribution_plot.r 程序运行失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        return True
        

def process_deresults(de_results_file, kns_df, run_type='deseq'):
    if de_results_file.endswith('xlsx'):
        de_df = load_table(de_results_file, dtype=str)
        readcounts_file = de_results_file.replace('.xlsx', '') + '_readCounts.matrix'
    else:
        de_df = load_table(de_results_file, sep='\t', dtype=str)
        readcounts_file = de_results_file + '_readCounts.matrix'
    
    # 只有在 deseq 模式下才尝试加载 readCounts.matrix 文件
    if run_type == 'deseq' and os.path.exists(readcounts_file):
        # readCounts.matrix 为制表符分隔，需显式指定 sep='\t'
        de_reads_df = load_table(readcounts_file, sep='\t', dtype=str)
        # 统一首列列名，确保与 DE 结果的键一致
        first_col = de_reads_df.columns.tolist()[0]
        if first_col != 'GeneID':
            de_reads_df.rename(columns={first_col: 'GeneID'}, inplace=True)
        de_reads_df.columns = [de_reads_df.columns.tolist()[0]] + [x + '_raw_reads' for x in de_reads_df.columns.tolist()[1:]]
        de_df = pd.merge(left=de_df, right=de_reads_df, on='GeneID', how='left')
    else:
        logger.info(f'跳过 readCounts.matrix 文件加载 (运行类型: {run_type}, 文件存在: {os.path.exists(readcounts_file) if run_type == "deseq" else "N/A"})')
    
    # 排序 down，up，NOsig。down 的 FC 值从小到大，up 的 FC 值从大到小
    # 先分三份，再合并
    de_df_down = de_df[de_df['regulation'] == 'Down'].copy()
    de_df_down.sort_values(by='FC', ascending=True, inplace=True)
    de_df_up = de_df[de_df['regulation'] == 'Up'].copy()
    de_df_up.sort_values(by='FC', ascending=False, inplace=True)
    de_df_nosig = de_df[de_df['regulation'] == 'NoSignificant'].copy()
    # 合并
    de_df = pd.concat([de_df_down, de_df_up, de_df_nosig])
    de_df = pd.merge(left=de_df, right=kns_df, on='GeneID', how='left')
    return de_df


def main():
    args = parse_input()
    
    # 运行 R 脚本
    if args.degvalue:
        if args.run_type == 'deseq':
            logger.info(f'正在执行 deseq 分析, {args.degvalue}')
            rep = deseq(args.fpkm, args.reads, args.samples, args.compare, args.filter_col, args.filter_value, args.degvalue, args.output_dir)
        else:
            logger.info(f'正在执行 limma_deseq 分析, {args.degvalue}')
            rep = deseq_limma(args.fpkm, args.samples, args.compare, args.filter_col, args.filter_value, args.degvalue, args.output_dir)
        if not rep:
            logger.critical(f'R 脚本运行失败')
            sys.exit(1)
        logger.info(f'对 DEG_summary 画图')
        deg_summary_df = load_table(os.path.join(args.output_dir, 'DEG_analysis_results', 'DEG_summary.txt'), comment='#', skipinitialspace=True)
        deg_summary_plot(deg_summary_df, os.path.join(args.output_dir, 'DEG_analysis_results', 'DEG_summary_plot.jpeg'))
    else:
        logger.info(f'未输入 degvalue，不执行 R 分析，将针对现有分析结果进行处理')
    if args.genego and args.keggclean:
        logger.info('正在执行 deg enrich')
        enrich_dir = os.path.join(args.output_dir, 'Pathway_enrichment_analysis')
        os.makedirs(enrich_dir, exist_ok=True)
        deg_enrich(
            compare = args.compare,
            degdata_dir = os.path.join(args.output_dir, 'DEG_analysis_results'),
            genego_file = args.genego,
            keggclean_file = args.keggclean,
            outputdir = enrich_dir
        )
        logger.info('正在执行 deg distribution 画图')
        deg_enrich_distribution(enrich_dir)
        deg_enrich_data_merge(
            os.path.join(args.output_dir, 'Pathway_enrichment_analysis', 'Pathway_enrichment_raw_data'),
            args.compare,
            os.path.join(args.output_dir, 'Pathway_enrichment_analysis', 'Pathway_enrichment_raw_data', 'DEG_enrichment_significant_pathway_summary.xlsx')
        )
    
    de_results_output_dir = os.path.join(args.output_dir, 'DEG_analysis_results/Expression_data')

    kns_df = load_table(args.kns, dtype={"GeneID": str})
    for de_results_file in os.listdir(args.output_dir):
        if de_results_file.endswith('DE_results') or de_results_file.endswith('DE_results.xlsx'):
            logger.info(f'processing---{de_results_file}')
            de_df = process_deresults(de_results_file, kns_df, args.run_type)
            deg_filename = os.path.basename(de_results_file).replace('DE_results', 'DEG_data.txt')
            deg_data_file = os.path.join(de_results_output_dir, deg_filename)
            write_output_df(de_df, deg_data_file, index=False)
        

    for up_down_id_file in os.listdir("DEG_analysis_results"):
        up_down_id_file = os.path.join("DEG_analysis_results", up_down_id_file)
        if up_down_id_file.endswith('Down_ID.txt') or up_down_id_file.endswith('Up_ID.txt'):
            up_down_id_df = load_table(up_down_id_file, dtype={"GeneID": str}, header=None, names=['GeneID'])
            result_df = pd.merge(left=up_down_id_df, right=kns_df, on='GeneID', how='left')
            write_output_df(result_df, up_down_id_file.replace('.txt', '_def.txt'), index=False)

    png_file = 'deg_line_plot.jpeg'
    count_summary_file = 'all_deg_count_summary.txt'
    
    if os.path.exists(png_file):
        os.rename(png_file, os.path.join('Expression_data_evaluation', png_file))
    if os.path.exists(count_summary_file):
        os.rename(count_summary_file, os.path.join('Expression_data_evaluation', count_summary_file))
    
    logger.success('Done!')


if __name__ == '__main__':
    main()


