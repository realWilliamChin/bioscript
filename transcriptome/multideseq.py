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

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools'))
from genedf_add_expression_and_def import add_kns_def
from check_SampDesAndCompInfo import check_sample_comp
from transcriptome_enrich import transcriptome_enrich


def parse_input():
    parser = argparse.ArgumentParser(description='输入 kegg, nr, swiss file 的路径')
    parser.add_argument('-w', '--workdir', type=str, help='输入工作目录，默认当前目录', default='./')
    parser.add_argument('-k', '--kegg', type=str, help='kegg_gene_def file')
    parser.add_argument('-n', '--nr', type=str, help='nr_gene_def file')
    parser.add_argument('-s', '--swiss', type=str, help='swiss_gene_def file')
    parser.add_argument('--kns', type=str, help='输入 kns_def.txt，则不用输入上面的三个')
    parser.add_argument('--degvalue', type=float, help='deg value FC 值小于多少')
    parser.add_argument('--genego', type=str, help='gene_go swiss 注释出来的文件')
    parser.add_argument('--keggclean', type=str, help='KEGG_clean.txt kegg 注释出来的文件')
    parser.add_argument('--degiddir', type=str, default='./DEG_analysis_results',
                        help='默认 DEG_analysis_results 目录，读取所有 ID.txt 文件')
    
    parser.add_argument('--samples', type=str, default='samples_described.txt', help='默认 samples_described.txt')
    parser.add_argument('--compare', type=str, default='compare_info.txt', help='默认 compare_info.txt')
    parser.add_argument('--fpkm', type=str, default='fpkm_matrix_filtered.txt', help='默认 fpkm_matrix_filtered.txt')
    parser.add_argument('--reads', type=str, default='reads_matrix_filtered.txt', help='默认 reads_matrix_filtered.txt')
    parser.add_argument('--filter-type', dest='filter_type', type=str, choices=['pvalue', 'padj'], default='padj', help='默认 padj')
    parser.add_argument('--filter-value', dest='filter_value', type=float, default=0.05, help='默认 0.05')

    args = parser.parse_args()
    
    args.workdir = os.path.abspath(args.workdir)

    return args


def transcriptome_r_deseq(work_dir, fpkm_file, reads_file, samples_file, compare_file, 
                          filter_type, filter_value, deg_value, output_dir):
    os.chdir(work_dir)
    # samples_file = 'samples_described.txt'
    # compare_file = 'compare_info.txt'
    # fpkm_file = 'fpkm_matrix_filtered.txt'
    # reads_file = 'reads_matrix_filtered.txt'
    # logger.info(f"检测组间比较 R 脚本准备文件是否缺失")
    # for f in [samples_file, compare_file, fpkm_file, reads_file]:
    #     if f not in os.listdir():
    #         logger.critical(f'无法执行组间比较脚本，{f} 文件不存在')
    #         sys.exit(1)
            
    logger.info(f'检查 {samples_file} 和 {compare_file} 文件是否符合要求')
    rep = check_sample_comp(samples_file, compare_file)
    if rep == 1:
        logger.critical(f'样本描述文件 {samples_file} 和比较文件 {compare_file} 不符合要求')
        sys.exit(1)
        
    logger.info(f'运行 multiple_samples_DESeq2.r 脚本中，fc 值为 {deg_value}')
    cmd = f"/opt/biosoft/R-4.2.2/bin/Rscript /home/colddata/qinqiang/script/transcriptome/multiple_samples_DESeq2.r \
        --degvalue {deg_value} \
        --fpkm {fpkm_file} \
        --reads {reads_file} \
        --samples {samples_file} \
        --compare {compare_file} \
        --filtertype {filter_type} \
        --filtervalue {filter_value} \
        --degvalue {deg_value} \
        --outputdir {output_dir}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"multiple_samples_DESeq2 程序运行失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        return True


# def transcriptome_enrich(work_dir, gene_go_file, kegg_clean_file, data_dir, compare_file):
#     os.chdir(work_dir)
#     logger.info(f'运行 enrich.r 脚本中...')
#     cmd = f"/opt/biosoft/R-4.2.2/bin/Rscript /home/colddata/qinqiang/script/transcriptome/enrich.r \
#         --genego {gene_go_file} \
#         --keggclean {kegg_clean_file} \
#         --compare {compare_file}"
#     ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     if ret.returncode != 0:
#         logger.error(f"enrich.r 程序运行失败")
#         logger.error(f"标准输出：{ret.stdout.decode()}")
#         logger.error(f"标准错误: {ret.stderr.decode()}")
#         return False
#     else:
#         return True


def transcriptome_enrich_distribution(work_dir):
    os.chdir(work_dir)
    logger.info(f'运行 GO_enrich_distribution 脚本中')
    cmd = f"/opt/biosoft/R-4.2.2/bin/Rscript /home/colddata/qinqiang/script/transcriptome/GO_enrich_distribution.r"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        logger.error(f"GO_enrich_distribution.r 程序运行失败")
        logger.error(f"标准输出：{ret.stdout.decode()}")
        logger.error(f"标准错误: {ret.stderr.decode()}")
        return False
    else:
        return True
        

def process_deresults(de_results_file, kegg_file, nr_file, swiss_file, kns_file):
    if de_results_file.endswith('xlsx'):
        de_df = pd.read_excel(de_results_file, engine='openpyxl', dtype=str)
        de_reads_df = pd.read_csv(de_results_file.replace('.xlsx', '') + '_readCounts.matrix', sep='\t', dtype=str)
    else:
        de_df = pd.read_csv(de_results_file, sep='\t', dtype={"GeneID": str})
        de_reads_df = pd.read_csv(de_results_file + '_readCounts.matrix', sep='\t', dtype={"GeneID": str})
    de_reads_df.columns = [de_reads_df.columns.tolist()[0]] + [x + '_raw_reads' for x in de_reads_df.columns.tolist()[1:]]
    de_df = pd.merge(left=de_df, right=de_reads_df, on=de_reads_df.columns.tolist()[0], how='left')
    # 排序 down，up，NOsig。down 的 FC 值从小到大，up 的 FC 值从大到小
    # 先分三份，再合并
    de_df_down = de_df[de_df['regulation'] == 'Down'].copy()
    de_df_down.sort_values(by='FC', ascending=True, inplace=True)
    de_df_up = de_df[de_df['regulation'] == 'Up'].copy()
    de_df_up.sort_values(by='FC', ascending=False, inplace=True)
    de_df_nosig = de_df[de_df['regulation'] == 'NoSignificant'].copy()
    # 合并
    de_df = pd.concat([de_df_down, de_df_up, de_df_nosig])
    if kns_file or kegg_file or nr_file or swiss_file:
        # 添加注释
        de_df = add_kns_def(de_df, kegg_file, nr_file, swiss_file, kns_file)
    return de_df
    


def main():
    args = parse_input()
    
    # 运行 R 脚本
    if args.degvalue:
        rep = transcriptome_r_deseq(args.workdir, args.fpkm, args.reads, args.samples, args.compare, 
                                    args.filter_type, args.filter_value, args.degvalue, args.workdir)
        if not rep:
            logger.critical(f'R 脚本运行失败')
            sys.exit(1)
    if args.degvalue and args.genego and args.keggclean:
        enrich_dir = os.path.join(args.workdir, 'Pathway_enrichment_analysis')
        transcriptome_enrich(
            compare = args.compare,
            degdata_dir = args.degiddir,
            genego_file = args.genego,
            keggclean_file = args.keggclean,
            outputdir = enrich_dir
        )
        transcriptome_enrich_distribution(enrich_dir)
        os.chdir(args.workdir)
    
    de_results_output_dir = os.path.join(args.workdir, 'DEG_analysis_results/Expression_data')

    for de_results_file in os.listdir(args.workdir):
        if de_results_file.endswith('DE_results') or de_results_file.endswith('DE_results.xlsx'):
            logger.info(f'processing---{de_results_file}')
            de_df = process_deresults(de_results_file, args.kegg, args.nr, args.swiss, args.kns)
            deg_filename = os.path.basename(de_results_file).replace('DE_results', 'DEG_data.txt')
            deg_data_file = os.path.join(de_results_output_dir, deg_filename)
            de_df.to_csv(deg_data_file, sep='\t', index=False)
        

    for up_down_id_file in os.listdir("DEG_analysis_results"):
        up_down_id_file = os.path.join("DEG_analysis_results", up_down_id_file)
        if up_down_id_file.endswith('Down_ID.txt') or up_down_id_file.endswith('Up_ID.txt'):
            up_down_id_df = pd.read_csv(up_down_id_file, sep='\t', names=['GeneID'], dtype={"GeneID": str})
            result_df = add_kns_def(up_down_id_df, args.kegg, args.nr, args.swiss, args.kns)
            result_df.to_csv(up_down_id_file.replace('.txt', '_def.txt'), sep='\t', index=False)

    png_file = 'deg_line_plot.jpeg'
    count_summary_file = 'all_deg_count_summary.txt'
    
    if os.path.exists(png_file):
        os.rename(png_file, os.path.join('Expression_data_evaluation', png_file))
    if os.path.exists(count_summary_file):
        os.rename(count_summary_file, os.path.join('Expression_data_evaluation', count_summary_file))

    logger.success('Done!')


if __name__ == '__main__':
    main()


