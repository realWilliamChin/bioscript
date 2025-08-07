import os, sys
import time
import pandas as pd
import argparse
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/Rscript/')
from Rscript import draw_multigroup_heatmap
sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df
from data_check import df_drop_row_sum_eq_zero

def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='KEGG_Pathway_ID 的一个 list 文件')
    p.add_argument('-k', '--keggclean', help='KEGG 注释解析出来的 KEGG_clean 文件')
    p.add_argument('-f', '--fpkmmatrix', default='fpkm_matrix_filtered.txt',
                   help='fpkm_matrix')
    p.add_argument('-s', '--samples', default='samples_described.txt',
                   help='samples_described.txt 样本描述文件')
    p.add_argument('-o', '--outputdir', default='./KEGG_pathway_heatmap',
                   help='所有 heatmap 输出文件夹，默认 KEGG_pathway_heatmap')
    
    args = p.parse_args()
    
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
    
    return args


def each_ko_gene_heatmap(kegg_id_list, kegg_clean_df, fpkm_matrix_df, samples_df, output_dir='./'):
    """针对一些 KEGG_ID 的相关基因，每一个 KEGG_ID 画出一个 heatmap 图

    Args:
        kegg_id_list (list): 张老师给的 KEGG_ID list
        kegg_clean_file (str): kegg 注释出来的 KEGG_clean.txt，只使用第一和第二列，GeneID 和 KEGG_Pathway
        fpkm_matrix_df (DataFrame): 第一列是 GeneID，其他列是 smaples_df 的 sample 列
        samples_df (DataFrame): pandas DataFrame 矩阵，包含两列 sample 和 group
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # samples_df = pd.read_csv(samples_file, sep='\t', usecols=[0, 1])
    samples_df = samples_df[['sample', 'group']]
    samples_list = ['GeneID'] + samples_df['sample'].values.tolist()

    # 对 expression fpkm matrix 文件只保留 fpkm 值
    fpkm_matrix_df['GeneID'] = fpkm_matrix_df['GeneID'].astype(str)
    
    # kegg_pathway_df = kegg_clean_df['GeneID', 'KEGG_Pathway']
    kegg_clean_df['KEGG_ID'] = kegg_clean_df['KEGG_Pathway'].str.split(':').str[0]
    kegg_clean_df = kegg_clean_df.drop(columns=['KEGG_Pathway'])
    kegg_clean_df = kegg_clean_df.drop_duplicates(subset=['GeneID'])
        
    # kegg_id_list = gene_df['KEGG_ID'].values.tolist()
    for each_kegg_id in kegg_id_list:
        # kegg_pic_dir = os.path.join(crt_kegg_id_dir, 'KEGG_pathway_heatmap')
        each_kegg_id_df = kegg_clean_df[kegg_clean_df['KEGG_ID'] == each_kegg_id]
        logger.info(f'尝试对 {each_kegg_id} 的相关基因画 heatmap 图，数量为 {each_kegg_id_df.shape[0]}')

        # 添加 fpkm
        each_kegg_id_gene_fpkm_df = pd.merge(each_kegg_id_df, fpkm_matrix_df, on='GeneID', how='left')
        each_kegg_id_gene_fpkm_df.dropna(how='any', axis=0, inplace=True) # 去除掉为空的行
        each_kegg_id_gene_fpkm_df = df_drop_row_sum_eq_zero(each_kegg_id_gene_fpkm_df)
        if each_kegg_id_gene_fpkm_df.shape[0] == 0:
            logger.error(f'{each_kegg_id} 相关基因表达量为空，跳过执行 multigroup heatmap')
            continue
        
        # kegg 相关的 id 小于 3 个就跳过 (2024_06_14:张老师：从 10 改为 3)
        if each_kegg_id_gene_fpkm_df.shape[0] < 3:
            logger.warning(f'{each_kegg_id} 相关基因数量小于 3 个，不对此画 heatmap 图')
            continue
        
        # multigroup heatmap 输入文件
        kegg_id_gene_ko_heatmap_filename = os.path.join(output_dir, f'{each_kegg_id}_ko_gene_heatmap.xlsx')
        with pd.ExcelWriter(kegg_id_gene_ko_heatmap_filename, engine='openpyxl') as writer:
            each_kegg_id_gene_fpkm_df[samples_list].to_excel(writer, sheet_name='Sheet1', index=False)
            samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
    
        # heatmap 画图
        kegg_id_heatmap_filename = os.path.join(output_dir, f"{each_kegg_id}_ko_gene_heatmap.jpeg")
        heatmap_result = draw_multigroup_heatmap(
            kegg_id_gene_ko_heatmap_filename,
            kegg_id_heatmap_filename,
            other_args='--cluster-rows'
            )
        if not heatmap_result:
            logger.error(f'{each_kegg_id} draw_multigroup_heatmap 结果失败，跳过执行 multigroup heatmap')


def main():
    args = parse_input()
    
    ko_df = load_table(args.input, dtype=str)
    ko_list = ko_df['KEGG_Pathway_ID'].values.tolist()
    
    kegg_clean_df = load_table(args.keggclean, usecols=[0, 1], dtype=str,
                               header=None, names=['GeneID', 'KEGG_Pathway'])
    fpkm_matrix_df = load_table(args.fpkmmatrix, dtype={'GeneID': str})
    if len(fpkm_matrix_df.columns) == 2:
        logger.error('fpkm_matrix 文件只有两列，无法进行 heatmap 绘制')
        sys.exit(1)
    samples_df = load_table(args.samples, dtype=str, usecols=[0, 1])
    
    each_ko_gene_heatmap(ko_list, kegg_clean_df, fpkm_matrix_df, samples_df, args.outputdir)


if __name__ == '__main__':
    time1 = time.time()
    main()
    time2 = time.time()
    logger.info('Time used: ' + str(time2-time1))