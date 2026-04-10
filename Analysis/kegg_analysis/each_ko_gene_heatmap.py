import os, sys
import time
import pandas as pd
import argparse
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/Rscript/')
from Rscript import smart_heatmap
sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
from load_input import load_table, write_output_df
from data_check import df_drop_row_sum_eq_zero, str_replace_illegal_folder_chars

def parse_input():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help='KEGG_pathway_ID 的一个 list 文件')
    p.add_argument('-k', '--keggclean', help='KEGG 注释解析出来的 KEGG_clean 文件')
    p.add_argument('-f', '--fpkmmatrix', default='fpkm_matrix_filtered.txt',
                   help='fpkm_matrix')
    p.add_argument('-s', '--samples', default='samples_described.txt',
                   help='samples_described.txt 样本描述文件')
    p.add_argument('--genesymbol', help='使用 GeneSymbol 作为索引画热图，人和大鼠小鼠通常使用，输入包含 GeneID GeneSymbol 两列的文件')
    p.add_argument('-o', '--outputdir', default='00_Each_Target_KEGG_pathway_heatmap',
                   help='所有 heatmap 输出文件夹，默认 00_Each_Target_KEGG_pathway_heatmap')
    
    args = p.parse_args()
    
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
    
    return args


def each_ko_gene_heatmap(input_ko_df, kegg_clean_df, fpkm_matrix_df, samples_df, output_dir='00_Each_Target_KEGG_pathway_heatmap', genesymbol_file=None):
    """针对一些 KEGG_ID 的相关基因，每一个 KEGG_ID 画出一个 heatmap 图

    Args:
        input_ko_df (DataFrame): 输入的 KEGG_pathway_ID 的 DataFrame，包含 KEGG_pathway_ID 和 KEGG_pathway_def 两列
        kegg_clean_file (str): kegg 注释出来的 KEGG_clean.txt，只使用第一和第二列，GeneID 和 KEGG_pathway
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
    
    # 如果提供了 genesymbol 参数，读取 GeneID 和 GeneSymbol 映射文件
    genesymbol_df = None
    if genesymbol_file:
        genesymbol_df = load_table(genesymbol_file)
        # 只保留 GeneID 和 GeneSymbol 两列
        genesymbol_df = genesymbol_df[['GeneID', 'GeneSymbol']].copy()
        # 删除空值
        genesymbol_df.dropna(subset=['GeneID', 'GeneSymbol'], inplace=True)
        logger.info(f'读取到 {genesymbol_df.shape[0]} 个 GeneID-GeneSymbol 映射')
    
    # kegg_pathway_df = kegg_clean_df['GeneID', 'KEGG_pathway']
    kegg_clean_df['KEGG_pathway_ID'] = kegg_clean_df['KEGG_pathway'].str.split(':').str[0]
    kegg_clean_df = kegg_clean_df.drop(columns=['KEGG_pathway'])
    kegg_clean_df = kegg_clean_df.drop_duplicates(subset=['GeneID'])
        
    # kegg_id_list = gene_df['KEGG_ID'].values.tolist()
    for _, row in input_ko_df.iterrows():
        each_kegg_id = row['KEGG_pathway_ID']
        safe_kegg_id_def = each_kegg_id.replace(':', '_') + '_' + str_replace_illegal_folder_chars(str(row['KEGG_pathway_def']))
        
        # kegg_pic_dir = os.path.join(crt_kegg_id_dir, 'KEGG_pathway_heatmap')
        each_kegg_id_df = kegg_clean_df[kegg_clean_df['KEGG_pathway_ID'] == each_kegg_id]
        logger.info(f'尝试对 {safe_kegg_id_def} 的相关基因画 heatmap 图，数量为 {each_kegg_id_df.shape[0]}')

        # 添加 fpkm
        each_kegg_id_gene_fpkm_df = pd.merge(each_kegg_id_df, fpkm_matrix_df, on='GeneID', how='left')
        each_kegg_id_gene_fpkm_df.dropna(how='any', axis=0, inplace=True) # 去除掉为空的行
        each_kegg_id_gene_fpkm_df = df_drop_row_sum_eq_zero(each_kegg_id_gene_fpkm_df)
        if each_kegg_id_gene_fpkm_df.shape[0] == 0:
            logger.error(f'{safe_kegg_id_def} 相关基因表达量为空，跳过执行 multigroup heatmap')
            continue
        
        # kegg 相关的 id 小于 3 个就跳过 (2024_06_14:张老师：从 10 改为 3)
        if each_kegg_id_gene_fpkm_df.shape[0] < 3:
            logger.warning(f'{safe_kegg_id_def} 相关基因数量小于 3 个，不对此画 heatmap 图')
            continue
        
        # 如果提供了 genesymbol_df，则替换 GeneID 为 GeneSymbol
        if genesymbol_df is not None:
            # 合并 GeneSymbol 映射
            each_kegg_id_gene_fpkm_df = pd.merge(each_kegg_id_gene_fpkm_df, genesymbol_df, on='GeneID', how='inner')
            # 将 GeneID 列替换为 GeneSymbol
            each_kegg_id_gene_fpkm_df['GeneID'] = each_kegg_id_gene_fpkm_df['GeneSymbol']
            each_kegg_id_gene_fpkm_df.drop(columns=['GeneSymbol'], inplace=True)
            each_kegg_id_gene_fpkm_df.drop_duplicates(subset=['GeneID'], inplace=True)  # 换成 GeneSymbol 可能会出现重复，进行去重处理，可能后续需要修改去重方式
            if each_kegg_id_gene_fpkm_df.shape[0] < 3:
                logger.warning(f'{safe_kegg_id_def} 相关 genesymbol 数量小于 3，跳过画 Heatmap 图')
                continue
        
        # multigroup heatmap 输入文件
        kegg_pathway_id_gene_heatmap_filename = os.path.join(output_dir, f'{safe_kegg_id_def}_gene_heatmap.xlsx')
        with pd.ExcelWriter(kegg_pathway_id_gene_heatmap_filename, engine='openpyxl') as writer:
            each_kegg_id_gene_fpkm_df[samples_list].to_excel(writer, sheet_name='Sheet1', index=False)
            samples_df.to_excel(writer, sheet_name='Sheet2', index=False)
    
        # heatmap 画图
        kegg_pathway_id_heatmap_filename = os.path.join(output_dir, os.path.basename(kegg_pathway_id_gene_heatmap_filename).replace('.xlsx', '.jpeg'))
        heatmap_result = smart_heatmap(
            kegg_pathway_id_gene_heatmap_filename,
            kegg_pathway_id_heatmap_filename,
            annotation_col=2,
            cluster_rows=True,
            cluster_cols=False,
            scale="row"
        )
        if not heatmap_result:
            logger.error(f'{safe_kegg_id_def} draw_multigroup_heatmap 程序失败')


def main():
    args = parse_input()
    
    ko_df = load_table(args.input, dtype=str)
    
    kegg_clean_df = load_table(args.keggclean, usecols=[0, 1], dtype=str,
                               header=None, names=['GeneID', 'KEGG_pathway'])
    fpkm_matrix_df = load_table(args.fpkmmatrix, dtype={'GeneID': str})
    if len(fpkm_matrix_df.columns) == 2:
        logger.error('fpkm_matrix 文件只有两列，无法进行 heatmap 绘制')
        sys.exit(1)
    samples_df = load_table(args.samples, dtype=str, usecols=[0, 1])
    
    each_ko_gene_heatmap(ko_df, kegg_clean_df, fpkm_matrix_df, samples_df, args.outputdir, args.genesymbol)


if __name__ == '__main__':
    time1 = time.time()
    main()
    time2 = time.time()
    logger.info('Time used: ' + str(time2-time1))