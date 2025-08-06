#!/home/train/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# Created Time  : 2023/6/6 20:32
# Author        : WilliamGoGo
import os, sys
import pandas as pd
import argparse
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df
from Blast.blast_drop_dup import deduplicate_blast_results


def parse_input():
    parser = argparse.ArgumentParser(description="biogrid, 使用 cds fasta 文件")
    parser.add_argument('-f', '--fasta', help='输入要比对的 cds_fasta 文件')
    parser.add_argument('-d', '--database', choices=['plant', 'animal', 'animal_human', 'animal_mouse', 'fungus', 'insect'],
                        help='输入比对的数据库, animal --> animal_human')
    parser.add_argument('-c', '--cpu', help='输入线程数，默认 20 线程', default=20)
    parser.add_argument('-p', '--prefix', required=True, help='输入输出文件的前缀')
    
    # 创建自定义参数组
    custom_group = parser.add_argument_group('自定义 biogrid 参考库，仍然需要输入 --fasta --cpu --prefix')
    custom_group.add_argument('--custom-database', dest='custom_database', help='自定义数据库')
    custom_group.add_argument('--custom-ref', dest='custom_ref', help='自定义输入参考文件')
    
    # 只进行解析
    parse_biogrid_group = parser.add_argument_group('直接对 blast 结果解析')
    parse_biogrid_group.add_argument('--biogrid', help='输入 blast 比对结果文件')
    parse_biogrid_group.add_argument('--parse-ref', dest='parse_ref', help='输入参考文件')
    
    args = parser.parse_args()
    
    # 添加已做好的 biogrid 数据库
    biogrid_database = '/home/colddata/qinqiang/02_Biogrid'
    if not args.database:
        pass
    elif args.database.lower() == 'plant':
        specie_name = 'Arabidopsis_thaliana'
        plant_dir = os.path.join(biogrid_database, 'Plant_Biogrid_Arabidopsis_thalinan')
        args.database = os.path.join(plant_dir, '00_Database', specie_name)
        args.parse_ref = os.path.join(plant_dir, 'BIOGRID-Arabidopsis_thaliana-embl.txt')
    elif args.database.lower() in ['animal_human', 'animal']:
        specie_name = 'Homo_sapiens'
        animal_dir = os.path.join(biogrid_database, 'Animal_Biogrid_Homo_sapiens')
        args.database = os.path.join(animal_dir, '00_Database', specie_name)
        args.parse_ref = os.path.join(animal_dir, 'BIOGRID-Homo_sapiens-embl.txt')
    elif args.database.lower() == 'fungus':
        specie_name = 'Saccharomyces_cerevisiae'
        fungus_dir = os.path.join(biogrid_database, 'Fungus_Biogrid_Saccharomyces_cerevisiae')
        args.database = os.path.join(fungus_dir, '00_Database', specie_name)
        args.parse_ref = os.path.join(fungus_dir, 'BIOGRID-Saccharomyces_cerevisiae_embl_2.txt')
    elif args.database.lower() == 'insect':
        specie_name = 'Drosophila_melanogaster'
        drosaphila_dir = os.path.join(biogrid_database, 'Insect_Biogrid_Drosophila_melanogaster')
        args.database = os.path.join(drosaphila_dir, '00_Database', specie_name)
        args.parse_ref = os.path.join(drosaphila_dir, 'BIOGRID-Drosophila_melanogaster-embl.txt')
    elif args.database.lower() == 'animal_mouse':
        specie_name = 'Mus_musculus'
        mus_dir = os.path.join(biogrid_database, 'Animal_Biogrid_Mus_musculus')
        args.database = os.path.join(mus_dir, '00_Database', specie_name)
        args.parse_ref = os.path.join(mus_dir, 'BIOGRID-Mus_musculus-embl.txt')
    else:
        pass

    return args, specie_name


def exec_blast(fasta_file, num_threads, database, prefix):
    """ 执行 blast 比对，生成 blast 结果

    Args:
        fasta_file (str): 输入 cds fasta 文件
        num_threads (int): 线程数量
        database (str): database 路径
        prefix (str): 输出前缀

    Returns:
        str: 结果文件名
    """
    blast_file_name = prefix + '.blast'
    pep_file_name = prefix + '_pep.fasta'
    # cds 转成 pep，后续需要清理一下序列
    os.system(f'seqkit translate {fasta_file} > {pep_file_name}')
    # pep 文件清理 _frame=1
    os.system(f"sed 's/_frame=1//g' -i {pep_file_name}")
    
    blast_command = f'/opt/biosoft/ncbi-blast-2.9.0+/bin/blastp \
        -db {database} \
        -query {pep_file_name} \
        -out {blast_file_name} \
        -evalue 1e-5 \
        -num_threads {num_threads} \
        -outfmt "6 qacc sacc qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle"'
        
    os.system(blast_command)
    
    return blast_file_name


def biogrid(in_blast_file, ref_file, out_file):
    dic = {}
    with open(in_blast_file, 'r') as f:
        num1 = 0
        for each_line1 in f:
            if num1 == 0:
                num1 = 1
                title1 = each_line1
            else:
                value = each_line1.split('\t')[0]
                key = each_line1.strip().split('\t')[1]
                if key not in dic:
                    dic[key] = value
                else:
                    dic[key] += ';' + value

    # Write results to file
    f2 = open(out_file, 'w')
    f2.write(
        'BioGRID Interaction ID\tEmbl_A\tEmbl_B\tThroughput\tSource_Database\n')
    with open(ref_file, 'r') as f1:
        num = 0
        for each_line in f1:
            if num == 0:
                title_part = each_line.split('\t')
                num = 1
            else:
                gene1 = each_line.strip().split('\t')[1]
                gene2 = each_line.strip().split('\t')[2]
                id = each_line.strip().split('\t')[0]
                other = each_line.strip().split('\t')[3]
                source_database = each_line.strip().split('\t')[4]
                if gene1 in dic:
                    if ';' in dic[gene1]:
                        geneAs = dic[gene1].split(';')
                        for each in geneAs:
                            geneA = each
                            if gene2 in dic:
                                if ';' in dic[gene2]:
                                    geneBs = dic[gene2].split(';')
                                    for each in geneBs:
                                        geneB = each
                                        content = id + '\t' + geneA + '\t' + geneB + '\t' + other + '\t' + source_database + '\n'
                                        f2.write(content)
                                else:
                                    geneB = dic[gene2]
                                    content = id + '\t' + geneA + '\t' + geneB + '\t' + other + '\t' + source_database + '\n'
                                    f2.write(content)
                    else:
                        geneA = dic[gene1]
                        if gene2 in dic:
                            if ';' in dic[gene2]:
                                geneBs = dic[gene2].split(';')
                                for each in geneBs:
                                    geneB = each
                                    content = id + '\t' + geneA + '\t' + geneB + '\t' + other + '\t' + source_database + '\n'
                                    f2.write(content)
                            else:
                                geneB = dic[gene2]
                                content = id + '\t' + geneA + '\t' + geneB + '\t' + other + '\t' + source_database + '\n'
#                                 f2.write(content)

# 改成 pandas 大概了解上面代码的逻辑
# def biogrid(in_blast_file, ref_file, out_file):
#     """
#     使用 Pandas 处理 BioGRID 数据，将基因名映射到 EMBL ID
    
#     Args:
#         in_blast_file: BLAST 结果文件路径
#         ref_file: 参考文件路径
#         out_file: 输出文件路径
#     """
#     blast_df = pd.read_csv(in_blast_file, sep='\t', header=0)
#     # 构建基因名到 EMBL ID 的映射字典
#     gene_to_embl = blast_df.groupby(blast_df.columns[1])[blast_df.columns[0]].agg(lambda x: ';'.join(x)).to_dict()
    
#     # 读取参考文件
#     ref_df = pd.read_csv(ref_file, sep='\t', header=0)
    
#     # 创建结果列表
#     results = []
    
#     # 处理每一行数据
#     for _, row in ref_df.iterrows():
#         gene1 = row[1]
#         gene2 = row[2]
#         id = row[0]
#         other = row[3]
#         source_database = row[4]
        
#         # 获取基因1的 EMBL ID
#         embl1_list = gene_to_embl.get(gene1, [])
#         if isinstance(embl1_list, str):
#             embl1_list = embl1_list.split(';')
            
#         # 获取基因2的 EMBL ID
#         embl2_list = gene_to_embl.get(gene2, [])
#         if isinstance(embl2_list, str):
#             embl2_list = embl2_list.split(';')
            
#         # 生成所有可能的组合
#         for embl1 in embl1_list:
#             for embl2 in embl2_list:
#                 results.append({
#                     'BioGRID Interaction ID': id,
#                     'Embl_A': embl1,
#                     'Embl_B': embl2,
#                     'Throughput': other,
#                     'Source_Database': source_database
#                 })
    
#     # 创建结果 DataFrame 并保存
#     result_df = pd.DataFrame(results)
#     result_df.to_csv(out_file, sep='\t', index=False)


def main():
    args, specie_name = parse_input()
    biogrid_output = os.path.join(
        os.path.dirname(args.prefix),
        f'Biogrid_PPI_relation_from_{specie_name}.txt'
    )
    if args.biogrid:
        biogrid(args.biogrid, args.parse_ref, biogrid_output)
        return
    # 比对
    exec_blast(args.fasta, args.cpu, args.database, args.prefix)
    
    # 处理 blast 去重
    blast_names = ['qacc','sacc','qcovhsp','ppos','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','stitle']
    blast_df = load_table(args.prefix + '.blast', header=None, names=blast_names, low_memory=False)
    uniq_blast_df = deduplicate_blast_results(blast_df, 'qacc', 'bitscore', 'max')
    uniq_blast_file = args.prefix + '_uniq.blast'
    write_output_df(uniq_blast_df, uniq_blast_file, index=False)

    # uniq_blast 进行 biogrid
    biogrid(uniq_blast_file, args.parse_ref, biogrid_output)
    
    logger.success(f'Done!')


if __name__ == '__main__':
    main()
