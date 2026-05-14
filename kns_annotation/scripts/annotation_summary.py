#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import argparse
import os
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(description='生成注释汇总统计信息')
    parser.add_argument('anno_fasta', help='注释的fasta文件')
    parser.add_argument('--swiss_goid', default=None, help='Swiss GO ID文件，默认自动查找swiss/*swiss_idNo_def.txt')
    parser.add_argument('--swiss_gene_def', default=None, help='Swiss基因定义文件，默认自动查找swiss/*swiss_gene_def.txt')
    parser.add_argument('--kegg_gene_def', default=None, help='KEGG基因定义文件，默认自动查找kegg/*KEGG_gene_def.txt')
    parser.add_argument('--nr_gene_def', default=None, help='NR基因定义文件，默认自动查找nr/*nr_gene_def.txt')
    parser.add_argument('--output_pic', default='annotation_summary.jpeg',help='输出图片路径')
    parser.add_argument('--output_summary_file', default='annotation_summary.txt', help='输出汇总文件路径')
    
    args = parser.parse_args()

    # 自动查找swiss_goid和swiss_gene_def的默认路径
    if args.swiss_goid is None:
        swiss_goid_candidates = glob.glob('swiss/*swiss_idNo_def.txt')
        if len(swiss_goid_candidates) == 1:
            args.swiss_goid = swiss_goid_candidates[0]
        elif not swiss_goid_candidates:
            parser.error("未找到默认的 swiss/*swiss_idNo_def.txt 文件，请使用 --swiss_goid 指定。")
        else:
            parser.error("有多个匹配的swiss/*swiss_idNo_def.txt文件，请使用 --swiss_goid 指定一个。")
    if args.swiss_gene_def is None:
        swiss_gene_def_candidates = glob.glob('swiss/*swiss_gene_def.txt')
        if len(swiss_gene_def_candidates) == 1:
            args.swiss_gene_def = swiss_gene_def_candidates[0]
        elif not swiss_gene_def_candidates:
            parser.error("未找到默认的 swiss/*swiss_gene_def.txt 文件，请使用 --swiss_gene_def 指定。")
        else:
            parser.error("有多个匹配的swiss/*swiss_gene_def.txt文件，请使用 --swiss_gene_def 指定一个。")
    if args.kegg_gene_def is None:
        kegg_gene_def_candidates = glob.glob('kegg/*KEGG_gene_def.txt')
        if len(kegg_gene_def_candidates) == 1:
            args.kegg_gene_def = kegg_gene_def_candidates[0]
        elif not kegg_gene_def_candidates:
            parser.error("未找到默认的 kegg/*KEGG_gene_def.txt 文件，请使用 --kegg_gene_def 指定。")
        else:
            parser.error("有多个匹配的kegg/*KEGG_gene_def.txt文件，请使用 --kegg_gene_def 指定一个。")
    if args.nr_gene_def is None:
        nr_gene_def_candidates = glob.glob('nr/*nr_gene_def.txt')
        if len(nr_gene_def_candidates) == 1:
            args.nr_gene_def = nr_gene_def_candidates[0]
        elif not nr_gene_def_candidates:
            parser.error("未找到默认的 nr/*nr_gene_def.txt 文件，请使用 --nr_gene_def 指定。")
        else:
            parser.error("有多个匹配的nr/*nr_gene_def.txt文件，请使用 --nr_gene_def 指定一个。")
    return args


def extract_first_column(input_file, output_file):
    """从输入文件中提取第一列，跳过第一行（与bash的tail -n +2 | cut -f 1行为一致）"""
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        lines = f_in.readlines()
        # 跳过第一行（标题行）
        for line in lines[1:]:
            # 提取第一列（以制表符分割，保持与cut -f 1一致的行为）
            # 注意：不strip整行，因为cut会保留原始行的分割
            parts = line.rstrip('\n\r').split('\t')
            first_col = parts[0] if parts else ''
            f_out.write(first_col + '\n')


def count_fasta_sequences(fasta_file):
    """统计fasta文件中的序列数量（统计'>'的数量）"""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def count_file_lines(file_path):
    """统计文件行数（与wc -l行为一致，包括空行）"""
    count = 0
    with open(file_path, 'r') as f:
        for line in f:
            count += 1
    return count


def draw_venn_diagram(go_id, swiss_id, kegg_id, nr_id, output_pic):
    """调用venn_plot.py绘制Venn图"""
    venn_script = '/home/colddata/qinqiang/script/Plot/Venn/venn_plot.py'
    
    cmd = [
        sys.executable,
        venn_script,
        '-i', go_id, swiss_id, kegg_id, nr_id,
        '--pic-name', output_pic
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"错误: 绘制Venn图失败", file=sys.stderr)
        print(f"标准错误: {result.stderr}", file=sys.stderr)
        print(f"标准输出: {result.stdout}", file=sys.stderr)
        sys.exit(1)


def generate_summary_file(output_file, total_id_count, kegg_count, go_count, nr_count, swiss_count):
    """生成汇总文件"""
    with open(output_file, 'w') as f:
        f.write("Database\tTotal\tAnnotated\tPercentage\n")
        
        # 计算百分比（注意：bash脚本中减1是因为wc -l会包含标题行，但我们已经跳过了标题行）
        # 为了保持与bash脚本完全一致的行为，我们也减1
        kegg_annotated = kegg_count - 1
        go_annotated = go_count - 1
        nr_annotated = nr_count - 1
        swiss_annotated = swiss_count - 1
        
        # 计算百分比（整数除法，与bash脚本一致）
        kegg_percent = (kegg_annotated * 100) // total_id_count if total_id_count > 0 else 0
        go_percent = (go_annotated * 100) // total_id_count if total_id_count > 0 else 0
        nr_percent = (nr_annotated * 100) // total_id_count if total_id_count > 0 else 0
        swiss_percent = (swiss_annotated * 100) // total_id_count if total_id_count > 0 else 0
        
        f.write(f"KEGG\t{total_id_count}\t{kegg_annotated}\t{kegg_percent}%\n")
        f.write(f"GO\t{total_id_count}\t{go_annotated}\t{go_percent}%\n")
        f.write(f"NR\t{total_id_count}\t{nr_annotated}\t{nr_percent}%\n")
        f.write(f"Swiss Protein\t{total_id_count}\t{swiss_annotated}\t{swiss_percent}%\n")


def main():
    args = parse_args()
    
    # 创建临时文件（在当前目录，与bash脚本一致）
    go_id = "GO.txt"
    swiss_id = "Swiss.txt"
    kegg_id = "KEGG.txt"
    nr_id = "NR.txt"
    
    try:
        # 处理输入文件
        extract_first_column(args.swiss_goid, go_id)
        extract_first_column(args.swiss_gene_def, swiss_id)
        extract_first_column(args.kegg_gene_def, kegg_id)
        extract_first_column(args.nr_gene_def, nr_id)
        
        # 绘制Venn图
        draw_venn_diagram(go_id, swiss_id, kegg_id, nr_id, args.output_pic)
        
        # 计算统计信息
        total_id_count = count_fasta_sequences(args.anno_fasta)
        kegg_count = count_file_lines(kegg_id)
        go_count = count_file_lines(go_id)
        nr_count = count_file_lines(nr_id)
        swiss_count = count_file_lines(swiss_id)
        
        # 生成汇总文件
        generate_summary_file(
            args.output_summary_file,
            total_id_count,
            kegg_count,
            go_count,
            nr_count,
            swiss_count
        )
        
    finally:
        # 清理临时文件
        for temp_file in [go_id, swiss_id, kegg_id, nr_id]:
            if os.path.exists(temp_file):
                os.remove(temp_file)


if __name__ == "__main__":
    main()

