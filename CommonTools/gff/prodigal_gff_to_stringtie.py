#!/usr/bin/env python3
"""
将 Prodigal GFF 文件转换为 StringTie 兼容的 GFF3 格式
"""

import argparse
import sys
import os
from collections import defaultdict

def parse_prodigal_gff(input_file):
    """解析 Prodigal GFF 文件，返回基因信息列表"""
    genes = []
    
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # 跳过注释行和空行
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            if len(fields) < 9:
                print(f"警告: 第 {line_num} 行字段数不足，跳过", file=sys.stderr)
                continue
                
            seqname, source, feature, start, end, score, strand, frame, attributes = fields
            
            # 只处理 CDS 特征
            if feature != 'CDS':
                continue
                
            # 解析属性字段
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key.strip()] = value.strip()
            
            gene_id = attr_dict.get('ID')
            if not gene_id:
                print(f"警告: 第 {line_num} 行缺少 ID 属性，跳过", file=sys.stderr)
                continue
                
            genes.append({
                'seqname': seqname,
                'source': source,
                'start': int(start),
                'end': int(end),
                'score': score,
                'strand': strand,
                'frame': frame,
                'gene_id': gene_id,
                'attributes': attr_dict
            })
    
    return genes

def convert_to_gff3(genes, output_file):
    """将基因信息转换为 GFF3 格式并写入文件"""
    with open(output_file, 'w') as f:
        # 写入 GFF3 文件头
        f.write('##gff-version 3\n')
        f.write(f'# Converted from Prodigal GFF to StringTie-compatible GFF3\n')
        f.write(f'# Total genes: {len(genes)}\n')
        
        for gene_info in genes:
            seqname = gene_info['seqname']
            source = 'Prodigal'  # 统一 source
            start = gene_info['start']
            end = gene_info['end']
            score = gene_info['score']
            strand = gene_info['strand']
            frame = gene_info['frame']
            gene_id = gene_info['gene_id']
            transcript_id = f"{gene_id}_t1"
            
            # 1. gene 条目
            gene_attrs = f"ID=gene_{gene_id};Name={gene_id};biotype=protein_coding"
            f.write(f'{seqname}\t{source}\tgene\t{start}\t{end}\t{score}\t{strand}\t.\t{gene_attrs}\n')
            
            # 2. mRNA 条目
            mrna_attrs = f"ID={transcript_id};Parent=gene_{gene_id};Name={transcript_id};biotype=protein_coding"
            f.write(f'{seqname}\t{source}\tmRNA\t{start}\t{end}\t{score}\t{strand}\t.\t{mrna_attrs}\n')
            
            # 3. exon 条目
            exon_attrs = f"ID=exon_{transcript_id}_1;Parent={transcript_id};Name=exon_{transcript_id}_1"
            f.write(f'{seqname}\t{source}\texon\t{start}\t{end}\t{score}\t{strand}\t.\t{exon_attrs}\n')
            
            # 4. CDS 条目（保留原始 frame 信息）
            cds_attrs = f"ID=cds_{transcript_id};Parent={transcript_id};Name=cds_{transcript_id}"
            f.write(f'{seqname}\t{source}\tCDS\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{cds_attrs}\n')


def check_sequence_names(gff_file, bam_file):
    """检查 GFF 和 BAM 文件的序列名是否匹配"""
    import subprocess
    
    print("\n=== 检查序列名匹配 ===")
    
    # 提取 GFF 中的序列名
    gff_seqs = set()
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 1:
                gff_seqs.add(fields[0])
    
    # 提取 BAM 中的序列名
    bam_seqs = set()
    try:
        result = subprocess.run(['samtools', 'view', '-H', bam_file], 
                              capture_output=True, text=True, check=True)
        for line in result.stdout.split('\n'):
            if line.startswith('@SQ'):
                parts = line.split('\t')
                for part in parts:
                    if part.startswith('SN:'):
                        bam_seqs.add(part[3:])
                        break
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"无法读取 BAM 文件头: {e}", file=sys.stderr)
        return
    
    print(f"GFF 序列数量: {len(gff_seqs)}")
    print(f"BAM 序列数量: {len(bam_seqs)}")
    
    common_seqs = gff_seqs & bam_seqs
    print(f"共同序列数量: {len(common_seqs)}")
    
    if len(common_seqs) == 0:
        print("错误: 没有共同的序列名!", file=sys.stderr)
        print("GFF 序列名示例:", list(gff_seqs)[:5])
        print("BAM 序列名示例:", list(bam_seqs)[:5])
    else:
        print("序列名匹配检查通过")

def main():
    parser = argparse.ArgumentParser(description='将 Prodigal GFF 转换为 StringTie 兼容的 GFF3 格式')
    parser.add_argument('-i', '--input', required=True, help='输入 Prodigal GFF 文件')
    parser.add_argument('-o', '--output', required=True, help='输出 GFF3 文件')
    parser.add_argument('-b', '--bam', help='BAM 文件（用于检查序列名匹配）')
    parser.add_argument('--check-only', action='store_true', help='仅检查序列名匹配，不转换')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"错误: 输入文件不存在: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if not args.check_only:
        print(f"正在转换 {args.input} -> {args.output}")
        
        # 解析和转换
        genes = parse_prodigal_gff(args.input)
        print(f"找到 {len(genes)} 个基因")
        
        if len(genes) == 0:
            print("错误: 没有找到有效的基因记录", file=sys.stderr)
            sys.exit(1)
            
        convert_to_gff3(genes, args.output)
        print(f"转换完成: {args.output}")
    
    # 检查序列名匹配
    if args.bam:
        if args.check_only and os.path.exists(args.output):
            check_sequence_names(args.output, args.bam)
        elif not args.check_only:
            check_sequence_names(args.output, args.bam)
    elif args.check_only:
        print("错误: --check-only 需要指定 --bam 参数", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()