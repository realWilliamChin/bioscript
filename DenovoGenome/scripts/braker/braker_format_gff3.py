#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2024/07/01 16:34
# Author        : William GoGo
# Description   : 处理 braker 注释出来的 gff3 文件
# Version       : 1.0
import os, sys
import argparse
import pandas as pd
import openpyxl
import subprocess
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df


def parse_input():
    parser = argparse.ArgumentParser(description="处理 braker 注释出来的 gff3 文件")
    parser.add_argument("-i", "--input_file", required=True, help="输入 braker 注释出来的 gff3 文件")
    parser.add_argument("-o", "--output_file", required=True, help="输出新的 gff3 文件")
    parser.add_argument('-p', '--prefix', required=True, help='修改的前缀')
    args = parser.parse_args()
    return args


def changeid_gff3(gff3_file, output_file, prefix):
    if prefix is None or str(prefix) == "":
        raise ValueError("prefix 不能为空，请使用 -p 指定非空前缀，例如 -p Avena_sativa_")

    perl_script = '/opt/biosoft/geta-2.4.3/bin/gff3_clear.pl'
    cmd = [perl_script, '--prefix', str(prefix), str(gff3_file)]

    # 直接将标准输出写入目标文件，避免使用 shell 重定向
    with open(output_file, 'w') as out_f:
        rep = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if rep.returncode != 0:
        logger.error(f"命令执行失败：{' '.join(cmd)}")
        logger.error(f'标准错误：{rep.stderr}')
        raise RuntimeError('gff3_clear.pl 运行失败，请检查上面的错误信息')

    # 校验输出文件是否为空
    try:
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            logger.warning(f'Perl 输出为空：{output_file}，尝试使用 Python 兼容方案进行转换')
            return False
    except OSError as e:
        logger.error(f'检查输出文件失败：{e}')
        return False

    logger.success('gff3_clear.pl 运行完成！')
    return True


def python_clear_gff3(gff3_file, output_file, prefix):
    """兼容处理：当 gff3_clear.pl 不输出时，基于 mRNA/feature 重构 gene 并重命名。

    逻辑：
    - 读取 GFF3 全部行，按第 3 列类型区分
    - 从 mRNA 的 Parent 作为 gene 的 ID，综合该 gene 下所有 feature 的坐标、链向，构造 gene 行
    - 生成新的 ID：gene 为 <prefix><num>，mRNA/特征追加 .tX / .exonY / .CDSZ 等
    - 输出 gene、mRNA、features，保持坐标排序
    """
    df = pd.read_csv(gff3_file, sep='\t', header=None, comment='#', dtype=str, na_filter=False)
    if df.shape[1] < 9:
        raise RuntimeError('输入 GFF3 列数不足 9，无法解析')

    df = df.iloc[:, :9]
    df.columns = ['seqid','source','type','start','end','score','strand','phase','attr']

    mrna_df = df[df['type'] == 'mRNA'].copy()
    feat_df = df[df['type'] != 'mRNA'].copy()

    # 解析 ID/Parent，安全处理可能的None值
    mrna_df['mRNA_id'] = mrna_df['attr'].astype(str).str.extract(r'ID=([^;\s]+)')[0].fillna('')
    mrna_df['gene_sym'] = mrna_df['attr'].astype(str).str.extract(r'Parent=([^;\s]+)')[0].fillna('')
    
    feat_df['parent'] = feat_df['attr'].astype(str).str.extract(r'Parent=([^;\s]+)')[0].fillna('')

    # 基于 gene_sym 分组，构造 gene 层级
    gene_list = []
    out_lines = []
    gene_num = 0

    # 建立索引: mRNA -> features
    mrna_to_feats = {}
    for m_id in mrna_df['mRNA_id'].dropna().unique():
        mrna_to_feats[m_id] = feat_df[feat_df['parent'] == m_id].copy()

    for gene_sym, sub_mrna in mrna_df.groupby('gene_sym'):
        # gene 级别坐标由其 mRNA 的坐标范围决定
        if sub_mrna.empty:
            continue
        gene_num += 1
        gene_name = f"{prefix}{str(gene_num).zfill( len(str(mrna_df['gene_sym'].nunique())) )}"

        # 汇总 gene 坐标与链向
        # 安全地转换数据类型，处理可能的非数字值
        try:
            start_min = pd.to_numeric(sub_mrna['start'], errors='coerce').dropna().astype(int).min()
            end_max = pd.to_numeric(sub_mrna['end'], errors='coerce').dropna().astype(int).max()
        except (ValueError, TypeError) as e:
            logger.error(f"坐标转换失败: {e}")
            logger.error(f"start列示例: {sub_mrna['start'].head()}")
            logger.error(f"end列示例: {sub_mrna['end'].head()}")
            raise RuntimeError(f"无法处理坐标数据，请检查GFF3文件格式")
        
        strand = sub_mrna['strand'].iloc[0] if not sub_mrna.empty else '.'
        seqid = sub_mrna['seqid'].iloc[0] if not sub_mrna.empty else 'unknown'
        source = sub_mrna['source'].iloc[0] if not sub_mrna.empty else 'unknown'
        score = sub_mrna['score'].iloc[0] if not sub_mrna.empty else '.'
        phase = sub_mrna['phase'].iloc[0] if not sub_mrna.empty else '.'

        # gene 行
        gene_attr = f"ID={gene_name};Name={gene_name}"
        out_lines.append([seqid, source, 'gene', start_min, end_max, score, strand, phase, gene_attr])

        # mRNA 及其 features
        # 安全排序，处理可能的非数字值
        try:
            sorted_sub_mrna = sub_mrna.sort_values(['start','end'])
        except Exception as e:
            logger.warning(f"mRNA排序失败，使用原始顺序: {e}")
            sorted_sub_mrna = sub_mrna
        
        for t_idx, (_, mrow) in enumerate(sorted_sub_mrna.iterrows(), start=1):
            mrna_attr_old = mrow['attr']
            mrna_attr = []
            # 安全处理属性字符串，处理可能的None或非字符串值
            if mrna_attr_old and isinstance(mrna_attr_old, str):
                # 从原有属性中剔除 ID/Name/Parent
                for part in [p for p in mrna_attr_old.replace(';', '; ').split(';') if p.strip()]:
                    if part.startswith(('ID=','Name=','Parent=')):
                        continue
                    mrna_attr.append(part.strip())
            else:
                logger.warning(f"跳过无效的mRNA属性: {mrna_attr_old}")
                mrna_attr = []
            mrna_attr = ';'.join([f"ID={gene_name}.t{t_idx}", f"Parent={gene_name}"] + mrna_attr).rstrip(';')
            out_lines.append([mrow['seqid'], mrow['source'], 'mRNA', mrow['start'], mrow['end'], mrow['score'], mrow['strand'], mrow['phase'], mrna_attr])

            # features
            feats = mrna_to_feats.get(mrow['mRNA_id'], pd.DataFrame())
            if feats.empty:
                continue
            exon_num = cds_num = utr5_num = utr3_num = 0
            # 排序：根据链向选择不同排序，尽量贴合 Perl 实现
            try:
                if mrow['strand'] == '+':
                    feats = feats.sort_values(by=['start','end','type'], key=lambda s: s.map(lambda x: str(x)))
                else:
                    feats = feats.sort_values(by=['end','start','type'], ascending=[False, False, False], key=lambda s: s.map(lambda x: str(x)))
            except Exception as e:
                logger.warning(f"排序失败，使用默认排序: {e}")
                feats = feats.sort_values(by=['start','end'])
            for _, frow in feats.iterrows():
                last_attr = []
                # 安全处理特征属性字符串
                if frow['attr'] and isinstance(frow['attr'], str):
                    for part in [p for p in frow['attr'].replace(';', '; ').split(';') if p.strip()]:
                        if part.startswith(('ID=','Name=','Parent=')):
                            continue
                        last_attr.append(part.strip())
                else:
                    logger.warning(f"跳过无效的特征属性: {frow['attr']}")
                    last_attr = []
                if frow['type'] == 'five_prime_UTR':
                    utr5_num += 1
                    new_id = f"ID={gene_name}.t{t_idx}.utr5p{utr5_num}"
                elif frow['type'] == 'three_prime_UTR':
                    utr3_num += 1
                    new_id = f"ID={gene_name}.t{t_idx}.utr3p{utr3_num}"
                elif frow['type'] == 'exon':
                    exon_num += 1
                    new_id = f"ID={gene_name}.t{t_idx}.exon{exon_num}"
                elif frow['type'] == 'CDS':
                    cds_num += 1
                    new_id = f"ID={gene_name}.t{t_idx}.CDS{cds_num}"
                else:
                    new_id = f"ID={gene_name}.t{t_idx}.{frow['type']}"
                attr = ';'.join([new_id, f"Parent={gene_name}.t{t_idx}"] + (last_attr if last_attr else [])).rstrip(';')
                out_lines.append([frow['seqid'], frow['source'], frow['type'], frow['start'], frow['end'], frow['score'], frow['strand'], frow['phase'], attr])

    # 输出
    out_df = pd.DataFrame(out_lines, columns=['seqid','source','type','start','end','score','strand','phase','attr'])
    # 类型顺序 gene -> mRNA -> 其他
    type_order = {'gene': 0, 'mRNA': 1}
    out_df['type_order'] = out_df['type'].map(lambda t: type_order.get(t, 2))
    
    # 安全排序，处理可能的非数字值
    try:
        out_df = out_df.sort_values(by=['seqid','start','end','type_order'])
    except Exception as e:
        logger.warning(f"最终排序失败，使用默认排序: {e}")
        out_df = out_df.sort_values(by=['seqid','type_order'])
    
    out_df = out_df.drop(columns=['type_order'])
    out_df.to_csv(output_file, sep='\t', header=False, index=False)
    logger.success(f'Python 兼容方案完成，已写入：{output_file}')


def parse_braker3_gff3(gff3_file, output_gff3_file, prefix):
    source_gff_df = load_table(gff3_file, header=None, usecols=[2, 8])
    source_gff_df.columns = ['Type', 'GeneID']
    logger.debug(f'读取原始 gff {source_gff_df}')
    source_gff_df = source_gff_df[source_gff_df['Type'] == 'CDS']
    source_gff_df['GeneID'] = source_gff_df['GeneID'].astype(str).str.split('Parent=').str[1].str.split(';').str[0].fillna('')
    source_gff_df = source_gff_df.drop_duplicates().reset_index()
    logger.debug(f'过滤之后 gff {source_gff_df}')
    source_gff_df = source_gff_df['GeneID']
    
    ok = changeid_gff3(gff3_file, output_gff3_file, prefix)
    if not ok:
        python_clear_gff3(gff3_file, output_gff3_file, prefix)
    
    # 如果 Perl 输出是空文件，此处会报 EmptyDataError，提前在 changeid_gff3 中已做校验
    gff_df = load_table(output_gff3_file, header=None, usecols=[2, 8])
    gff_df.columns = ['Type', 'TargetGeneID']
    gff_df = gff_df[gff_df['Type'] == 'CDS']
    gff_df['TargetGeneID'] = gff_df['TargetGeneID'].astype(str).str.split('Parent=').str[1].str.split(';').str[0].fillna('')
    gff_df = gff_df.drop_duplicates().reset_index()
    print(gff_df)
    gff_df = gff_df['TargetGeneID']
    
    gene_id_replace_df = pd.concat([source_gff_df, gff_df], axis=1)
    gene_id_replace_df.fillna('NA', inplace=True)
    gene_id_replace_df.to_csv('gene_id_change.txt', sep='\t', index=False)


def main():
    args = parse_input()
    parse_braker3_gff3(args.input_file, args.output_file, args.prefix)


if __name__ == '__main__':
    main()