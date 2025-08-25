# -*- coding: utf-8 -*-
# Created Time  : 2024/11/04 15:47
# Author        : William GoGo
import os, sys
import pandas as pd
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/CommonTools'))
from Fasta.get_sequence_from_list import get_seq_from_idlist


def process_braker_output(gff3_file, braker_new_gff, braker_aa, braker_cds, braker_new_aa, braker_new_cds, output_dir="./"):
    """处理 braker 注释出来的 gff3 文件，对 gff3 文件进行更名
    对 cds 和 蛋白 fasta 的序列名称也进行更改

    Args:
        gff3_file (str): 通常是 braker3 注释出来的 gff3 文件
        prefix (前缀): 通常是 specie 例如 Gallus_gallus

    Returns:
        tuple: output_gff3_file, output_aa_fasta, output_cds_fasta, geneid_replace_fname
    """
    # cmd = f"/opt/biosoft/geta-2.4.3/bin/gff3_clear.pl --prefix {prefix}_ {gff3_file} > aa; mv aa {output_gff3_file}"
    # run_cmd(cmd, "gff3 更名")

    source_gff_df = pd.read_csv(gff3_file, sep="\t", header=None, usecols=[2, 8])
    source_gff_df.columns = ["Type", "GeneID"]
    print(source_gff_df)
    source_gff_df = source_gff_df[source_gff_df["Type"] == "CDS"]
    source_gff_df["GeneID"] = (
        source_gff_df["GeneID"].str.split("Parent=").str[1].str.split(";").str[0]
    )
    source_gff_df = source_gff_df.drop_duplicates().reset_index()
    print(source_gff_df)
    source_gff_df = source_gff_df["GeneID"]

    gff3_columns = [
        "seqid", "source", "type", "start", "end",
        "score", "strand", "phase", "attribute",
    ]
    gff_df = pd.read_csv(braker_new_gff, sep="\t", header=None, names=gff3_columns)
    gff_df.to_csv(braker_new_gff, sep="\t", index=False)
    # gff_df = pd.read_csv(output_gff3_file, sep='\t', header=None, usecols=[2, 8])
    gff_df = gff_df.drop(
        columns=["seqid", "source", "start", "end", "score", "strand", "phase"]
    )
    gff_df.columns = ["Type", "TargetGeneID"]
    gff_df = gff_df[gff_df["Type"] == "CDS"]
    gff_df["TargetGeneID"] = (
        gff_df["TargetGeneID"].str.split("Parent=").str[1].str.split(";").str[0]
    )
    gff_df = gff_df.drop_duplicates().reset_index()
    gff_df = gff_df["TargetGeneID"]

    gene_id_replace_df = pd.concat([source_gff_df, gff_df], axis=1)
    gene_id_replace_df.fillna("NA", inplace=True)

    geneid_replace_fname = os.path.join(output_dir, "GeneIDReplace.txt")
    print(gene_id_replace_df)
    gene_id_replace_df.to_csv(geneid_replace_fname, sep="\t", index=False)

    get_seq_from_idlist(gene_id_replace_df, braker_aa, "on", braker_new_aa)
    get_seq_from_idlist(gene_id_replace_df, braker_cds, "on", braker_new_cds)


gff3_file = snakemake.input.braker_gff
braker_new_gff = snakemake.input.braker_new_gff[0]
braker_aa = snakemake.input.braker_aa
braker_cds = snakemake.input.braker_cds
braker_new_aa = snakemake.output.braker_new_aa[0]
braker_new_cds = snakemake.output.braker_new_cds[0]
output_dir = snakemake.params.output_dir
print(type(braker_new_aa), braker_new_aa)
print(type(braker_new_cds), braker_new_cds)
process_braker_output(gff3_file, braker_new_gff, braker_aa, braker_cds, braker_new_aa, braker_new_cds, output_dir)