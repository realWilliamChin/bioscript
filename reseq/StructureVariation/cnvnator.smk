#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2026/06/01 10:34
# Author        : William GoGo

import os

# -------------------------- 配置部分 --------------------------
# 通过 snakemake --config 传入，或在 config.yaml 中指定
# 例: snakemake -s cnvnator.smk --config reference=/path/ref.fa chromosomes="chr1 chr2"
configfile: "cnvnator_config.yaml" if os.path.exists("cnvnator_config.yaml") else "/dev/null"

# 软件路径配置
GATK = config.get("gatk", "/home/data/opt/biosoft/gatk-4.4.0.0/gatk")
CNVNATOR = config.get("cnvnator", "cnvnator")
CNVNATOR2VCF = config.get("cnvnator2vcf", "cnvnator2VCF.pl")

# 参考基因组相关配置
reference = config.get("reference", "")  # 参考基因组 (必填)
mask_fasta = config.get("mask_fasta", "")  # mask 基因组 (必填)
mask_fasta_dir = config.get("mask_fasta_dir", "Mask_fasta")  # 每个染色体 fasta 文件目录

# 输入 SNP/INDEL VCF 前缀 (合并前)
snp_vcf_in = config.get("snp_vcf", "merged_SNPs_filterPASSED.vcf")
indel_vcf_in = config.get("indel_vcf", "merged_INDELs_filterPASSED.vcf")
merged_vcf_out = config.get("merged_vcf", "merged_SNPs_InDels_passd.vcf")
root_basename = config.get("root_basename", "cnvnator")  # *.root 文件前缀

# 样本列表
SAMPLES = [x.strip() for x in open(config.get("samples_lst", "samples.lst"))]

# 目标染色体列表 (空字符串表示全染色体)
CHROMOSOMES = config.get("chromosomes", "")

# 参数校验
if not reference:
    raise ValueError("必须通过 --config reference=... 指定参考基因组路径")
if not mask_fasta:
    raise ValueError("必须通过 --config mask_fasta=... 指定 mask 基因组路径")

# -------------------------- 规则定义 --------------------------
rule all:
    input:
        expand("{sample}/{sample}_CNV.vcf", sample=SAMPLES)

# 1. 合并SNP和INDEL的VCF文件
rule merge_snp_indel_vcf:
    input:
        snp_vcf = snp_vcf_in,
        indel_vcf = indel_vcf_in
    output:
        merged_vcf = merged_vcf_out
    threads: 4
    resources:
        mem_mb = 16000
    shell:
        """
        {GATK} MergeVcfs \
            -I {input.snp_vcf} \
            -I {input.indel_vcf} \
            -O {output.merged_vcf}
        """

# 2. 按样本拆分合并后的VCF文件
rule select_sample_variants:
    input:
        merged_vcf = merged_vcf_out,
        ref = reference
    output:
        sample_vcf = "{sample}/SNPs_InDels_passd.vcf"
    threads: 4
    resources:
        mem_mb = 128000
    shell:
        """
        mkdir -p $(dirname {output.sample_vcf})
        {GATK} --java-options -Xmx128G SelectVariants \
            -R {input.ref} \
            -V {input.merged_vcf} \
            -sn {wildcards.sample} \
            -O {output.sample_vcf}
        """

# 3. CNVnator全流程分析
rule run_cnvnator_analysis:
    input:
        bam = "../../02_Mapping/{sample}.final.bam",
        sample_vcf = "{sample}/SNPs_InDels_passd.vcf",
        mask_fa = mask_fasta
    output:
        cnv_result = "{sample}/{sample}_CNV.txt",
        root_file = "{sample}/" + root_basename + ".root"
    threads: 8
    resources:
        mem_mb = 64000
    shell:
        """
        cd {wildcards.sample}
        ROOT_FILE={root_basename}.root

        # 从BAM文件生成root文件
        {CNVNATOR} -root $ROOT_FILE -tree {input.bam} -chrom {CHROMOSOMES}

        # Mask基因组重复区域
        {CNVNATOR} -root $ROOT_FILE -mask {input.mask_fa}

        # 构建500bp窗口的histogram
        {CNVNATOR} -root $ROOT_FILE -his 500 -chrom {CHROMOSOMES} -fasta {input.mask_fa}

        # 导入样本VCF变异信息
        {CNVNATOR} -root $ROOT_FILE -vcf {input.sample_vcf}

        # 统计计算
        {CNVNATOR} -root $ROOT_FILE -stat 500

        # 分区处理
        {CNVNATOR} -root $ROOT_FILE -partition 500

        # BAF计算
        {CNVNATOR} -root $ROOT_FILE -baf 10000

        # CNV Calling并输出结果
        echo -e 'CNV_type\\tcoordinates\\tCNV_size\\tnormalized_RD\\te-val1\\te-val2\\te-val3\\te-val4\\tq0' > {wildcards.sample}_CNV.txt
        {CNVNATOR} -root $ROOT_FILE -call 500 >> {wildcards.sample}_CNV.txt
        """

# 4. 将CNV结果转换为VCF格式
rule convert_cnvnator_to_vcf:
    input:
        cnv_result = "{sample}/{sample}_CNV.txt"
    output:
        cnv_vcf = "{sample}/{sample}_CNV.vcf"
    shell:
        """
        {CNVNATOR2VCF} \
            -prefix {wildcards.sample} \
            -reference {reference} \
            {input.cnv_result} {mask_fasta_dir} \
            > {output.cnv_vcf}
        """