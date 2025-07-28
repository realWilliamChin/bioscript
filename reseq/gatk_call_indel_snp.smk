#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/24 16:34
# Author        : William GoGo

import os

# 配置参数
# configfile: "config.yaml"
# snakemake -s /home/colddata/qinqiang/script/reseq/gatk_call_indel_snp.smk --config input=Father.gvcf.gz reference="../00_seq_dict/Homo_sapiens.GRCh38.dna.toplevel.fa" prefix="Father" -c 20  &

# 读取配置参数
input_file = config.get("input")
reference = config.get("reference")
sample_ploidy = config.get("ploidy", 2)
prefix = config.get("prefix", "gatk")
outdir = config.get("outdir", os.getcwd())

# 确保输出目录存在
os.makedirs(outdir, exist_ok=True)

# 定义最终输出文件
final_outputs = [
    os.path.join(outdir, f"{prefix}_INDEL_filtered_PASSED.vcf"),
    os.path.join(outdir, f"{prefix}_SNP_filtered_PASSED.vcf")
]

localrules: activate_env

rule all:
    input:
        final_outputs

rule activate_env:
    shell:
        """
        source /opt/biosoft/gatk-4.4.0.0/env.sh
        """


rule bam2gvcf:
    input:
    output:
    shell:
        r'''
        gatk --java-options -Xmx128G HaplotypeCaller \
            -R 00_seq_dict/Homo_sapiens.GRCh38.dna.toplevel.fa \
            -I ../01_Bwa_alignment/Father.sorted.markdup.bam \
            -O Father.gvcf.gz \
            --emit-ref-confidence GVCF \
            -stand-call-conf 10 \
            --sample-ploidy 2 \
            --max-genotype-count 4 &
        '''


rule gvcf2vcf:
    input:
        gvcf = input_file,
        ref = reference
    output:
        vcf = os.path.join(outdir, f"{prefix}.vcf")
    log:
        os.path.join(outdir, "logs", "gvcf2vcf.log")
    shell:
        """
        mkdir -p {outdir}/logs
        /opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" GenotypeGVCFs \
            -R {input.ref} \
            -V {input.gvcf} \
            -O {output.vcf} \
            2> {log}
        """


rule call_snp_indel_from_vcf:
    input:
        vcf = os.path.join(outdir, f"{prefix}.vcf"),
        ref = reference
    output:
        indel_vcf = os.path.join(outdir, f"{prefix}_INDEL.vcf"),
        snp_vcf = os.path.join(outdir, f"{prefix}_SNP.vcf")
    log:
        indel_log = os.path.join(outdir, "logs", "select_indel.log"),
        snp_log = os.path.join(outdir, "logs", "select_snp.log")
    shell:
        """
        /opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -select-type INDEL \
            -O {output.indel_vcf} \
            2> {log.indel_log}
        
        /opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -select-type SNP \
            -O {output.snp_vcf} \
            2> {log.snp_log}
        """


rule filter_snp:
    input:
        snp_vcf = os.path.join(outdir, f"{prefix}_SNP.vcf"),
        ref = reference
    output:
        snp_filter_vcf = os.path.join(outdir, f"{prefix}_SNP_filtered.vcf")
    log:
        os.path.join(outdir, "logs", "filter_snp.log")
    shell:
        """
        /opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" VariantFiltration \
            -R {input.ref} \
            -V {input.snp_vcf} \
            -filter "MQ < 30.0" --filter-name "MQ_filter_SNP" \
            -filter "QD < 2.0" --filter-name "QD_filter_SNP" \
            -filter "FS > 60.0" --filter-name "FS_filter_SNP" \
            -filter "SOR > 3.0" --filter-name "SOR_filter_SNP" \
            -O {output.snp_filter_vcf} \
            2> {log}
        """


rule filter_indel:
    input:
        indel_vcf = os.path.join(outdir, f"{prefix}_INDEL.vcf"),
        ref = reference
    output:
        indel_filter_vcf = os.path.join(outdir, f"{prefix}_INDEL_filtered.vcf")
    log:
        os.path.join(outdir, "logs", "filter_indel.log")
    shell:
        """
        /opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" VariantFiltration \
            -R {input.ref} \
            -V {input.indel_vcf} \
            -filter "MQ < 30.0" --filter-name "MQ_filter_INDEL" \
            -filter "SOR > 10.0" --filter-name "SOR_filter_INDEL" \
            -filter "QD < 2.0" --filter-name "QD_filter_INDEL" \
            -filter "FS > 200.0" --filter-name "FS_filter_INDEL" \
            -O {output.indel_filter_vcf} \
            2> {log}
        """


rule filter_passed_snp:
    input:
        snp_filter_vcf = os.path.join(outdir, f"{prefix}_SNP_filtered.vcf")
    output:
        snp_passed_vcf = os.path.join(outdir, f"{prefix}_SNP_filtered_PASSED.vcf")
    log:
        os.path.join(outdir, "logs", "filter_passed_snp.log")
    shell:
        """
        grep -E '^#|PASS' {input.snp_filter_vcf} > {output.snp_passed_vcf} 2> {log}
        """


rule filter_passed_indel:
    input:
        indel_filter_vcf = os.path.join(outdir, f"{prefix}_INDEL_filtered.vcf")
    output:
        indel_passed_vcf = os.path.join(outdir, f"{prefix}_INDEL_filtered_PASSED.vcf")
    log:
        os.path.join(outdir, "logs", "filter_passed_indel.log")
    shell:
        """
        grep -E '^#|PASS' {input.indel_filter_vcf} > {output.indel_passed_vcf} 2> {log}
        """