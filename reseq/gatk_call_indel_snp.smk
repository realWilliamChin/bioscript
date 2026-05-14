#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/24 16:34
# Author        : William GoGo

import os

# 配置参数
# configfile: "config.yaml"
# 运行示例1：从gvcf文件开始
# snakemake -s /home/colddata/qinqiang/script/reseq/gatk_call_indel_snp.smk \
#   --config input=Father.g.vcf.gz \
#            reference="../00_seq_dict/Homo_sapiens.GRCh38.dna.toplevel.fa" \
#            prefix="Father" \
#            ploidy=2 \
#            max_genotype_count=4 \
#   -c 20 &
#
# 运行示例2：从bam文件开始（全基因组模式，自动运行bam2gvcf步骤）
# snakemake -s /home/colddata/qinqiang/script/reseq/gatk_call_indel_snp.smk \
#   --config input_bam=02_bam/BSA_Sg-T_Mu_sorted_markdup.bam \
#            reference="../00_index/Luffa_aegyptiaca_genome.fna" \
#            prefix="BSA_Sg-T_Mu" \
#            ploidy=2 \
#            max_genotype_count=4 \
#            java_mem=128G \
#   -c 20 &
#
# 运行示例3：按所有染色体拆分并行运行（自动获取参考基因组所有染色体）
# snakemake -s /home/colddata/qinqiang/script/reseq/gatk_call_indel_snp.smk \
#   --config input_bam=02_bam/BSA_Sg-T_Mu_sorted_markdup.bam \
#            reference="../00_index/Luffa_aegyptiaca_genome.fna" \
#            prefix="BSA_Sg-T_Mu" \
#            ploidy=2 \
#            split_by_chr=true \
#            java_mem=64G \
#            hc_extra_args="--genotyping-mode DISCOVERY -stand-call-conf 30" \
#            delete_chr_tmp=false \
#   -c 100 # 同时跑多个染色体任务，根据CPU核心数调整
#
# 运行示例4：只跑指定染色体
# snakemake -s /home/colddata/qinqiang/script/reseq/gatk_call_indel_snp.smk \
#   --config input_bam=02_bam/BSA_Sg-T_Mu_sorted_markdup.bam \
#            reference="../00_index/Luffa_aegyptiaca_genome.fna" \
#            prefix="BSA_Sg-T_Mu" \
#            intervals="chr1,chr2,chr3,chr4" \
#            merge_chr_output=true \
#   -c 32
#
# 运行示例5：批量处理文件夹下所有bam文件（自动拆分染色体并行）
# snakemake -s /home/colddata/qinqiang/script/reseq/gatk_call_indel_snp.smk \
#   --config input_folder="./02_bam/" \
#            reference="../00_index/Luffa_aegyptiaca_genome.fna" \
#            split_by_chr=true \
#            java_mem=64G \
#            bam_suffix=".sorted.markdup.bam" \
#            outdir_by_sample=true \
#   -c 200
#
# 运行示例6：从原始vcf文件开始，直接进行SNP/INDEL分离和过滤
# snakemake -s /home/colddata/qinqiang/script/reseq/gatk_call_indel_snp.smk \
#   --config input_vcf=sample.raw.vcf \
#            reference="../00_seq_dict/Homo_sapiens.GRCh38.dna.toplevel.fa" \
#            prefix="sample" \
#   -c 10 &

# 读取配置参数
input_file = config.get("input") # 输入gvcf文件
input_vcf = config.get("input_vcf") # 输入原始vcf文件，直接进入变异过滤流程
input_bam = config.get("input_bam") # 输入bam文件
input_folder = config.get("input_folder", None) # 批量bam文件夹路径，优先级最高
reference = config.get("reference")
sample_ploidy = config.get("ploidy", 2)
max_genotype_count = config.get("max_genotype_count", 4)
prefix = config.get("prefix", "gatk") # 单样本模式时的前缀
outdir = config.get("outdir", os.getcwd())

# 批量处理参数
bam_suffix = config.get("bam_suffix", ".bam") # bam文件后缀
outdir_by_sample = config.get("outdir_by_sample", True) # 是否按样本名创建子目录

# 新增功能参数
split_by_chr = config.get("split_by_chr", False) # 是否按染色体拆分运行
intervals = config.get("intervals", None) # 自定义运行的区间/染色体列表，逗号分隔，优先级高于split_by_chr
java_mem = config.get("java_mem", "128G") # HaplotypeCaller Java堆内存
hc_extra_args = config.get("hc_extra_args", "") # 额外传递给HaplotypeCaller的参数
merge_chr_output = config.get("merge_chr_output", True) # 拆分运行后是否合并所有染色体结果
delete_chr_tmp = config.get("delete_chr_tmp", False) # 合并后是否删除单个染色体的临时gvcf文件

# 处理样本列表
sample_list = []
sample_bam_map = {}
sample_outdir_map = {}
gvcf_path_map = {}

if input_folder:
    # 批量模式：扫描文件夹下所有bam文件
    if not os.path.isdir(input_folder):
        raise ValueError(f"input_folder {input_folder} 不是有效的文件夹路径")
    # 扫描所有符合后缀的文件
    for file in os.listdir(input_folder):
        if file.endswith(bam_suffix) and not (file.endswith('.bai') or file.endswith('.bam.bai')):
            sample_name = file[:-len(bam_suffix)]
            sample_list.append(sample_name)
            sample_bam_map[sample_name] = os.path.join(input_folder, file)
            # 确定样本输出目录
            if outdir_by_sample:
                sample_outdir = os.path.join(outdir, sample_name)
            else:
                sample_outdir = outdir
            sample_outdir_map[sample_name] = sample_outdir
            gvcf_path_map[sample_name] = os.path.join(sample_outdir, f"{sample_name}.HC.g.vcf.gz")
            os.makedirs(sample_outdir, exist_ok=True)
    if not sample_list:
        raise ValueError(f"在 {input_folder} 下没有找到后缀为 {bam_suffix} 的bam文件")
else:
    # 单样本模式
    sample_list = [prefix]
    sample_bam_map[prefix] = input_bam
    sample_outdir_map[prefix] = outdir
    if input_bam:
        gvcf_path_map[prefix] = os.path.join(outdir, f"{prefix}.HC.g.vcf.gz")
    elif input_vcf:
        # 直接输入vcf的情况，跳过gvcf转vcf步骤
        vcf_path_map = {prefix: input_vcf}
    else:
        gvcf_path_map[prefix] = input_file

# 获取染色体列表
chr_list = []
if (input_bam or input_folder) and (split_by_chr or intervals):
    if intervals:
        # 用户指定了自定义区间
        chr_list = [x.strip() for x in intervals.split(",") if x.strip()]
    else:
        # 自动从参考基因组fai文件获取所有染色体
        fai_file = reference + ".fai"
        if not os.path.exists(fai_file):
            # 自动创建fai索引
            shell("samtools faidx " + reference)
        # 读取fai第一列获取所有序列名
        with open(fai_file, "r") as f:
            chr_list = [line.split("\t")[0].strip() for line in f if line.strip()]

# 检查输入参数
input_count = sum([1 for x in [input_file, input_vcf, input_bam, input_folder] if x])
if input_count > 1:
    raise ValueError("input(gvcf)、input_vcf(原始vcf)、input_bam(单个bam)、input_folder(bam文件夹)只能提供其中一个")
if input_count == 0:
    raise ValueError("必须提供input(gvcf)、input_vcf(原始vcf)、input_bam(单个bam)、input_folder(bam文件夹)其中之一")

# 定义最终输出文件
final_outputs = []
for sample in sample_list:
    sample_outdir = sample_outdir_map[sample]
    final_outputs.extend([
        os.path.join(sample_outdir, f"{sample}_INDEL_filtered_PASSED.vcf"),
        os.path.join(sample_outdir, f"{sample}_SNP_filtered_PASSED.vcf")
    ])

localrules: activate_env

rule all:
    input:
        final_outputs

rule activate_env:
    shell:
        """
        source /home/data/opt/biosoft/gatk-4.4.0.0/env.sh
        """

# 只有在输入是bam文件或bam文件夹时才定义bam2gvcf相关规则
if input_bam or input_folder:
    if split_by_chr or intervals:
        # 按染色体拆分运行模式（支持多样本）
        rule bam2gvcf_per_chr:
            input:
                bam = lambda wildcards: sample_bam_map[wildcards.sample],
                ref = reference,
                fai = reference + ".fai"
            output:
                gvcf = "{sample_outdir}/{sample}.{chr}.HC.g.vcf.gz",
                tbi = "{sample_outdir}/{sample}.{chr}.HC.g.vcf.gz.tbi"
            log:
                "{sample_outdir}/logs/bam2gvcf_{sample}_{chr}.log"
            threads: 1
            shell:
                r'''
                mkdir -p $(dirname {log})
                /home/data/opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx{java_mem}" HaplotypeCaller \
                    -R {input.ref} \
                    -I {input.bam} \
                    -O {output.gvcf} \
                    --emit-ref-confidence GVCF \
                    --sample-ploidy {sample_ploidy} \
                    --max-genotype-count {max_genotype_count} \
                    --intervals {wildcards.chr} \
                    {hc_extra_args} \
                    2> {log}
                '''
        
        # 合并所有染色体结果规则（支持多样本）
        rule merge_chr_gvcf:
            input:
                lambda wildcards: expand("{outdir}/{sample}.{chr}.HC.g.vcf.gz", outdir=sample_outdir_map[wildcards.sample], sample=wildcards.sample, chr=chr_list)
            output:
                gvcf = "{sample_outdir}/{sample}.HC.g.vcf.gz",
                tbi = "{sample_outdir}/{sample}.HC.g.vcf.gz.tbi"
            log:
                "{sample_outdir}/logs/merge_chr_gvcf_{sample}.log"
            shell:
                r'''
                mkdir -p $(dirname {log})
                # 合并所有染色体gvcf
                INPUT_ARGS=""
                for GVCF in {input:q}; do
                    INPUT_ARGS="$INPUT_ARGS -I $GVCF"
                done

                /home/data/opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" GatherVcfs \
                    $INPUT_ARGS \
                    -O {output.gvcf} \
                    2>> {log}

                # 为合并后的gvcf文件创建索引
                tabix -p vcf {output.gvcf} \
                    2>> {log}

                # 可选删除临时文件
                if [ "{delete_chr_tmp}" = "True" ]; then
                    for GVCF in {input:q}; do
                        rm -f $GVCF ${{GVCF}}.tbi
                    done
                fi
                '''
    else:
        # 全基因组非拆分模式（支持多样本）
        rule bam2gvcf:
            input:
                bam = lambda wildcards: sample_bam_map[wildcards.sample],
                ref = reference
            output:
                gvcf = "{sample_outdir}/{sample}.HC.g.vcf.gz"
            log:
                "{sample_outdir}/logs/bam2gvcf_{sample}.log"
            shell:
                r'''
                mkdir -p $(dirname {log})
                /home/data/opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx{java_mem}" HaplotypeCaller \
                    -R {input.ref} \
                    -I {input.bam} \
                    -O {output.gvcf} \
                    --emit-ref-confidence GVCF \
                    --sample-ploidy {sample_ploidy} \
                    --max-genotype-count {max_genotype_count} \
                    {hc_extra_args} \
                    2> {log}
                '''


rule gvcf2vcf:
    input:
        gvcf = lambda wildcards: gvcf_path_map[wildcards.sample],
        ref = reference
    output:
        vcf = "{sample_outdir}/{sample}.vcf"
    log:
        "{sample_outdir}/logs/gvcf2vcf_{sample}.log"
    shell:
        """
        mkdir -p $(dirname {log})
        /home/data/opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" GenotypeGVCFs \
            -R {input.ref} \
            -V {input.gvcf} \
            -O {output.vcf} \
            2> {log}
        """


rule call_snp_indel_from_vcf:
    input:
        vcf = lambda wildcards: vcf_path_map[wildcards.sample] if 'vcf_path_map' in globals() else f"{sample_outdir_map[wildcards.sample]}/{wildcards.sample}.vcf",
        ref = reference
    output:
        indel_vcf = "{sample_outdir}/{sample}_INDEL.vcf",
        snp_vcf = "{sample_outdir}/{sample}_SNP.vcf"
    log:
        indel_log = "{sample_outdir}/logs/select_indel_{sample}.log",
        snp_log = "{sample_outdir}/logs/select_snp_{sample}.log"
    shell:
        """
        /home/data/opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -select-type INDEL \
            -O {output.indel_vcf} \
            2> {log.indel_log}
        
        /home/data/opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -select-type SNP \
            -O {output.snp_vcf} \
            2> {log.snp_log}
        """


rule filter_snp:
    input:
        snp_vcf = "{sample_outdir}/{sample}_SNP.vcf",
        ref = reference
    output:
        snp_filter_vcf = "{sample_outdir}/{sample}_SNP_filtered.vcf"
    log:
        "{sample_outdir}/logs/filter_snp_{sample}.log"
    shell:
        """
        /home/data/opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" VariantFiltration \
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
        indel_vcf = "{sample_outdir}/{sample}_INDEL.vcf",
        ref = reference
    output:
        indel_filter_vcf = "{sample_outdir}/{sample}_INDEL_filtered.vcf"
    log:
        "{sample_outdir}/logs/filter_indel_{sample}.log"
    shell:
        """
        /home/data/opt/biosoft/gatk-4.4.0.0/gatk --java-options "-Xmx200g" VariantFiltration \
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
        snp_filter_vcf = "{sample_outdir}/{sample}_SNP_filtered.vcf"
    output:
        snp_passed_vcf = "{sample_outdir}/{sample}_SNP_filtered_PASSED.vcf"
    log:
        "{sample_outdir}/logs/filter_passed_snp_{sample}.log"
    shell:
        """
        grep -E '^#|PASS' {input.snp_filter_vcf} > {output.snp_passed_vcf} 2> {log}
        """


rule filter_passed_indel:
    input:
        indel_filter_vcf = "{sample_outdir}/{sample}_INDEL_filtered.vcf"
    output:
        indel_passed_vcf = "{sample_outdir}/{sample}_INDEL_filtered_PASSED.vcf"
    log:
        "{sample_outdir}/logs/filter_passed_indel_{sample}.log"
    shell:
        """
        grep -E '^#|PASS' {input.indel_filter_vcf} > {output.indel_passed_vcf} 2> {log}
        """