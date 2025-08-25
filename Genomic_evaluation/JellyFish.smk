import pandas as pd
import os
import sys

sys.path.append('/home/colddata/qinqiang/script/CommonTools')
from load_input import load_table, write_output_df

# 读取参数
samples_file = config.get("samples", "samples_described.txt")
kmer = config.get("kmer", 21)
genomescope_version = config.get("genomescope_version", 1)  # Genomescope, (2倍体，v1) (多倍体, v2)
threads_num = int(config.get("threads", 20))
outdir = config.get("outdir", "results")
ploidy = config.get("ploidy", 2)  # 默认二倍体
max_kmer_cov = config.get("max_kmer_cov", 10000)  # 最大k-mer覆盖度

def get_samples():
    df = load_table(samples_file)
    return df.to_dict("records")

SAMPLES = get_samples()

rule all:
    input:
        expand("{outdir}/{sample}_{kmer}mer_out.histo",
               zip,
               outdir=[outdir]*len(SAMPLES),
               sample=[s["sample"] for s in SAMPLES],
               kmer=[kmer]*len(SAMPLES)),
        expand("{outdir}/{sample}_{kmer}mer_genomescope_v{genomescope_version}",
               zip,
               outdir=[outdir]*len(SAMPLES),
               sample=[s["sample"] for s in SAMPLES],
               kmer=[kmer]*len(SAMPLES),
               genomescope_version=[genomescope_version]*len(SAMPLES))

rule jellyfish_count:
    input:
        r1=lambda wildcards: next(s["R1"] for s in SAMPLES if s["sample"] == wildcards.sample),
        r2=lambda wildcards: next(s["R2"] for s in SAMPLES if s["sample"] == wildcards.sample)
    output:
        jf="{outdir}/{sample}_{kmer}mer_out"
    threads: threads_num
    resources:
        mem_mb=20000
    shell:
        r'''
        jellyfish count -m {kmer} -s 20G -t {threads} -o {output.jf} -C \
            <(zcat {input.r1}) <(zcat {input.r2})
        '''

rule jellyfish_histo:
    input:
        jf="{outdir}/{sample}_{kmer}mer_out"
    output:
        histo="{outdir}/{sample}_{kmer}mer_out.histo"
    threads: threads_num
    resources:
        mem_mb=20000
    shell:
        r'''
        jellyfish histo -t {threads} -o {output.histo} {input.jf}
        '''


# 根据版本自动选择GenomeScope命令
rule genomescope:
    input:
        histo="{outdir}/{sample}_{kmer}mer_out.histo"
    output:
        directory("{outdir}/{sample}_{kmer}mer_genomescope_v{genomescope_version}")
    threads: 1
    resources:
        mem_mb=8000
    run:
        if genomescope_version == 1:
            shell(f"""
                mkdir -p {output.directory}
                Rscript /home/colddata/qinqiang/script/Rscript/genomescope_v1.R \
                    {input.histo} \
                    {kmer} \
                    {ploidy} \
                    {output.directory} \
                    {max_kmer_cov}
            """)
        elif genomescope_version == 2:
            shell(f"""
                mkdir -p {output.directory}
                Rscript /home/colddata/qinqiang/script/Rscript/genomescope.R \
                    -i {input.histo} \
                    -k {kmer} \
                    -o {output.directory}
            """)
        else:
            raise ValueError(f"Unsupported GenomeScope version: {genomescope_version}")