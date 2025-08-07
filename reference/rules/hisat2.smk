# --------------------------
# 比对流程
# --------------------------
rule hisat2_alignment:
    input:
        read1 = lambda wildcards: f'{DATA_DIR}/{samples_df.loc[samples_df["sample"] == wildcards.sample, "R1"].iloc[0]}',
        read2 = lambda wildcards: f'{DATA_DIR}/{samples_df.loc[samples_df["sample"] == wildcards.sample, "R2"].iloc[0]}'
    output:
        mapping_results = f'{MAPPING_DIR}/{{sample}}_mapping_stat.txt',
        bam_results = f'{MAPPING_DIR}/{{sample}}.bam'
    threads: config['threads']
    params:
        ref_index = config['alignment']['ref_index']
    shell:
        """
        hisat2 -x {params.ref_index} \
            -p {threads} \
            --fr \
            -I 200 \
            -X 400 \
            --min-intronlen 20 \
            --max-intronlen 4000 \
            -1 {input.read1} \
            -2 {input.read2} \
            2> {output.mapping_results} | \
            samtools sort \
            --threads {threads} \
            -O BAM \
            -o \
            - > {output.bam_results}
        tail -n 1 {output.mapping_results}
        """

rule stringtie:
    input:
        bam_file = f'{MAPPING_DIR}/{{sample}}.bam',
        gff_file = config['alignment']['gff_file']
    output:
        fpkm_file = f'{EXPRESSION_DIR}/fpkm/{{sample}}_fpkm.txt',
        gtf_file = f'{EXPRESSION_DIR}/ballgown/{{sample}}/{{sample}}.gtf'
    shell:
        """
        stringtie -e -B \
            -G {input.gff_file} \
            -A {output.fpkm_file} \
            -o {output.gtf_file} \
            {input.bam_file}
        """


rule mapping_summary:
    input:
        mapping_results = expand(f'{MAPPING_DIR}/{{sample}}_mapping_stat.txt', sample=samples_list),
        samples_file = f'{SAMPLES_FILE}'
    output:
        mapping_summary_file = f'{MAPPING_DIR}/alignment_report.txt'
    shell:
        """
        python /home/colddata/qinqiang/script/CommonTools/alignment/alignment_mapping_summary.py \
            -i {MAPPING_DIR} \
            -s {input.samples_file} \
            -o {output.mapping_summary_file}
        """