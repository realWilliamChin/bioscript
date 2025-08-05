# --------------------------
# 比对流程
# --------------------------
rule rsem_prepare:
    input:
        trinity_fasta = f'{TRINITY_DIR}/{SPECIE_NAME}_unigene.fasta'
    output:
        fa = f'{REFERENCE_DIR}/{SPECIE_NAME}.transcripts.fa'
        # dynamic(f'{REFERENCE_DIR}/{SPECIE_NAME}.*')
    params:
        bowtie2 = '--bowtie2',
        reference_dir = f'{REFERENCE_DIR}/{SPECIE_NAME}'
    shell:
        """
        rsem-prepare-reference {params.bowtie2} \
            {input.trinity_fasta} \
            {params.reference_dir}
        """


rule rsem_calculate_expression:
    input:
        read1 = lambda wildcards: f'{DATA_DIR}/{samples_df.loc[samples_df["sample"] == wildcards.sample, "R1"].iloc[0]}',
        read2 = lambda wildcards: f'{DATA_DIR}/{samples_df.loc[samples_df["sample"] == wildcards.sample, "R2"].iloc[0]}',
        fa = f'{REFERENCE_DIR}/{SPECIE_NAME}.transcripts.fa'
    output:
        mapping_results = f'{MAPPING_DIR}/{{sample}}_mapping_stat.txt',
        genes_results = f'{MAPPING_DIR}/{{sample}}.genes.results',
        isoforms_results = f'{MAPPING_DIR}/{{sample}}.isoforms.results',
        transcript_bam = f'{MAPPING_DIR}/{{sample}}.transcript.bam'
    params:
        reference = f'{REFERENCE_DIR}/{SPECIE_NAME}',  # RSEM 参考基因组
    threads: config['rsem']['threads']
    log:
        f'{LOG_DIR}/rsem/{{sample}}_mapping.log'
    shell:
        """
        rsem-calculate-expression \
            -p {threads} \
            --bowtie2 \
            --strandedness reverse \
            --paired-end \
            {input.read1} {input.read2} \
            {params.reference} {MAPPING_DIR}/{wildcards.sample} \
            2> {output.mapping_results} \
            1> {log}
        """


rule mapping_summary:
    input:
        mapping_results = expand(f'{MAPPING_DIR}/{{sample}}_mapping_stat.txt', sample=samples_list),
        samples_file = f'{SAMPLES_FILE}'
    output:
        mapping_summary_file = f'{MAPPING_DIR}/mapping_summary.txt'
    shell:
        """
        python /home/colddata/qinqiang/script/noreference/scripts/rsem_mapping_summary.py \
            -i {MAPPING_DIR} \
            -s {input.samples_file} \
            -o {output.mapping_summary_file}
        """