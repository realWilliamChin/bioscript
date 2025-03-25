# --------------------------
# Expression Data 处理流程
# --------------------------

# 提出 fpkm
rule get_samples_fpkm:
    input:
        fpkm_file = expand(f'{EXPRESSION_DIR}/fpkm/{{sample}}_fpkm.txt', sample=samples_list),
        gtf_file = expand(f'{EXPRESSION_DIR}/ballgown/{{sample}}/{{sample}}.gtf', sample=samples_list)
    output:
        gene_matrix_file = f'{EXPRESSION_DIR}/gene_fpkm_matrix.csv',
        transcript_file = f'{EXPRESSION_DIR}/transcript_fpkm_matrix.csv'
    shell:
        """
        python2 /home/colddata/qinqiang/script/reference/scripts/getFPKM.py \
            -i {EXPRESSION_DIR}/ballgown \
            -g {output.gene_matrix_file} \
            -t {output.transcript_file}
        """

# 提出 raw reads
rule get_samples_reads:
    input:
        fpkm_file = expand(f'{EXPRESSION_DIR}/fpkm/{{sample}}_fpkm.txt', sample=samples_list),
        gtf_file = expand(f'{EXPRESSION_DIR}/ballgown/{{sample}}/{{sample}}.gtf', sample=samples_list)
    output:
        gene_matrix_file = f'{EXPRESSION_DIR}/gene_count_matrix.csv',
        transcript_file = f'{EXPRESSION_DIR}/transcript_count_matrix.csv'
    shell:
        """
        python2 /home/colddata/qinqiang/script/reference/scripts/prepDE.py \
            -i {EXPRESSION_DIR}/ballgown \
            -g {output.gene_matrix_file} \
            -t {output.transcript_file}
        """


# count filtered
rule count_filtered:
    input:
        gene_count_matrix_file = f'{EXPRESSION_DIR}/gene_count_matrix.csv',
        samples_file = f'{SAMPLES_FILE}'
    output:
        gene_count_file = f'{EXPRESSION_DIR}/gene_count_matrix.txt',
        gene_count_filtered_file = f'{EXPRESSION_DIR}/reads_matrix_filtered.txt'
    params:
        count_threshold = config['count_threshold']
    shell:
        """
        python /home/colddata/qinqiang/script/reference/scripts/count_filtered.py \
            -n {params.count_threshold} \
            -i {input.gene_count_matrix_file} \
            -o {output.gene_count_filtered_file} \
            -s {input.samples_file}
        """

# 计算 FPKM 并过滤
rule fpkm_filtered:
    input:
        gene_count_filtered_file = f'{EXPRESSION_DIR}/reads_matrix_filtered.txt',
        gene_fpkm_matrix_file = f'{EXPRESSION_DIR}/gene_fpkm_matrix.csv',
        samples_file = f'{SAMPLES_FILE}'
    output:
        gene_fpkm_file = f'{EXPRESSION_DIR}/gene_fpkm_matrix.txt',
        gene_fpkm_filtered_file = f'{EXPRESSION_DIR}/fpkm_matrix_filtered.txt'
    shell:
        """
        python /home/colddata/qinqiang/script/reference/scripts/fpkm_filtered.py \
            -r {input.gene_count_filtered_file} \
            -i {input.gene_fpkm_matrix_file} \
            -o {output.gene_fpkm_filtered_file} \
            -s {input.samples_file}
        """


# 合并 fpkm 和 reads 加上基因定义
rule merge_fpkm_reads:
    input:
        gene_count_file = f'{EXPRESSION_DIR}/gene_count_matrix.txt',
        gene_fpkm_file = f'{EXPRESSION_DIR}/gene_fpkm_matrix.txt',
        gene_count_filtered_file = f'{EXPRESSION_DIR}/reads_matrix_filtered.txt',
        gene_fpkm_filtered_file = f'{EXPRESSION_DIR}/fpkm_matrix_filtered.txt',
        kns_file = config['paths']['kns_file']
    output:
        expression_matrix = f'{EXPRESSION_DIR}/fpkm_reads_matraix_gene_def.txt',
        expression_filtered_matrix = f'{EXPRESSION_DIR}/fpkm_reads_matraix_filtered_gene_def.txt'
    shell:
        """
        python /home/colddata/qinqiang/script/transcriptome/merge_fpkm_reads_matrix.py \
            -f {input.gene_fpkm_file} \
            -r {input.gene_count_file} \
            --kns {input.kns_file} \
            -o {output.expression_matrix}
        python /home/colddata/qinqiang/script/transcriptome/merge_fpkm_reads_matrix.py \
            -f {input.gene_fpkm_filtered_file} \
            -r {input.gene_count_filtered_file} \
            --kns {input.kns_file} \
            -o {output.expression_filtered_matrix}
        """
