# --------------------------
# Expression Data 处理流程
# --------------------------

# 提出 raw reads
rule get_samples_reads:
    input:
        map_result = expand(f'{MAPPING_DIR}/{{sample}}.genes.results', sample=samples_list)
    output:
        gene_count_matrix_file = f'{EXPRESSION_DIR}/gene_count_matrix.txt'
    shell:
        """
        python /home/colddata/qinqiang/script/noreference/scripts/all_sample_raw_reads.py \
            -i {MAPPING_DIR} \
            -o {output.gene_count_matrix_file}
        """

# 提出 fpkm
rule get_samples_fpkm:
    input:
        map_result = expand(f'{MAPPING_DIR}/{{sample}}.genes.results', sample=samples_list)
    output:
        gene_fpkm_matrix_file = f'{EXPRESSION_DIR}/gene_fpkm_matrix.txt'
    shell:
        """
        python /home/colddata/qinqiang/script/noreference/scripts/all_sample_fpkm.py \
            -i {MAPPING_DIR} \
            -o {output.gene_fpkm_matrix_file}
        """

# count filtered
rule count_filtered:
    input:
        gene_count_matrix_file = f'{EXPRESSION_DIR}/gene_count_matrix.txt'
    output:
        gene_count_filtered_file = f'{EXPRESSION_DIR}/reads_matrix_filtered.txt'
    params:
        count_threshold = config['count_threshold']
    shell:
        """
        python /home/colddata/qinqiang/script/noreference/scripts/count_filtered.py \
            -n {params.count_threshold} \
            -i {input.gene_count_matrix_file} \
            -o {output.gene_count_filtered_file}
        """

# 计算 FPKM 并过滤
rule fpkm_filtered:
    input:
        gene_count_filtered_file = f'{EXPRESSION_DIR}/reads_matrix_filtered.txt',
        gene_fpkm_matrix_file = f'{EXPRESSION_DIR}/gene_fpkm_matrix.txt'
    output:
        gene_fpkm_filtered_file = f'{EXPRESSION_DIR}/fpkm_matrix_filtered.txt'
    shell:
        """
        python /home/colddata/qinqiang/script/noreference/scripts/fpkm_filtered.py \
            -r {input.gene_count_filtered_file} \
            -i {input.gene_fpkm_matrix_file} \
            -o {output.gene_fpkm_filtered_file}
        """


# 合并 fpkm 和 reads 加上基因定义，需要先运行 kns_annotation 流程，输出 kns_gene_def.txt 文件
# rule merge_fpkm_reads:
#     input:
#         gene_count_file = f'{EXPRESSION_DIR}/gene_count_matrix.txt',
#         gene_fpkm_file = f'{EXPRESSION_DIR}/gene_fpkm_matrix.txt',
#         gene_count_filtered_file = f'{EXPRESSION_DIR}/reads_matrix_filtered.txt',
#         gene_fpkm_filtered_file = f'{EXPRESSION_DIR}/fpkm_matrix_filtered.txt',
#         kns_file = config['paths']['kns_file']
#     output:
#         expression_matrix = f'{EXPRESSION_DIR}/fpkm_reads_matraix_gene_def.txt',
#         expression_filtered_matrix = f'{EXPRESSION_DIR}/fpkm_reads_matraix_filtered_gene_def.txt'
#     shell:
#         """
#         python /home/colddata/qinqiang/script/transcriptome/merge_fpkm_reads_matrix.py \
#             -f {input.gene_fpkm_file} \
#             -r {input.gene_count_file} \
#             --kns {input.kns_file} \
#             -o {output.expression_matrix}
#         python /home/colddata/qinqiang/script/transcriptome/merge_fpkm_reads_matrix.py \
#             -f {input.gene_fpkm_filtered_file} \
#             -r {input.gene_count_filtered_file} \
#             --kns {input.kns_file} \
#             -o {output.expression_filtered_matrix}
#         """