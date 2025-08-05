# --------------------------
# Trinity 拼接流程
# --------------------------
rule trinity_prep:
    input:
        samples_file = f'{SAMPLES_FILE}',
        data_dir = f'{DATA_DIR}',
        samples_check_done = f'{WORK_DIR}/logs/samples_check.done'
    output:
        samples_trinity_file = f'{WORK_DIR}/samples_trinity.txt'
    params:
        run_type = config['trinity']['run_type']
    shell:
        """
        python /home/colddata/qinqiang/script/noreference/scripts/trinity_samples_file.py \
        -s {input.samples_file} \
        -d {input.data_dir} \
        -o {output.samples_trinity_file} \
        -t {params.run_type}
        """


rule run_trinity:
    input:
        samples_trinity_file = f'{WORK_DIR}/samples_trinity.txt'
    output:
        trinity_fasta = f'{TRINITY_DIR}/Trinity.fasta'
    params:
        output_dir = TRINITY_DIR
    threads: workflow.cores
    shell:
        """
        Trinity --seqType fq \
            --max_memory 500G \
            --no_salmon \
            --no_version_check \
            --samples_file {input.samples_trinity_file} \
            --output {params.output_dir} \
            --CPU {threads} \
            --SS_lib_type RF \
            --normalize_reads \
            --min_contig_length 500
        """

# --------------------------
# 提取最长转录本
# --------------------------
rule extract_longest_isoforms:
    input:
        trinity_fasta = f'{TRINITY_DIR}/Trinity.fasta'
    output:
        longest_fasta = f'{TRINITY_DIR}/unigene_longest.fasta'
    shell:
        """
        perl /home/colddata/qinqiang/script/noreference/scripts/extract_longest_isoforms_from_TrinityFasta.pl \
            {input.trinity_fasta} \
            > {output.longest_fasta}
        """

# --------------------------
# CD-HIT 去冗余
# --------------------------
rule run_cdhit:
    input:
        longest_fasta = f'{TRINITY_DIR}/unigene_longest.fasta'
    output:
        cdhit_fasta = f'{TRINITY_DIR}/{SPECIE_NAME}_unigene.fasta'
    params:
        identity_threshold = 0.95  # CD-HIT 相似度阈值
    shell:
        """
        cd-hit-est -i {input.longest_fasta} \
            -o {output.cdhit_fasta} \
            -c {params.identity_threshold}
        """


rule assemble_trinitystats_report:
    input:
        unigene_fasta = f'{TRINITY_DIR}/{SPECIE_NAME}_unigene.fasta'
    output:
        report_file = f'{TRINITY_DIR}/trinity_stats.txt'
    shell:
        """
        perl /opt/biosoft/Trinity-v2.8.5/util/TrinityStats.pl \
        {input.unigene_fasta} > {output.report_file}
        """


rule assemble_seqkitstats_report:
    input:
        unigene_fasta = f'{TRINITY_DIR}/{SPECIE_NAME}_unigene.fasta'
    output:
        report_file = f'{TRINITY_DIR}/seqkit_stats.txt'
    shell:
        """
        seqkit stats {input.unigene_fasta} > {output.report_file}
        """

rule assemble_report:
    input:
        trinity_stats_file = f'{TRINITY_DIR}/trinity_stats.txt',
        seqkit_stats_file = f'{TRINITY_DIR}/seqkit_stats.txt'
    output:
        report_file = f'{TRINITY_DIR}/assemble_report.txt'
    run:
        # 从 trinity_stats.txt 中提取数据
        with open(input.trinity_stats_file, 'r') as f:
            trinity_content = f.read()

        total_genes = trinity_content.split("Total trinity 'genes':")[1].split()[0].strip()
        total_transcripts = trinity_content.split('Total trinity transcripts:')[1].split()[0].strip()
        percent_gc = trinity_content.split('Percent GC:')[1].split()[0].strip()
        total_assembled_bases = trinity_content.split('Total assembled bases:')[1].split()[0].strip()
        n50 = trinity_content.split('Contig N50:')[1].split()[0].strip()
        average_length = trinity_content.split('Average contig:')[1].split()[0].strip()

        # 从 seqkit_stats.txt 中提取数据
        with open(input.seqkit_stats_file, 'r') as f:
            seqkit_content = f.readlines()[1]  # 读取第二行（数据行）

        seqkit_data = seqkit_content.strip().split()
        num_seqs = seqkit_data[3].replace(',', '')
        min_len = seqkit_data[5].replace(',', '')
        max_len = seqkit_data[7].replace(',', '')

        # 写入输出文件
        with open(output.report_file, 'w') as f:
            f.write(f'Total_genes\t{total_genes}\n')
            f.write(f'Total_transcripts\t{total_transcripts}\n')
            f.write(f'Percent_GC\t{percent_gc}\n')
            f.write(f'Total_assembled_bases\t{total_assembled_bases}\n')
            f.write(f'N50\t{n50}\n')
            f.write(f'Average_length\t{average_length}\n')
            f.write(f'Num_seqs\t{num_seqs}\n')
            f.write(f'Min_len\t{min_len}\n')
            f.write(f'Max_len\t{max_len}\n')
