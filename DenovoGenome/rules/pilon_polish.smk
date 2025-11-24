# configfile: "/home/colddata/qinqiang/work/2024_06_25_denovegenome_test/Test2/denovo_genome_config.yaml"

############################################
# 参数预定义（可通过 configfile 或 --config 覆盖）
############################################

# 参考基因组与测序数据
INPUTFASTA           = config['inputfasta']        # 例如：01.Assemble/assembly.fasta
READ1              = config['f1']                # 例如：sample_1.fq.gz
READ2              = config['f2']                # 例如：sample_2.fq.gz

# 物种及线程数
SPECIE             = config['specie']
THREADS            = config['threads']

# 路径参数
BWA_INDEX_DIR              = "bwa_index"
BWA_INDEX_PREFIX           = BWA_INDEX_DIR + "/" + SPECIE

BWA_SORTED_BAM             = "bwa_alignment_sorted.bam"
BWA_MARKDUP_BAM            = "bwa_alignment_sorted_markdup.bam"
BWA_FILTER_BAM             = "bwa_alignment_filter.bam"
BWA_FILTER_SORTED_BAM      = "bwa_alignment_filter_sorted.bam"
BWA_FILTER_SORTED_BAM_BAI  = BWA_FILTER_SORTED_BAM + ".bai"

PILON_OUTPUT_PREFIX        = "pilon_polished"
PILON_POLISHED_FASTA       = PILON_OUTPUT_PREFIX + ".fasta"

# Pilon 相关参数（如需也可放入 config 中再在此读取）
PILON_JAVA_MEMORY          = "-Xmx128G"
PILON_JAR                  = "/home/data/opt/biosoft/pilon/pilon-1.23.jar"
PILON_FIX_OPTIONS          = "--fix snps,indels"

rule all:
    input:
        PILON_POLISHED_FASTA

# 建立索引（参考序列由运行时 --config 提供，如：--config inputfasta=01.Assemble/assembly.fasta）
rule bwa_index:
    input:
        INPUTFASTA
    output:
        expand("{prefix}.{ext}", prefix=BWA_INDEX_PREFIX, ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        prefix=BWA_INDEX_PREFIX
    conda: 'base'
    message:
        "执行 bwa index 命令中\n"
        "输入文件: {input}\n"
        "输出文件: {output}"
    shell:
        'bwa index -p {params.prefix} {input}'

# 比对（reads 文件由运行时 --config 提供，如：--config f1=xxx.fq.gz f2=yyy.fq.gz）
# 这里显式声明依赖 BWA 索引文件，从而保证先执行 rule bwa_index 再执行比对
rule bwa_alignment:
    input:
        f1 = READ1,
        f2 = READ2,
        idx = expand("{prefix}.{ext}", prefix=BWA_INDEX_PREFIX, ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        BWA_SORTED_BAM
    params:
        prefix = BWA_INDEX_PREFIX
    threads: THREADS
    conda: 'base'
    message:
        "执行 bwa mem 命令中\n"
        "输入文件: {input}\n"
        "输出文件: {output}"
    shell:
        "bwa mem -t {threads} {params.prefix} {input.f1} {input.f2} | samtools sort -@ {threads} -O bam -o {output}"

# 标记重复
rule sambamba_markdup:
    input:
        BWA_SORTED_BAM
    output:
        BWA_MARKDUP_BAM
    threads: THREADS
    conda: "base"
    message:
        "执行 sambamba markdup 命令中\n"
        "输入文件: {input}\n"
        "输出文件: {output}"
    shell:
        "sambamba markdup -t {threads} {input} {output}"

# 过滤出高质量比对的 read
rule samtools_filter:
    input:
        BWA_MARKDUP_BAM
    output:
        BWA_FILTER_BAM
    threads: THREADS
    conda: "base"
    message:
        "执行 samtools view 过滤出高质量比对的 reads\n"
        "输入文件: {input}\n"
        "输出文件: {output}"
    shell:
        "samtools view -t {threads} -b -@ {threads} -q 30 {input} > {output}"

# sort
rule samtools_sort_and_index:
    input:
        BWA_FILTER_BAM
    output:
        BWA_FILTER_SORTED_BAM,
        BWA_FILTER_SORTED_BAM_BAI
    threads: THREADS
    conda: "base"
    message:
        "执行 samtools sort 命令中\n"
        "输入文件: {input}\n"
        "输出文件: {output}"
    shell:
        "samtools sort -@ {threads} -o {output[0]} {input} && "
        "samtools index -@ {threads} {output[0]} {output[1]}"
        
# polish
rule pilon_polish:
    input:
        INPUTFASTA,
        BWA_FILTER_SORTED_BAM,
        BWA_FILTER_SORTED_BAM_BAI
    output:
        PILON_POLISHED_FASTA
    conda: 'base'
    message:
        "执行 Pilon Polish 命令中\n"
        "输入文件: {input}\n"
        "输出文件: {output} "
    params:
        output_prefix=PILON_OUTPUT_PREFIX,
        java_memory=PILON_JAVA_MEMORY,
        pilon_jar=PILON_JAR,
        fix_options=PILON_FIX_OPTIONS
    threads: THREADS
    shell:
        'java {params.java_memory} -jar {params.pilon_jar} '
        '--genome {input[0]} '
        '--frags {input[1]} '
        '{params.fix_options} '
        '--threads {threads} '
        '--output {params.output_prefix}'

