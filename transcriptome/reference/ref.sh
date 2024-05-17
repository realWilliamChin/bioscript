#!/bin/bash
# set -e -o pipefail


### 比对
hisat2_alignment() {
    samplename_list=($(tail -n+2 ${samples_described_f} | grep -v '^$' | cut -f 2))
    filename_1_list=($(tail -n+2 ${samples_described_f} | grep -v '^$' | cut -f 3))
    filename_2_list=($(tail -n+2 ${samples_described_f} | grep -v '^$' | cut -f 4))
    sample_count=${#filename_1_list[@]}
    log INFO "即将比对的样本数量是 $sample_count"

    # 检查 filename_1_list 和 filename_2_list 的文件是否存在
    log INFO "检查 ${samples_described_f} 中的文件是否存在"
    for ((i=0; i<sample_count; i++)); do
        file_1="${cleandata_d}/${filename_1_list[i]}"
        file_2="${cleandata_d}/${filename_2_list[i]}"

        if [[ ! -e "$file_1" ]]; then
            log ERROR "未找到文件: $file_1"
            ls -l "$file_1"  # 添加调试信息
            exit 1
        fi

        if [[ ! -e "$file_2" ]]; then
            log ERROR "未找到文件: $file_2"
            ls -l "$file_2"  # 添加调试信息
            exit 1
        fi
    done
    log INFO "检查完成，无错误"

    # 检查比对的参考基因组的库是否存在
    if [[ ! -f ${database_d}/${ref_specie}.1.ht2 ]]; then
        log ERROR "未找到 ${database_d}/${ref_specie} 相关库文件"
        exit 1
    fi

    # 创建 ${bam_d} 和 ${mapping_d} 目录
    if [[ ! -d ${bam_d} ]]; then
        mkdir ${bam_d}
    fi
    if [[ ! -d ${mapping_d} ]]; then
        mkdir ${mapping_d}
    fi

    for ((i=0; i<sample_count; i++)); do
        log INFO "正在比对第${i}个，使用线程数 ${num_threads}，${filename_1_list[i]} ${filename_2_list[i]} hisat2 alignment ... "
        hisat2 -x ${database_d}/${ref_specie} \
            -p ${num_threads} \
            -I 200 -X 400 --fr \
            --min-intronlen 20 --max-intronlen 4000 \
            -1 ${cleandata_d}/${filename_1_list[i]} \
            -2 ${cleandata_d}/${filename_2_list[i]} \
            2> ${mapping_d}/${samplename_list[i]}_mapping.txt | \
            samtools sort --threads ${num_threads} -O BAM -o - > ${bam_d}/${samplename_list[i]}.bam

        tail -n 1 ${mapping_d}/${samplename_list[i]}_mapping.txt
    done

    log INFO "hisat2 mapping summary"
    cd ${mapping_d} || exit
    python ${script}/hisat2_mapping_summary.py -s ${samples_described_f}
    cat ${mapping_d}/mapping_summary.txt
    cd ${work_dir} || exit
}

bowtie2_alignment() {
    # 创建 ${bam_d} 和 ${mapping_d} 目录
    if [[ ! -d ${bam_d} ]]; then
        mkdir ${bam_d}
    fi
    if [[ ! -d ${mapping_d} ]]; then
        mkdir ${mapping_d}
    fi

    samplename_list=($(tail -n+2 ${samples_described_f} | grep -v '^$' | cut -f 2))
    filename_1_list=($(tail -n+2 ${samples_described_f} | grep -v '^$' | cut -f 3))
    filename_2_list=($(tail -n+2 ${samples_described_f} | grep -v '^$' | cut -f 4))
    sample_count=${#filename_1_list[@]}
    log INFO "即将比对的样本数量是 $sample_count"

    # 检查 filename_1_list 和 filename_2_list 的文件是否存在
    log INFO "检查 ${samples_described_f} 中的文件是否存在"
    for ((i=0; i<sample_count; i++)); do
        if [[ ! -e "${cleandata_d}/${filename_1_list[i]}" ]]; then
            log ERROR "未找到 ${cleandata_d}/${filename_1_list[i]} 文件"
            exit 1
        fi
        if [[ ! -e "${cleandata_d}/${filename_2_list[i]}" ]]; then
            log ERROR "未找到 ${cleandata_d}/${filename_2_list[i]} 文件"
            exit 1
        fi
    done
    log INFO "检查完成，无错误"

    # 检查比对的参考基因组的库是否存在
    if [[ ! -f ${database_d}/${ref_specie}.1.bt2 ]]; then
        log ERROR "未找到 ${database_d}/${ref_specie} 相关库文件"
        exit 1
    fi

    for ((i=0; i<sample_count; i++)); do
        # 创建单独的 样本 bam 目录，存放没比对上的 fastq 和比对上的 fastq 文件
        if [ ! -d ${bam_d}/${samplename_list[i]} ]; then
            mkdir ${bam_d}/${samplename_list[i]}
        fi

        log INFO "正在比对第${i}个，使用线程数 ${num_threads}，${filename_1_list[i]} ${filename_2_list[i]} bowtie2 alignment ... "
        bowtie2 -x ${database_d}/${ref_specie} \
            -p ${num_threads} \
            -1 ${cleandata_d}/${filename_1_list[i]} \
            -2 ${cleandata_d}/${filename_2_list[i]} \
            --un-conc-gz ${bam_d}/${samplename_list[i]}/ \
            --al-conc-gz ${bam_d}/${samplename_list[i]}/ \
            2> ${mapping_d}/${samplename_list[i]}_mapping.txt \
            | samtools view -@ ${num_threads} -bu - \
            | samtools sort -@ ${num_threads} -o ${bam_d}/${samplename_list[i]}.bam

        tail -n 1 ${mapping_d}/${samplename_list[i]}_mapping.txt
    done
    cd ${mapping_d} || exit
    python ${script}/hisat2_mapping_summary.py -s ${samples_described_f}
    cd ${work_dir} || exit
}

### stringtie
exec_stringtie() {
    # 执行前检查
    # 检查 gff 文件是否存在
    gff_file=$(find ${database_d} -maxdepth 1 -type f -name '*.g*f*' | head -n 1)
    if [[ ! -f ${gff_file} ]]; then
        log ERROR "未找到基因注释文件，gff 或 gtf，优先使用 gtf 文件"
        exit 1
    else
        log INFO "基因注释文件使用 ${gff_file}"
    fi

    log INFO "执行 stringtie 步骤"
    samplename_list=($(tail -n+2 ${samples_described_f} | grep -v '^$' | cut -f 2))
    sample_count=${#samplename_list[@]}

    # 循环 samplename_list 进行 stringtie 处理
    stringtie_psnum_cmd="ps -ef | grep stringtie | grep ${work_dir} | wc -l"
    cut -f 2 ${samples_described_f} | grep -v '^$' | tail -n+2 | while read samplename; do
        # stringtie 单线程，并行跑 
        # 检测当前运行多少 stringtie, 超过配置中的线程数则等待
        stringtie_psnum=$(echo $stringtie_psnum_cmd | bash)
        while [[ $stringtie_psnum -ge $num_threads ]]; do
            sleep 1
            stringtie_psnum=$(echo $stringtie_psnum_cmd | bash)
        done
        log INFO "stringtie processing ${samplename}"
        echo "nohup stringtie -e -B \
            -G ${gff_file} \
            -A ${bam_d}/fpkm/${samplename}_fpkm.txt \
            -o ${bam_d}/ballgown/${samplename}/${samplename}.gtf \
            ${bam_d}/${samplename}.bam > ${log_d}/stringtie_${ref_specie}_${samplename}.log 2>&1 &" | bash
    done

    # 检测 stringtie 还在运行，则等待 stringtie 运行完成
    stringtie_psnum=$(echo $stringtie_psnum_cmd | bash)
    while [[ $stringtie_psnum -gt 0 ]]; do
        sleep 1
        stringtie_psnum=$(echo $stringtie_psnum_cmd | bash)
    done

    log INFO "stringtie 处理完成"
}

### salmon
exec_salmon() {
    
}


###
process_fpkm_reads() {
    cd ${bam_d} || exit
    log INFO "正在合并每个样本的 fpkm 和 reads 文件"
    conda deactivate
    python2 ${script}/getFPKM.py -i ballgown
    python2 ${script}/getTPM.py -i ballgown
    prepDE.py
    conda activate base
    log INFO "对 fpkm 和 reads 文件进行过滤"
    python ${script}/count_filtered.py
    python ${script}/fpkm_filtered.py

    # 列名按照 ${samples_described_f} 重新排序
    log INFO "对 fpkm_matrix 和 reads_matrix 重新排序"
    for file in fpkm_matrix_filtered.txt reads_matrix_filtered.txt gene_count_matrix.txt gene_fpkm_matrix.txt; do
        python ${script}/reorder_genetable_with_samplesdes.py \
            -s ${samples_described_f} \
            -f ${file} \
            -o ${file}
    done
    # 合并 fpkm 和 reads 矩阵
    log INFO "正在合并 fpkm 和 reads 矩阵"

    if [[ ! -f ${annotation_d}/${ref_specie}_kns_gene_def.txt ]]; then
        log ERROR "未找到  ${annotation_d}/${ref_specie}_kns_gene_def.txt 文件，无法进行合并"
        return 1
    fi

    python ${script}/merge_fpkm_reads_matrix.py \
        -f fpkm_matrix_filtered.txt \
        -r reads_matrix_filtered.txt \
        --kns ${annotation_d}/${ref_specie}_kns_gene_def.txt \
        -o fpkm_and_reads_matrix_filtered_data_def.txt
    python ${script}/merge_fpkm_reads_matrix.py \
        -f gene_fpkm_matrix.txt \
        -r gene_count_matrix.txt \
        --kns ${annotation_d}/${ref_specie}_kns_gene_def.txt \
        -o fpkm_and_reads_matrix_data_def.txt
    cd ${work_dir} || exit
}


### jiaofu
jiaofu_prepare() {
    mkdir -p ${jiaofu}/00_Background_materials \
        ${jiaofu}/03_PPI_analysis_KEGG_pathways \
        ${jiaofu}/02_DEG_analysis/Analysis_Funrich_report/Funrich_DEG_analysis_report \
        ${jiaofu}/02_DEG_analysis/Analysis_Funrich_report/DEG_Barplot_graphs \
        ${jiaofu}/02_DEG_analysis/Analysis_Funrich_report/DEG_Bubble_graphs
    
    log INFO "copy files to 00_Reference_genome_annotation_files"
    mycp ${annotation_d}/*_gene_def.txt ${jiaofu}/00_Reference_genome_annotation_files/
    # TODO: copy cds fasta 文件
    
    log INFO "copy files to 01_Original_expression_data"
    mycp ${bam_d}/*_data_def.txt ${jiaofu}/01_Original_expression_data

}

exec_all() {
    exec_fastqc
    hisat2_alignment
    exec_stringtie
    process_fpkm_reads
    exec_multi_deseq
    jiaofu_prepare
}

# 目录变量
script=/home/colddata/qinqiang/ProjectScript/02_Reference
work_dir=$(pwd)
log_d=${work_dir}/log
database_d=${work_dir}/00_Database
annotation_d=${work_dir}/01_Annotation
cleandata_d=${work_dir}/02_Cleandata
bam_d=${work_dir}/03_Bam
mapping_d=${work_dir}/04_Mapping
multideseq_d=${work_dir}/05_Multi_DESeq
jiaofu=${work_dir}/jiaofu_prepare

samples_described_f=${work_dir}/samples_described.txt
compare_info_f=${work_dir}/compare_info.txt

# 激活配置
source $1
source ${script}/transcriptome.sh


log INFO "
####################################
run:\t\t\t${run}
ref_specie:\t\t${ref_specie}
specie:\t\t\t${specie}
threads:\t\t${num_threads}
max_memory:\t\t${max_memory}
rlog_number:\t\t${rlog_number}

word_dir:\t\t${work_dir}
database:\t\t${database_d}
annotation:\t\t${annotation_d}
cleandata:\t\t${cleandata_d}
bam:\t\t\t${bam_d}
mapping:\t\t${mapping_d}
multideseq:\t\t${multideseq_d}
jiaofu:\t\t\t${jiaofu}

samples_described_f:\t${samples_described_f}
compare_info_f:\t\t${compare_info_f}
####################################
"
# 切换到 base 环境下，python 程序都是在 base 环境下编写的
source /home/train/miniconda3/bin/activate base
if [[ ! -d ${log_d} ]]; then
    mkdir ${log_d}
fi
# 预检测工作目录是否包含必须文件支持脚本运行
# pre_check

# 执行流程
if [[ "$(declare -p run 2>/dev/null)" =~ "declare -a" ]]; then
    for item in "${run[@]}"; do
        $item
    done
else
    $run
fi
