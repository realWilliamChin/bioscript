#!/bin/bash
# set -e -o pipefail


### 比对
hisat2_alignment() {
    gff_file=$(find ${database_d} -maxdepth 1 -type f -name '*.g*f*' | head -n 1)
    python /home/colddata/qinqiang/script/transcriptome/reference/hisat2_alignment.py \
        hisat2 \   # 这里是 alignment 类型，选择 hisat2 或者 bowtie2
        -s ${samples_described_f} \
        --gff ${gff_file}
        --md ${mapping_d} \
        --cd ${cleandata_d} \
        --bd ${bam_d} \
        --ref ${database_d}/${ref_specie} \
        --cpu ${num_threads}
        
    cat ${mapping_d}/mapping_summary.txt
    cd ${work_dir} || exit
}

bowtie2_alignment() {
    python /home/colddata/qinqiang/script/transcriptome/reference/hisat2_alignment.py \
        bowtie2 \   # 这里是 alignment 类型，选择 hisat2 或者 bowtie2
        -s ${samples_described_f} \
        --md ${mapping_d} \
        --cd ${cleandata_d} \
        --bd ${bam_d} \
        --ref ${database_d}/${ref_specie} \
        --cpu ${num_threads}
        
    cat ${mapping_d}/mapping_summary.txt
    cd ${work_dir} || exit
}


# TODO: 没写完继续写
### salmon
# exec_salmon() {
    
# }


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
    jiaofu_prepare
}

# 目录变量
script_d=/home/colddata/qinqiang/ProjectScript/02_Reference
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
        eval "$item"
    done
else
    eval "$run"
fi