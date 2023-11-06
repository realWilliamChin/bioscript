#!/bin/bash
set -e -o pipefail


### 比对
hisat_compare() {
    mkdir ${bam} ${mapping} > /dev/null 2>&1
    samplename_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f 2))
    filename_1_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f 3))
    filename_2_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f 4))
    sample_count=${#filename_1_list[@]}

    for ((i=0; i<sample_count; i++)); do
        log INFO "${filename_1_list[i]} ${filename_2_list[i]} hisat compare, ${num_threads}"
        hisat2 -x 00_Database/${specie} \
        -p ${num_threads} \
        -I 200 -X 400 --fr \
        --min-intronlen 20 --max-intronlen 4000 \
        -1 ${cleandata}/${filename_1_list[i]} \
        -2 ${cleandata}/${filename_2_list[i]} \
        2> ${mapping}/${samplename_list[i]}_mapping.txt | \
        samtools sort --threads ${num_threads} -O BAM -o ${bam}/${samplename_list[i]}.bam

        tail -n 1 ${mapping}/${samplename_list[i]}_mapping.txt
    done
    cd ${mapping} || exit
    python ${script}/hisat2_mapping_summary.py
    cd ${work_dir} || exit
}

### stringtie
exec_stringtie() {
    log INFO "执行 stringtie 步骤"
    mkdir ${log}/bam >/dev/null 2>&1
    gff_file=$(ls ${database} | grep "\.g.f")
    log INFO "stringtie gff 或 gtf 文件使用 ${gff_file}"

    for bam_file in $(ls ${bam} | grep ".bam"); do
        # stringtie 单线程，并行跑 
        # 检测当前运行多少 stringtie, 超过配置中的线程数则等待
        stringtie_psnum=$(ps -ef | grep stringtie | wc -l)
        try_count=1
        while [ $try_count -lt 65500 ]; do
            stringtie_psnum=$(ps -ef | grep stringtie | grep -v "grep" | wc -l)
            if [[ $stringtie_psnum -lt $num_threads ]]; then
                break
            else
                try_count=$(expr $try_count + 1)
                sleep 1
                continue
            fi
        done

        log INFO "stringtie processing ${bam_file}"
        bam_name=$(basename ${bam_file} .bam)
        echo "nohup stringtie -e -B \
            -G ${database}/${gff_file} \
            -A ${bam}/fpkm/${bam_name}_fpkm.txt \
            -o ${bam}/ballgown/${bam_name}/${bam_name}.gtf \
            ${bam}/${bam_file} > ${log}/bam/${bam_name}.log 2>&1 &" | bash
    done
    
    # 检测运行结束
    try_count=1
    while [[ $try_count -lt 65500 ]]; do
        stringtie_psnum=$(ps -ef | grep stringtie | grep -v "grep" | wc -l)
        if [[ $stringtie_psnum -eq 0 ]]; then
            log INFO "stringtie 已执行完成"
            return 0
        else
            try_count=$(expr $try_count + 1)
            sleep 1
        fi
    done
}


###
process_fpkm_reads() {
    cd ${bam} || exit
    log INFO "正在生成 fpkm 和 reads 文件"
    conda deactivate
    python2 ${script}/getFPKM.py -i ballgown
    python2 ${script}/getTPM.py -i ballgown
    prepDE.py
    conda activate base
    log INFO "对 fpkm 和 reads 文件进行过滤"
    python ${script}/count_filtered.py
    python ${script}/fpkm_filtered.py

    # 列名按照 samples_described.txt 重新排序
    log INFO "对 fpkm_matrix 和 reads_matrix 重新排序"
    python ${script}/realignment_fpkm_columns.py -r \
        -s ${work_dir}/samples_described.txt \
        -f reads_matrix_filtered.txt
    python ${script}/realignment_fpkm_columns.py -r \
        -s ${work_dir}/samples_described.txt \
        -f fpkm_matrix_filtered.txt
    python ${script}/realignment_fpkm_columns.py -r \
        -s ${work_dir}/samples_described.txt \
        -f gene_count_matrix.txt
    python ${script}/realignment_fpkm_columns.py -r \
        -s ${work_dir}/samples_described.txt \
        -f gene_fpkm_matrix.txt
    # 合并 fpkm 和 reads 矩阵
    log INFO "正在合并 fpkm 和 reads 矩阵"

    if [[ ! -f ${annotation}/${specie}_kns_gene_def.txt ]]; then
        log ERROR "未找到  ${annotation}/${specie}_kns_gene_def.txt 文件，无法进行合并"
        return 1
    fi

    python ${script}/merge_fpkm_reads_matrix.py \
        -f fpkm_matrix_filtered.txt \
        -r reads_matrix_filtered.txt \
        --kns ${annotation}/${specie}_kns_gene_def.txt \
        -o fpkm_and_reads_matrix_filtered_data_def.txt
    python ${script}/merge_fpkm_reads_matrix.py \
        -f gene_fpkm_matrix.txt \
        -r gene_count_matrix.txt \
        --kns ${annotation}/${specie}_kns_gene_def.txt \
        -o fpkm_and_reads_data_def.txt
    cd ${work_dir} || exit
}


### jiaofu
jiaofu_prepare() {
    mkdir -p ${jiaofu}/00_Background_materials \
        ${jiaofu}/03_PPI_analysis_KEGG_pathways \
    
    log INFO "copy files to 00_Reference_genome_annotation_files"
    mycp ${annotation}/*_gene_def.txt ${jiaofu}/00_Reference_genome_annotation_files/
    # TODO: copy cds fasta 文件
    
    log INFO "copy files to 01_Original_expression_data"
    mycp ${bam}/*_data_def.txt ${jiaofu}/01_Original_expression_data

    # 检测多个 compare_info 创建多个组间比较交付目录
    comp_file_list=($(find ${work_dir} -maxdepth 1 -type f -name 'compare_info*' | sort))
    if [ ${#comp_file_list[@]} -eq 0 ]; then
        log ERROR "没有找到 compare_info_{group}.txt 或 compare_info_default 文件，不执行此步骤"
        return 1
    elif [ ${#comp_file_list[@]} -gt 1 ]; then
        # 循环所有 compare_info.txt
        for comp_file in "${comp_file_list[@]}"; do
            compare_group=$(echo "${comp_file}" | sed 's/.*compare_info//' | sed 's/\.txt//')
            log INFO "copy ${compare_group} files to 02_DEG_analysis"
            mycp ${multideseq}/${compare_group}/DEG_summary.txt     ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/
            mycp ${multideseq}/${compare_group}/*Down_ID*.txt        ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/
            mycp ${multideseq}/${compare_group}/*Up_ID*.txt          ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/
            mycp ${multideseq}/${compare_group}/*_DEG_data.txt      ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/Expression_data
            mycp ${multideseq}/${compare_group}/*vs*_heatmap.jpeg   ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/Expression_data_graphs
            mycp ${multideseq}/${compare_group}/*vs*_volcano.jpeg   ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/Expression_data_graphs

            log INFO "复制 ${compare_group} 的 5 个图到 Expression_data_evaluation"
            mycp ${multideseq}/${compare_group}/all_gene_heatmap.jpeg   ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/
            mycp ${multideseq}/${compare_group}/correlation.png         ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/
            mycp ${multideseq}/${compare_group}/fpkm_boxplot.jpeg       ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/
            mycp ${multideseq}/${compare_group}/fpkm_density.jpeg       ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/
            mycp ${multideseq}/${compare_group}/PCA.jpeg                ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/ 
        done
    elif [ ${#comp_file_list[@]} -eq 1 ]; then
        log INFO "copy files to 02_DEG_analysis"
        mycp ${multideseq}/DEG_summary.txt      ${jiaofu}/02_DEG_analysis/Analysis/
        mycp ${multideseq}/*Down_ID*.txt         ${jiaofu}/02_DEG_analysis/Analysis/
        mycp ${multideseq}/*Up_ID*.txt           ${jiaofu}/02_DEG_analysis/Analysis/
        mycp ${multideseq}/*_DEG_data.txt       ${jiaofu}/02_DEG_analysis/Analysis/Expression_data
        mycp ${multideseq}/*vs*_heatmap.jpeg    ${jiaofu}/02_DEG_analysis/Analysis/Expression_data_graphs
        mycp ${multideseq}/*vs*_volcano.jpeg    ${jiaofu}/02_DEG_analysis/Analysis/Expression_data_graphs

        log INFO "复制 5 个图到 Expression_data_evaluation"
        mycp ${multideseq}/all_gene_heatmap.jpeg    ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
        mycp ${multideseq}/correlation.png          ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
        mycp ${multideseq}/fpkm_boxplot.jpeg        ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
        mycp ${multideseq}/fpkm_density.jpeg        ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
        mycp ${multideseq}/PCA.jpeg                 ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/ 
    fi

}

exec_all() {
    exec_fastqc
    hisat_compare
    exec_stringtie
    process_fpkm_reads
    exec_multi_deseq
    jiaofu_prepare
}

# 目录变量
script=/home/colddata/qinqiang/ProjectScript/02_Reference
work_dir=$(pwd)
log=${work_dir}/log
database=${work_dir}/00_Database
annotation=${work_dir}/01_Annotation
cleandata=${work_dir}/02_Cleandata
bam=${work_dir}/03_Bam
mapping=${work_dir}/04_Mapping
multideseq=${work_dir}/05_Multi_DESeq
biogrid=${work_dir}/06_Biogrid
jiaofu=${work_dir}/jiaofu_prepare

# 激活配置
source ${script}/transcriptome.sh
source $1

log INFO "
####################################
run:${run} 
specie:${specie}
threads:${num_threads}
max_memory:${max_memory}
specie_type:${specie_type}
rlog_number:${rlog_number}
####################################"
# 切换到 base 环境下，python 程序都是在 base 环境下编写的
source /home/train/miniconda3/bin/activate base

# 预检测工作目录是否包含必须文件支持脚本运行
#pre_check

if [[ -f ${work_dir}/compare_info.txt && "$run" == "0" ]]; then
    check_compareinfo_and_samplesdescribed
fi

# 执行流程
# 判断是否是列表，如果是列表，则按照列表内循环运行，列表内不能含有 0
if [ "${#run[@]}" -gt 1 ]; then
    for item in "${run[@]}"; do
        $item
    done
else
    $run
fi
