#!/bin/bash
set -e -o pipefail


### 比对
hisat2_alignment() {
    # 创建 ${bam_d} 和 ${mapping_d} 目录
    if [[ ! -d ${bam_d} ]]; then
        mkdir ${bam_d}
    fi
    if [[ ! -d ${mapping_d} ]]; then
        mkdir ${mapping_d}
    fi

    samplename_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f 2))
    filename_1_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f 3))
    filename_2_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f 4))
    sample_count=${#filename_1_list[@]}
    log INFO "即将比对的样本数量是 $sample_count"

    # 检查 filename_1_list 和 filename_2_list 的文件是否存在
    log INFO "检查 samples_described.txt 中的文件是否存在"
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

    for ((i=0; i<sample_count; i++)); do
        log INFO "正在比对第${i}个，使用线程数 ${num_threads}，${filename_1_list[i]} ${filename_2_list[i]} hisat alignment ... "
        hisat2 -x 00_Database/${ref_specie} \
            -p ${num_threads} \
            -I 200 -X 400 --fr \
            --min-intronlen 20 --max-intronlen 4000 \
            -1 ${cleandata_d}/${filename_1_list[i]} \
            -2 ${cleandata_d}/${filename_2_list[i]} \
            2> ${mapping_d}/${samplename_list[i]}_mapping.txt | \
            samtools sort --threads ${num_threads} -O BAM -o - > ${bam_d}/${samplename_list[i]}.bam

        tail -n 1 ${mapping_d}/${samplename_list[i]}_mapping.txt
    done
    cd ${mapping_d} || exit
    python ${script}/hisat2_mapping_summary.py
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

    samplename_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f 2))
    filename_1_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f 3))
    filename_2_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f))
    sample_count=${#filename_1_list[@]}
    log INFO "即将比对的样本数量是 $sample_count"

    # 检查 filename_1_list 和 filename_2_list 的文件是否存在
    log INFO "检查 samples_described.txt 中的文件是否存在"
    file_not_fount=()
    for ((i=0; i<sample_count; i++)); do
        if [[ ! -e "${cleandata_d}/${filename_1_list[i]}" ]]; then
            log ERROR "未找到 ${cleandata_d}/${filename_1_list[i]} 文件"
            file_not_fount+=("${filename_1_list[i]}")
        fi
        if [[ ! -e "${cleandata_d}/${filename_2_list[i]}" ]]; then
            log ERROR "未找到 ${cleandata_d}/${filename_2_list[i]} 文件"
            file_not_fount+=("${filename_2_list[i]}")
        fi
    done

    if [[ ${#file_not_fount[@]} -gt 0 ]]; then
        log ERROR "未找到 ${file_not_fount[@]} 文件"
        exit 1
    fi
    log INFO "检查完成，无错误"

    for ((i=0; i<sample_count; i++)); do
        # 创建单独的 样本 bam 目录，存放没比对上的 fastq 和比对上的 fastq 文件
        if [ ! -d ${bam_d}/${samplename_list[i]} ]; then
            mkdir ${bam_d}/${samplename_list[i]}
        fi
        log INFO "${filename_1_list[i]} ${filename_2_list[i]} bowtie2 alignment ... 使用线程数 ${num_threads}"
        bowtie2 -x ${database_d}/${ref_specie} \
            -p ${num_threads} \
            -1 ${cleandata_d}/${filename_1_list[i]} \
            -2 ${cleandata_d}/${filename_2_list[i]} \
            --un-conc-gz ${bam_d}/${samplename_list[i]}/ \
            --al-conc-gz ${bam_d}/${samplename_list[i]}/ \
            2> ${mapping_d}/${samplename_list[i]}_mapping.txt > ${bam_d}/${samplename_list[i]}.bam

        tail -n 1 ${mapping_d}/${samplename_list[i]}_mapping.txt
    done
    cd ${mapping_d} || exit
    python ${script}/hisat2_mapping_summary.py
    cd ${work_dir} || exit
}

### stringtie
exec_stringtie() {
    log INFO "执行 stringtie 步骤"
    filename_1_list=($(tail -n+2 samples_described.txt | grep -v '^$' | cut -f 3))
    sample_count=${#filename_1_list[@]}

    gff_file=$(ls ${database_d} | grep "\.g.f")
    log INFO "stringtie gff 或 gtf 文件使用 ${gff_file}"

    for bam_file in $(ls ${bam_d} | grep ".bam"); do
        # stringtie 单线程，并行跑 
        # 检测当前运行多少 stringtie, 超过配置中的线程数则等待
        stringtie_psnum=$(ps -ef | grep stringtie | wc -l)
        while [[ $stringtie_psnum -ge $num_threads ]]; do
            sleep 1
            stringtie_psnum=$(ps -ef | grep stringtie | wc -l)
        done

        log INFO "stringtie processing ${bam_file}"
        bam_name=$(basename ${bam_file} .bam)
        echo "nohup stringtie -e -B \
            -G ${database_d}/${gff_file} \
            -A ${bam_d}/fpkm/${bam_name}_fpkm.txt \
            -o ${bam_d}/ballgown/${bam_name}/${bam_name}.gtf \
            ${bam_d}/${bam_file} > ${log_d}/stringtie_${bam_name}.log 2>&1 &" | bash
    done

    sleep 60
    # 检测 stringtie 是否全部生成 fpkm，如果没有全部生成，并检测 stringtie 是否还在运行
    # 如果 stringtie 还在运行，则等待 stringtie 运行完成
    stringtie_fpkm_num=$(ls ${bam_d}/fpkm | grep "_fpkm.txt" | wc -l)
    while [[ $stringtie_fpkm_num -lt ${#filename_1_list[@]} ]]; do
        sleep 1
        stringtie_fpkm_num=$(ls ${bam_d}/fpkm | grep "_fpkm.txt" | wc -l)
        stringtie_psnum=$(ps -ef | grep stringtie | wc -l)
        if [[ $stringtie_psnum -ne $sample_count ]]; then
            log ERROR "stringtie 进程全部结束，但是 fpkm 文件没有全部生成"
            return 1
        fi
    done

    log INFO "stringtie 处理完成"
}


###
process_fpkm_reads() {
    cd ${bam_d} || exit
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
    python ${script}/reorder_genetable_with_samplesdes.py \
        -s ${samples_described_f} \
        -f reads_matrix_filtered.txt \
        -o reads_matrix_filtered.txt
    python ${script}/reorder_genetable_with_samplesdes.py \
        -s ${samples_described_f} \
        -f fpkm_matrix_filtered.txt \
        -o fpkm_matrix_filtered.txt
    python ${script}/reorder_genetable_with_samplesdes.py \
        -s ${samples_described_f} \
        -f gene_count_matrix.txt \
        -o gene_count_matrix.txt
    python ${script}/reorder_genetable_with_samplesdes.py \
        -s ${samples_described_f} \
        -f gene_fpkm_matrix.txt \
        -o gene_fpkm_matrix.txt
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
        -o fpkm_and_reads_data_def.txt
    cd ${work_dir} || exit
}


### jiaofu
jiaofu_prepare() {
    mkdir -p ${jiaofu}/00_Background_materials \
        ${jiaofu}/03_PPI_analysis_KEGG_pathways \
    
    log INFO "copy files to 00_Reference_genome_annotation_files"
    mycp ${annotation_d}/*_gene_def.txt ${jiaofu}/00_Reference_genome_annotation_files/
    # TODO: copy cds fasta 文件
    
    log INFO "copy files to 01_Original_expression_data"
    mycp ${bam_d}/*_data_def.txt ${jiaofu}/01_Original_expression_data

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
            mycp ${multideseq_d}/${compare_group}/DEG_summary.txt     ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/
            mycp ${multideseq_d}/${compare_group}/*Down_ID*.txt        ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/
            mycp ${multideseq_d}/${compare_group}/*Up_ID*.txt          ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/
            mycp ${multideseq_d}/${compare_group}/*_DEG_data.txt      ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/Expression_data
            mycp ${multideseq_d}/${compare_group}/*vs*_heatmap.jpeg   ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/Expression_data_graphs
            mycp ${multideseq_d}/${compare_group}/*vs*_volcano.jpeg   ${jiaofu}/02_DEG_analysis${compare_group}/Analysis/Expression_data_graphs

            log INFO "复制 ${compare_group} 的 5 个图到 Expression_data_evaluation"
            mycp ${multideseq_d}/${compare_group}/all_gene_heatmap.jpeg   ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/
            mycp ${multideseq_d}/${compare_group}/correlation.png         ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/
            mycp ${multideseq_d}/${compare_group}/fpkm_boxplot.jpeg       ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/
            mycp ${multideseq_d}/${compare_group}/fpkm_density.jpeg       ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/
            mycp ${multideseq_d}/${compare_group}/PCA.jpeg                ${jiaofu}/02_DEG_analysis${compare_group}/Expression_data_evaluation/ 
        done
    elif [ ${#comp_file_list[@]} -eq 1 ]; then
        log INFO "copy files to 02_DEG_analysis"
        mycp ${multideseq_d}/DEG_summary.txt      ${jiaofu}/02_DEG_analysis/Analysis/
        mycp ${multideseq_d}/*Down_ID*.txt         ${jiaofu}/02_DEG_analysis/Analysis/
        mycp ${multideseq_d}/*Up_ID*.txt           ${jiaofu}/02_DEG_analysis/Analysis/
        mycp ${multideseq_d}/*_DEG_data.txt       ${jiaofu}/02_DEG_analysis/Analysis/Expression_data
        mycp ${multideseq_d}/*vs*_heatmap.jpeg    ${jiaofu}/02_DEG_analysis/Analysis/Expression_data_graphs
        mycp ${multideseq_d}/*vs*_volcano.jpeg    ${jiaofu}/02_DEG_analysis/Analysis/Expression_data_graphs

        log INFO "复制 5 个图到 Expression_data_evaluation"
        mycp ${multideseq_d}/all_gene_heatmap.jpeg    ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
        mycp ${multideseq_d}/correlation.png          ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
        mycp ${multideseq_d}/fpkm_boxplot.jpeg        ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
        mycp ${multideseq_d}/fpkm_density.jpeg        ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
        mycp ${multideseq_d}/PCA.jpeg                 ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/ 
    fi

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
database_d:\t\t${database_d}
annotation_d:\t\t${annotation_d}
cleandata_d:\t\t${cleandata_d}
bam_d:\t\t\t${bam_d}
mapping_d:\t\t${mapping_d}
multideseq:\t\t${multideseq_d}
jiaofu:\t\t\t${jiaofu}
####################################"
# 切换到 base 环境下，python 程序都是在 base 环境下编写的
source /home/train/miniconda3/bin/activate base

# 预检测工作目录是否包含必须文件支持脚本运行
#pre_check

# 执行流程
# 判断是否是列表，如果是列表，则按照列表内循环运行，列表内不能含有 0
if [ "${#run[@]}" -gt 1 ]; then
    for item in "${run[@]}"; do
        $item
    done
else
    $run
fi
