#!/bin/bash
# set -e -o pipefail

GREEN="\e[32m"  # 绿色
YELLOW="\e[33m" # 黄色
RED="\e[31m"    # 红色
RESET="\e[0m"   # 重置颜色

# 定义日志级别
INFO="[INFO]"
WARNING="[WARNING]"
ERROR="[ERROR]"

# 定义日志函数
log() {
  local level="$1"
  local message="$2"
  local timestamp="$(date +'%Y-%m-%d_%H:%M:%S')"
  local colored_message="[$timestamp $level]:$message"
  
  # 根据级别添加颜色
  case "$level" in
    INFO)
      colored_message="${GREEN}${colored_message}${RESET}"
      ;;
    WARNING)
      colored_message="${YELLOW}${colored_message}${RESET}"
      ;;
    ERROR)
      colored_message="${RED}${colored_message}${RESET}"
      ;;
    *)
      # 默认为INFO级别
      colored_message="${GREEN}${colored_message}${RESET}"
      ;;
  esac
  
  echo -e "$colored_message"
}

mycp() {
    # 检查参数数量是否小于2
    if [ "$#" -lt 2 ]; then
        log WARNING "用法: $0 源文件... 目的地目录"
        return 1
    fi

    # 获取最后一个参数作为目的地目录
    dest="${!#}"

    if [ ! -d $dest ]; then
        mkdir -p $dest
    fi

    # 复制源文件到目的地
    for ((i = 1; i < $#; i++)); do
        src="${!i}"
        if [ -e "$src" ]; then
            cp -r "$src" "$dest"
            log INFO "复制 '$src' 到 '$dest'"
        else
            log WARNING "源文件 '$src' 不存在"
        fi
    done
}

# 检查 conda 当前环境
check_conda_env() {
    conda_env=$(conda env list | grep "*" | cut -d " " -f 1)
    log INFO $conda_env
}

# 检查 compare_info.txt 和 samples_described.txt 中名称是否有不匹配的地方
check_compareinfo_and_samplesdescribed() {
    # 如果 compare info 和 samples described 文件不存在则跳过检查
    log INFO "检查 compare_info.txt 和 samples_described.txt 中..."
    if [[ ! -f ${work_dir}/compare_info.txt ]] || [[ ! -f ${work_dir}/samples_described.txt ]]; then
        log ERROR "compare_info.txt 或 samples_described.txt 文件不存在或没有放在工作目录"
        return 1
    fi
    python ${script}/check_SampDesAndCompInfo.py
    if [[ $? -ne 0 ]]; then
        log ERROR "compare_info 与 samples_described 不匹配，程序退出"
        exit 1
    fi
}

### fastqc
exec_fastqc() {
    mkdir ${pinjiedata}/fastqc
    log INFO "正在后台生成 fastqc 报告"
    echo "nohup fastqc ${pinjiedata}/*.fq.* \
    -o ${work_dir}/fastqc > ${log}/fastqc.log 2>&1 & " | bash
}


### 比对
hisat_compare() {
    mkdir ${bam} ${mapping} > /dev/null 2>&1
    samplename_list=($(cut -f 2 samples_described.txt | grep -v sample))
    filename_1_list=($(cut -f 3 samples_described.txt | grep -v filename))
    filename_2_list=($(cut -f 4 samples_described.txt | grep -v filename))
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
    python2 ${script}/getFPKM.py -i ballgown
    python2 ${script}/getTPM.py -i ballgown
    prepDE.py
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

### multi deseq
### 需要创建一个 compare.txt 和 samples_described.txt
exec_multi_deseq() {
    mkdir ${multideseq} >/dev/null 2>&1
    cd ${multideseq} || return 1
    cp ${bam}/reads_matrix_filtered.txt ${multideseq}/
    cp ${bam}/fpkm_matrix_filtered.txt ${multideseq}/
    cut -f 1,2 ${work_dir}/samples_described.txt > ${multideseq}/samples_described.txt
    cp ${work_dir}/compare_info.txt ${multideseq}/

    log INFO "正在执行 multi_deseq 流程，Rlog number 为 ${rlog_number}"
    Rscript ${script}/multiple_samples_DESeq2.r ${rlog_number}
    
    cd ${multideseq} || return 1
    python ${script}/de_results_add_def.py \
        --kns ${annotation}/${specie}_kns_gene_def.txt \
    cat DEG_summary.txt
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

    log INFO "copy files to 02_DEG_analysis"
    mycp ${multideseq}/DEG_summary.txt ${jiaofu}/02_DEG_analysis/Analysis/
    mycp ${multideseq}/*Down_ID.txt ${jiaofu}/02_DEG_analysis/Analysis/
    mycp ${multideseq}/*Up_ID.txt ${jiaofu}/02_DEG_analysis/Analysis/
    mycp ${multideseq}/*_DEG_data.txt ${jiaofu}/02_DEG_analysis/Analysis/Expression_data
    mycp ${multideseq}/*vs*_heatmap.jpeg ${jiaofu}/02_DEG_analysis/Analysis/Expression_data_graphs
    mycp ${multideseq}/*vs*_volcano.jpeg ${jiaofu}/02_DEG_analysis/Analysis/Expression_data_graphs

    log INFO "复制 5 个图到 Expression_data_evaluation"
    mycp ${multideseq}/all_gene_heatmap.jpeg ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
    mycp ${multideseq}/correlation.png ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
    mycp ${multideseq}/fpkm_boxplot.jpeg ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
    mycp ${multideseq}/fpkm_density.jpeg ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/
    mycp ${multideseq}/PCA.jpeg ${jiaofu}/02_DEG_analysis/Expression_data_evaluation/ 

}

run_program() {
    case "$1" in
    0)
        exec_fastqc
        hisat_compare
        exec_stringtie
        process_fpkm_reads
        exec_multi_deseq
        jiaofu_prepare
        ;;
    1)
        exec_fastqc
        ;;
    2)
        hisat_compare
        ;;
    3)
        exec_stringtie
        ;;
    4)
        process_fpkm_reads
        ;;
    5)
        exec_multi_deseq
        ;;
    6)
        jiaofu_prepare
        ;;
    *)
        log ERROR "无效的选项：$1"
        exit 0
        ;;
    esac
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
mkdir ${log} > /dev/null 2>&1

# 如果已经有了组间比较文件，则首先检查组间比较文件和样本文件是否正确，以便后续不出错
if [[ -f ${work_dir}/compare_info.txt && $run -eq 0 ]]; then
    check_compareinfo_and_samplesdescribed
fi

# 执行流程
if [ ${#run[@]}  -gt 1 ]; then
    run_program $run
else
    for item in "${run[@]}"; do
        run_program $item
    done
fi
