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

gen_allgeneid() {
    log INFO "正在生成 ${specie}_all_gene_id.txt"
    cd ${reference_genome_fs} || exit
    python ${script}/process_gff.py \
        -p ${specie} \
        -g ${gene_description_file}
    if [[ -f ${reference_genome_fs}/${specie}_all_gene_id.txt ]]; then
        log INFO "${specie}_all_gene_id.txt 已生成"
        return 0
    else
        log ERROR "${specie}_all_geneid.txt 未生成"
        return 1
    fi
}

exec_hisat2_build() {
    log INFO "正在进行 hisat2-build 建库"

    if [[ ! -d ${hisat_database} ]]; then
        mkdir ${hisat_database}
    fi
    
    hisat2-build ${all_gene_sequence_file} ${hisat_database}/${specie} -p ${num_threads}

    if [[ $? -eq 0 ]]; then
        log INFO "建库完成"
    fi

}

### annotation
## Swiss
exec_swiss() {
    log INFO "执行 Annotation - Swiss 步骤"
    if [[ ! -d ${prep_fs}/swiss ]]; then
        mkdir -p ${prep_fs}/swiss
    fi

    python /home/colddata/qinqiang/script/transcriptome/annotation/swiss.py \
        -b ${prep_fs}/swiss/${specie}_swiss.blast \
        -p ${prep_fs}/swiss/${specie} \
        --fasta ${cds_gene_sequence_file} \
        -t ${num_threads}
}

## nr
exec_nr() {
    log INFO "执行 Annotation - nr 步骤"
    if [[ ! -d ${prep_fs}/nr/temp ]]; then
        mkdir -p ${prep_fs}/nr/temp
    fi

    python /home/colddata/qinqiang/script/transcriptome/annotation/nr.py \
        -b ${prep_fs}/nr/${specie}_nr.blast \
        --basicinfo ${reference_genome_fs}/${specie}_gene_basicinfo.txt \
        -p ${prep_fs}/nr/${specie} \
        --fasta ${cds_gene_sequence_file} \
        -t ${num_threads}
}

### kegg
exec_kegg() {
    log INFO "[KEGG]执行 kegg 注释步骤"
    annotation_d=${prep_fs}/kegg
    if [[ ! -d ${prep_fs}/kegg ]]; then
        mkdir -p ${prep_fs}/kegg
    fi
    
    gene_count=$(grep -c '>' ${cds_gene_sequence_file})
    split_count=$((gene_count / 50000))
    python ${script_d}/kegg_annotation.py \
        -f ${cds_gene_sequence_file} \
        -k ${annotation_d}/${specie}.keg \
        -l "${kegg_org}" \
        -s $split_count \
        --allid ${reference_genome_fs}/${specie}_all_gene_id.txt
}

merge_kns_def() {
    if [[ ! -f ${reference_genome_fs}/${specie}_all_gene_id.txt ]]; then
        log ERROR "未找到 ${prep_fs}/${specie}_all_gene_id.txt 文件，无法进行合并"
        return 1
    fi
    kegg_gene_def=$(realpath -s  ${prep_fs}/kegg/*KEGG_gene_def.txt)
    nr_gene_def=$(realpath -s ${prep_fs}/nr/*nr_gene_def.txt)
    swiss_gene_def=$(realpath -s ${prep_fs}/swiss/*swiss_gene_def.txt)
    log INFO "正在合并 kegg swiss nr 基因注释，使用${kegg_gene_def}, ${nr_gene_def}, ${swiss_gene_def}"
    python /home/colddata/qinqiang/script/transcriptome/genedf_add_expression_and_def.py \
        -k ${kegg_gene_def} \
        -n ${nr_gene_def} \
        -s ${swiss_gene_def} \
        -i ${reference_genome_fs}/${specie}_all_gene_id.txt \
        --input-header 0 \
        -o ${prep_fs}/${specie}_kns_gene_def.txt
    if [[ ! -f ${prep_fs}/${specie}_kns_gene_def.txt || $? -ne 0 ]]; then
        log ERROR "${prep_fs}/${specie}_kns_gene_def.txt 合并失败"
        return 1
    else
        log INFO "合并成功，${prep_fs}/${specie}_kns_gene_def.txt"
        return 0
    fi
}

# 执行 Biogrid 流程
exec_biogrid() {
    log INFO "正在执行 Biogrid 步骤"
    mkdir ${biogrid_fs} > /dev/null 2>&1
    cd ${biogrid_fs} || exit
    python ${script}/biogrid.py \
        -f ${cds_gene_sequence_file} \
        -d ${specie_type} \
        -p ${specie} \
        -c ${num_threads}
    if [[ -f ${biogrid_fs}/Biogrid_PPI_relation_from_${specie}.txt || $? -eq 0 ]]; then
        log INFO "${biogrid_fs}/Biogrid_PPI_relation_from_${specie}.txt 文件已生成"
        return 0
    else
        log ERROR "Biogrid 错误，文件未生成"
        return 1
    fi
}

exec_gmt() {
    log INFO "正在生成 .GMT 文件"
    mkdir ${gsea_gmt_fs} > /dev/null 2>&1
    cd ${gsea_gmt_fs} || exit
    mycp ${prep_fs}/swiss/*GO* ${gsea_gmt_fs}/
    mycp ${prep_fs}/kegg/*KEGG.txt ${gsea_gmt_fs}/
    mycp ${prep_fs}/kegg/*_tier* ${gsea_gmt_fs}/
    tmp_a=$(ls | grep KEGG.txt)
    tmp_b=$(ls | grep tier2 | tr '2' '3')
    mv $tmp_a $tmp_b
    python ${script}/gsea.py
    if [[ $? -eq 0 ]]; then
        cd ${gsea_gmt_fs} || exit
        rm *txt
        cd ${work_dir}
        log INFO "GMT 文件生成结束"
    fi
}


cp_files() {
    if [[ ! -d ${gene_annotation_fs} ]] ; then
        mkdir ${gene_annotation_fs}
    fi
    if [[ ! -d ${funrich_def_fs} ]] ; then
        mkdir ${funrich_def_fs}
    fi
    mycp ${prep_fs}/nr/*def.txt ${gene_annotation_fs}
    mycp ${prep_fs}/swiss/*_gene_def.txt ${gene_annotation_fs}
    mycp ${prep_fs}/kegg/*_gene_def.txt ${gene_annotation_fs}

    mycp ${prep_fs}/kegg/*KEGG.txt ${funrich_def_fs}
    mycp ${prep_fs}/swiss/*GO* ${funrich_def_fs}
    mycp ${reference_genome_fs}/*_all_gene_id.txt ${funrich_def_fs}
    mycp ${prep_fs}/kegg/*shortname.txt ${funrich_def_fs}
    mycp ${biogrid_fs}/Biogrid_PPI* ${funrich_def_fs}
}

exec_all() {
    gen_allgeneid
    exec_hisat2_build
    exec_kegg &
    exec_swiss
    exec_nr
    merge_kns_def
    exec_biogrid
    exec_gmt
    cp_files
}

exec_annotation() {
    exec_kegg &
    exec_swiss
    exec_nr
    merge_kns_def
    exec_biogrid
    exec_gmt
}

# 目录变量
script=/home/colddata/qinqiang/ProjectScript/02_Reference
work_dir=$(pwd)
hisat_database=${work_dir}/00_Database
gene_annotation_fs=${work_dir}/01_Gene_annotation_files
funrich_def_fs=${work_dir}/02_Funrich_def_files
gsea_gmt_fs=${work_dir}/03_GSEA_GMT_files
prep_fs=${work_dir}/Prep_files
biogrid_fs=${prep_fs}/Biogrid
reference_genome_fs=${work_dir}/Reference_genome_data
log=${work_dir}/log

# 激活配置
source $1

log INFO "
####################################
run:${run} 
specie:${specie}
threads:${num_threads}
specie_type:${specie_type}
kegg_org:${kegg_org}
gene_description_file:${gene_description}
all_gene_sequence_file:${all_gene_sequnce}
cds_gene_sequence_file:${cds_gene_sequnce}
####################################"
# 切换到 base 环境下，python 程序都是在 base 环境下编写的
source /home/train/miniconda3/bin/activate base
gene_description_file=${reference_genome_fs}/${gene_description}
all_gene_sequence_file=${reference_genome_fs}/${all_gene_sequnce}
cds_gene_sequence_file=${reference_genome_fs}/${cds_gene_sequnce}

# 创建一个 log 目录
if [[ ! -d ${log} ]]; then
    mkdir ${log}
fi

# 执行流程
if [[ "$(declare -p run 2>/dev/null)" =~ "declare -a" ]]; then
    for item in "${run[@]}"; do
        $item
    done
else
    $run
fi
