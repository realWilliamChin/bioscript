#!/bin/bash
set -e -o pipefail

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
    python ${script}/process_gff.py -p ${specie}
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
    cd ${hisat_database} || exit
    genomic_fna=$(ls ${reference_genome_fs} | grep "genomic.fna")
    hisat2-build ${reference_genome_fs}/${genomic_fna} ${specie} -p ${num_threads}
    if [[ $? -eq 0 ]]; then
        log INFO "建库完成"
    fi
}

### annotation
## Swiss
exec_swiss() {
    mkdir -p ${prep_fs}/swiss >/dev/null 2>&1
    log INFO "执行 Annotation - Swiss 步骤"
    /opt/biosoft/ncbi-blast-2.9.0+/bin/blastx \
        -db /home/data/ref_data/Linux_centos_databases/2019_Unprot_databases/swissprot \
        -query ${reference_genome_fs}/${specie}_cds.fasta \
        -out ${prep_fs}/swiss/${specie}_swiss.blast \
        -max_target_seqs 20 \
        -evalue 1e-5 \
        -num_threads ${num_threads} \
        -outfmt "6 qacc sacc pident qcovs qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    # 如果 swiss 注释成功，则对 swiss blast 结果进行处理
    if [[ -f ${prep_fs}/swiss/${specie}_swiss.blast ]]; then
        cd ${prep_fs}/swiss
        python ${script}/swiss.py
        cd ${work_dir} || exit
    else
        log ERROR "swiss blast 失败"
    fi
}

## nr
exec_nr() {
    log INFO "执行 Annotation - nr 步骤"
    mkdir -p ${prep_fs}/nr/temp >/dev/null 2>&1
    diamond blastx --db /home/data/ref_data/db/diamond_nr/diamond_nr \
        --threads ${num_threads} \
        --query ${reference_genome_fs}/${specie}_cds.fasta \
        --out ${prep_fs}/nr/${specie}_nr.blast \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --id 30 \
        --block-size 20.0 \
        --tmpdir ${prep_fs}/nr/temp \
        --index-chunks 1
    # 如果 nr 执行成功，则对 nr 结果进行处理
    if [[ -f ${prep_fs}/nr/${specie}_nr.blast || $? -eq 0 ]]; then
        cd ${prep_fs}/nr || exit
        gff_file=$(ls ${reference_genome_fs} | grep ".g.f")
        python ${script}/nr.py && python ${script}/nr_gff.py -g ${reference_genome_fs}/${gff_file}
        if [[ $? -eq 0 ]]; then
            cd ${work_dir} || exit
            log INFO "nr 注释完成"
            return 0
        else
            log ERROR "nr 注释失败"
            return 1
        fi
    else
        log ERROR "nr 注释失败"
        return 1
    fi
}

### kegg
exec_kegg() {
    mkdir -p ${prep_fs}/kegg >/dev/null 2>&1
    log INFO "[KEGG]执行 Annotation - kegg 步骤"
    # 判断 $specie_cds.fasta 是否大于 50000 条，如果大于 50000 条，则需要分批生成 keg 文件
    if [[ $(grep -c '>' ${reference_genome_fs}/${specie}_cds.fasta) -gt 50000 ]]; then
        log INFO "[KEGG]unigene.fasta 文件大于 50000 条，切割进行注释"
        split -l 100000 ${reference_genome_fs}/${specie}_cds.fasta ${prep_fs}/kegg/${specie}_unigene.fasta_
        for i in $(ls ${prep_fs}/kegg/ | grep "${specie}_cds.fasta_"); do
            log INFO "[KEGG]正在生成 ${i} 的 keg 文件"
            python ${script}/kegg_prep_fs.py \
                -f ${prep_fs}/kegg/${i} \
                -o ${prep_fs}/kegg/${i}_keg \
                -l "${kegg_org}" >>${prep_fs}/kegg/kegg.log
            # 检查文件是否生成
            if [[ -f ${prep_fs}/kegg/${i}_keg ]]; then
                cat ${prep_fs}/kegg/${i}_keg >>${prep_fs}/kegg/${specie}.keg
            elif [[ ! -f ${prep_fs}/kegg/${i}_keg ]]; then
                log ERROR "[KEGG]${i} 的 keg 文件生成失败"
            fi
        done
    else
        log INFO "[KEGG]unigene.fasta 文件小于 50000 条，直接进行注释"
        python ${script}/kegg_annotation.py \
            -f ${reference_genome_fs}/${specie}_cds.fasta \
            -o ${prep_fs}/kegg/${specie}.keg \
            -l "${kegg_org}"
    fi

    # 拿着 $specie_cds.fasta 去 kegg 网站生成 keg 文件，再做下面的东西 -t plant/animal 动物或植物
    cd ${prep_fs}/kegg || exit
    if [[ -f ${reference_genome_fs}/${specie}_all_gene_id.txt ]]; then
        python ${script}/kegg.py -t ${specie_type} -i ${reference_genome_fs}/${specie}_all_gene_id.txt
    else
        log WARNING "[KEGG]未检测到 ${reference_genome_fs}/${specie}_all_gene_id.txt 文件，将不会生成 ${specie}_shortname.txt"
        python ${script}/kegg.py -t ${specie_type}
    fi
    cd ${work_dir} || exit
}

merge_kns_def() {
    kegg_gene_def=$(realpath -s  ${prep_fs}/kegg/*KEGG_gene_def.txt)
    nr_gene_def=$(realpath -s ${prep_fs}/nr/*nr_gene_def.txt)
    swiss_gene_def=$(realpath -s ${prep_fs}/swiss/*swiss_gene_def.txt)
    log INFO "正在合并 kegg swiss nr 基因注释，使用${kegg_gene_def}, ${nr_gene_def}, ${swiss_gene_def}"
    python ${script}/kns_def_merge.py \
        -k ${kegg_gene_def} \
        -n ${nr_gene_def} \
        -s ${swiss_gene_def} \
        -i ${reference_genome_fs}/${specie}_all_gene_id.txt \
        -o ${prep_fs}/${specie}_kns_gene_def.txt
    if [[ ! -f ${prep_fs}/${specie}_kns_gene_def.txt || $? -ne 0 ]]; then
        log INFO "${prep_fs}/${specie}_kns_gene_def.txt 合并失败"
        return 1
    else
        log INFO "合并成功，${prep_fs}/${specie}_kegg_gene_def.txt"
        return 0
    fi
}

# 执行 Biogrid 流程
exec_biogrid() {
    log INFO "正在执行 Biogrid 步骤"
    mkdir ${biogrid_fs} > /dev/null 2>&1
    cd ${biogrid_fs} || exit
    python ${script}/biogrid.py \
        -f ${reference_genome_fs}/${specie}_cds.fasta \
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

run_program() {
    case "$1" in
    0)
        gen_allgeneid
        #    exec_kegg &
        exec_swiss
        exec_nr
        exec_biogrid
        exec_hisat2_build
        ;;
    1)
        exec_nr
        ;;
    2)
        exec_biogrid
        exec_gmt
        ;;
    3)
        merge_kns_def
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
hisat_database=${work_dir}/00_Database
gene_annotation_fs=${work_dir}/01_Gene_annotation_files
funrich_def_fs=${work_dir}/02_Funrich_def_files
gsea_gmt_fs=${work_dir}/03_GSEA_GMT_files
biogrid_fs=${work_dir}/Biogrid
prep_fs=${work_dir}/Prep_files
reference_genome_fs=${work_dir}/Reference_genome_data

# 激活配置
source $1

log INFO "
####################################
run:${run} 
specie:${specie}
threads:${num_threads}
specie_type:${specie_type}
kegg_org:${kegg_org}
####################################"
# 切换到 base 环境下，python 程序都是在 base 环境下编写的
source /home/train/miniconda3/bin/activate base

# 执行流程
# 判断是否是列表，如果是列表，则按照列表内循环运行，列表内不能含有 0
if [ "${#run[@]}" -gt 1 ]; then
    for item in "${run[@]}"; do
        run_program $item
    done
else
    run_program $run
fi
