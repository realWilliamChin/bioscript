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
    mkdir ${work_dir}/fastqc > /dev/null 2>&1
    log INFO "正在后台生成 fastqc 报告"
    echo "nohup fastqc ${pinjiedata}/*.fq* \
    -o ${work_dir}/fastqc > ${log}/fastqc.log 2>&1 & " | bash
}

### Trinity pinjie 步骤
exec_pinjie() {
    python ${script}/trinity_samples_file.py \
    -i ${pinjiedata} \
    -s ${work_dir}/samples_described.txt \
    -o samples_trinity.txt \
    -t ${trinity_type} > /dev/null
    cat samples_trinity.txt
    log WARNING "请确认 samples_trinity.txt，回车继续（如需修改，重新开个终端修改完再回车，勿重新启动程序）"
    read -p ""

    mkdir ${assemble_trinity} > /dev/null 2>&1
    # 生成 trinity.fasta 文件
    log INFO "正在执行 Trinity 拼接流程"
    Trinity --seqType fq \
        --max_memory ${max_memory}G \
        --no_salmon \
        --no_version_check \
        --samples_file ${specie}_samples_trinity.txt \
        --output ${assemble_trinity} \
        --CPU ${num_threads} \
        --SS_lib_type RF \
        --normalize_reads \
        --min_contig_length 500 \
        && log INFO "Trinity 结束"
    if [[ ! -f ${assemble_trinity}/Trinity.fasta ]]; then
        log ERROR "Trinity 拼接错误，未生成 Trinity.fasta 文件，程序退出"
        exit 1
    fi

    log INFO "已生成 Trinity.fasta"
    log INFO "正在执行 extract_longest_isoforms_from_TrinityFasta.pl"
    perl /home/train/trainingStuff/bin/extract_longest_isoforms_from_TrinityFasta.pl \
        ${assemble_trinity}/Trinity.fasta \
        > ${assemble_trinity}/unigene_longest.fasta
    if [[ ! -f ${assemble_trinity}/unigene_longest.fasta ]]; then
        log ERROR "生成 unigene_longest.fasta 错误，程序退出"
        exit 1
    fi

    log INFO "cd-hit-est 正在执行"
    cd-hit-est -i ${assemble_trinity}/unigene_longest.fasta \
        -o ${assemble_trinity}/${specie}_unigene.fasta \
        -c 0.95
    if [[ ! -f ${assemble_trinity}/${specie}_unigene.fasta ]]; then
        log ERROR "cd-hit-est 程序错误，未生成 ${specie}_unigene.fasta，程序退出"
        exit 1
    fi
    log INFO "cd-hit-est 已完成，已生成 ${specie}_unigene.fasta"

    log INFO "正在生成拼接报告 assemble_stat.txt"
    /opt/biosoft/Trinity-v2.8.5/util/TrinityStats.pl \
        ${assemble_trinity}/${specie}_unigene.fasta \
        > ${assemble_trinity}/assemble_stat.txt && log INFO "ok!"
    seqkit stats ${assemble_trinity}/${specie}_unigene.fasta >> ${assemble_trinity}/assemble_stat.txt
    grep '>' ${assemble_trinity}/${specie}_unigene.fasta | cut -d ' ' -f 1 | tr -d '>' > ${specie}_all_gene_id.txt
    log INFO "pinjie 流程已完成"
    cd ${assemble_trinity}
    assemble_report
    cd ${work_dir} || exit
}
# assemble_stat.txt for report
assemble_report() {
    # 针对 TrinityStats.pl 和 seqkit stats 生成的 assemble_stat.txt 进行处理，可直接插入报告使用
    Total_sequence_num=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $4}' | tr -d ',')
    Total_sequence_bases=$(grep 'Total assembled bases' assemble_stat.txt | head -n 1 | awk -F':' '{print $2}' | tr -d ' ')
    Percent_GC=$(grep 'Percent GC' assemble_stat.txt | awk -F':' '{print $2}' | tr -d ' ')
    Largest_transcript=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $8}' | tr -d ',')
    Smallest_transcript=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $6}' | tr -d ',')
    Average_length=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $7}' | tr -d ',')
    N50=$(grep N50 assemble_stat.txt | head -n 1 | awk -F':' '{print $2}' | tr -d ' ')
    echo -en "Total_sequence_num\t${Total_sequence_num}\n" >> assemble_stat_report.txt
    echo -en "Total_sequence_bases\t${Total_sequence_bases}\n" >> assemble_stat_report.txt
    echo -en "Percent_GC\t${Percent_GC}\n" >> assemble_stat_report.txt
    echo -en "Largest_transcript\t${Largest_transcript}\n" >> assemble_stat_report.txt
    echo -en "Smallest_transcript\t${Smallest_transcript}\n" >> assemble_stat_report.txt
    echo -en "Average_length\t${Average_length}\n" >> assemble_stat_report.txt
    echo -en "N50\t${N50}" >> assemble_stat_report.txt
}

### Annotation
## Swiss
exec_swiss() {
    mkdir ${annotation} > /dev/null 2>&1
    log INFO "执行 Annotation - Swiss 步骤"
    /opt/biosoft/ncbi-blast-2.9.0+/bin/blastx \
        -db /home/data/ref_data/Linux_centos_databases/2019_Unprot_databases/swissprot \
        -query ${assemble_trinity}/${specie}_unigene.fasta \
        -out ${annotation}/${specie}_unigene_swiss.blast \
        -max_target_seqs 20 \
        -evalue 1e-5 \
        -num_threads ${num_threads} \
        -outfmt "6 qacc sacc pident qcovs qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle" 
    # 如果 swiss 注释成功，则对 swiss blast 结果进行处理
    if [[ -f ${annotation}/${specie}_unigene_swiss.blast ]]; then
        cd ${annotation}
        python ${script}/swiss.py
        cd ${work_dir} || exit
    else
        log ERROR "swiss blast 失败"
    fi
}

## nr
exec_nr() {
    log INFO "执行 Annotation - nr 步骤"
    mkdir -p ${annotation}/temp > /dev/null 2>&1
    diamond blastx --db /home/data/ref_data/db/diamond_nr/diamond_nr \
        --threads ${num_threads} \
        --query ${assemble_trinity}/${specie}_unigene.fasta \
        --out ${annotation}/${specie}_unigene_nr.blast \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --id 30 \
        --block-size 20.0 \
        --tmpdir ${annotation}/temp \
        --index-chunks 1
    # 如果 nr 执行成功，则对 nr 结果进行处理
    if [[ -f ${annotation}/${specie}_unigene_nr.blast ]]; then
        cd ${annotation} || exit
        python ${script}/nr.py
        cd ${work_dir} || exit
        log INFO "nr 注释完成"
    else
        log ERROR "nr 注释失败"
    fi
}

### cog
# emappey.py youhuma 45分钟，运行完需要 conda deactivate
exec_cog() {
    log INFO "正在执行 cog 步骤"
    mkdir ${annotation} > /dev/null 2>&1
    cd ${annotation} || exit
    log INFO "正在切换到 python27 conda 环境"
    source /home/train/miniconda3/bin/activate python27
    check_conda_env
    emapper.py -i ${assemble_trinity}/${specie}_unigene.fasta \
        --output ${specie}_unigene \
        --cpu ${num_threads} \
        --data_dir /home/data/xuezhen/download/eggNOG/eggnog.db \
        -m diamond --translate
    log INFO "正在切换到 base conda 环境"
    source /home/train/miniconda3/bin/activate base
    check_conda_env
    cd ${work_dir} || exit
}

### kegg
exec_kegg() {
    log INFO "[KEGG]执行 kegg 注释步骤"
    # 判断 $specie_unigene.fasta 是否大于 50000 条，如果大于 50000 条，则需要分批生成 keg 文件
    if [[ $(grep -c '>' ${assemble_trinity}/${specie}_unigene.fasta) -gt 50000 ]]; then
        log INFO "[KEGG]unigene.fasta 文件大于 50000 条，切割进行注释"
        split -l 100000 ${assemble_trinity}/${specie}_unigene.fasta ${annotation}/${specie}_unigene.fasta_
        for i in $(ls ${annotation} | grep "${specie}_unigene.fasta_"); do
            log INFO "[KEGG]正在生成 ${i} 的 keg 文件"
            python ${script}/kegg_annotation.py \
            -f ${annotation}/${i} \
            -o ${annotation}/${i}_keg \
            -l "${kegg_org}" >> ${log}/kegg.log
            # 检查文件是否生成
            if [[ -f ${annotation}/${i}_keg ]]; then
                cat ${annotation}/${i}_keg >> ${annotation}/${specie}_unigene.keg
            elif [[ ! -f ${annotation}/${i}_keg ]]; then
                log ERROR "[KEGG]${i} 的 keg 文件生成失败"
            fi
        done
    else
        log INFO "[KEGG]unigene.fasta 文件小于 50000 条，直接进行注释"
        python ${script}/kegg_annotation.py \
            -f ${assemble_trinity}/${specie}_unigene.fasta \
            -o ${annotation}/${specie}_unigene.keg \
            -l "${kegg_org}"
    fi

    # 拿着 $specie_unigene.fasta 去 kegg 网站生成 keg 文件，再做下面的东西 -t plant/animal 动物或植物
    cd ${annotation} || exit
    if [[ -f ${work_dir}/${specie}_all_gene_id.txt ]]; then
        python ${script}/kegg.py -t ${specie_type} -i ${work_dir}/${specie}_all_gene_id.txt
    else
        log WARNING "[KEGG]未检测到 ${work_dir}/${specie}_all_gene_id.txt 文件，将不会生成 ${specie}_shortname.txt"
        python ${script}/kegg.py -t ${specie_type}
    fi
    cd ${work_dir} || exit
}

### transdecoder
transdecoder() {
    log INFO "执行 transdecoder 步骤"
    cd ${annotation} || exit
    mkdir transdecoder >/dev/null 2>&1
    cd transdecoder || exit
    TransDecoder.LongOrfs -t ${assemble_trinity}/${specie}_unigene.fasta
    TransDecoder.Predict -t ${assemble_trinity}/${specie}_unigene.fasta
    cd ${work_dir} || exit
}


merge_kns_def() {
    kegg_gene_def=$(realpath -s  ${annotation}/*KEGG_gene_def.txt)
    nr_gene_def=$(realpath -s ${annotation}/*nr_gene_def.txt)
    swiss_gene_def=$(realpath -s ${annotation}/*swiss_gene_def.txt)
    python ${script}/kns_def_merge.py \
        -k ${kegg_gene_def}
        -n ${nr_gene_def}
        -s ${swiss_gene_def}
        -i ${work_dir}/${specie}_all_gene_id.txt
        -o ${annotation}/${specie}_kns_gene_def.txt
    if [[ ! -f ${annotation}/${specie}_kns_gene_def.txt ]]; then
        log INFO "${annotation}/${specie}_kns_gene_def.txt 合并失败"
    fi
}


exec_annotation_report() {
    log INFO "正在准备 annotation_report.r 画图所需文件"
    cd ${annotation} || return 1
    mkdir annotation_report > /dev/null 2>&1
    cp ${work_dir}/${specie}_all_gene_id.txt annotation_report/all_gene_id.txt
    cp ${specie}_unigene_KEGG_clean.txt annotation_report/KEGG_clean.txt
    cut -f 1 ${specie}_unigene_swiss_idNo_def.txt | sort -u >  annotation_report/GO_ID.list
    cut -f 1 ${specie}_unigene_KEGG_clean.txt | sort -u > annotation_report/KEGG_ID.list
    cut -f 1 ${specie}_unigene.emapper.seed_orthologs | grep -v '^#' | sort -u > annotation_report/COG_ID.list
    tail -n +2 ${specie}_unigene_swiss_gene_def.txt | cut -f 1 | sort -u > annotation_report/Swiss_ID.list
    tail -n +2 ${specie}_unigene_nr_uniq.blast | cut -f 1 | sort -u > annotation_report/NR_ID.list
    tail -n +2 ${specie}_unigene_nr_uniq.blast | cut -f 3  > annotation_report/identity.txt
    tail -n +2 ${specie}_unigene_nr_uniq.blast | cut -f 11 > annotation_report/evalue.txt
    tail -n +2 tiannanxing_unigene_nr_uniq.blast |  awk -F '[][]' '{print $(NF-1)}' | sort | uniq -c | sort -nr |\
        awk -v OFS='\t' '{print $2, $1}' | head > annotation_report/species_count.txt
    python ${script}/cog_count.py -c ${specie}_unigene.emapper.annotations -o annotation_report
    cp *GO*ID.txt annotation_report/
    
    cd ${annotation}/annotation_report/ || return 1
    mmv "${specie}_unigene_*" "#1"

    log INFO "检查 annotation_report.r 画图所需文件是否全部生成"
    # annotation_report.r 画图所需的全部文件
    file_list=(
        "all_gene_id.txt" "KEGG_clean.txt" "GO_ID.list" "KEGG_ID.list" "COG_ID.list"
        "Swiss_ID.list" "NR_ID.list" "identity.txt" "evalue.txt" "species_count.txt"
        "COG_count.txt" "swiss_GO_BP_ID.txt" "swiss_GO_CC_ID.txt" swiss_GO_MF_ID.txt
        )

    annotation_report_flag=0
    for file in $(ls); do
        # 检查文件是否在列表中
        if [[ " ${file_list[@]} " =~ " $file " ]]; then
            continue
        else
            log ERROR "$file 不存在."
            annotation_report_flag=1
        fi
    done

    if [ $annotation_report_flag -eq 0 ]; then
        log INFO "正在执行 annotation_report.r"
        Rscript ${script}/annotation_report.r
        cd ${work_dir}
        return 0
    else
        log ERROR "annotation_report.r 缺少所需文件，将不会执行"
        cd ${work_dir}
        return 1
    fi
}

# 执行 Biogrid 流程
exec_biogrid() {
    log INFO "正在执行 Biogrid 步骤"
    mkdir ${biogrid}
    cd ${biogrid} || exit
    ${script}/biogrid.py \
        -f ${assemble_trinity}/${specie}_unigene.fasta \
        -d ${specie_type} \
        -p ${specie} \
        -c ${num_threads}
    if [[ -f ${biogrid}/Biogrid_PPI_relation_from_${specie}.txt ]]; then
        log INFO "${biogrid}/Biogrid_PPI_relation_from_${specie}.txt 文件已生成"
    else
        log ERROR "Biogrid 错误，文件未生成"
        return 1
    fi
}


exec_annotation() {
    mkdir ${annotation} > /dev/null 2>&1
    cp ${assemble_trinity}/${specie}_unigene.fasta ${annotation}
    exec_kegg &
    exec_swiss
    exec_nr
    exec_cog
    transdecoder
    merge_kns_def
    exec_biogrid
    exec_annotation_report
}

### rsem
exec_rsem() {
    mkdir ${reference}
    # 建库
    cd ${reference} || exit
    log INFO "正在建库 ${assemble_trinity}/${specie}_unigene.fasta ${specie}"
    rsem-prepare-reference --bowtie2 ${assemble_trinity}/${specie}_unigene.fasta ${specie}
    cd ${work_dir} || exit
    mkdir ${mapping}
    # 读取 samples_described.txt 循环比对
    tail -n +2 ${work_dir}/samples_described.txt | grep -v '^$' | while read line; do
        sample=$(echo $line | awk '{print $2}')
        sample1=$sample$(ls ${pinjiedata} | grep ${sample} | awk -F"${sample}" '{print $2}' | grep '1')
        sample2=$sample$(ls ${pinjiedata} | grep ${sample} | awk -F"${sample}" '{print $2}' | grep '2')
        log INFO "正在执行 $sample 的 rsem 流程, $sample1 $sample2 正在比对"
        rsem-calculate-expression -p ${num_threads} \
            --bowtie2 \
            --strandedness reverse \
            --paired-end \
            ${pinjiedata}/${sample1} ${pinjiedata}/${sample2} \
            ${reference}/${specie} ${mapping}/${sample} \
            2>${mapping}/${sample}_mapping.txt \
            1>${mapping}/${sample}.log
        log INFO "$sample 的 rsem 比对流程已完成"
    done
}

process_fpkm_reads(){
    log INFO "处理 fpkm 和 reads 矩阵中 ..."
    # 使用 base 的 python3 环境运行 python 脚本
    cd ${mapping} || exit
    # 提出 rawreads
    python ${script}/all_sample_raw_reads.py
    # 提出 fpkm
    python ${script}/all_sample_fpkm.py
    # 过滤掉 rawreads 小于 50 的 gene
    python ${script}/count_filtered.py
    # 提出过滤后的基因的 fpkm 值
    python ${script}/fpkm_filtered.py

    log INFO "正在对 fpkm 和 reads 矩阵的列重新排序"
    # 列名按照 samples_described.txt 重新排序
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
    cp ${mapping}/reads_matrix_filtered.txt ${multideseq}/
    cp ${mapping}/fpkm_matrix_filtered.txt ${multideseq}/
    cut -f 1,2 ${work_dir}/samples_described.txt > ${multideseq}/samples_described.txt
    cp ${work_dir}/compare_info.txt ${multideseq}/

    log INFO "正在执行 multi_deseq 流程，Rlog number 为 ${rlog_number}"
    Rscript ${script}/multiple_samples_DESeq2.r ${rlog_number}

    cd ${multideseq} || return 1
    python ${script}/de_results_add_def.py \
        --kns ${annotation}/${specie}_kns_gene_def.txt
    cat DEG_summary.txt
    cd ${work_dir} || exit
}


### 整理交付目录的文件，比较麻烦
jiaofu_prepare() {

    mkdir -p ${jiaofu}/00_Background_materials \
        ${jiaofu}/03_PPI_analysis_KEGG_pathways \

    log INFO "copy files to 00_Funrich_software_def_files"
    mycp ${annotation}/*GO*ID.txt ${jiaofu}/00_Funrich_software_def_files/
    mycp ${annotation}/*KEGG.txt ${jiaofu}/00_Funrich_software_def_files/
    mycp ${annotation}/*shortname.txt ${jiaofu}/00_Funrich_software_def_files/
    mycp ${work_dir}/${specie}_all_gene_id.txt ${jiaofu}/00_Funrich_software_def_files/
    mycp ${biogrid}/Biogrid_PPI* ${jiaofu}/00_Funrich_software_def_files/

    log INFO "copy files to 00_Reference_genome_annotation_files"
    mycp ${annotation}/*_gene_def.txt ${jiaofu}/00_Reference_genome_annotation_files/
    mycp ${assemble_trinity}/${specie}_unigene.fasta ${jiaofu}/00_Reference_genome_annotation_files/

    log INFO "copy files to 00_测序质控文件"
    mycp ${pinjiedata}/fastqc/*.zip ${jiaofu}/00_测序质控文件

    log INFO "copy files to 01_Original_expression_data"
    mycp ${mapping}/*_data_def.txt ${jiaofu}/01_Original_expression_data

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

    # Prep_files
    log INFO "\ncopy files to Prep_files"
    mycp ${assemble_trinity}/assemble_stat_report.txt ${jiaofu}/Prep_files

    # wego
    # log INFO "\ncopy files to wego"
    # cp ${annotation}/*GO*ID.txt "$wego"
    # cp ${annotation}/*idNo_def.txt "$wego"/all_geneid_GO.txt
    # cp ${multideseq}/*Down_ID.txt "$wego"
    # cp ${multideseq}/*Up_ID.txt "$wego"
    # cd ${wego} || exit
    # python ${script}/wego.py
    # rm *Down_ID.txt *Up_ID.txt *GO*ID.txt
    cd ${word_dir} || exit
}


run_program() {
    case "$1" in
    0)
        log INFO "程序开始执行"
        # check_source_data
        exec_fastqc
        exec_pinjie
        exec_annotation
        exec_rsem
        exec_multi_deseq
        ;;
    1)
        exec_fastqc
        ;;
    2)
        exec_pinjie
        ;;
    3)
        exec_annotation
        ;;
    3.1)
        exec_swiss
        ;;
    3.2)
        exec_nr
        ;;
    3.3)
        exec_cog
        ;;
    3.4)
        exec_kegg &
        ;;
    3.5)
        transdecoder
        ;;
    3.6)
        exec_biogrid
        ;;
    3.7)
        exec_annotation_report
        ;;
    4)
        exec_rsem
        ;;
    5)
        process_fpkm_reads
        ;;
    6)
        check_compareinfo_and_samplesdescribed
        exec_multi_deseq
        ;;
    7)
        jiaofu_prepare
        ;;
    *)
        echo "无效的选项：$1"
        exit 0
        ;;
    esac
}

# 目录变量
script=/home/colddata/qinqiang/ProjectScript/01_NoReference
work_dir=$(pwd)
log=${work_dir}/log
rawdata=${work_dir}/01_Rawdata
pinjiedata=${work_dir}/02_Pinjiedata
assemble_trinity=${work_dir}/03_Assemble_trinity
annotation=${work_dir}/04_Annotation
reference=${work_dir}/05_Reference
mapping=${work_dir}/06_Mapping_mapping
multideseq=${work_dir}/07_Multi_DESeq
biogrid=${work_dir}/08_Biogrid
jiaofu=${work_dir}/jiaofu_prepare

# 激活配置文件
source $1

log INFO "
####################################
程序开始执行：初始变量如下
run:${run} 
specie:${specie}
num_threads:${num_threads}
max_memory:${max_memory}
trinity_type:${trinity_type}
specie_type:${specie_type}
kegg_org:${kegg_org}
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
# 判断是否是列表，如果是列表，则按照列表内循环运行，列表内不能含有 0
if [ ${#run[@]} -gt 1 ]; then
    run_program $run
else
    for item in "${run[@]}"; do
        run_program $item
    done
fi
