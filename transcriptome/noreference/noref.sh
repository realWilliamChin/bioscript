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
pre_check() {
    # 如果 compare info 和 samples described 文件不存在则跳过检查
    log INFO "检查 compare_info.txt 和 samples_described.txt 中..."
    if [[ ! -f ${work_dir}/compare_info.txt ]] || [[ ! -f ${work_dir}/samples_described.txt ]]; then
        log ERROR "compare_info.txt 或 samples_described.txt 文件不存在或没有放在工作目录"
        return 1
    fi
    python ${script_d}/check_SampDesAndCompInfo.py
    if [[ $? -ne 0 ]]; then
        log ERROR "compare_info 与 samples_described 不匹配，程序退出"
        exit 1
    fi
}

### fastqc
exec_fastqc() {
    if [[ ! -d ${work_dir}/fastqc ]]; then
        mkdir ${work_dir}/fastqc
    fi
    
    log INFO "正在后台生成 fastqc 报告"
    echo "nohup fastqc ${cleandata_d}/*.fq.* \
    -o ${work_dir}/fastqc > ${log_d}/fastqc.log 2>&1 & " | bash
}

### Trinity pinjie 步骤
exec_pinjie() {
    if [[ ! -d ${cleandata_d} ]]; then
        log ERROR "未找到 ${cleandata_d} 目录，无法开始拼接"
        exit 1
    fi

    if [[ -f ${work_dir}/${assemble_trinity_d}/samples_trinity.txt ]]; then
        log INFO "检测到 samples_trinity.txt 文件"
        read -p "是否使用已存在的 samples_trinity.txt 文件进行拼接？[y/n]" answer
        
    fi
    
    if [[ -f ${work_dir}/samples_described.txt ]]; then
        python ${script_d}/trinity_samples_file.py \
        -i ${cleandata_d} \
        -s ${work_dir}/samples_described.txt \
        -o samples_trinity.txt \
        -t ${trinity_type} >/dev/null
        cat samples_trinity.txt
        log WARNING "请确认 samples_trinity.txt, 回车继续（如需修改，重新开个终端修改完再回车，勿重新启动程序）"
        read -p ""
    else
        log ERROR "未找到 samples_described.txt 文件，无法开始拼接"
        exit 1
    fi
    
    if [[ ! -d ${assemble_trinity_d} ]]; then
        mkdir ${assemble_trinity_d}
    fi
    
    # 生成 trinity.fasta 文件
    log INFO "正在执行 Trinity 拼接流程"
    Trinity --seqType fq \
        --max_memory ${max_memory}G \
        --no_salmon \
        --no_version_check \
        --samples_file samples_trinity.txt \
        --output ${assemble_trinity_d} \
        --CPU ${num_threads} \
        --SS_lib_type RF \
        --normalize_reads \
        --min_contig_length 500
    
    log INFO "Trinity 结束"
    
    if [[ ! -f ${assemble_trinity_d}/Trinity.fasta ]]; then
        log ERROR "Trinity 拼接错误，未生成 Trinity.fasta 文件，程序退出"
        exit 1
    fi
    
    log INFO "已生成 Trinity.fasta"
    log INFO "正在执行 extract_longest_isoforms_from_TrinityFasta.pl"
    perl /home/train/trainingStuff/bin/extract_longest_isoforms_from_TrinityFasta.pl \
    ${assemble_trinity_d}/Trinity.fasta \
    >${assemble_trinity_d}/unigene_longest.fasta
    if [[ ! -f ${assemble_trinity_d}/unigene_longest.fasta ]]; then
        log ERROR "生成 unigene_longest.fasta 错误，程序退出"
        exit 1
    fi
    
    log INFO "cd-hit-est 正在执行"
    cd-hit-est -i ${assemble_trinity_d}/unigene_longest.fasta \
    -o ${assemble_trinity_d}/${specie}_unigene.fasta \
    -c 0.95
    if [[ ! -f ${assemble_trinity_d}/${specie}_unigene.fasta ]]; then
        log ERROR "cd-hit-est 程序错误，未生成 ${specie}_unigene.fasta，程序退出"
        exit 1
    fi
    log INFO "cd-hit-est 已完成，已生成 ${specie}_unigene.fasta"
    
    log INFO "正在生成拼接报告 assemble_stat.txt"
    /opt/biosoft/Trinity-v2.8.5/util/TrinityStats.pl \
    ${assemble_trinity_d}/${specie}_unigene.fasta \
    >${assemble_trinity_d}/assemble_stat.txt && log INFO "ok!"
    seqkit stats ${assemble_trinity_d}/${specie}_unigene.fasta >>${assemble_trinity_d}/assemble_stat.txt
    grep '>' ${assemble_trinity_d}/${specie}_unigene.fasta | cut -d ' ' -f 1 | tr -d '>' >${specie}_all_gene_id.txt
    log INFO "pinjie 流程已完成"
    cd ${assemble_trinity_d}
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
    echo -en "Total_sequence_num\t${Total_sequence_num}\n" >>assemble_stat_report.txt
    echo -en "Total_sequence_bases\t${Total_sequence_bases}\n" >>assemble_stat_report.txt
    echo -en "Percent_GC\t${Percent_GC}\n" >>assemble_stat_report.txt
    echo -en "Largest_transcript\t${Largest_transcript}\n" >>assemble_stat_report.txt
    echo -en "Smallest_transcript\t${Smallest_transcript}\n" >>assemble_stat_report.txt
    echo -en "Average_length\t${Average_length}\n" >>assemble_stat_report.txt
    echo -en "N50\t${N50}" >> assemble_stat_report.txt
}

### transdecoder
transdecoder() {
    log INFO "正在执行 transdecoder 步骤"
    cd ${annotation_d} || exit
    if [[ ! -d transdecoder ]]; then
        mkdir transdecoder
    fi
    cd transdecoder || exit
    TransDecoder.LongOrfs -t ${assemble_trinity_d}/${specie}_unigene.fasta
    TransDecoder.Predict -t ${assemble_trinity_d}/${specie}_unigene.fasta
    cd ${work_dir} || exit
}


### Annotation
## Swiss
exec_swiss() {
    if [[ ! -d ${annotation_d} ]]; then
        mkdir ${annotation_d}
    fi
    log INFO "执行 Annotation - Swiss 步骤"
    /opt/biosoft/ncbi-blast-2.9.0+/bin/blastx \
        -db /home/data/ref_data/Linux_centos_databases/2019_Unprot_databases/swissprot \
        -query ${assemble_trinity_d}/${specie}_unigene.fasta \
        -out ${annotation_d}/${specie}_swiss.blast \
        -max_target_seqs 20 \
        -evalue 1e-5 \
        -num_threads ${num_threads} \
        -outfmt "6 qacc sacc pident qcovs qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    # 如果 swiss 注释成功，则对 swiss blast 结果进行处理
    if [[ -f ${annotation_d}/${specie}_swiss.blast ]]; then
        cd ${annotation_d}
        python ${script_d}/swiss.py
        cd ${work_dir} || exit
    else
        log ERROR "swiss blast 失败"
    fi
}

## nr
exec_nr() {
    log INFO "执行 Annotation - nr 步骤"
    if [[ ! -d ${annotation_d}/temp ]]; then
        mkdir ${annotation_d}/temp
    fi
    diamond blastx --db /home/data/ref_data/db/diamond_nr/diamond_nr \
        --threads ${num_threads} \
        --query ${assemble_trinity_d}/${specie}_unigene.fasta \
        --out ${annotation_d}/${specie}_nr.blast \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --id 30 \
        --block-size 20.0 \
        --tmpdir ${annotation_d}/temp \
        --index-chunks 1
    # 如果 nr 执行成功，则对 nr 结果进行处理
    if [[ -f ${annotation_d}/${specie}_nr.blast ]]; then
        cd ${annotation_d} || exit
        python ${script_d}/nr.py
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
    if [[ ! -d ${annotation_d} ]]; then
        mkdir ${annotation_d}
    fi
    cd ${annotation_d} || exit
    log INFO "正在切换到 eggnog2 conda 环境"
    source /home/train/miniconda3/bin/activate eggnog2
    check_conda_env
    emapper.py -i ${assemble_trinity_d}/${specie}_unigene.fasta \
        --output ${specie} \
        --cpu ${num_threads} \
        --data_dir /opt/biosoft/eggnog5.0.0 \
        --dmnd_db /opt/biosoft/eggnog5.0.0/eggnog_proteins.dmnd \
        -m diamond \
        --translate \
        --itype genome \
        --genepred search

    log INFO "正在切换到 base conda 环境"
    source /home/train/miniconda3/bin/activate base
    check_conda_env
    cd ${work_dir} || exit
}

### kegg
exec_kegg() {
    if [[ ! -d ${annotation_d} ]]; then
        mkdir ${annotation_d}
    fi
    
    log INFO "[KEGG]执行 kegg 注释步骤"
    # 判断 $specie${annotation_d}.fasta 是否大于 50000 条，如果大于 50000 条，则需要分批生成 keg 文件
    if [[ $(grep -c '>' ${assemble_trinity_d}/${specie}_unigene.fasta) -gt 50000 ]]; then
        log INFO "[KEGG]unigene.fasta 文件大于 50000 条，切割进行注释"
        split -l 100000 ${assemble_trinity_d}/${specie}_unigene.fasta ${annotation_d}/${specie}_unigene.fasta_
        for i in $(ls ${annotation_d} | grep "${specie}_unigene.fasta_"); do
            log INFO "[KEGG]正在生成 ${i} 的 keg 文件"
            python ${script_d}/kegg_annotation.py \
                -f ${annotation_d}/${i} \
                -o ${annotation_d}/${i}_keg \
                -l "${kegg_org}" >>${log_d}/kegg.log
            # 检查文件是否生成
            if [[ -f ${annotation_d}/${i}_keg ]]; then
                cat ${annotation_d}/${i}_keg >>${annotation_d}/${specie}.keg
            elif [[ ! -f ${annotation_d}/${i}_keg ]]; then
                log ERROR "[KEGG]${i} 的 keg 文件生成失败"
            fi
        done
    else
        log INFO "[KEGG]unigene.fasta 文件小于 50000 条，直接进行注释"
        python ${script_d}/kegg_annotation.py \
            -f ${assemble_trinity_d}/${specie}_unigene.fasta \
            -o ${annotation_d}/${specie}.keg \
            -l "${kegg_org}"
    fi
    
    # 拿着 $specie_unigene.fasta 去 kegg 网站生成 keg 文件，再做下面的东西 -t plant/animal 动物或植物
    cd ${annotation_d} || exit
    if [[ -f ${work_dir}/${specie}_all_gene_id.txt ]]; then
        python ${script_d}/kegg.py -t ${specie_type} -i ${work_dir}/${specie}_all_gene_id.txt
    else
        log WARNING "[KEGG]未检测到 ${work_dir}/${specie}_all_gene_id.txt 文件，将不会生成 ${specie}_shortname.txt"
        python ${script_d}/kegg.py -t ${specie_type}
    fi
    cd ${work_dir} || exit
}

merge_kns_def() {
    kegg_gene_def=$(realpath -s ${annotation_d}/*KEGG_gene_def.txt)
    nr_gene_def=$(realpath -s ${annotation_d}/*nr_gene_def.txt)
    swiss_gene_def=$(realpath -s ${annotation_d}/*swiss_gene_def.txt)
    python ${script_d}/kns_def_merge.py \
        -k ${kegg_gene_def} \
        -n ${nr_gene_def} \
        -s ${swiss_gene_def} \
        -i ${work_dir}/${specie}_all_gene_id.txt \
        -o ${annotation_d}/${specie}_kns_gene_def.txt
    if [[ ! -f ${annotation_d}/${specie}_kns_gene_def.txt ]]; then
        log INFO "${annotation_d}/${specie}_kns_gene_def.txt 合并失败"
    fi
}

exec_annotation_report() {
    log INFO "正在准备 annotation_report.r 画图所需文件"
    cd ${annotation_d} || return 1
    if [[ ! -d annotation_report ]]; then
        mkdir annotation_report
    fi
    cp ${work_dir}/${specie}_all_gene_id.txt annotation_report/all_gene_id.txt
    cp ${specie}_KEGG_clean.txt annotation_report/KEGG_clean.txt
    cut -f 1 ${specie}_swiss_idNo_def.txt | sort -u > annotation_report/GO_ID.list
    cut -f 1 ${specie}_KEGG_clean.txt | sort -u > annotation_report/KEGG_ID.list
    cut -f 1 ${specie}.emapper.seed_orthologs | grep -v '^#' | sort -u > annotation_report/COG_ID.list
    tail -n +2 ${specie}_swiss_gene_def.txt | cut -f 1 | sort -u > annotation_report/Swiss_ID.list
    tail -n +2 ${specie}_nr_uniq.blast | cut -f 1 | sort -u > annotation_report/NR_ID.list
    tail -n +2 ${specie}_nr_uniq.blast | cut -f 3 > annotation_report/identity.txt
    tail -n +2 ${specie}_nr_uniq.blast | cut -f 11 > annotation_report/evalue.txt
    tail -n +2 ${specie}_nr_uniq.blast \
        | awk -F '[][]' '{print $(NF-1)}' | sort | uniq -c | sort -nr \
        | awk '{ column1 = substr($0, 1, 7); column2 = substr($0, 8); gsub(/^[[:space:]]+|[[:space:]]+$/, "", column1); gsub(/^[[:space:]]+|[[:space:]]+$/, "", column2); print column2"\t"column1 }' \
        | head > annotation_report/species_count.txt
    python ${script_d}/cog_count.py -c ${specie}.emapper.annotations -o annotation_report
    cp *GO*ID.txt annotation_report/
    
    cd ${annotation_d}/annotation_report/ || return 1
    mmv "${specie}_*" "#1"
    
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
    log INFO "检查完成"
    
    if [ $annotation_report_flag -eq 0 ]; then
        log INFO "正在执行 annotation_report.r"
        Rscript ${script_d}/annotation_report.r
        cd ${work_dir}
        return 0
    else
        log ERROR "annotation_report.r 缺少所需文件，将不会执行"
        cd ${work_dir}
        return 1
    fi
}


### multi deseq
### 需要创建一个 compare.txt 和 samples_described.txt
exec_multi_deseq() {
    if [[ ! -d ${multideseq_d} ]]; then
        mkdir ${multideseq_d}
    fi
    
    # 查看工作目录有多少个 compare_info_{group}.txt 文件
    comp_file_list=($(find ${work_dir} -maxdepth 1 -type f -name 'compare_info*' | sort))
    if [ ${#comp_file_list[@]} -eq 0 ]; then
        log ERROR "没有找到 compare_info_{group}.txt 或 compare_info.txt 文件，不执行此步骤"
        return 1
    fi
    
    for comp_file in "${comp_file_list[@]}"; do
        # 从文件名中提取组名
        comp_group=$(echo "${comp_file}" | sed 's/.*compare_info_//' | sed 's/\.txt//')
        
        # 检查文件名是否为 'compare_info.txt'
        if [ "${comp_file##*/}" = 'compare_info.txt' ]; then
            comp_group="default"
        fi
        
        # 设置目录路径
        comp_multideseq="${multideseq_d}/${comp_group}_${rlog_number}"
        
        if [ ! -d "${comp_multideseq}" ]; then
            mkdir -p "${comp_multideseq}"
        else
            log ERROR "${comp_multideseq} 目录已存在，请手动删除后再执行此步骤"
            continue
        fi
        
        # 复制所需文件
        cd ${comp_multideseq} || return 1
        mycp ${mapping_d}/*matrix_filtered.txt ${comp_multideseq}
        cp ${comp_file} ${comp_multideseq}/compare_info.txt
        cut -f 1,2 ${samples_described_f} > ${comp_multideseq}/samples_described.txt
        
        # 挑出指定样本
        log INFO "从 compare_info 中的组名中挑出 sampels_described.txt、reads 和 fpkm 文件的指定样本"
        python ${script_d}/filter_samples_from_comp.py
        python ${script_d}/check_SampDesAndCompInfo.py
        python ${script_d}/reorder_genetable_with_samplesdes.py \
            -f fpkm_matrix_filtered.txt \
            -o fpkm_matrix_filtered.txt
        python ${script_d}/reorder_genetable_with_samplesdes.py \
            -f reads_matrix_filtered.txt \
            -o reads_matrix_filtered.txt
        
        # 执行流程
        log INFO "正在执行 ${comp_group} multi_deseq 流程，Rlog number 为 ${rlog_number}"
        Rscript "${script_d}/multiple_samples_DESeq2.r" "${rlog_number}"
        log INFO "正在给 multideseq 程序生成的文件添加定义"
        kns_def_file=$(find ${annotation_d} -maxdepth 1 -type f -name '*kns_gene_def.txt')
        python "${script_d}/de_results_add_def.py" \
        --kns "${kns_def_file}"
        mv *_DEG_data.txt Analysis/Expression_data/
        cat Analysis/DEG_summary.txt
    done
    cd "${work_dir}"
}


# 执行 Biogrid 流程
exec_biogrid() {
    log INFO "正在执行 Biogrid 步骤"
    if [[ ! -d ${biogrid_d} ]]; then
        mkdir ${biogrid_d}
    fi
    cd ${biogrid_d} || exit
    ${script_d}/biogrid.py \
    -f ${assemble_trinity_d}/${specie}_unigene.fasta \
    -d ${specie_type} \
    -p ${specie} \
    -c ${num_threads}
    if [[ -f ${biogrid_d}/Biogrid_PPI_relation_from_${specie}.txt ]]; then
        log INFO "${biogrid_d}/Biogrid_PPI_relation_from_${specie}.txt 文件已生成"
    else
        log ERROR "Biogrid 错误，文件未生成"
        return 1
    fi
}

### rsem
exec_rsem() {
    if [[ ! -d ${reference_d} ]]; then
        mkdir ${reference_d}
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
        file_1="${cleandata_d}/${filename_1_list[i]}"
        file_2="${cleandata_d}/${filename_2_list[i]}"
        
        if [[ ! -e "$file_1" ]]; then
            log ERROR "未找到文件: $file_1"
            ls -l "$file_1"
            exit 1
        fi
        
        if [[ ! -e "$file_2" ]]; then
            log ERROR "未找到文件: $file_2"
            ls -l "$file_2"
            exit 1
        fi
    done
    log INFO "检查完成，无错误"
    
    # 建库
    if [[ ! -f ${reference_d}/${specie}.transcripts.fa ]]; then
        cd ${reference_d} || exit
        log INFO "正在建库 ${assemble_trinity_d}/${specie}_unigene.fasta ${specie}"
        rsem-prepare-reference --bowtie2 ${assemble_trinity_d}/${specie}_unigene.fasta ${specie}
        cd ${work_dir} || exit
    else
        log INFO "检测到已经建库，跳过建库"
    fi
    
    for ((i=0; i<sample_count; i++)); do
        log INFO "正在比对第${i}个，使用线程数 ${num_threads}，${filename_1_list[i]} ${filename_2_list[i]} ... "
        rsem-calculate-expression -p ${num_threads} \
            --bowtie2 \
            --strandedness reverse \
            --paired-end \
            ${cleandata_d}/${filename_1_list[i]} ${cleandata_d}/${filename_2_list[i]} \
            ${reference_d}/${specie} ${mapping_d}/${samplename_list[i]} \
            2>${mapping_d}/${samplename_list[i]}_mapping.txt \
            1>${mapping_d}/${samplename_list[i]}.log
        
        if [[ $? -ne 0 ]]; then
            log ERROR "比对失败"
        else
            log INFO "${samplename_list[i]} 的 rsem 比对流程已完成"
            tail -n 1 ${mapping_d}/${samplename_list[i]}_mapping.txt
        fi
        
    done
    
    # 读取 samples_described.txt 循环比对
    # tail -n +2 ${work_dir}/samples_described.txt | grep -v '^$' | while read line; do
    #     # sample=$(echo $line | awk '{print $2}')
    #     # sample1=$sample$(ls ${cleandata_d} | grep ${sample} | awk -F"${sample}" '{print $2}' | grep '1')
    #     # sample2=$sample$(ls ${cleandata_d} | grep ${sample} | awk -F"${sample}" '{print $2}' | grep '2')
    #     sample=$(echo $line | cut -f 2)
    #     sample_file1=$(echo $line | cut -f 3)
    #     sample_file2=$(echo $line | cut -f 4)
    #     log INFO "正在执行样本名： $sample 的 rsem 流程 \n$sample_file1 $sample_file2 正在比对"
    #     rsem-calculate-expression -p ${num_threads} \
    #         --bowtie2 \
    #         --strandedness reverse \
    #         --paired-end \
    #         ${cleandata_d}/${filename_1_list[i]} ${cleandata_d}/${filename_2_list[i]} \
    #         ${reference_d}/${specie} ${mapping_d}/${sample} \
    #         2>${mapping_d}/${sample}_mapping.txt \
    #         1>${mapping_d}/${sample}.log
    #     log INFO "$sample 的 rsem 比对流程已完成"
    # done
}

rsem_mapping_summary() {
    echo -e "Sample\tAllReads\tSingleMappedReads\tPairedMappedReads\tUnmappedReads\tAlignRatio(%)" > Mapping_summary.txt
    for i in $(ls *.txt); do
        Sample=$(echo ${i} | cut -d '_' -f 1)
        AllReads=$(head -n 1 ${i} | cut -d ' ' -f 1)
        Sing=$(head -n 11 ${i} | tail -n 1 | cut -d ' ' -f 1)
        tmp_Paire=$(head -n 10 ${i} | tail -n 1 | cut -d ' ' -f 1)
        Paire=$((tmp_Paire / 2))
        UnmappedReads=$((AllReads - Sing - Paire))
        Ratio=$(echo "scale=4; (\(${Sing} + ${Paire}\)) / ${AllReads}" | bc)
        echo -e "${Sample}\t${AllReads}\t${Sing}\t${Paire}\t${UnmappedReads}\t${Ratio}" >> Mapping_summary.txt
    done
}

process_fpkm_reads() {
    log INFO "处理 fpkm 和 reads 矩阵中 ..."
    # 使用 base 的 python3 环境运行 python 脚本
    cd ${mapping_d} || exit
    # 提出 rawreads
    python ${script_d}/all_sample_raw_reads.py
    # 提出 fpkm
    python ${script_d}/all_sample_fpkm.py
    # 过滤掉 rawreads 小于 50 的 gene
    python ${script_d}/count_filtered.py
    # 提出过滤后的基因的 fpkm 值
    python ${script_d}/fpkm_filtered.py
    
    # 列名按照 ${samples_described_f} 重新排序
    for file in fpkm_matrix_filtered.txt reads_matrix_filtered.txt gene_count_matrix.txt gene_fpkm_matrix.txt; do
        log INFO "对 ${file} 重新排序"
        python ${script_d}/reorder_genetable_with_samplesdes.py \
        -s ${samples_described_f} \
        -f ${file} \
        -o ${file}
    done
    
    # 合并 fpkm 和 reads 矩阵
    log INFO "正在合并 fpkm 和 reads 矩阵"
    
    if [[ ! -f ${annotation_d}/${specie}_kns_gene_def.txt ]]; then
        log ERROR "未找到  ${annotation_d}/${specie}_kns_gene_def.txt 文件，无法进行合并"
        return 1
    fi
    
    python ${script_d}/merge_fpkm_reads_matrix.py \
    -f fpkm_matrix_filtered.txt \
    -r reads_matrix_filtered.txt \
    --kns ${annotation_d}/${specie}_kns_gene_def.txt \
    -o fpkm_and_reads_matrix_filtered_data_def.txt
    python ${script_d}/merge_fpkm_reads_matrix.py \
    -f gene_fpkm_matrix.txt \
    -r gene_count_matrix.txt \
    --kns ${annotation_d}/${specie}_kns_gene_def.txt \
    -o fpkm_and_reads_matrix_data_def.txt
    
    cd ${work_dir} || exit
}


exec_annotation() {
    if [[ ! -d ${annotation_d} ]]; then
        mkdir ${annotation_d}
    fi
    if [[ ! -f ${assemble_trinity_d}/${specie}_unigene.fasta ]]; then
        log ERROR "未找到 ${assemble_trinity_d}/${specie}_unigene.fasta 文件，无法开始注释"
        exit 1
    fi
    cp ${assemble_trinity_d}/${specie}_unigene.fasta ${annotation_d}
    exec_kegg &
    exec_swiss
    exec_nr
    exec_cog
    transdecoder
    merge_kns_def
    exec_biogrid
    exec_annotation_report
}

exec_alignment() {
    exec_rsem
    rsem_mapping_summary
    process_fpkm_reads
}

exec_all() {
    pre_check
    exec_pinjie
    exec_annotation
    exec_alignment
}


### 整理交付目录的文件，比较麻烦
jiaofu_prepare() {
    
    mkdir -p ${jiaofu}/00_Background_materials \
    ${jiaofu}/03_PPI_analysis_KEGG_pathways
    
    log INFO "copy files to 00_Funrich_software_def_files"
    mycp ${annotation_d}/*GO*ID.txt ${jiaofu}/00_Funrich_software_def_files/
    mycp ${annotation_d}/*KEGG.txt ${jiaofu}/00_Funrich_software_def_files/
    mycp ${annotation_d}/*shortname.txt ${jiaofu}/00_Funrich_software_def_files/
    mycp ${work_dir}/${specie}_all_gene_id.txt ${jiaofu}/00_Funrich_software_def_files/
    mycp ${biogrid_d}/Biogrid_PPI* ${jiaofu}/00_Funrich_software_def_files/
    
    log INFO "copy files to 00_Reference_genome_annotation_files"
    mycp ${annotation_d}/*_gene_def.txt ${jiaofu}/00_Reference_genome_annotation_files/
    mycp ${assemble_trinity_d}/${specie}_unigene.fasta ${jiaofu}/00_Reference_genome_annotation_files/
    
    # 这个不用了
    # log INFO "copy files to 00_测序质控文件"
    # mycp ${cleandata_d}/fastqc/*.zip ${jiaofu}/00_测序质控文件
    
    log INFO "copy files to 01_Original_expression_data"
    mycp ${mapping_d}/*_data_def.txt ${jiaofu}/01_Original_expression_data
    
    # TODO: 这段儿需要更新
    # 检测多个 compare_info 创建多个组间比较交付目录
    # compare_count=$(ls -i compare_info*.txt | wc -l)
    # if [ $compare_count -gt 1 ]; then
    #     for compare_gourp in $(find ${multideseq_d} -maxdepth1 -type d | tail -n +2); do
    #         log INFO "copy ${compare_group} files to 02_DEG_analysis"
            
    #     done
    # else
    #     log INFO "copy files to 02_DEG_analysis"
    # fi
    
    # Prep_files
    log INFO "\ncopy files to Prep_files"
    mycp ${assemble_trinity_d}/assemble_stat_report.txt ${jiaofu}/Prep_files
    
    # wego
    # log INFO "\ncopy files to wego"
    # cp ${annotation_d}/*GO*ID.txt "$wego"
    # cp ${annotation_d}/*idNo_def.txt "$wego"/all_geneid_GO.txt
    # cp ${multideseq_d}/*Down_ID.txt "$wego"
    # cp ${multideseq_d}/*Up_ID.txt "$wego"
    # cd ${wego} || exit
    # python ${script_d}/wego.py
    # rm *Down_ID.txt *Up_ID.txt *GO*ID.txt
    cd ${word_dir} || exit
}

# 目录变量
script_d=/home/colddata/qinqiang/ProjectScript/01_NoReference
work_dir=$(pwd)
log_d=${work_dir}/log
rawdata_d=${work_dir}/01_Rawdata
cleandata_d=${work_dir}/02_Cleandata
assemble_trinity_d=${work_dir}/03_Assemble_trinity
annotation_d=${work_dir}/04_Annotation
reference_d=${work_dir}/05_Reference
mapping_d=${work_dir}/06_Mapping_mapping
multideseq_d=${work_dir}/07_Multi_DESeq
biogrid_d=${work_dir}/08_Biogrid
jiaofu=${work_dir}/jiaofu_prepare

samples_described_f=${work_dir}/samples_described.txt
compare_info_f=${work_dir}/compare_info.txt

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

word_dir:\t\t${work_dir}
annotation_d:\t\t${annotation_d}
cleandata_d:\t\t${cleandata_d}
bam_d:\t\t\t${bam_d}
mapping_d:\t\t${mapping_d}
multideseq:\t\t${multideseq_d}
jiaofu:\t\t\t${jiaofu}

samples_described_f:\t${samples_described_f}
compare_info_f:\t\t${compare_info_f}
####################################"

# 切换到 base 环境下，python 程序都是在 base 环境下编写的
source /home/train/miniconda3/bin/activate base
if [[ ! -d ${log_d} ]]; then
    mkdir ${log_d}
fi

# 执行流程
if [[ "$(declare -p run 2>/dev/null)" =~ "declare -a" ]]; then
    for item in "${run[@]}"; do
        if [[ $item == "exec_kegg" ]]; then
            exec_kegg &
            continue
        fi
        $item
    done
else
    $run
fi

