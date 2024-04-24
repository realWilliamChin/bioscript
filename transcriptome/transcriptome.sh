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
    python ${script}/check_SampDesAndCompInfo.py
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
    echo "nohup fastqc ${cleandata}/*.fq.* \
    -o ${work_dir}/fastqc > ${log_d}/fastqc.log 2>&1 & " | bash
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
        python ${script}/swiss.py
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
    if [[ ! -d ${annotation_d} ]]; then
        mkdir ${annotation_d}
    fi
    cd ${annotation_d} || exit
    log INFO "正在切换到 python27 conda 环境"
    source /home/train/miniconda3/bin/activate python27
    check_conda_env
    emapper.py -i ${assemble_trinity_d}/${specie}_unigene.fasta \
        --output ${specie} \
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
            python ${script}/kegg_annotation.py \
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
        python ${script}/kegg_annotation.py \
            -f ${assemble_trinity_d}/${specie}_unigene.fasta \
            -o ${annotation_d}/${specie}.keg \
            -l "${kegg_org}"
    fi

    # 拿着 $specie_unigene.fasta 去 kegg 网站生成 keg 文件，再做下面的东西 -t plant/animal 动物或植物
    cd ${annotation_d} || exit
    if [[ -f ${work_dir}/${specie}_all_gene_id.txt ]]; then
        python ${script}/kegg.py -t ${specie_type} -i ${work_dir}/${specie}_all_gene_id.txt
    else
        log WARNING "[KEGG]未检测到 ${work_dir}/${specie}_all_gene_id.txt 文件，将不会生成 ${specie}_shortname.txt"
        python ${script}/kegg.py -t ${specie_type}
    fi
    cd ${work_dir} || exit
}

merge_kns_def() {
    kegg_gene_def=$(realpath -s ${annotation_d}/*KEGG_gene_def.txt)
    nr_gene_def=$(realpath -s ${annotation_d}/*nr_gene_def.txt)
    swiss_gene_def=$(realpath -s ${annotation_d}/*swiss_gene_def.txt)
    python ${script}/kns_def_merge.py \
        -k ${kegg_gene_def} \
        -n ${nr_gene_def} \
        -s ${swiss_gene_def} \
        -i ${work_dir}/${specie}_all_gene_id.txt \
        -o ${annotation_d}/${specie}_kns_gene_def.txt
    if [[ ! -f ${annotation_d}/${specie}_kns_gene_def.txt ]]; then
        log INFO "${annotation_d}/${specie}_kns_gene_def.txt 合并失败"
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
        # else
        #     log WARNING "${comp_multideseq} 目录已存在"
        #     read -p "是否删除 ${comp_multideseq} 目录 Y/N" choice
        #     if ${choice}
        fi
        
        # 复制所需文件
        cd ${comp_multideseq} || return 1
        mycp "${bam_d}"/*matrix_filtered.txt ${comp_multideseq}
        cp ${comp_file} ${comp_multideseq}/compare_info.txt
        cut -f 1,2 ${samples_described_f} > ${comp_multideseq}/samples_described.txt

        # 挑出指定样本
        log INFO "从 compare_info 中的组名中挑出 sampels_described.txt、reads 和 fpkm 文件的指定样本"
        python ${script}/filter_samples_from_comp.py
        python ${script}/check_SampDesAndCompInfo.py
        python ${script}/reorder_genetable_with_samplesdes.py \
            -f fpkm_matrix_filtered.txt \
            -o fpkm_matrix_filtered.txt
        python ${script}/reorder_genetable_with_samplesdes.py \
            -f reads_matrix_filtered.txt \
            -o reads_matrix_filtered.txt

        # 执行流程
        log INFO "正在执行 ${comp_group} multi_deseq 流程，Rlog number 为 ${rlog_number}"
        /opt/biosoft/R-4.2.2/bin/Rscript "${script}/multiple_samples_DESeq2.r" "${rlog_number}"
        if [[ $? -ne 0 ]]; then
            log ERROR "multi_deseq 流程执行失败"
            return 1
        fi
        log INFO "正在给 multideseq 程序生成的文件添加定义"
        kns_def_file=$(find ${annotation_d} -maxdepth 1 -type f -name '*kns_gene_def.txt')
        python "${script}/de_results_add_def.py" \
            --kns "${kns_def_file}"
        mv *_DEG_data.txt Analysis/Expression_data/
        cat Analysis/DEG_summary.txt
    done
    cd "${work_dir}"
}

