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
    python /home/colddata/qinqiang/script/transcriptome/annotation/swiss.py \
        -b ${specie}_swiss.blast \
        -p ${specie} \
        --fasta ${annotation_fasta} \
        -t ${num_threads}
}

## nr
exec_nr() {
    log INFO "执行 Annotation - nr 步骤"
    python /home/colddata/qinqiang/script/transcriptome/annotation/nr.py \
        -b ${specie}_nr.blast \
        -p ${specie} \
        --fasta ${annotation_fasta} \
        -t ${num_threads}
}

### kegg
exec_kegg() {
    if [[ ! -d ${annotation_d} ]]; then
        mkdir ${annotation_d}
    fi
    
    log INFO "[KEGG]执行 kegg 注释步骤"
    # 判断 $specie${annotation_d}.fasta 是否大于 50000 条，如果大于 50000 条，则需要分批生成 keg 文件
    unigene_count=$(grep -c '>' ${assemble_trinity_d}/${specie}_unigene.fasta)
    split_count=$((unigene_count / 50000))
    python ${script_d}/kegg_annotation.py \
        -f ${assemble_trinity_d}/${specie}_unigene.fasta \
        -k ${annotation_d}/${specie}.keg \
        -l "${kegg_org}" \
        -s $split_count \
        --allid ${work_dir}/${specie}_all_gene_id.txt
}

merge_kns_def() {
    kegg_gene_def=$(realpath -s ${annotation_d}/*KEGG_gene_def.txt)
    nr_gene_def=$(realpath -s ${annotation_d}/*nr_gene_def.txt)
    swiss_gene_def=$(realpath -s ${annotation_d}/*swiss_gene_def.txt)
    python /home/colddata/qinqiang/script/transcriptome/genedf_add_expression_and_def.py \
        -k ${kegg_gene_def} \
        -n ${nr_gene_def} \
        -s ${swiss_gene_def} \
        -i ${work_dir}/${specie}_all_gene_id.txt \
        --input-header 0 \
        -o ${annotation_d}/${specie}_kns_gene_def.txt
    if [[ ! -f ${annotation_d}/${specie}_kns_gene_def.txt ]]; then
        log INFO "${annotation_d}/${specie}_kns_gene_def.txt 合并失败"
    fi
}

