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


exec_bwa() {
  # 1. 建立索引
  log INFO "正在使用 bwa 建立索引"
  bwa index ${genome_fasta_File} -p ${specie}

  samplename_list=($(tail -n +2 samples_described.txt))
  filename_1_list=($(tail -n +2 samples_described.txt))
  filename_2_list=($(tail -n +2 samples_described.txt))
  sample_count=${#filename_1_list[@]}

  for ((i=0; i<sample_count; i++)); do
    log INFO "${filename_1_list[i]} ${filename_2_list[i]} bwa compare, ${num_threads}"
    bwa mem \
      -t ${num_threads} \
      -M \
      -R "@RG\\tID:${samplename_list[i]}\\tSM:${samplename_list[i]}\\tLB:WES\\tPL:Illumina" \
      ../00_index/${specie} \
      ${filename_1_list[i]} \
      ${filename_2_list[i]} \
      | samtools sort -@ ${num_threads} -o > ${samplename_list[i]}.bam
}


exec_gatk() {

}