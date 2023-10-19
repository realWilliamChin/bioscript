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
  cd ${ref_d} || exit
  log INFO "正在使用 bwa 建立索引"
  bwa index ${genome_fasta_File} -p ${specie}
  cd ${work_dir} || exit

  samplename_list=($(tail -n +2 samples_described.txt))
  filename_1_list=($(tail -n +2 samples_described.txt))
  filename_2_list=($(tail -n +2 samples_described.txt))
  sample_count=${#filename_1_list[@]}

  for ((i=0; i<sample_count; i++)); do

    sample_name=${samplename_list[i]}
    file_1=${cleandata_d}/${filename_1_list[i]}
    file_2=${cleandata_d}/${filename_2_list[i]}

    log INFO "${file_1} ${file_2} bwa compare, ${num_threads}"
    bwa mem \
      -t ${num_threads} \
      -M \
      -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:WES\\tPL:Illumina" \
      ../${ref_d}/${specie} \
      $file_1 $file_2 \
      | samtools sort -@ ${num_threads} -o > ${mapping_d}/${sample_name}.bam
    
    if [[ -f ${mapping_d}/${sample_name} || $? -eq 0 ]]; then
      log INFO "${sample_name} bwa 比对完成. bam 文件已生成"
      return 0
    else
      log ERROR "${sample_name} bwa 比对失败，bam 文件未生成"
      return 1
  done 
}


get_vcf() {
  # 1. 创建元数据
  gatk CreateSequenceDictionary -R ${genome_fasta_file} -O ${genome_fasta_file}.dict
  cd ${work_dir} || exit

  # 2. 生成 bam 文件
  mkdir -p ${vcf_d}/01_vcfs >/dev/null 2>&1
  cd ${mapping_d} || exit
  samplename_list=($(tail -n +2 samples_described.txt))
  sample_count=${#filename_1_list[@]}

  for ((i=0; i<sample_count; i++)); do
    sample_name=${samplename_list[i]}

    gatk MarkDuplicates \
      -I ${sample_name}.bam \
      -O ${sample_name}_markdup.bam \
      -M ${sample_name}_markdup_metrics.txt
    samtools index ${sample_name}_markdup.bam

    # 3. 变异检测
    gatk --java-options -Xmx128G HaplotypeCaller \
      -R ${ref_d}/${genome_fasta_file} \
      -I ${sample_name}_markdup.bam \
      -O ${vcf_d}/01_vcfs/${sample_name}.gvcf.gz \
      --emit-ref-confidence GVCF \
      -stand-call-conf 10 \
      --sample-ploidy ${ploidy_num} \
      --max-genotype-count 8

  cd ${work_dir} || exit
}


call_snp_indel() {
  # 合并
  ls ${vcf_d}/01_vcfs/*.gvcf.gz > ${vcf_d}/gvcfs.list
  gatk CombineGVCFs \
    -R ${genome_fasta_file} \
    -V ${vcf_d}/gvcfs.list \
    -O ${vcf_d}/${specie}_combined_gvcf.gz
  
  # 转换为 vcf
  gatk GenotypeGVCFs \
    -R ${genome_fasta_file} \
    -V ${vcf_d}/${specie}_combined_gvcf.gz \
    -O ${vcf_d}/${specie}_combined.vcf

  # call snp 和 indel
  gatk SelectVariants \
    -R ${genome_fasta_file} \
    -V ${vcf_d}/${specie}_combined.vcf
    -select-type SNP \
    -O ${vcf_d}/${specie}_SNPs.vcf
  gatk SelectVariants \
    -R ${genome_fasta_file} \
    -V ${vcf_d}/${specie}_combined.vcf
    -select-type INDEL \
    -O ${vcf_d}/${specie}_INDELs.vcf

  # 过滤
  gatk VariantFiltration \
    -V ${vcf_d}/${specie}_SNPs.vcf \
    -filter "MQ < 30.0" --filter-name "MQ_filter_SNP" \
    -filter "QD < 2.0" --filter-name "QD_filter_SNP" \
    -O ${vcf_d}/${specie}_SNPs_filter.vcf \
  grep -E "^#|PASS" ${vcf_d}/${specie}_SNPs_filter.vcf > ${vcf_d}/${specie}_SNPs_filter_PASSED.vcf
  gatk VariantFiltration \
    -V ${vcf_d}/${specie}_INDELs.vcf \
    -filter "MQ < 30.0" --filter-name "MQ_filter_INDEL" \
    -filter "SOR > 10.0" --filter-name "SOR_filter_INDEL" \
    -filter "QD < 2.0" --filter-name "QD_filter_INDEL" \
    -filter "FS > 200.000" --filter-name "FS_filter_INDEL" \
    -O ${vcf_d}/${specie}_INDELs_filter.vcf \
  grep -E "^#|PASS" ${vcf_d}/${specie}_INDELs_filter.vcf > ${vcf_d}/${specie}_INDELs_filter_PASSED.vcf

}


cd ${ref_d} || exit
fasta_files=$(*.fna *.fasta)
if [[ ${#files[@]} -eq 0 ]]; then
  log ERROR "没有找到 .fna 或 .fasta 文件"
  exit 1
elif [[ ${#files[@]} -gt 1 ]]; then
  log ERROR "请确保 ${ref_d} 下只有一个 fasta 或 fna 文件"
  exit 1
else
  genome_fasta_file=${ref_d}/${files[0]}
fi

# 目录结构
script=/home/colddata/qinqiang/ProjectScript/04_Reseq
work_dir=$(pwd)
ref_d=${work_dir}/00_Reference
cleandata_d=${work_dir}/01_Cleandata
mapping_d=${work_dir}/02_Mapping
vcf_d=${work_dir}/03_Vcf

run_program() {
    case "$1" in
    0)
        ;;
    1)
        ;;
    2)
        ;;
    3)
        ;;
    *)
        log ERROR "无效的选项：$1"
        exit 0
        ;;
    esac
}

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
