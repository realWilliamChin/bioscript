#!/bin/bash

work_dir=$(pwd)
cleandata=${work_dir}/02_Cleandata
analysis=${work_dir}/03_Analysis

specie=$1
num_threads=$2

check_conda_env() {
    conda_env=$(conda env list | grep "*" | cut -d " " -f 1)
    echo $conda_env
}

### merge
# dir format 参考 youhuma /home/colddata/train/pll/rna-seq/youhuma/01_rawdata
# merge 很慢
merge() {
    cd 01_rawdata || exit
    for dir in *; do
        if [[ $dir == *"HM"* ]]; then
            # TODO sample 名字需要核实正确，后面重写
            file_name=${dir:16}

            echo -e "merge processing $file_name"

            # TODO L3_1 也需要改，重写
            $(zcat $dir"/"*L3_1* >$dir"/"$file_name"_clean_1.fq")
            $(zcat $dir"/"*L3_2* >$dir"/"$file_name"_clean_2.fq")
        fi
    done
    cd ../
}

### fastqc
# TODO: 需要判断是否是 gz 格式，有可能是 fq 或者是 fastq 格式
# 这个也很慢，在后台运行就行
fastqc() {
    cd 02_pinjiedata
    mkdir fastqc
    suffix=$(ls | head -n 1 | cut -d '.' -f 2)
    nohup fastqc *.$suffix -o fastqc > fastqc.log 2>&1 &
    cd ../
}
exec_fastqc() {
    echo -e "${GR}正在执行 fastqc 生成质检报告，已放入后台。${PLN}"
    fastqc
    for ((tc=1;tc<=86400;tc++)); do
        # 1 is running, 0 is not running
        status=$(ps -ef | grep fastqc | grep $(ls | head -n 10 | tail -n 1) | wc -l)
        # 1 is completed, 0 is not completed
        complete=$(cat fastqc.log | grep "Analysis complete for" | wc -l)
        if [[ $status -eq 0 ]] && [ $complete -eq 1 ]]; then
            echo -e "${GR}fastqc 报告已完成${PLN}"
            break
        elif [[ $status -eq 0 ]] && [[ $complete -eq 0 ]]; then
            echo -e "${RD}fastqc 可能出错了，不影响后续执行步骤，可以稍后再查看${PLN}"
            break
        else
            continue
        fi
    done &
}
### hisat 建库
hisat_build() {
    cd 00_database || exit
    fasta=$(ls ./*.fasta | wc -l)
    fasta_file=$(ls ./*.fasta)
    if [[ ${fasta} -eq 1 ]]; then
        hisat2-build "${fasta_file}" "$specie" -p 32
    fi
}

### 比对
hisat_compare() {
    mkdir 02_Cleandata/mapping 02_Cleandata/bam;
    for name in $(cut -f 2 samples_described.txt | grep -v sample); do
        echo "${name} hisat compare, ${num_threads}, ${specie}"
        hisat2 -x 00_Database/${specie} -p ${num_threads} -I 200 -X 400 --fr --min-intronlen 20 --max-intronlen 4000 \
        -1 ./02_Cleandata/${name}_1.clean.fq.gz -2 ./02_Cleandata/${name}_2.clean.fq.gz \
        2> ./02_Cleandata/mapping/${name}_mapping.txt | \
        samtools sort -@ ${num_threads} -O BAM -o ./02_Cleandata/bam/${name}.bam
    done
    python /home/data/command/Reference_transcriptome/mapping_summary_Rerf.py
}

### stringtie
exec_stringtie() {
    gff3_file=$(ls ./00_Database/ | grep gff3)
    for bam_file in $(ls 02_Cleandata/bam/ | grep '.bam'); do
        stringtie_psnum=$(ps -ef | grep stringtie | wc -l)
        if [ ${stringtie_psnum} lt ${}]
        bam_name=$(echo ${bam_file} | cut -d "." -f 1)
        echo "nohup stringtie -e -B \
        -G ./00_Database/${gff3_file} \
        -A ./02_Cleandata/bam/fpkm/${bam_name}_fpkm.txt \
        -o ./02_Cleandata/bam/ballgown/${bam_name}/${bam_name}.gtf \
        ./02_Cleandata/bam/${bam_name}.bam \
        > ${log}/${bam_name}_stringtie.log 2>&1 &" | bash
    done
}


###
fpkm_reads() {
    cd 05_Reference
    # 提出 rawreads
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/all_sample_raw_reads.py
    # 提出 fpkm
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/all_sample_fpkm.py
    # 过滤掉 rawreads 小于 50 的 gene
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/count_filtered.py
    # 提出过滤后的基因的 fpkm 值
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/fpkm_filtered.py
    cd ../
}

### multi deseq
multi_deseq() {
    cp ${cleandata}/bam/*filtered.txt ${analysis}/
    cp ${work_dir}/compare_info.txt ${analysis}/
    cp ${work_dir}/samples_described.txt ${analysis}/
    cd ${analysis}
    Rscript /home/data/command/Reference_transcriptome/multiple_samples_DESeq2_V7.r 1
    # 需要核对，数量不对得重新跑
    cd ${work_dir}
}


### Funrich

### Bubble and Bar graph

### WEGO





show_menu() {
    echo -e "
    ${GR}无参执行程序${PLN}
    ${GR}物种名:${PLN} $specie

    ${GR}0${PLN} 执行所有流程
    ————————————————
    ${GR}1${PLN} 执行 hisat_build 建库
    ${GR}2${PLN} 执行 hisat 比对
    ${GR}3${PLN} 执行 stringtie
    ${GR}4${PLN} 执行 pinjie result
    ${GR}5${PLN} 执行 multideseq
    ${GR}6${PLN} 执行 rsem
    ${GR}7${PLN} 执行 fpkm_reads
 "
    # show_status
    echo -e "${YL}"
    read -p "请输入选择(其他任意退出): " select
    echo -e "${PLN}"

    case "${select}" in
    0)
        echo -e "${GR}程序开始执行${PLN}"
        # check_source_data
        exec_fastqc
        exec_pinjie
        exec_annotation
        rsem
        fpkm_reads
        ;;
    1)
        merge
        ;;
    2)
        hisat_compare
        ;;
    3)
        exec_stringtie
        ;;
    4)
        pinjie_result
        ;;
    5)
        multi_deseq
        ;;
    6)
        exec_rsem
        ;;
    *)
        # LOGE "请输入正确的数字: "
        exit 0
        ;;
    esac
}
show_menu

