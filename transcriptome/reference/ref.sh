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
    samplename_list=($(cut -f 2 samples_described.txt | grep -v sample))
    filename_1_list=($(cut -f 3 samples_described.txt | grep -v filename))
    filename_2_list=($(cut -f 4 samples_described.txt | grep -v filename))
    len=${#filename_1_list[@]}

    for ((i=0; i<len; i++)); do
        echo "${filename_1_list[i]} ${filename_2_list[i]} hisat compare, ${num_threads}"
        hisat2 -x 00_Database/${specie} \
        -p ${num_threads} \
        -I 200 -X 400 --fr \
        --min-intronlen 20 --max-intronlen 4000 \
        -1 ./02_Cleandata/${filename_1_list[i]} \
        -2 ./02_Cleandata/${filename_2_list[i]} \
        2> ./02_Cleandata/mapping/${samplename_list[i]}_mapping.txt | \
        samtools sort -@ ${num_threads} -O BAM -o ./02_Cleandata/bam/${samplename_list[i]}.bam
    done

    python /home/data/command/Reference_transcriptome/mapping_summary_Rerf.py
}

### stringtie
exec_stringtie() {
    gff3_file=$(ls ./00_Database/ | grep gff3)
    # for bam_file in $(ls 02_Cleandata/bam/ | grep '.bam'); do
    #     stringtie_psnum=$(ps -ef | grep stringtie | wc -l)
    #     if [ ${stringtie_psnum} lt ${}]
    #         bam_name=$(echo ${bam_file} | cut -d "." -f 1)
    #         echo "nohup stringtie -e -B \
    #         -G ./00_Database/${gff3_file} \
    #         -A ./02_Cleandata/bam/fpkm/${bam_name}_fpkm.txt \
    #         -o ./02_Cleandata/bam/ballgown/${bam_name}/${bam_name}.gtf \
    #         ./02_Cleandata/bam/${bam_name}.bam \
    #         > ${log}/${bam_name}_stringtie.log 2>&1 &" | bash
    #     fi
    # done
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


max_memory=500
num_threads=60
rlog_number=1

run_program() {
    case "$1" in
    0)
        runall;;
    1)
        hisat_compare;;
    2)
        stringtie;;
    *)
        # LOGE "请输入正确的数字: "
    exit 0
    ;;
esac
}

help_info() {
    echo "帮助信息：
    ————————————————————————————————————————————————————————
    0 执行所有流程

    1 执行 hisat-compare（需指定 --specie, --threads, --max-memory）
    ————————————————————————————————————————————————————————
    使用方法：./noref.sh [参数]
    参数：
        -r, --run <run> 指定运行的流程 (指上面的数字，例如 -r 5)
        -s, --specie <specie> 指定物种
        -c, --threads <threads> 指定线程数
        -m, --max-memory <max_memory> 指定最大内存
        --rlog-number <rlog_number> multi deseq 的倍数
        -h 显示帮助信息
"
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r | --run)
            run="$2"
            echo "run: $2"
            shift 2
            ;;
        -s | --specie)
            specie="$2"
            echo "specie: $2"
            shift 2
            ;;
        -c | --threads)
            num_threads="$2"
            echo "num_threads: $2"
            shift 2
            ;;
        -m | --max-memory)
            max_memory="$2"
            echo "max_memory: $2"
            shift 2
            ;;
        --rlog-number)
            rlog_number="$2"
            echo "rlog_number: $2"
            shift 2
            ;;
        -h)
            help_info
            exit 0
            ;;
        *)
            echo "无效的选项：$1"
            exit 1
            ;;
    esac
done
if [[ -z $run ]]; then
    echo "未指定运行流程，使用 -h 查看帮助"
    exit 0
fi
run_program $run