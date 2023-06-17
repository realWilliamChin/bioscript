#!/bin/bash
# set -e -o pipefail

# 输入参数
specie=$1
num_threads=$2
max_memory=$3

work_dir=$(pwd)
rawdata=${work_dir}/01_Rawdata
pinjiedata=${work_dir}/02_Pinjiedata
assemble_trinity=${work_dir}/03_Assemble_trinity
annotation=${work_dir}/04_Annotation
reference=${work_dir}/05_Reference
mapping=${work_dir}/06_Mapping_mapping
multideseq=${work_dir}/07_Multi_deseq
funrich=${work_dir}/08_Funrich

unigene_fasta=${assemble_trinity}/${specie}_unigene.fasta

check_conda_env() {
    conda_env=$(conda env list | grep "*" | cut -d " " -f 1)
    echo $conda_env
}

### fastqc
exec_fastqc() {
    run "mkdir ${pinjiedata}/fastqc" "正在后台生成 fastqc 报告"
    # 检查 pinjiedata 目录下都是 fq 文件还是文件夹，如果是文件夹，则循环文件夹下执行 fastqc
    # 如果是文件，则直接执行 fastqc
    if [[ $(ls ${pinjiedata} | grep ".fq" | wc -l) -gt 1 ]]; then
        exec "nohup fastqc ${pinjiedata}/*.fq \
        -o ${pinjiedata}/fastqc > ${work_dir}/fastqc.log 2>&1 && echo "fastqc 已完成" & " | bash 
    elif [[ $(ls ${pinjiedata} | grep ".fq" | wc -l) -lt 1 ]]; then
        for i in $(ls ${pinjiedata}); do
            run "nohup fastqc ${pinjiedata}/${i}/*.fq \
            --o ${pinjiedata}/fastqc > ${work_dir}/fastqc.log 2>&1 && echo "fastqc 已完成" & " | bash 
        done
    fi
}

### pinjie 步骤
exec_pinjie() {
    # 生成 sample trinity 文件, 可能需要修改一下 sample trinity，太多了拼接太慢了
    source /home/train/miniconda3/bin/activate base
    python /home/colddata/qinqiang/script/noreference/trinity_samples_file.py
    cat ${work_dir}/samples_trinity.txt
    echo "##############################################\n"
    read -p "是否更改 samples_trinity.txt，回车继续"
    mv ${work_dir}/samples_trinity.txt ${work_dir}/${specie}_samples_trinity.txt

    mkdir ${assemble_trinity}
    # 生成 trinity.fasta 文件
    echo "正在执行pinjie流程"
    Trinity --seqType fq --max_memory ${max_memory}G \
    --no_salmon --no_version_check \
    --samples_file ${specie}_samples_trinity.txt \
    --output ${assemble_trinity} \
    --CPU ${num_threads} --SS_lib_type RF --normalize_reads \
    --min_contig_length 500
    wait
    # 1s 执行完成
    perl /home/train/trainingStuff/bin/extract_longest_isoforms_from_TrinityFasta.pl \
    ${assemble_trinity}/Trinity.fasta > ${assemble_trinity}/unigene_longest.fasta

    # 大概 1 分钟执行完成
    cd-hit-est -i ${assemble_trinity}/unigene_longest.fasta \
    -o ${unigene_fasta} -c 0.95 && echo "cd-hit-est 已完成，已生成 ${specie}_unigene.fasta"

    #TODO 后续可以把 assemble_stat.txt 选取需要的自动插入到报告中
    /opt/biosoft/Trinity-v2.8.5/util/TrinityStats.pl \
    ${unigene_fasta} > ${assemble_trinity}/assemble_stat.txt
    seqkit stats ${unigene_fasta} >> ${assemble_trinity}/assemble_stat.txt
    echo "pinjie流程已完成"
}

### Annotation
### Swiss
# Swiss 注释挺慢的，youhuma 运行这个运行了小 3 个小时，运行没有日志生成，直接就完成了
py_swiss() {
    cd ${annotation} || exit
    # TODO 上一条命令运行完才能运行下面一句。这个也挺慢的，二十来分钟吧
    echo "python /home/data/command/swiss/swiss_to_go_new_v4.py &" | bash
    cd ${work_dir} || exit
}
exec_swiss() {
    mkdir ${annotation}
    echo "执行 Annotation - Swiss 步骤"
    /opt/biosoft/ncbi-blast-2.9.0+/bin/blastx \
        -db /home/data/ref_data/Linux_centos_databases/2019_Unprot_databases/swissprot \
        -query ${unigene_fasta} -out ${annotation}/${specie}_unigene_swiss.blast \
        -max_target_seqs 20 -evalue 1e-5 -num_threads ${num_threads} \
        -outfmt "6 qacc sacc pident qcovs qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        && py_swiss; echo "swiss 注释完成"
}

## nr
# diamond blastx 非常非常慢，7个小时左右
py_nr() {
cd ${annotation} || exit
    python /home/data/command/nr/nr_v2.py
    cd ${work_dir} || exit
}
exec_nr() {
    echo "执行 Annotation - nr 步骤"
    mkdir -p ${annotation}/temp
    diamond blastx --db /home/data/ref_data/db/diamond_nr/diamond_nr \
        --query ${unigene_fasta} --out ${annotation}/${specie}_unigene_nr_diamond.blast \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --sensitive --max-target-seqs 1 --evalue 1e-5 --id 30 --block-size 20.0 \
        --tmpdir ${annotation}/temp --index-chunks 1 && echo 'nr 注释完成'
}

### cog
# emappey.py youhuma 45分钟，运行完需要 conda deactivate
cog() {
    mkdir ${annotation}
    cd ${annotation} || exit
    source /home/train/miniconda3/bin/activate python27
    echo "正在切换到 python27 conda 环境"
    check_conda_env
    nohup emapper.py -i $specie"_unigene.fasta" --output $specie"_unigene" \
        --cpu 32 --data_dir /home/data/xuezhen/download/eggNOG/eggnog.db \
        -m diamond --translate >cog.log 2>&1 &
    cd ../
}
exec_cog() {
    echo -e "${BL}执行 Annotation - cog 步骤"
    cog
    for ((tc=1;tc<=8640;tc++)); do
        cd ${annotation} || exit
        status=$(ps -ef | grep "emapper.py" | grep $specie | wc -l)
        generate_file_count=$(ls | grep $specie"_unigene.emapper" | wc -l)
        complete=$(cat cog.log | grep "Total time:" | wc -l)
        cd ../
        if [[ $status -eq 1 ]]; then
            for ((i=1;i<=10;i++)); do
                echo -ne "cog is running $i \r"
                sleep 1
            done
        elif [[ $status -eq 0 ]] && [[ $complete -ge 1 ]] && [[ $generate_file_count -ge 1 ]]; then
            echo -e "完成 cog 步骤"
            conda deactivate
            break
        elif [[ $status -eq 0 ]] && [[ $complete -eq 0 ]]; then
            echo -e "${RD}cog 可能出错了，不影响后续执行步骤，可以稍后再查看"
            break
        else
            continue
        fi
    done
}

### kegg
kegg() {
    # TODO: 后面继续写爬虫程序
    # 拿着 $specie_unigene.fasta 去 kegg 网站生成 keg 文件，再做下面的东西 -t plant/animal 动物或植物
    python /home/data/command/kegg/kegg_v5.py -t $4 -i ${all_gene_id}
}

### transdecoder
transdecoder() {
    # 也可以自动生成到报告里
    cd ${annotation} || exit
    mkdir transdecoder
    cd transdecoder || exit
    TransDecoder.LongOrfs -t ${unigene_fasta}
    # 先运行完上面那句，才能运行下面这句
    TransDecoder.Predict -t ${unigene_fasta}
    cd ${work_dir} || exit
}
exec_annotation() {
    mkdir ${annotation}
    cp ${unigene_fasta} ${annotation}
    exec_swiss
    exec_nr
    # exec_cog
    transdecoder
}

### rsem
rsem() {
    mkdir ${reference}
    # 建库
    cd ${reference} || exit
    rsem-prepare-reference --bowtie2 ${unigene_fasta} ${specie}
    cd ${work_dir}; mkdir ${mapping}
    # 读取 samples_described.txt 循环比对
    cat ${work_dir}/samples_described.txt | grep -v sample | while read line; do
        sample=$(echo $line | awk '{print $2}')
        echo "正在执行 $sample 的 rsem 流程"
        rsem-calculate-expression -p ${num_threads} \
        --bowtie2 --strandedness reverse --paired-end \
        ${pinjiedata}/${sample}_1.clean.fq.gz ${pinjiedata}/${sample}_2.clean.fq.gz \
        ${reference}/${specie} ${mapping}/${sample} \
        1>${mapping}/${sample}_mapping.txt \
        2>${mapping}/${sample}.log
    done
    cd ${mapping} || exit
    # 提出 rawreads
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/all_sample_raw_reads.py
    # 提出 fpkm
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/all_sample_fpkm.py
    # 过滤掉 rawreads 小于 50 的 gene
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/count_filtered.py
    # 提出过滤后的基因的 fpkm 值
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/fpkm_filtered.py
    cd ${work_dir} || exit
}

### multi deseq
### 需要创建一个 compare.txt 和 samples_described.txt
multi_deseq() {
    mkdir ${multideseq}
    cd ${multideseq} || exit
    cp ${mapping}/fpkm_matrix_filtered.txt ${multideseq}/
    cp ${mapping}/reads_matrix_filtered.txt ${multideseq}/
    cp ${work_dir}/samples_described.txt ${multideseq}/
    cp ${work_dir}/compare_info.txt ${multideseq}/
    Rscript /home/data/command/Reference_transcriptome/multiple_samples_DESeq2_V7.r 2
    cd ${work_dir}
}

### Funrich
funrich() {
    mkdir ${funrich}
    cd ${funrich} || exit
    cp ${annotation}/*GO*ID.txt ${funrich}/
    cp ${annotation}/*KEGG.txt ${funrich}/
    grep '>' ${unigene_fasta} | sed 's/>//g' > ${funrich}/all_gene_id.txt
    cp ${multideseq}/*_ID.txt ${funrich}/
    cd ${work_dir}
}

### Bubble and Bar graph


### WEGO





show_menu() {
    echo "
    无参执行程序
    物种名: $specie

    0 执行所有流程
    ————————————————
    1 执行 merge
    2 执行 fastqc
    3 执行 pinjie
    4 执行 pinjie result
    5 执行 Annotation
        5.1 执行 swiss
        5.2 执行 nr
        5.3 执行 cog
        5.4 执行 kegg
        5.5 执行 transdecoder
    6 执行 rsem
    7 执行 multi_deseq
    8 复制 funrich 所需文件到 ${funrich} 目录
 "
    # show_status
    read -p "请输入选择(其他任意退出): " select

    case "${select}" in
    0)
        echo "程序开始执行"
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
        exec_fastqc
        ;;
    3)
        exec_pinjie
        ;;
    4)
        pinjie_result
        ;;
    5)
        exec_annotation
        ;;
    5.1)
        exec_swiss
        ;;
    5.2)
        exec_nr
        ;;
    5.3)
        exec_cog
        ;;
    5.4)
        exec_kegg
        ;;
    5.5)
        transdecoder
        ;;
    6)
        rsem
        ;;
    7)
        multi_deseq
        ;;
    8)
        funrich
        ;;
    *)
        # LOGE "请输入正确的数字: "
        exit 0
        ;;
    esac
}
show_menu

